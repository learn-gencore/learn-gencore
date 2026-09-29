#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
tutorial_dir="$(cd "${script_dir}/.." && pwd)"
raw_dir="${tutorial_dir}/data/raw"
processed_dir="${tutorial_dir}/data/processed"
singlecell_csv="${raw_dir}/singlecell.csv"
fragments_file="${raw_dir}/fragments.tsv.gz"
barcodes_file="${processed_dir}/cell_barcodes.tsv"
metrics_file="${processed_dir}/fragment_metrics.tsv"
histogram_file="${processed_dir}/fragment_length_histogram.tsv"

mkdir -p "${processed_dir}"

if [[ ! -s "${singlecell_csv}" || ! -s "${fragments_file}" ]]; then
    echo "Missing singlecell.csv or fragments.tsv.gz. Run download_data.sh first." >&2
    exit 1
fi

awk -F',' 'NR > 1 && $10 == 1 {print $1}' "${singlecell_csv}" > "${barcodes_file}"

gzip -dc "${fragments_file}" | awk \
    -F'\t' \
    -v metrics_file="${metrics_file}" \
    -v histogram_file="${histogram_file}" '
    BEGIN { OFS = "\t" }
    NR == FNR {
        keep[$1] = 1
        next
    }
    /^#/ { next }
    ($4 in keep) {
        length_bp = $3 - $2
        total[$4] += 1
        if (length_bp < 147) {
            nfr[$4] += 1
        } else if (length_bp < 294) {
            mono[$4] += 1
        }
        if (length_bp >= 0 && length_bp <= 1000) {
            length_hist[length_bp] += 1
        }
    }
    END {
        print "barcode\tn_fragments\tnucleosome_free_fragments\tmononucleosomal_fragments\tnucleosome_signal" > metrics_file
        for (barcode in keep) {
            n_total = total[barcode] + 0
            n_nfr = nfr[barcode] + 0
            n_mono = mono[barcode] + 0
            ratio = n_nfr > 0 ? n_mono / n_nfr : "NA"
            print barcode, n_total, n_nfr, n_mono, ratio > metrics_file
        }

        print "fragment_length\tn_fragments" > histogram_file
        for (length_bp = 0; length_bp <= 1000; length_bp++) {
            print length_bp, length_hist[length_bp] + 0 > histogram_file
        }
    }
' "${barcodes_file}" -

echo "Wrote ${metrics_file}"
echo "Wrote ${histogram_file}"
