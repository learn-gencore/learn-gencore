#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
tutorial_dir="$(cd "${script_dir}/.." && pwd)"
raw_dir="${tutorial_dir}/data/raw"
base_url="https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_unsorted_3k/pbmc_unsorted_3k"

mkdir -p "${raw_dir}"

download() {
    local suffix="$1"
    local output_name="$2"
    if [[ -s "${raw_dir}/${output_name}" ]]; then
        echo "Already present: ${output_name}"
        return
    fi
    curl -L --fail --retry 3 \
        -o "${raw_dir}/${output_name}" \
        "${base_url}_${suffix}"
}

download "filtered_feature_bc_matrix.h5" "filtered_feature_bc_matrix.h5"
download "per_barcode_metrics.csv" "per_barcode_metrics.csv"
download "atac_peaks.bed" "atac_peaks.bed"
download "atac_peak_annotation.tsv" "atac_peak_annotation.tsv"
download "summary.csv" "summary.csv"
download "atac_fragments.tsv.gz" "atac_fragments.tsv.gz"
download "atac_fragments.tsv.gz.tbi" "atac_fragments.tsv.gz.tbi"

if command -v sha256sum >/dev/null 2>&1; then
    (cd "${raw_dir}" && sha256sum -c "${tutorial_dir}/checksums.sha256")
elif command -v shasum >/dev/null 2>&1; then
    (cd "${raw_dir}" && shasum -a 256 -c "${tutorial_dir}/checksums.sha256")
else
    python - "${tutorial_dir}/checksums.sha256" "${raw_dir}" <<'PY'
import hashlib
import sys
from pathlib import Path

checksum_path = Path(sys.argv[1])
raw_dir = Path(sys.argv[2])
failed = False
for line in checksum_path.read_text().splitlines():
    if not line.strip():
        continue
    expected, filename = line.split(maxsplit=1)
    path = raw_dir / filename.strip()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    actual = digest.hexdigest()
    status = "OK" if actual == expected else "FAILED"
    print(f"{filename}: {status}")
    failed |= actual != expected
raise SystemExit(1 if failed else 0)
PY
fi

echo "Downloaded the public 10x PBMC multiome tutorial inputs to ${raw_dir}"
