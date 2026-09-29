#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
tutorial_dir="$(cd "${script_dir}/.." && pwd)"
raw_dir="${tutorial_dir}/data/raw"
base_url="https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_X/10k_pbmc_ATACv2_nextgem_Chromium_X"

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

download "filtered_peak_bc_matrix.h5" "filtered_peak_bc_matrix.h5"
download "filtered_tf_bc_matrix.h5" "filtered_tf_bc_matrix.h5"
download "singlecell.csv" "singlecell.csv"
download "peak_annotation.tsv" "peak_annotation.tsv"
download "summary.csv" "summary.csv"
download "summary.json" "summary.json"
download "fragments.tsv.gz" "fragments.tsv.gz"
download "fragments.tsv.gz.tbi" "fragments.tsv.gz.tbi"

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

echo "Downloaded the public 10x PBMC ATAC tutorial inputs to ${raw_dir}"
