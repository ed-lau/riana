"""
Build ground_truth.csv for a calibration cell line.

The ground-truth signal is the nominal D2O mixing proportion per LC-MS run.
The SDRF samplesheet's `source name` column encodes it directly as
`<LINE>_<PROPORTION>` (e.g. `AC16_12.5`, `iPSC_87.5`), so parsing is just a
suffix split. The script joins that against `comment[file uri]` (JPOST URL)
and `comment[data file]` (raw filename) for provenance.

Output:  tests/data/calibration_d2o_mixing/<line>/ground_truth.csv
Columns: source_name, raw_filename, mzml_filename, riana_filename, nominal_proportion, raw_uri
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
LINES = {
    'ac16': REPO_ROOT / 'data' / 'calibration_ac16' / 'samplesheet_ac16_alpine.sdrf.tsv',
    'ipsc': REPO_ROOT / 'data' / 'calibration_ipsc' / 'samplesheet_ipsc_alpine.sdrf.tsv',
}
OUT_DIR = REPO_ROOT / 'tests' / 'data' / 'calibration_d2o_mixing'


def parse_sdrf(sdrf_path: Path) -> list[dict]:
    with sdrf_path.open() as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = list(reader)
    out = []
    for row in rows:
        source = row['source name']
        proportion = float(source.rsplit('_', 1)[1])
        raw_filename = row['comment[data file]']
        mzml_filename = raw_filename.rsplit('.raw', 1)[0] + '.mzML.gz'
        riana_filename = f'time{proportion:g}_riana.txt'
        out.append({
            'source_name': source,
            'raw_filename': raw_filename,
            'mzml_filename': mzml_filename,
            'riana_filename': riana_filename,
            'nominal_proportion': proportion,
            'raw_uri': row['comment[file uri]'],
        })
    out.sort(key=lambda r: r['nominal_proportion'])
    return out


def write_csv(rows: list[dict], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ['source_name', 'raw_filename', 'mzml_filename',
                  'riana_filename', 'nominal_proportion', 'raw_uri']
    with out_path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--line', choices=['ac16', 'ipsc', 'all'], default='all')
    args = parser.parse_args()

    targets = ['ac16', 'ipsc'] if args.line == 'all' else [args.line]
    for line in targets:
        sdrf = LINES[line]
        if not sdrf.exists():
            raise FileNotFoundError(f'SDRF not found: {sdrf}')
        rows = parse_sdrf(sdrf)
        out_path = OUT_DIR / line / 'ground_truth.csv'
        write_csv(rows, out_path)
        print(f'[{line}] wrote {len(rows)} rows to {out_path.relative_to(REPO_ROOT)}')


if __name__ == '__main__':
    main()
