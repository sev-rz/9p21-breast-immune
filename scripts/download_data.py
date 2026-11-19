"""Download the raw files listed in metadata/sources.tsv (GDC) and verify checksums.

Usage:  python scripts/download_data.py
Files already present with the right checksum are skipped. The cBioPortal files
(cbp_clinical.csv, cbp_patient.csv) are fetched by scripts/download_cbioportal.py.
"""
import csv, hashlib, os, subprocess, sys

def digest(path, algo):
    h = hashlib.new(algo)
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()

sha = {}
for line in open('metadata/checksums.sha256'):
    if line.strip():
        value, path = line.split(); sha[path] = value

bad = 0
for row in csv.DictReader(open('metadata/sources.tsv'), delimiter='\t'):
    url, dest = row['location'], row['local_file']
    if 'api.gdc.cancer.gov' not in url:
        continue  # cBioPortal studies are handled by download_cbioportal.py
    algo, want = ('sha256', sha.get(dest)) if dest in sha else ('md5', row['md5'])
    if os.path.exists(dest) and want not in (None, 'NA') and digest(dest, algo) == want:
        print('ok (already there)', dest); continue
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    print('downloading', dest, flush=True)
    # curl: the GDC resets connections from Python's default HTTP client
    subprocess.run(['curl', '-L', '--fail', '--silent', '--show-error', '--retry', '6', '--retry-all-errors', '--retry-delay', '5',
                    '-o', dest + '.part', url], check=True)
    os.replace(dest + '.part', dest)
    got = digest(dest, algo)
    if want not in (None, 'NA') and got != want:
        print(f'CHECKSUM MISMATCH {dest}: expected {want}, got {got}'); bad += 1
    else:
        print('ok', dest, got)
sys.exit(1 if bad else 0)
