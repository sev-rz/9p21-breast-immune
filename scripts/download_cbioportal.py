"""Download the two cBioPortal tables used by the pipeline and check their checksums.

Writes data/raw/cbp_clinical.csv (per-sample clinical values) and
data/raw/cbp_patient.csv (per-patient values, including the PAM50 SUBTYPE) for the
TCGA breast PanCanAtlas study. Same queries as the early helpers cbp.py and cbp2.py.
"""
import hashlib, json, os, sys, time, urllib.request
import pandas as pd

BASE, STUDY = 'https://www.cbioportal.org/api', 'brca_tcga_pan_can_atlas_2018'

def get(kind):
    url = f'{BASE}/studies/{STUDY}/clinical-data?clinicalDataType={kind}&projection=SUMMARY&pageSize=100000'
    for k in range(4):
        try:
            req = urllib.request.Request(url, headers={'Accept': 'application/json'})
            return pd.DataFrame(json.load(urllib.request.urlopen(req, timeout=300)))
        except Exception as e:
            print(f'  attempt {k + 1} failed: {e}', flush=True)
            if k == 3: raise
            time.sleep(5 * (k + 1))

os.makedirs('data/raw', exist_ok=True)
get('SAMPLE').pivot_table(index='sampleId', columns='clinicalAttributeId', values='value', aggfunc='first').to_csv('data/raw/cbp_clinical.csv')
get('PATIENT').pivot_table(index='patientId', columns='clinicalAttributeId', values='value', aggfunc='first').to_csv('data/raw/cbp_patient.csv')

want = {p: v for v, p in (l.split() for l in open('metadata/checksums.sha256') if l.strip())}
bad = 0
for path in ('data/raw/cbp_clinical.csv', 'data/raw/cbp_patient.csv'):
    got = hashlib.sha256(open(path, 'rb').read()).hexdigest()
    ok = got == want.get(path)
    print('ok' if ok else 'CHECKSUM DIFFERS (the live API may have changed since 2026-10-16)', path, got)
    bad += not ok
sys.exit(1 if bad else 0)
