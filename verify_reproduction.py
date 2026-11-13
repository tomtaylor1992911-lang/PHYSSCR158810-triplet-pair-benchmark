"""Check fresh output against preserved final-v2 numerical references."""
from pathlib import Path
import argparse,csv,json,math
HERE=Path(__file__).resolve().parent
KEYS=['tau_projected','tau_local','eta_tau','eta_R','tau_projected_scaled','tau_local_scaled','Gamma_tauB_projected','Gamma_tauB_local']

def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path);a=p.parse_args()
    refs={(r['case_id'],int(r['N'])):r for r in read(HERE/'reference_comparisons.csv')}
    checks=[]
    environment=json.loads((a.output/'environment.json').read_text(encoding='utf-8'))
    for row in read(a.output/'comparisons.csv'):
        key=(row['case_id'],int(row['N']));ref=refs[key]
        errors={k:abs(float(row[k])-float(ref[k]))/max(abs(float(ref[k])),1e-30) for k in KEYS}
        ok=all(math.isclose(float(row[k]),float(ref[k]),rel_tol=1e-9,abs_tol=1e-12) for k in KEYS)
        checks.append({'case':key[0],'N':key[1],'pass':ok,'relative_differences':errors})
    complete=len(checks)==environment['comparison_count'] and not (a.output/'failure.json').exists()
    result={'pass':bool(checks) and complete and all(x['pass'] for x in checks),'complete':complete,'expected_comparisons':environment['comparison_count'],'tolerance':{'relative':1e-9,'absolute':1e-12},'scope':'Floating-point reproduction of the reported population implementation, not validation of exact open-system dynamics.','checks':checks}
    (a.output/'reference_verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2));raise SystemExit(0 if result['pass'] else 1)
