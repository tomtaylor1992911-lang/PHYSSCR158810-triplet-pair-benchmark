"""Reproduce the declared single-pair comparison; default is V=6, N=96."""
from pathlib import Path
import argparse,csv,json,sys,platform
import numpy as np
import scipy
from solver import solve,write_json
from analysis import compare

HERE=Path(__file__).resolve().parent

def csvwrite(path,rows):
    if rows:
        with path.open('w',encoding='utf-8',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--suite',choices=['pilot','main','all'],default='pilot',help='pilot: two V6/N96 generators; main: four V at 96/192/384; all: main plus the two declared controls.')
    p.add_argument('--output',type=Path,default=None,help='New output directory; existing output is never overwritten.')
    args=p.parse_args()
    cases=json.loads((HERE/'parameters.json').read_text(encoding='utf-8'))['cases']
    if args.suite=='pilot':jobs=[(next(c for c in cases if c['case_id']=='MAIN_V6'),96)]
    else:jobs=[(c,N) for c in cases if args.suite=='all' or c['role']=='main' for N in c['grids']]
    output=args.output or HERE/'outputs'/args.suite
    output.mkdir(parents=True,exist_ok=False)
    write_json(output/'environment.json',{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'suite':args.suite,'comparison_count':len(jobs),'note':'Fresh reproduction; floating-point last digits may vary with BLAS/platform. No absolute machine paths are recorded.'})
    comparisons=[];budgets=[];identities=[]
    try:
        for case,N in jobs:
            loaded={}
            for kernel in ['projected','local_second_order']:
                run=output/f"{case['case_id']}__{kernel}__N{N}";run.mkdir()
                summary=solve(case,N,kernel,run)
                with np.load(run/'arrays.npz') as z:arrays={k:z[k] for k in z.files}
                loaded[kernel]=(summary,arrays)
            row,bins,identity=compare(case['case_id'],N,*loaded['projected'],*loaded['local_second_order'])
            comparisons.append(row);budgets.extend(bins);identities.append(identity)
            csvwrite(output/'comparisons.csv',comparisons);csvwrite(output/'signed_Q_budget.csv',budgets)
            write_json(output/'signed_Q_identities.json',identities)
            print(json.dumps({'case':case['case_id'],'N':N,'tau_projected':row['tau_projected'],'tau_local':row['tau_local'],'eta_tau':row['eta_tau'],'eta_R':row['eta_R']}),flush=True)
    except Exception as exc:
        write_json(output/'failure.json',{'type':type(exc).__name__,'message':str(exc),'completed_comparisons':len(comparisons)})
        raise
    return 0

if __name__=='__main__':raise SystemExit(main())
