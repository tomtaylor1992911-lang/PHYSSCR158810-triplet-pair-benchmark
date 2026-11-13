"""Plot saved response arrays and comparison errors; no simulation is run."""
from pathlib import Path
import argparse,csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analysis import response

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path);a=p.parse_args()
    with (a.output/'comparisons.csv').open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
    plt.rcParams.update({'font.family':'DejaVu Serif','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    for row in rows:
        case=row['case_id'];N=int(row['N']);tau=float(row['tau_projected']);fig,ax=plt.subplots(figsize=(5.4,3.6),layout='constrained')
        for kernel,color,style,label in [('projected','#0072B2','-','Full projection'),('local_second_order','#D55E00','--','Local, second order')]:
            folder=a.output/f'{case}__{kernel}__N{N}'
            summary=json.loads((folder/'summary.json').read_text(encoding='utf-8'))
            with np.load(folder/'arrays.npz') as z:arrays={k:z[k] for k in z.files}
            c=summary['case'];times=np.linspace(0,5*tau,501);scaled=times*(c['beta']*c['epsilon_zero_T_rms'])**2*c['W']
            ax.plot(scaled,response(arrays,times),color=color,ls=style,label=label)
        ax.set(xlabel=r'$s=(\beta\epsilon)^2Wt/\hbar$',ylabel=r'$R_v(t)$',title=f'{case}, N={N}; time error {100*float(row["eta_tau"]):.2f}%');ax.legend()
        for ext in ['png','pdf','svg']:fig.savefig(a.output/f'{case}_N{N}_response.{ext}',dpi=220)
        plt.close(fig)
    print('Saved standalone response plots using the same physical time axis for both vertices.')
