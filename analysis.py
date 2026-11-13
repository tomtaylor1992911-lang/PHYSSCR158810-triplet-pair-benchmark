"""Saved-generator analysis from the corrected final v2 implementation.

The error extremum search is numerical refinement, not a rigorous enclosure.
The signed-Q budget is an algebraic identity, not a channel-removal experiment.
"""
import numpy as np
from scipy.optimize import brentq

def response(z,t):
    mask=~z["zero_mask"]
    return z["weights"][mask]@np.exp(z["eigenvalues"][mask,None]*np.atleast_1d(t)[None,:])

def response_derivative(z,x,tau):
    mask=~z["zero_mask"]; lam=z["eigenvalues"][mask]*tau
    return float(np.sum(z["weights"][mask]*lam*np.exp(lam*x)))

def curve_error(zp,zl,tau):
    # Locate derivative sign changes, then refine extrema; retained grid checks
    # document the resolution of this smooth finite modal sum.
    xs=np.linspace(0,5,2001)
    deriv=np.array([response_derivative(zl,x,tau)-response_derivative(zp,x,tau) for x in xs])
    candidates=[0.,5.]
    for a,b,da,db in zip(xs[:-1],xs[1:],deriv[:-1],deriv[1:]):
        if da*db<0:
            candidates.append(brentq(lambda x:response_derivative(zl,x,tau)-response_derivative(zp,x,tau),a,b,xtol=1e-13))
    delta=response(zl,np.array(candidates)*tau)-response(zp,np.array(candidates)*tau)
    k=int(np.argmax(np.abs(delta)))
    dense=np.linspace(0,5,4001)
    gridmax=float(np.max(np.abs(response(zl,dense*tau)-response(zp,dense*tau))))
    return float(abs(delta[k])),float(candidates[k]),gridmax,len(candidates)-2

def potentials(z):
    q=z["v"]*np.sqrt(z["pi"]);mask=~z["zero_mask"]
    y=z["eigenvectors"][:,mask]@((z["eigenvectors"][:,mask].T@q)/(-z["eigenvalues"][mask]))
    return y/np.sqrt(z["pi"])

def compare(case,N,sp,zp,sl,zl):
    tp,tl=sp["tau_v"],sl["tau_v"]
    er,peak,gridmax,extrema=curve_error(zp,zl,tp)
    et=abs(tl/tp-1)
    cat="usable" if et<=.05 and er<=.05 else ("resolved_difference" if et>=.10 or er>=.10 else "gray")
    row={"case_id":case,"V":sp["V"],"N":N,"tau_projected":tp,"tau_local":tl,
         "tau_projected_scaled":sp["tau_v_in_common_scaled_units"],
         "tau_local_scaled":sl["tau_v_in_common_scaled_units"],
         "signed_tau_local_over_projected_minus1":tl/tp-1,"eta_tau":et,"eta_R":er,
         "R_error_peak_time_over_tau_projected":peak,
         "R_error_fine_grid_lower_estimate":gridmax,
         "R_error_optimizer_minus_grid":er-gridmax,"interior_extrema_found":extrema,
         "classification_before_grid_gate":cat,
         "Gamma_tauB_projected":sp["checks"]["Gamma_max_tau_B"],
         "Gamma_tauB_local":sl["checks"]["Gamma_max_tau_B"]}
    hp,hl=potentials(zp),potentials(zl)
    dW=zl["L"]-zp["L"];np.fill_diagonal(dW,0.)
    norm=float(np.sum(zp["pi"]*zp["v"]**2))
    contrib=-.5*dW*zp["pi"][None,:]*(hl[:,None]-hl[None,:])*(hp[:,None]-hp[None,:])/norm
    total=float(contrib.sum());expected=tl-tp
    # Use exact mesh indices for half-open bins so floating differences at
    # pi/3 and 2pi/3 cannot split one momentum-transfer shell between bins.
    directed=(np.arange(N)[:,None]-np.arange(N)[None,:])%N
    distance=np.minimum(directed,N-directed)
    masks=[6*distance<N,(6*distance>=N)&(3*distance<N),3*distance>=N]
    bins=[0,np.pi/3,2*np.pi/3,np.pi];rows=[]
    for b,(lo,hi) in enumerate(zip(bins[:-1],bins[1:])):
        mask=masks[b]
        rows.append({"case_id":case,"N":N,"Q_bin":b+1,"Q_min_over_pi":lo/np.pi,
                     "Q_max_over_pi":min(1.,hi/np.pi),"signed_delta_tau_contribution":float(contrib[mask].sum()),
                     "contribution_over_tau_projected":float(contrib[mask].sum()/tp)})
    attr={"case_id":case,"N":N,"delta_tau_actual":expected,"resolvent_sum":total,
          "identity_relative_residual":abs(total-expected)/max(abs(expected),tp*1e-15),
          "h_projected_backward_residual":float(np.max(np.abs(-zp["L"].T@hp-zp["v"])))/max(abs(zp["v"])),
          "h_local_backward_residual":float(np.max(np.abs(-zl["L"].T@hl-zl["v"])))/max(abs(zl["v"])),
          "interpretation":"Signed resolvent decomposition of existing generators; no new mechanism or new physical case."}
    return row,rows,attr
