"""Portable entry to the validated population solver; no computation at import.

The numerical solve body is retained from the archived implementation.
Only workspace release-file checks and source-snapshot I/O are externalized.
"""
import csv,json
import numpy as np
from scipy.linalg import eigh
from scipy.sparse.csgraph import connected_components
from model import energy,grid,vertex_matrix
from bath import Bath

def assert_gate(condition,message):
    if not condition:
        raise RuntimeError(message)

def write_json(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+"\n",encoding="utf-8")

def solve(case,N,kernel,run):
    bath=Bath(case["W"],epsilon=case["epsilon_zero_T_rms"],
              theta=case["T_over_W"],omega_ratio=case["Omega_over_W"],
              gamma_ratio=case["gamma_over_W"],cutoff_ratio=case["Ec_over_W"])
    K=grid(N); E=energy(K,case["V"])
    g=vertex_matrix(K,case["V"],kernel,case["beta"])
    delta_E=E[None,:]-E[:,None]
    spectrum=np.asarray(bath.spectrum(delta_E),float)
    rates=g*g*spectrum/N
    np.fill_diagonal(rates,0.)
    assert_gate(np.isfinite(rates).all(),"Nonfinite rates")
    rate_min=float(rates.min())
    assert_gate(rate_min>=0,"Negative rate: no clipping allowed")
    Gamma=rates.sum(axis=0)
    rate_scale=float(Gamma.max())
    assert_gate(rate_scale>0,"No population process: tau is nonfinite; stop finite-tau workflow")
    L=rates.copy()
    np.fill_diagonal(L,-Gamma)
    pi=np.exp(-(E-E.min())/case["T"]);pi/=pi.sum()
    v=(case["W"]/2)*np.sin(K)
    p0=pi*(1+.20*np.sin(K))
    rootpi=np.sqrt(pi)
    A=L*rootpi[None,:]/rootpi[:,None]
    flux=L*pi[None,:]
    checks={
        "offdiagonal_min": rate_min,
        "Gamma_max": rate_scale, "Gamma_min":float(Gamma.min()),
        "column_sum_abs":float(np.max(np.abs(L.sum(axis=0)))),
        "column_sum_scaled":float(np.max(np.abs(L.sum(axis=0)))/rate_scale),
        "stationary_abs":float(np.max(np.abs(L@pi))),
        "stationary_scaled":float(np.max(np.abs(L@pi))/(rate_scale*pi.max())),
        "detailed_balance_abs":float(np.max(np.abs(flux-flux.T))),
        "detailed_balance_scaled":float(np.max(np.abs(flux-flux.T))/(rate_scale*pi.max())),
        "A_symmetry_scaled":float(np.max(np.abs(A-A.T))/rate_scale),
        "vertex_Hermiticity_scaled":float(np.max(np.abs(g-g.T))/(4*abs(case["beta"])/case["V"])),
        "initial_mass_error":float(abs(p0.sum()-1)),
        "initial_probability_min":float(p0.min()),
        "bath_tau_B":float(bath.tau_B),"bath_smax":float(bath.smax),
        "bath_u_rms":float(bath.u_rms),"bath_normalization":float(bath.normalization),
        "Gamma_max_tau_B":float(rate_scale*bath.tau_B),
        "global_gap_lower_bound":case["global_gap_lower_bound"],
        "cutoff_margin":case["cutoff_margin"],
        "rms_margin_ratio":float(4*abs(case["beta"])*bath.u_rms/case["cutoff_margin"]),
    }
    write_json(run/"generator_checks_pre_eigensolve.json",checks)
    for key in ["column_sum_scaled","stationary_scaled","detailed_balance_scaled",
                "A_symmetry_scaled","vertex_Hermiticity_scaled"]:
        assert_gate(checks[key]<=1e-12,f"{key} exceeds frozen tolerance")
    assert_gate(checks["Gamma_max_tau_B"]<.01,"Weak-memory scale gate failed")
    assert_gate(checks["rms_margin_ratio"]<.05,"RMS/margin scale gate failed")
    assert_gate(case["cutoff_margin"]>0,"Single-quantum band protection failed")
    # Every strictly positive edge is kept. No arbitrary rate pruning.
    ncomp,labels=connected_components(rates>0,directed=True,connection="strong")
    checks["strong_components"]=int(ncomp)
    scaled_values,U=eigh(A/rate_scale,check_finite=True,driver="evr")
    zero_tol=100*N*np.finfo(float).eps
    is_zero=np.abs(scaled_values)<=zero_tol
    assert_gate(float(scaled_values.max())<=zero_tol,"Positive generator eigenvalue")
    assert_gate(int(is_zero.sum())==ncomp,"Spectral zero-mode count and graph differ; investigate before tau")
    q=v*rootpi; normq=float(q@q)
    c=U.T@q
    weights=c*c/normq
    zero_weight=float(weights[is_zero].sum())
    checks.update({
        "zero_mode_relative_tolerance":float(zero_tol),
        "zero_mode_count":int(is_zero.sum()),
        "velocity_zero_mode_weight":zero_weight,
        "eigen_residual_scaled":float(np.linalg.norm((A/rate_scale)@U-U*scaled_values[None,:],ord=np.inf)),
        "eigenvector_orthogonality_inf":float(np.linalg.norm(U.T@U-np.eye(N),ord=np.inf)),
        "largest_raw_eigenvalue":float(scaled_values.max()*rate_scale),
        "smallest_nonzero_decay_over_Gamma":float(-scaled_values[~is_zero].max()),
    })
    assert_gate(zero_weight<=1e-20,"Velocity overlaps a nondecaying mode; integral tau is nonfinite")
    assert_gate(checks["eigen_residual_scaled"]<=1e-10,"Eigensystem residual too large")
    # The unique numerical stationary eigenvalue is set to exactly zero only
    # after the conservation, graph, and zero-overlap checks above. L is unchanged.
    values=scaled_values*rate_scale
    values[is_zero]=0.
    nonzero=~is_zero
    decay=-values[nonzero]
    tau=float(np.sum(weights[nonzero]/decay))
    assert_gate(np.isfinite(tau) and tau>0,"Invalid integral relaxation time")
    def response(times):
        return weights[nonzero]@np.exp(values[nonzero,None]*np.atleast_1d(times)[None,:])
    times=np.linspace(0,5*tau,2001)
    R=response(times)
    ptimes=tau*np.array([0,.01,.1,.5,1,2,5,10,20])
    dc=U.T@((p0-pi)/rootpi)
    populations=pi[:,None]+rootpi[:,None]*(U@(dc[:,None]*np.exp(values[:,None]*ptimes[None,:])))
    pop_R=(v@populations)/(v@p0)
    checks.update({
        "population_mass_max_error":float(np.max(np.abs(populations.sum(axis=0)-1))),
        "population_min":float(populations.min()),
        "population_response_max_abs_error":float(np.max(np.abs(pop_R-response(ptimes)))),
        "R_at_zero_error":float(abs(R[0]-1)),
        "R_nonincreasing_max_violation":float(max(0,np.diff(R).max())),
        "tau_tail_beyond_5tau_fraction":float(np.sum(weights[nonzero]*np.exp(-decay*5*tau)/decay)/tau),
        "tau_tail_beyond_20tau_fraction":float(np.sum(weights[nonzero]*np.exp(-decay*20*tau)/decay)/tau),
    })
    assert_gate(checks["population_mass_max_error"]<=1e-10,"Population mass drift")
    assert_gate(checks["population_min"]>=-1e-12,"Negative population")
    assert_gate(checks["population_response_max_abs_error"]<=1e-10,"Population/correlation mismatch")
    summary={
        "case_id":case["case_id"],"V":case["V"],"N":N,"kernel":kernel,
        "pass":True,"tau_v":tau,
        "tau_v_in_common_scaled_units":tau*(case["beta"]*case["epsilon_zero_T_rms"])**2*case["W"],
        "checks":checks,"case":case,
        "physical_scope":"single-band second-order Born-Markov population dynamics of the frozen effective coordinate bath",
        "not_validation_of_full_nonMarkov_dynamics":True,
    }
    np.savez_compressed(run/"arrays.npz",K=K,E=E,pi=pi,v=v,p0=p0,
                        g=g,spectrum=spectrum,Gamma=Gamma,L=L,A=A,
                        eigenvalues=values,eigenvalues_raw=scaled_values*rate_scale,
                        eigenvectors=U,zero_mask=is_zero,weights=weights,
                        R_times=times,R_values=R,p_times=ptimes,p_values=populations)
    with (run/"response.csv").open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f);w.writerow(["case_id","kernel","N","time_hbar_over_t0","R_v"])
        w.writerows((case["case_id"],kernel,N,float(t),float(r)) for t,r in zip(times,R))
    write_json(run/"generator_checks.json",checks)
    write_json(run/"summary.json",summary)
    return summary
