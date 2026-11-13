"""Independent two-fermion translation-block eigensolver.

Geometric wavefunctions are not used as eigenstates. The nearest-neighbour
Bloch state fixes the phase only. Defaults reproduce the V=6 reference.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix,csr_matrix,diags
import csv
T=1.0
BETA=0.30
ABS_VERTEX_TOL=1e-9
SCALED_VERTEX_TOL=1e-8
NUMERICAL_TOL=1e-12
POINTS=(('uniform',1/3,0),('zone_boundary',0,1),('wrap',5/6,1/2),('local_node',1/3,-2/3))

def max_abs_sparse(matrix):
    return float(np.max(np.abs(matrix.data))) if matrix.nnz else 0.0

def hop(state,src,dst):
    if not (state>>src)&1 or (state>>dst)&1:
        return None
    sign=(-1)**((state & ((1<<src)-1)).bit_count())
    reduced=state^(1<<src)
    sign*=(-1)**((reduced & ((1<<dst)-1)).bit_count())
    return reduced|(1<<dst),sign

def build_ring(N):
    pairs=[(i,j) for i in range(N) for j in range(i+1,N)]
    states=[(1<<i)|(1<<j) for i,j in pairs]
    index={s:i for i,s in enumerate(states)}
    rows=[]; cols=[]; bonds=[]; signs=[]
    adjacent=[]
    urows=[]; usigns=[]
    for col,(state,pair) in enumerate(zip(states,pairs)):
        adjacent.append(sum(bool((state>>j)&1 and (state>>((j+1)%N))&1) for j in range(N)))
        shifted=[(j+1)%N for j in pair]
        urows.append(index[(1<<shifted[0])|(1<<shifted[1])])
        usigns.append(1 if shifted[0]<shifted[1] else -1)
        for j in range(N):
            k=(j+1)%N
            for src,dst in ((j,k),(k,j)):
                result=hop(state,src,dst)
                if result is not None:
                    dest,sign=result
                    rows.append(index[dest]); cols.append(col)
                    bonds.append(j); signs.append(sign)
    dim=len(states)
    rows=np.array(rows); cols=np.array(cols)
    bonds=np.array(bonds); signs=np.array(signs)
    Hkin=coo_matrix((-T*signs,(rows,cols)),shape=(dim,dim)).tocsr()
    U=coo_matrix((usigns,(urows,np.arange(dim))),shape=(dim,dim)).tocsr()
    seen=set(); orbits=[]
    for seed in range(dim):
        if seed in seen:
            continue
        orbit=[]; accumulated=1; pos=seed
        while pos not in seen:
            seen.add(pos)
            orbit.append((pos,accumulated))
            accumulated*=usigns[pos]
            pos=urows[pos]
        if pos!=seed:
            raise AssertionError('Translation orbit did not close on its seed')
        orbits.append((orbit,accumulated))
    return dict(N=N,dim=dim,pairs=pairs,states=states,index=index,
        rows=rows,cols=cols,bonds=bonds,signs=signs,
        adjacent=np.array(adjacent),Hkin=Hkin,U=U,orbits=orbits)

def canonical_m(m,N):
    return int((m+N//2)%N-N//2)

def momentum_basis(ring,m):
    N=ring['N']; K=2*np.pi*m/N
    rows=[]; cols=[]; data=[]; block_dim=0
    for orbit,cycle_sign in ring['orbits']:
        length=len(orbit)
        if abs(np.exp(1j*K*length)*cycle_sign-1)>1e-10:
            continue
        for step,(index,sign) in enumerate(orbit):
            rows.append(index); cols.append(block_dim)
            data.append(sign*np.exp(1j*K*step)/np.sqrt(length))
        block_dim+=1
    B=coo_matrix((data,(rows,cols)),shape=(ring['dim'],block_dim)).tocsr()
    return B

def bond_anchor(ring,K):
    """Nearest-neighbor Bloch reference only fixes phase; no geometric tail."""
    N=ring['N']; anchor=np.zeros(ring['dim'],dtype=complex)
    for j in range(N):
        k=(j+1)%N
        sign=1 if j<k else -1
        anchor[ring['index'][(1<<j)|(1<<k)]]+=sign*np.exp(1j*K*j)/np.sqrt(N)
    return anchor

def analytic_energy(K,V):
    return -V-4*T*T*np.cos(K/2)**2/V

def analytic_vertex(K,Q,V):
    z=2*T*np.cos(K/2)/V
    zp=2*T*np.cos((K+Q)/2)/V
    if not (abs(z)<1 and abs(zp)<1):
        raise ValueError('Analytic endpoint is not a bound state')
    w=z*zp
    norm=np.sqrt((1-z*z)*(1-zp*zp))
    num=zp*(np.cos((K-Q)/2)-w*np.cos(K/2))
    num+=z*(np.cos(K/2+Q)-w*np.cos((K+Q)/2))
    return -2*T*BETA*norm*num/(1-2*w*np.cos(Q/2)+w*w)

def vertex_operator(ring,Q):
    # Periodic cell gauge, with no sqrt(N) or bath amplitude included.
    values=-T*BETA*ring['signs']*np.exp(1j*Q*ring['bonds'])
    return coo_matrix((values,(ring['rows'],ring['cols'])),shape=(ring['dim'],ring['dim'])).tocsr()

def solve_state(ring,H,m,V):
    N=ring['N']; K=2*np.pi*m/N
    B=momentum_basis(ring,m)
    gram=(B.conj().T@B).toarray()
    Hblock=(B.conj().T@H@B).toarray()
    # Actual eigensolve on the matrix obtained by projection of the full Fock H.
    eigenvalues,eigenvectors=np.linalg.eigh(Hblock)
    coeff=eigenvectors[:,0]
    psi=np.asarray(B@coeff).ravel()
    anchor=bond_anchor(ring,K)
    overlap=np.vdot(anchor,psi)
    if abs(overlap)<1e-8:
        raise AssertionError('Nearest-neighbor phase anchor is singular')
    phase=np.exp(-1j*np.angle(overlap))
    psi*=phase; coeff*=phase
    overlap=np.vdot(anchor,psi)
    energy=float(eigenvalues[0])
    h_res=float(np.linalg.norm(H@psi-energy*psi))
    trans_res=float(np.linalg.norm(ring['U']@psi-np.exp(-1j*K)*psi))
    block_leak=H@B-B@csr_matrix(Hblock)
    record=dict(N=N,V=V,m=m,K=K,fock_dimension=ring['dim'],block_dimension=B.shape[1],
        solver='numpy.linalg.eigh on full-Fock translation block',
        energy=energy,analytic_energy=float(analytic_energy(K,V)),
        energy_abs_error=float(abs(energy-analytic_energy(K,V))),
        first_excited_energy=float(eigenvalues[1]),finite_block_gap=float(eigenvalues[1]-energy),
        continuum_lower_edge=float(-4*T*abs(np.cos(K/2))),
        norm_error=float(abs(np.vdot(psi,psi)-1)),eigen_residual=h_res,
        translation_residual=trans_res,
        block_hermiticity_error=float(np.max(np.abs(Hblock-Hblock.conj().T))),
        basis_orthonormality_error=float(np.max(np.abs(gram-np.eye(B.shape[1])))),
        block_invariance_max_error=max_abs_sparse(block_leak),
        phase_anchor_real=float(overlap.real),phase_anchor_imag=float(overlap.imag),
        geometric_state_used=False)
    record['numerical_pass']=bool(
        record['norm_error']<=NUMERICAL_TOL and h_res<=NUMERICAL_TOL*T
        and trans_res<=NUMERICAL_TOL and record['block_hermiticity_error']<=NUMERICAL_TOL*T
        and record['basis_orthonormality_error']<=NUMERICAL_TOL
        and record['block_invariance_max_error']<=NUMERICAL_TOL*T
        and energy<record['continuum_lower_edge'] and overlap.real>0)
    return psi,record,eigenvalues

def write_csv(path,rows):
    if not rows:
        return
    with path.open('w',encoding='utf-8-sig',newline='') as file:
        writer=csv.DictWriter(file,fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

def dump(path,value):
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2),encoding='utf-8')

def run_validation(Vs,sizes,output):
    output.mkdir(parents=True,exist_ok=False)
    states=[];points=[]
    for N in sizes:
        ring=build_ring(N)
        for V in Vs:
            H=ring['Hkin']+diags(-V*ring['adjacent'],format='csr')
            if max_abs_sparse(H-H.conj().T)>NUMERICAL_TOL or max_abs_sparse(H@ring['U']-ring['U']@H)>NUMERICAL_TOL:
                raise AssertionError('Hamiltonian Hermiticity/translation check failed')
            cache={}
            def get_state(m):
                m=canonical_m(m,N)
                if m not in cache:
                    psi,state,eigs=solve_state(ring,H,m,V)
                    cache[m]=(psi,state);states.append(state)
                return cache[m]
            for name,k_pi,q_pi in POINTS:
                mi=int(round(k_pi*N/2));qi=int(round(q_pi*N/2))
                if abs(mi-k_pi*N/2)>1e-10 or abs(qi-q_pi*N/2)>1e-10:
                    raise AssertionError('Diagnostic momentum is not on the ring grid')
                mf=canonical_m(mi+qi,N)
                pi,si=get_state(mi);pf,sf=get_state(mf)
                K=2*np.pi*mi/N;Q=2*np.pi*qi/N
                D=vertex_operator(ring,Q);Dm=vertex_operator(ring,-Q);Dp=vertex_operator(ring,Q+2*np.pi)
                actual=np.vdot(pf,D@pi);reverse=np.vdot(pi,Dm@pf);periodic=np.vdot(pf,Dp@pi)
                expected=analytic_vertex(K,Q,V);scale=4*abs(BETA)*T*T/V
                absolute=float(abs(actual-expected));scaled=absolute/scale
                herm=float(abs(actual-np.conj(reverse)))/scale
                periodic_error=float(abs(actual-periodic))/scale
                numerical=bool(si['numerical_pass'] and sf['numerical_pass'] and herm<=NUMERICAL_TOL and periodic_error<=NUMERICAL_TOL and max_abs_sparse(D.conj().T-Dm)<=NUMERICAL_TOL)
                comparison=bool(absolute<=ABS_VERTEX_TOL*T and scaled<=SCALED_VERTEX_TOL)
                points.append(dict(N=N,V=V,point=name,g_fock_real=float(actual.real),g_fock_imag=float(actual.imag),g_infinite=float(expected),g_abs_error=absolute,g_scaled_error=scaled,numerical_pass=numerical,vertex_comparison_pass=comparison,point_pass=bool(numerical and comparison)))
    finest=[p for p in points if p['N']==max(sizes)]
    passed=all(s['numerical_pass'] for s in states) and all(p['point_pass'] for p in finest)
    result={'gate_pass':passed,'eigensolves':len(states),'matrix_elements':len(points),'sizes':sizes,'V':Vs,'failures':[p for p in points if not p['point_pass']],'finest_max_abs_error':max(p['g_abs_error'] for p in finest),'note':'The finest-ring gate does not erase coarser finite-size failures; no bath or dynamics claim follows.'}
    write_csv(output/'eigenstates.csv',states);write_csv(output/'matrix_elements.csv',points);dump(output/'summary.json',result)
    print(json.dumps(result,indent=2))
    return 0 if passed else 2

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--all',action='store_true',help='Use all four reported V values; default is V=6.')
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parent/'outputs'/'finite_ring')
    args=parser.parse_args()
    raise SystemExit(run_validation([5,6,8,12] if args.all else [6],[24,36,48],args.output))
