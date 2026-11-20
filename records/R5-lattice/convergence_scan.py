"""Fixed-density, first-crossed one-axis scans and static RPA charge diagnostic.

Each process owns its output directory. Grid arguments are N:Nw pairs. A warm
lambda path starts at (2,.03), then descends for each phonon frequency; cold
starts and fixed-temperature lambda continuation are checked at every root.
"""
import json
import os
import sys
import time
import numpy as np
from scipy.optimize import brentq
from scipy.special import expit
from scipy.sparse.linalg import LinearOperator, eigs
import lattice

out = os.environ.get('SM_OUT', os.getcwd())
os.makedirs(out, exist_ok=True)
lattice.OUT = os.path.join(out, 'progress.jsonl')
fill = .6


def emit(**kwargs):
    lattice.log(**kwargs)


def rpa_charge(solver, T):
    n, nw = solver.n, solver.nw
    w, u, h, delta = solver.history[-1][1]
    x, _, _, t = lattice.setup(n, nw, solver.lam, solver.om, T)
    xi = x + 1 - fill - delta
    position = .3*n
    lower = int(np.floor(position))
    fraction = position-lower
    f = expit(-xi/t)
    g = 1/(1j*u-(xi[None,:]+h))
    g0 = 1/(1j*w[:,None]-xi[None,:])
    def bubble(shift):
        xq = np.roll(xi, -shift)
        fq = np.roll(f, -shift)
        d = xq-xi
        free = np.empty_like(d)
        mask = np.abs(d)>1e-12
        free[mask] = (f[mask]-fq[mask])/d[mask]
        free[~mask] = f[~mask]*(1-f[~mask])/t
        gq = np.roll(g, -shift, axis=1)
        g0q = np.roll(g0, -shift, axis=1)
        return float(2*np.mean(free)-2*t*np.sum(np.mean((g*gq-g0*g0q).real,axis=1)))
    chi_lower, chi_upper = bubble(lower), bubble(lower+1)
    chi0 = (1-fraction)*chi_lower + fraction*chi_upper
    eigenvalue = 2*solver.lam*chi0
    emit(kind='rpa_charge', n=n, nw=nw, lam=solver.lam, T=T, q_target=.6*np.pi,
         q_lower=2*np.pi*lower/n, q_upper=2*np.pi*(lower+1)/n,
         chi0=chi0, chi0_lower=chi_lower, chi0_upper=chi_upper,
         ladder_eigenvalue=eigenvalue,
         method='static RPA, crossed-normal bubble, linearly interpolated q; exact shifted-free tail; total spin density')
    return eigenvalue


def run_grid(n, nw, resume_omega=None):
    lattice.Solver.cache.clear()
    previous = None
    omega_anchor = None
    for omega in (.03, .06, .12):
        if resume_omega is not None and omega > resume_omega:
            break
        if omega_anchor is not None:
            previous = omega_anchor
        for lam in (2., 1., .5):
            if resume_omega is not None and omega < resume_omega and lam != 2.:
                continue
            if time.time() >= lattice.DEADLINE:
                emit(kind='scan_deadline',n=n,nw=nw,lam=lam,omega=omega)
                return
            # Seed the first lambda at a new omega from the previously converged
            # state at the same grid. Its frequency grid remains valid at its T.
            if previous is not None:
                oldT, state, old_lam = previous
                # Solver accepts only seeds whose tagged temperature is at
                # least the requested T. Its stored Matsubara grid still
                # belongs to oldT, so it shifts u to the new grid itself.
                seedT = max(oldT, 1.8*lattice.C*omega*np.sqrt(lam)/lattice.KB)
                key_lam = old_lam if old_lam != lam else old_lam + 1e-9
                lattice.Solver.cache[(n,nw,omega,True,fill,key_lam)] = [(seedT,state)]
            solver = lattice.Solver(n,nw,lam,omega,True,fill)
            try:
                solver.solve()
                rootT, state = solver.history[-1]
                previous = (rootT, state, lam)
                if lam == 2.:
                    omega_anchor = previous
                if n == 514 and nw == 64 and lam == 2. and omega == .03:
                    rpa_charge(solver, rootT)
                # Independent cold normal solve at the warm onset. If it differs,
                # find its own onset, retaining the discontinuity guard.
                lattice.Solver.cache.clear()
                cold = lattice.Solver(n,nw,lam,omega,True,fill)
                cold_err = cold.evaluate(rootT)
                emit(kind='cold_at_warm_root',n=n,nw=nw,lam=lam,omega=omega,
                     T=rootT,eigenvalue_error=cold_err,
                     mu_shift=float(cold.history[-1][1][3]))
                if abs(cold_err)>1e-4:
                    # Every temperature is a fresh cold start, not a warm scan.
                    def cold_eval(T):
                        lattice.Solver.cache.clear()
                        return lattice.Solver(n,nw,lam,omega,True,fill).evaluate(T)
                    hi=1.8*lattice.C*omega*np.sqrt(lam)/lattice.KB
                    if cold_eval(hi)>=0: raise RuntimeError('cold upper bracket unstable')
                    for _ in range(30):
                        lo=hi/1.25
                        if cold_eval(lo)>0: break
                        hi=lo
                    else: raise RuntimeError('cold lower bracket not found')
                    cr=brentq(cold_eval,lo,hi,xtol=.005,rtol=1e-5)
                    ce=cold_eval(cr)
                    if abs(ce)>1e-4: raise RuntimeError(f'cold branch discontinuity: eigenvalue-1={ce:.6g}')
                    emit(kind='cold_root',n=n,nw=nw,lam=lam,omega=omega,T=cr,eigenvalue_error=ce)
                # Walk lambda from the already established lambda=2 normal
                # solution at fixed T. Where unavailable, mark as unchecked.
                # The next normal state is the previous warm solution, and each
                # incremental solve is independently residual checked.
                if lam == 2.:
                    emit(kind='lambda_check',n=n,nw=nw,lam=lam,omega=omega,T=rootT,
                         eigenvalue_error=0.,status='anchor')
                else:
                    lattice.Solver.cache.clear()
                    if omega_anchor is None:
                        raise RuntimeError('lambda=2 anchor unavailable')
                    anchorT, anchor_state, _ = omega_anchor
                    lattice.Solver.cache[(n,nw,omega,True,fill,2.000000001)] = [
                        (max(anchorT, rootT), anchor_state)]
                    start = lattice.Solver(n,nw,2.,omega,True,fill)
                    start.evaluate(rootT)
                    for l in np.linspace(2.,lam, int(np.ceil((2.-lam)/.25))+1)[1:]:
                        cont = lattice.Solver(n,nw,float(l),omega,True,fill)
                        err=cont.evaluate(rootT)
                    emit(kind='lambda_check',n=n,nw=nw,lam=lam,omega=omega,T=rootT,
                         eigenvalue_error=err,mu_shift=float(cont.history[-1][1][3]),status='completed')
            except Exception as exc:
                emit(kind='scan_failure',n=n,nw=nw,lam=lam,omega=omega,error=repr(exc))
            finally:
                lattice.Solver.cache.clear()


def check_anchor(n, nw, omega, T):
    """Compare the lambda=2 cold state with continuation from lambda=1."""
    lattice.Solver.cache.clear()
    cold = lattice.Solver(n, nw, 2., omega, True, fill)
    cold_error = cold.evaluate(T)
    cold_mu = float(cold.history[-1][1][3])
    lattice.Solver.cache.clear()
    start = lattice.Solver(n, nw, 1., omega, True, fill)
    start.evaluate(T)
    for lam in (1.25, 1.5, 1.75, 2.):
        continued = lattice.Solver(n, nw, lam, omega, True, fill)
        continued_error = continued.evaluate(T)
    emit(kind='lambda_anchor_check', n=n, nw=nw, omega=omega, T=T,
         cold_eigenvalue_error=cold_error,
         continued_eigenvalue_error=continued_error,
         eigenvalue_difference=continued_error-cold_error,
         mu_difference=float(continued.history[-1][1][3])-cold_mu)


def state_metrics(solver, T):
    w, u, h, delta = solver.history[-1][1]
    x, _, _, t = lattice.setup(solver.n, solver.nw, solver.lam, solver.om, T)
    density = lattice.density(u, h, x+1-fill, w, t, solver.nw, delta)
    inverse = (-np.arange(solver.n)) % solver.n
    return dict(density=float(density), mu_shift=float(delta),
                inversion_u=float(np.max(np.abs(u-u[:,inverse]))),
                inversion_h=float(np.max(np.abs(h-h[:,inverse]))),
                frequency_u=float(np.max(np.abs(u+u[::-1]))),
                frequency_h=float(np.max(np.abs(h-h[::-1]))),
                min_positive_u_over_w=float(np.min(u[solver.nw:]/w[solver.nw:,None])),
                h_span=float(np.max(h)-np.min(h)))


def pairing_symmetry(solver, T):
    """Parity of the unrestricted leading linearized pairing eigenvector."""
    _, u, h, delta = solver.history[-1][1]
    x, _, v, t = lattice.setup(solver.n,solver.nw,solver.lam,solver.om,T)
    apply=lattice.gap(u,h,x+1-fill-delta,v,t)
    shape=u.shape
    op=LinearOperator((u.size,u.size),
                      matvec=lambda y:apply(y.reshape(shape)).ravel(),dtype=float)
    values, vectors=eigs(op,k=1,which='LR',tol=2e-5,maxiter=100,
                         v0=np.ones(u.size))
    vec=vectors[:,0].reshape(shape)
    inverse=(-np.arange(solver.n)) % solver.n
    norm=np.linalg.norm(vec)
    return dict(eigenvalue_real=float(values[0].real),
                eigenvalue_imag=float(values[0].imag),
                inversion_even_defect=float(np.linalg.norm(vec-vec[:,inverse])/norm),
                inversion_odd_defect=float(np.linalg.norm(vec+vec[:,inverse])/norm),
                frequency_even_defect=float(np.linalg.norm(vec-vec[::-1])/norm),
                frequency_odd_defect=float(np.linalg.norm(vec+vec[::-1])/norm))


def branch_profile(n, nw, T, symmetry=False):
    """Compare descending lambda continuation with cold normal states."""
    omega = .03
    lattice.Solver.cache.clear()
    anchor = lattice.Solver(n,nw,2.,omega,True,fill)
    # Reproduce the warm lambda=2 anchor before lowering T, as in run_grid.
    anchor.evaluate(78.79365349128435 if nw == 64 else 78.82136095105994)
    anchor.evaluate(T)
    for lam in (1.75,1.5,1.25,1.):
        cont = lattice.Solver(n,nw,lam,omega,True,fill)
        ce = cont.evaluate(T)
        cm = state_metrics(cont,T)
        _, cu, ch, _ = cont.history[-1][1]
        saved_cache = dict(lattice.Solver.cache)
        lattice.Solver.cache.clear()
        cold = lattice.Solver(n,nw,lam,omega,True,fill)
        try:
            ke = cold.evaluate(T)
            km = state_metrics(cold,T)
            _, ku, kh, _ = cold.history[-1][1]
            diff = max(float(np.max(np.abs(cu-ku))),
                       float(np.max(np.abs(ch-kh))))
            emit(kind='branch_profile',n=n,nw=nw,lam=lam,omega=omega,T=T,
                 continued_eigenvalue_error=ce,cold_eigenvalue_error=ke,
                 continued=cm,cold=km,max_state_difference=diff)
            if symmetry and lam == 1.:
                emit(kind='pairing_symmetry',n=n,nw=nw,lam=lam,omega=omega,T=T,
                     continued=pairing_symmetry(cont,T),cold=pairing_symmetry(cold,T))
        except Exception as exc:
            emit(kind='branch_profile_failure',n=n,nw=nw,lam=lam,omega=omega,T=T,
                 continued_eigenvalue_error=ce,continued=cm,error=repr(exc))
        finally:
            lattice.Solver.cache.clear()
            lattice.Solver.cache.update(saved_cache)


def continuation_root(n, nw):
    """Temperature root of the lambda=2 to 1 normal continuation branch."""
    omega, lam = .03, 1.
    anchor = lattice.Solver(n,nw,2.,omega,True,fill)
    def f(T):
        lattice.Solver.cache.clear()
        anchor.evaluate(T)
        for coupling in (1.75,1.5,1.25,1.):
            state = lattice.Solver(n,nw,coupling,omega,True,fill)
            err = state.evaluate(T)
        emit(kind='continuation_probe',n=n,nw=nw,lam=lam,omega=omega,
             T=T,eigenvalue_error=err,state=state_metrics(state,T))
        return err
    try:
        ceiling=lattice.C*omega/lattice.KB
        hi=1.8*ceiling
        fhi=f(hi)
        if fhi>=0:
            raise RuntimeError(f'upper bracket unstable: eigenvalue-1={fhi:.6g}')
        at_ceiling=f(ceiling)
        lo=ceiling
        flo=at_ceiling
        for _ in range(20):
            if flo>0:
                break
            hi, fhi = lo, flo
            lo /= 1.25
            flo=f(lo)
        else:
            raise RuntimeError('unstable lower bracket not found')
        root=brentq(f,lo,hi,xtol=.005,rtol=1e-5)
        err=f(root)
        if abs(err)>1e-4:
            raise RuntimeError(f'continuation branch discontinuity: eigenvalue-1={err:.6g}, T={root:.9g}')
        emit(kind='continuation_root',n=n,nw=nw,lam=lam,omega=omega,T=root,
             eigenvalue_error=err,ceiling=ceiling,ratio=root/ceiling,
             eigenvalue_error_at_ceiling=at_ceiling)
    except Exception as exc:
        emit(kind='continuation_failure',n=n,nw=nw,lam=lam,omega=omega,error=repr(exc))


def charge_onset(n, nw, lam):
    """Track the homogeneous normal branch upward and bracket the RPA pole."""
    Tpair={2.:78.793653491,1.:40.84878815096565}[lam]
    omega=.03
    try:
        lattice.Solver.cache.clear()
        solver=lattice.Solver(n,nw,lam,omega,True,fill)
        pair_error=solver.evaluate(Tpair)
        value=rpa_charge(solver,Tpair)
        emit(kind='charge_onset_pair',n=n,nw=nw,lam=lam,T=Tpair,
             pair_eigenvalue_error=pair_error,charge_eigenvalue=value)
        if value <= 1:
            raise RuntimeError(f'charge ladder below pole at pairing reference: {value:.9g}')
        lo=Tpair
        seed=solver.history[-1][1]
        def probe(T, initial):
            lattice.Solver.cache.clear()
            state=lattice.Solver(n,nw,lam,omega,True,fill)
            state.history=[(T,initial)]
            state.evaluate(T)
            return rpa_charge(state,T),state.history[-1][1]
        for _ in range(20):
            hi=2*lo
            high_value,high_state=probe(hi,seed)
            if high_value<=1:
                break
            lo,seed=hi,high_state
        else:
            raise RuntimeError('charge pole not bracketed')
        root=brentq(lambda T:probe(T,seed)[0]-1,lo,hi,xtol=.05,rtol=1e-5)
        residual,_=probe(root,seed)
        if abs(residual-1)>1e-4:
            raise RuntimeError(f'charge branch discontinuity: eigenvalue-1={residual-1:.6g}, T={root:.9g}')
        emit(kind='charge_onset',n=n,nw=nw,lam=lam,T=root,T_pair=Tpair,
             ratio=root/Tpair,charge_eigenvalue=residual)
    except Exception as exc:
        emit(kind='charge_onset_failure',n=n,nw=nw,lam=lam,error=repr(exc))


def high_branch_root(n, nw):
    """Heat the large-eigenvalue normal state, then guard its pairing root."""
    T0=40.84878815096565
    omega=.03
    try:
        lattice.Solver.cache.clear()
        anchor=lattice.Solver(n,nw,2.,omega,True,fill)
        anchor.evaluate(78.82136095105994)
        anchor.evaluate(T0)
        for lam in (1.75,1.5,1.25,1.):
            state=lattice.Solver(n,nw,lam,omega,True,fill)
            err=state.evaluate(T0)
        emit(kind='high_branch_probe',n=n,nw=nw,T=T0,eigenvalue_error=err,
             state=state_metrics(state,T0))
        if err < 1.:
            raise RuntimeError(f'high branch not reproduced: eigenvalue-1={err:.6g}')
        ceiling=lattice.C*omega/lattice.KB
        previous_T=T0
        previous_state=state.history[-1][1]
        def evaluate_from(T, seed):
            lattice.Solver.cache.clear()
            solver=lattice.Solver(n,nw,1.,omega,True,fill)
            solver.history=[(T,seed)]
            value=solver.evaluate(T)
            emit(kind='high_branch_probe',n=n,nw=nw,T=T,eigenvalue_error=value,
                 state=state_metrics(solver,T))
            return value, solver.history[-1][1]
        while previous_T < 1.8*ceiling:
            T=min(previous_T+3.,1.8*ceiling)
            if previous_T < ceiling < T:
                T=ceiling
            value, new_state=evaluate_from(T,previous_state)
            if value <= 0.:
                seed=previous_state
                root=brentq(lambda t:evaluate_from(t,seed)[0],
                            previous_T,T,xtol=.005,rtol=1e-5)
                residual,_=evaluate_from(root,seed)
                if abs(residual)>1e-4:
                    raise RuntimeError(f'high branch discontinuity: eigenvalue-1={residual:.6g}, T={root:.9g}')
                emit(kind='high_branch_root',n=n,nw=nw,T=root,
                     eigenvalue_error=residual,ceiling=ceiling,ratio=root/ceiling)
                return
            previous_T,previous_state=T,new_state
        raise RuntimeError('no pairing crossing below 1.8 times ceiling')
    except Exception as exc:
        emit(kind='high_branch_failure',n=n,nw=nw,error=repr(exc))


if __name__ == '__main__':
    if sys.argv[1] == 'chargeonset':
        n,nw=int(sys.argv[2]),int(sys.argv[3])
        for arg in sys.argv[4:]:
            charge_onset(n,nw,float(arg))
    elif sys.argv[1] == 'charge':
        n, nw, T = int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
        solver = lattice.Solver(n,nw,2.,.03,True,fill)
        error = solver.evaluate(T)
        emit(kind='charge_pairing_check',n=n,nw=nw,T=T,eigenvalue_error=error)
        rpa_charge(solver, T)
    elif sys.argv[1] == 'anchorcheck':
        for spec in sys.argv[2:]:
            n, nw, omega, T = spec.split(':')
            check_anchor(int(n), int(nw), float(omega), float(T))
    elif sys.argv[1] in ('profile', 'symprofile'):
        for spec in sys.argv[2:]:
            n, nw, T = spec.split(':')
            branch_profile(int(n), int(nw), float(T), sys.argv[1]=='symprofile')
    elif sys.argv[1] == 'branchroot':
        n, nw = map(int, sys.argv[2:4])
        continuation_root(n,nw)
    elif sys.argv[1] == 'highbranch':
        n, nw = map(int, sys.argv[2:4])
        high_branch_root(n,nw)
    elif sys.argv[1] == 'resume':
        n, nw = map(int, sys.argv[2:4])
        run_grid(n,nw,resume_omega=float(sys.argv[4]))
    else:
        for arg in sys.argv[1:]:
            n,nw=map(int,arg.split(':'))
            run_grid(n,nw)
