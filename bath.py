"""Frozen open internal-bond coordinate bath, natural units t0=a=hbar=kB=1.

u_j=epsilon X_j, integral f(E)dE=1, hence epsilon is zero-temperature
RMS displacement. spectrum() contains 2*pi*epsilon**2, but NO electronic
coupling beta, t0, vertex, or Fourier measure 1/N.

Shape-only quadrature and FFT correlation integrals are computed lazily once
per (theta,omega_ratio,gamma_ratio,cutoff_ratio), then cached in memory by
functools.lru_cache. No constructor/import writes files. Default tau_B uses
N_positive=8192, zero-padding factor=128, dimensionless time window=1600;
audit_bath.py verifies frequency/time-grid convergence and the omitted tail.
For epsilon=0, tau_B denotes its continuous epsilon->0 shape limit.
"""
from functools import lru_cache
import numpy as np
from scipy.integrate import quad, simpson
from scipy.fft import fft
from scipy.optimize import brentq, minimize_scalar

__all__ = ["Bath"]

class _Shape:
    def __init__(self, theta, omega, gamma, cutoff):
        self.theta, self.omega, self.gamma, self.cutoff = theta, omega, gamma, cutoff
        self.raw_norm, self.norm_error = quad(self.raw, 0., cutoff, epsabs=2e-13,
                                             epsrel=2e-13, points=[omega] if omega<cutoff else None)
        self.slope0 = gamma / omega**4 / self.raw_norm
        self.c0, self.c0_error = quad(self.thermal_weight, 0., cutoff, epsabs=2e-13,
                                    epsrel=2e-13, points=[omega] if omega<cutoff else None)
        self.normalization = quad(self.density,0.,cutoff,epsabs=2e-13,epsrel=2e-13)[0]
        self.peak_x, self.peak_value, self.stationary_points = self.optimize_peak()
        self.tau_shape, self.memory_metadata = self.memory_integral(8192,128,1600.)

    def raw(self, x):
        a=np.asarray(x,dtype=float); y=np.zeros_like(a)
        m=(a>0)&(a<self.cutoff); z=a[m]
        y[m]=self.gamma*z/((z*z-self.omega*self.omega)**2+self.gamma**2*z*z)*(1-(z/self.cutoff)**2)**4
        return float(y) if y.ndim==0 else y

    def density(self,x):
        return self.raw(x)/self.raw_norm

    def positive_weight(self,x):
        a=np.asarray(x,dtype=float); y=np.zeros_like(a)
        m=(a>0)&(a<self.cutoff); z=a[m]
        y[m]=self.density(z)/(-np.expm1(-z/self.theta))
        y[a==0]=self.theta*self.slope0
        return float(y) if y.ndim==0 else y

    def thermal_weight(self,x):
        a=np.asarray(x,dtype=float)
        p=self.positive_weight(a)
        y=p*(1+np.exp(-a/self.theta))
        return float(y) if y.ndim==0 else y

    def log_derivative_positive(self,x):
        x=float(x); q=x/self.theta
        if q<1e-3:
            thermal=1/(2*self.theta)-x/(12*self.theta**2)+x**3/(720*self.theta**4)-x**5/(30240*self.theta**6)
        else:
            thermal=1/x-1/(self.theta*np.expm1(q))
        d=(x*x-self.omega*self.omega)**2+self.gamma*self.gamma*x*x
        dp=4*x*(x*x-self.omega*self.omega)+2*self.gamma*self.gamma*x
        return thermal-dp/d-8*x/(self.cutoff*self.cutoff-x*x)

    def optimize_peak(self):
        grid=np.linspace(self.cutoff*1e-10,self.cutoff*(1-1e-10),2049)
        der=np.array([self.log_derivative_positive(x) for x in grid])
        roots=[]
        for a,b,da,db in zip(grid[:-1],grid[1:],der[:-1],der[1:]):
            if da*db<0:
                roots.append(brentq(self.log_derivative_positive,a,b,xtol=5e-15,rtol=1e-14))
        candidates=[0.,self.cutoff]+roots
        # Separate bounded optimization provides an additional optimizer check.
        opt=minimize_scalar(lambda x:-self.positive_weight(x),bounds=(0,self.cutoff),method='bounded',options={'xatol':1e-14})
        candidates.append(float(opt.x))
        values=np.array([self.positive_weight(x) for x in candidates]); k=int(np.argmax(values))
        self.bounded_peak={'x':float(opt.x),'value':float(-opt.fun),'success':bool(opt.success)}
        return candidates[k],float(values[k]),roots

    def memory_curve(self,n_positive=8192,pad_factor=128,max_time=1600.):
        n_positive=int(n_positive);pad_factor=int(pad_factor)
        dx=self.cutoff/n_positive
        nfft=2*n_positive*pad_factor
        x=np.arange(n_positive+1)*dx
        plus=self.positive_weight(x);minus=plus*np.exp(-x/self.theta)
        arr=np.zeros(nfft,dtype=float)
        arr[:n_positive+1]=plus
        arr[-n_positive:]=minus[1:][::-1]
        corr=fft(arr)*dx
        dt=2*np.pi/(nfft*dx)
        last=int(np.floor(max_time/dt))
        if last>=nfft//2: raise ValueError('Correlation window reaches FFT recurrence/alias region.')
        times=np.arange(last+1)*dt
        cc=corr[:last+1].copy()
        if times[-1]<max_time:
            fraction=(max_time-times[-1])/dt
            cc=np.append(cc,corr[last]*(1-fraction)+corr[last+1]*fraction)
            times=np.append(times,max_time)
        return times,cc,{'n_positive':n_positive,'pad_factor':pad_factor,'nfft':nfft,'dx':dx,'dt':dt,'max_time':max_time,'recurrence_period':2*np.pi/dx}

    def memory_integral(self,n_positive=8192,pad_factor=128,max_time=1600.):
        times,cc,meta=self.memory_curve(n_positive,pad_factor,max_time)
        val=simpson(np.abs(cc)/self.c0,x=times)
        meta.update(tau_shape_truncated=float(val),correlation_at_zero_fft=float(cc[0].real),final_abs_normalized=float(abs(cc[-1])/self.c0))
        return float(val),meta

    def correlation_quad(self,s):
        re,er=quad(self.thermal_weight,0,self.cutoff,weight='cos',wvar=float(s),epsabs=2e-12,epsrel=2e-12,limit=200)
        im,ei=quad(self.density,0,self.cutoff,weight='sin',wvar=float(s),epsabs=2e-12,epsrel=2e-12,limit=200)
        return re-1j*im,er+ei

@lru_cache(maxsize=16)
def _get_shape(theta,omega,gamma,cutoff):
    return _Shape(theta,omega,gamma,cutoff)

class Bath:
    """spectrum(delta_E) returns unsymmetrized S_u(delta_E); positive E is energy accepted by bath."""
    def __init__(self,W,epsilon=.002,theta=.25,omega_ratio=.5,gamma_ratio=.25,cutoff_ratio=1.1):
        values=[W,epsilon,theta,omega_ratio,gamma_ratio,cutoff_ratio]
        if not np.all(np.isfinite(values)) or W<=0 or epsilon<0 or min(theta,omega_ratio,gamma_ratio,cutoff_ratio)<=0:
            raise ValueError('W, theta and spectral ratios must be positive finite; epsilon nonnegative finite.')
        self.W=float(W);self.epsilon=float(epsilon);self.theta=float(theta)
        self.omega_ratio=float(omega_ratio);self.gamma_ratio=float(gamma_ratio);self.cutoff_ratio=float(cutoff_ratio)
        self.T=self.theta*self.W;self.Omega=self.omega_ratio*self.W;self.gamma=self.gamma_ratio*self.W;self.Ec=self.cutoff_ratio*self.W
        self._shape=_get_shape(self.theta,self.omega_ratio,self.gamma_ratio,self.cutoff_ratio)
        self.tau_B=self._shape.tau_shape/self.W
        self.u_rms=self.epsilon*np.sqrt(self._shape.c0)
        self.smax=2*np.pi*self.epsilon**2/self.W*self._shape.peak_value
        self.smax_energy=self._shape.peak_x*self.W
        self.normalization=self._shape.normalization
        self.raw_normalization=self._shape.raw_norm
        self.memory_cache_source='In-memory lru_cache of direct shape quadrature/FFT; no disk lookup; see memory_metadata.'
        self.memory_metadata=dict(self._shape.memory_metadata)
        self.s0=2*np.pi*self.epsilon**2/self.W*self.theta*self._shape.slope0

    def spectrum(self,delta_E):
        a=np.asarray(delta_E,dtype=float)
        if not np.all(np.isfinite(a)):raise ValueError('delta_E must be finite.')
        x=np.abs(a)/self.W
        values=2*np.pi*self.epsilon**2/self.W*self._shape.positive_weight(x)
        values=np.where(a<0,values*np.exp(-x/self.theta),values)
        return float(values) if values.ndim==0 else values

    def density(self,E):
        """Positive-frequency f(E), normalized to one; no epsilon or 2*pi."""
        return self._shape.density(np.asarray(E,dtype=float)/self.W)/self.W

    def correlation(self,time):
        """Independent oscillatory quadrature of C_u(time); intended for a few check points."""
        a=np.asarray(time,dtype=float)
        vals=np.array([self._shape.correlation_quad(float(s)*self.W)[0] for s in a.ravel()])*self.epsilon**2
        vals=vals.reshape(a.shape)
        return complex(vals) if vals.ndim==0 else vals
