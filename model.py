"""Frozen step-8 electronic model. No output or dynamics at import."""
import numpy as np

def energy(K, V, t=1.0):
    return -V - 4*t*t/V*np.cos(np.asarray(K)/2)**2

def g_projected(K, Q, V, beta=.30, t=1.0):
    """Exact bound-state projection, compact V9 formula, periodic cell gauge."""
    K, Q = np.broadcast_arrays(np.asarray(K, float), np.asarray(Q, float))
    c = 2*t/V
    C, Cp = np.cos(K/2), np.cos((K+Q)/2)
    s = np.cos(Q) + np.cos(Q/2)*np.cos(K+Q/2)
    n = np.sqrt((1-c*c*C*C)*(1-c*c*Cp*Cp))
    den = 1-2*c*c*C*Cp*np.cos(Q/2)+c**4*C*C*Cp*Cp
    return -4*t*t*beta/V*n*(s-2*c*c*C*C*Cp*Cp)/den

def g_local(K, Q, V, beta=.30, t=1.0):
    return -4*t*t*beta/V*(np.cos(Q)+np.cos(Q/2)*np.cos(K+Q/2))

def g_series(K, Q, V, terms=32, beta=.30, t=1.0):
    """Paired real-space geometric sum, not used to build rates."""
    K, Q = np.broadcast_arrays(np.asarray(K, float), np.asarray(Q, float))
    z, zp = 2*t/V*np.cos(K/2), 2*t/V*np.cos((K+Q)/2)
    w, n = z*zp, np.sqrt((1-z*z)*(1-zp*zp))
    value = np.zeros_like(K)
    for r in range(1, terms+1):
        value += w**(r-1)*(zp*np.cos(K/2-Q*r/2)+z*np.cos((K+Q)/2+Q*r/2))
    tail = 2*abs(t*beta)*n*(np.abs(z)+np.abs(zp))*np.abs(w)**terms/(1-np.abs(w))
    return -2*t*beta*n*value, tail

def grid(N):
    return -np.pi+2*np.pi*np.arange(N)/N

def vertex_matrix(K, V, kernel, beta=.30):
    initial = K[None, :]
    transfer = K[:, None]-initial
    fn = {"projected": g_projected, "local_second_order": g_local}[kernel]
    return fn(initial, transfer, V, beta)
