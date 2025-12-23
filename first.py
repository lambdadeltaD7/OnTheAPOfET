import numpy as np
from scipy.integrate import nquad
from scipy import stats
import time


def compute_integrals(f, g, d2_g):
    integrals = dict()

    print("computing J1...")
    integrals["J1"] = nquad(
        lambda x,y: g(x - y) * f(x) * f(y),
        ranges=[(-np.inf, +np.inf),
                (-np.inf, +np.inf)])[0]
    print("J1 =", integrals["J1"])
    
    print("computing J2...")
    integrals["J2"] = nquad(
        lambda x,y: (g(x - y) ** 2) * f(x) * f(y),
        ranges=[(-np.inf, +np.inf),
                (-np.inf, +np.inf)])[0]
    print("J2 =", integrals["J2"])
    
    t0= time.time()
    R = 30
    print("computing J3...")
    integrals["J3"] = nquad(
        lambda x,y,z: g(x - y) * g(x - z) * f(x) * f(y) * f(z),
        ranges=[(-R, +R),
                (-R, +R),
                (-R, +R)], opts={"epsabs":1e-3, "epsrel":1e-3})[0]
    # integrals["J3"] = 0.763368 # for normal
    # integrals["J3"] = 6.881056 # for cauchy
    print("J3 =", integrals["J3"])
    print("elapsed time for J3:", time.time()-t0)
    
    
    print(f"R = {R}")
    print()
    print("computing J1_star...")
    integrals["J1_star"] = 0.5 * nquad(
        lambda x,y: d2_g(x - y) * f(x) * f(y),
        ranges=[(-R, +R),
                (-R, +R)], opts={"epsabs":1e-4, "epsrel":1e-4})[0]
    print("J1_star =", integrals["J1_star"])

    print("computing J2_star...")
    integrals["J2_star"] = 0.5 * nquad(
        lambda x,y: (y**2 - 0.5 * (x-y)**2) * d2_g(x - y) * f(x) * f(y),
        ranges=[(-R, +R),
                (-R, +R)], opts={"epsabs":1e-4, "epsrel":1e-4})[0]
    print("J2_star =", integrals["J2_star"])


    return integrals

def g(x):
    return abs(x)

def d2_g(x):
    return 2 * (1 - x**2) / (1 + x**2)**2


integrals = compute_integrals(lambda x: stats.t.pdf(x,df=7), g, d2_g)
print("integrals for Student distribution")
print(integrals)



