from first import *
from scipy import stats
from tqdm import tqdm

def g(x):
    return np.log(1 + x**2)

def d2_g(x):
    return 2 * (1 - x**2) / ((1 + x**2)**2) 
 
templates = {
    "normal": lambda n,h1,h2: stats.norm(
                                        loc = -h1 / (h2 + n**0.5),
                                        scale= 1 / (1 + h2/(n**0.5)) 
                                        ),

    "cauchy": lambda n,h1,h2: stats.cauchy(
                                        loc = -h1 / (h2 + n**0.5),
                                        scale= 1 / (1 + h2/(n**0.5))
                                        ),
    "student": lambda n,h1,h2: stats.t(
                                    loc = -h1 / (h2 + n**0.5),
                                    scale= 1 / (1 + h2/(n**0.5))
                                    )
}


def experiment_step(g, d2_g, template,
                   h1=0, h2=2.1, alpha=0.05,
                   N=1000, M=17, sample_sizes=[100, 400, 900],
                   integrals1=None, crit_vals=None):
    """
    params:
    ---
    ``h1`` - shift diff (h1 >= 0)

    ``h2`` - scale diff

    ``alpha`` - level of significance 

    ``N`` - number of repeats computing empirical power

    ``M`` - number of permutations for computing critical value

    ``sample_sizes`` - sizes of samples to compute empirical power

    return:
    ---------
    ``empirical_powers`` and ``asymptotic_power``
    """
    if integrals1 is None or crit_vals is None:
        raise Exception("integrals1 of crit_vals is None")

    if h1 < 0:
        raise Exception("h1 must be geq than 0")

    d1 = template(1, 0, 0)
    f = d1.pdf

    print("computing integrals2...")
    integrals = compute_integrals2(f, d2_g, h1, h2, integrals=integrals1)

    b1 = np.sqrt( np.abs( integrals["J1_star"] ) )
    b2 = np.sqrt( np.abs( integrals["J2_star"] ) )
    
    a = (integrals["J2"] + integrals["J1"] ** 2 \
         - 2 * integrals["J3"]) ** 0.25
    b = np.sqrt( b1 ** 2 + b2 ** 2)

    emp_powers = []
    asp_power = compute_asymptotic_power(alpha, b, a)

    print("\ncomputing empirical_powers...")
    for n in sample_sizes:
        print(f"n={n}")
        d2_n = template(n, h1, h2)
        print("computing e_pow for n*T_n...")
        e_pow = compute_empirical_power(n, N, crit_vals[n], d1, d2_n, g)
        emp_powers.append(e_pow)

    return emp_powers, asp_power



def run_experiment(g, d2_g, template, h1=0,
                   h2_vals=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5],
                   alpha=0.05, N=1000, M=700,
                   sample_sizes=[100, 400, 900, 1600]):
    
    print("precomputing integrals1...")
    d1 = template(1, 0, 0)
    f = d1.pdf
    integrals1 = compute_integrals1(f, g)

    print("-"*100)

    print("precomputing critical_values for different n...")
    crit_vals = dict()
    for n in tqdm(sample_sizes):
        print(f"n={n}")
        crit_val = compute_crit_val(n, M, alpha, template, g)
        crit_vals[n] = crit_val

    print("-"*100)

    for h2 in h2_vals:
        print('experiment step with h2 =', h2, '\n')
        template = templates["normal"]
        e_pow, a_pow = experiment_step(g=g, d2_g=d2_g, template=template,
                                    h1=h1, h2=h2, alpha=alpha, N=N,
                                    M = M, sample_sizes=sample_sizes,
                                    integrals1=integrals1,
                                    crit_vals=crit_vals)
        print('\n')
        print(e_pow)
        print(a_pow)
        print("="*100)
        print("="*100)
        print()




def test1():
    template = templates["cauchy"]
    h1_vals = np.linspace(0,5,3)
    h2_vals = np.linspace(0,5,3)
    n_vals = [100,400,900]
    x_vals = [-1, -0.5, 0, 0.5, 1]
    d1 = template(1,0,0)
    for h1 in h1_vals:
        for h2 in h2_vals:
            for n in n_vals:
                d2 = template(n, h1, h2)
                for x in x_vals:
                    print( abs( d1.cdf( x*(1+h2/n**0.5) + h1/n**0.5 ) - d2.cdf(x) ) )

def test2():
    template = templates["cauchy"]
    d1 = template(1,0,0)

    n = 100
    M = 700
    alpha = 0.05
    N = 1000

    c = compute_crit_val(n, M, alpha, template, g)

    cnt = 0

    for _ in range(N):
        X = d1.rvs(n)
        Y = d1.rvs(n)
        cnt += int(n * compute_etest(g,X,Y) >= c)

    print(cnt / N)

def test3():
    template = templates["normal"]
    d1 = template(1,0,0)
    n=900
    h1=0
    
    N = 10000
    M = 10000
    cv = compute_crit_val(n, M, 0.05, template, g)


    for h2 in range(1,11):
        d2_n = template(n,h1,h2)
        e_pow = compute_empirical_power(n, N, cv, d1, d2_n, g)
        print(h2, e_pow)
    

def test4():
    X = [14, -10, 1, 2, 3, 4, 4, 2, 2, 3]
    Y = [1, 20, 9, 9, 9, 9, 1, 3, 8, 8]
    print(compute_etest2(g, X, Y))

def test5():
    X = [ 2, 2, 3, 11, 12, 17, 6, 20, 14]
    print(np.quantile(X, 0.39))

def test6():
    print(compute_asymptotic_power(0.14, 0.7, 0.1))


def test7():
    n = 200
    M = 700
    alpha = 0.01
    template = templates["normal"]
    print(compute_crit_val(n, M, alpha, template, g))

def test8():
    n = 100
    N = 10000
    cv = 51
    template = templates["normal"]
    d1 = stats.norm(0, 1)
    d2 = stats.norm(0.2, 5)
    ep = compute_empirical_power(n, N, cv, d1, d2, g)
    print(ep)

def test():
    # test1()
    # test2()
    # test3()
    # test4()
    # test5()
    # test6()
    # test7()
    test8()




def main():
    template = templates["cauchy"]

    run_experiment(g, d2_g, template, h1=0,
                   h2_vals=range(1, 11),
                   alpha=0.05, N=1000, M=700,
                   sample_sizes=[100, 400, 900, 1600, 2000])
    


# test()
# main()

f = templates["student"](1, 0, 0).pdf
integrals = compute_integrals(f, g, d2_g)
print(integrals)

# from scipy.stats import anderson
# import numpy as np

# # Пример использования
# data = [1, 2, 3, 4, 5, 5, 5, 5]
# result = anderson(data, dist='norm')
# print(f"Статистика: {result.statistic}")
# print(f"Критические значения: {result.critical_values}")
# print(f"Уровни значимости: {result.significance_level}")