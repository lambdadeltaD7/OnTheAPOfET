#include "boost/random/cauchy_distribution.hpp"
#include "boost/random/normal_distribution.hpp"
#include <boost/random.hpp>
#include <boost/math/distributions.hpp>
#include <iostream>
#include "first.h"
#include <random>
#include <chrono>
#include "omp.h"

boost::random::mt19937 gen(1337);

double g(double x){
    return log(1 + x*x);
}



void ex_tmp_st(std::vector<double> h1_vals, double h2, double alpha, int N, int M,
            int df, std::vector<int>sample_sizes, std::vector<double> integrals)
{
    std::cout<< "ET_" << "Student" << "\n";
    run_experiment_st(h1_vals, h2, alpha, N, M, df, sample_sizes, integrals,
         [](double* X, double *Y, int n ){return compute_etest(g, X, Y, n);}, true, gen);
    std::cout << "================================\n";
    std::cout << "================================\n\n";

    // std::cout<< "HT_" << name << "\n";
    // run_experiment(h1_vals, h2, alpha, N, M, sample_sizes, get_dist, integrals,
    // [](double* X, double *Y, int n ){return compute_etest(g2, X, Y, n);}, false, gen);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "WMW_" << name << "\n";
    // run_experiment(h1_vals, h2, alpha, N, M, sample_sizes, get_dist, integrals, compute_wmw, false, gen);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "AD_" << name << "\n";
    // run_experiment(h1_vals, h2, alpha, N, M, sample_sizes, get_dist, integrals, compute_ad, false, gen);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "KS_" << name << "\n";
    // run_experiment(h1_vals, h2, alpha, N, M, sample_sizes, get_dist, integrals, compute_ks,  false, gen);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "CM_" << name << "\n";
    // run_experiment(h1_vals, h2, alpha, N, M, sample_sizes, get_dist, integrals, compute_cm,  false, gen);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";
}


template <typename T>
void ex_tmp(std::vector<double> h1_vals, double h2, double alpha, int N, int M,
            std::vector<int>sample_sizes, std::vector<double> integrals, T get_dist, const char name[])
{
    // std::cout<< "ET_" << name << "\n";
    // run_experiment(h1_vals, h2, alpha, N, M, sample_sizes, get_dist, integrals,
    //      [](double* X, double *Y, int n ){return compute_etest(g, X, Y, n);}, true, gen);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "HT_" << name << "\n";
    // run_experiment(h1_vals, h2, alpha, N, M, sample_sizes, get_dist, integrals,
    // [](double* X, double *Y, int n ){return compute_etest(g2, X, Y, n);}, false, gen);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "WMW_" << name << "\n";
    // run_experiment(h1_vals, h2, alpha, N, M, sample_sizes, get_dist, integrals, compute_wmw, false, gen);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "AD_" << name << "\n";
    // run_experiment(h1_vals, h2, alpha, N, M, sample_sizes, get_dist, integrals, compute_ad, false, gen);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "KS_" << name << "\n";
    // run_experiment(h1_vals, h2, alpha, N, M, sample_sizes, get_dist, integrals, compute_ks,  false, gen);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "CM_" << name << "\n";
    // run_experiment(h1_vals, h2, alpha, N, M, sample_sizes, get_dist, integrals, compute_cm,  false, gen);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";
}


template <typename T>
void ex_tmp(double h1, std::vector<double> h2_vals,  double alpha, int N, int M,
            std::vector<int>sample_sizes, std::vector<double> integrals, T get_dist, const char name[])
{

    std::cout<< "ET_" << name << "\n";
    run_experiment(  h1, h2_vals, alpha, N, M, sample_sizes, get_dist, integrals,
         [](double* X, double *Y, int n ){return compute_etest(g, X, Y, n);}, true, gen);
    std::cout << "================================\n";
    std::cout << "================================\n\n";

    // std::cout<< "HT_" << name << "\n";
    // run_experiment(h1, h2_vals, alpha, N, M, sample_sizes, get_dist, integrals,
    // [](double* X, double *Y, int n ){return compute_etest(g2, X, Y, n);}, false, gen);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "WMW_" << name << "\n";
    // run_experiment(h1, h2_vals, alpha, N, M, sample_sizes, get_dist, integrals, compute_wmw, false, gen);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "AD_" << name << "\n";
    // run_experiment(h1, h2_vals, alpha, N, M, sample_sizes, get_dist, integrals, compute_ad, false, gen);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "KS_" << name << "\n";
    // run_experiment(h1, h2_vals, alpha, N, M, sample_sizes, get_dist, integrals, compute_ks,  false, gen);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";

    // std::cout<< "CM_" << name << "\n";
    // run_experiment(h1, h2_vals, alpha, N, M, sample_sizes, get_dist, integrals, compute_cm,  false, gen);
    // std::cout << "================================\n";
    // std::cout << "================================\n\n";
}


// normal with h1=0
void ex1(){
    std::vector<double>  h2_vals{0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5};
    double h1 = 0.0;
    double alpha = 0.05;
    std::vector<int> sample_sizes{100};
    int N = 5000;
    int M = 5000;
    std::vector<double> integrals = read_integrals(DistributionType::NORMAL);
    ex_tmp(h1, h2_vals, alpha, N, M, sample_sizes, integrals, get_normal, "NORMAL");
}


// normal with h2=0
void ex2(){
    std::vector<double>  h1_vals{1, 2, 3, 4, 5};
    double h2 = 0.0;
    double alpha = 0.05;
    std::vector<int> sample_sizes{30, 70, 100, 200, 300};
    int N = 500;
    int M = 500;
    std::vector<double> integrals = read_integrals(DistributionType::NORMAL);
    ex_tmp(h1_vals, h2, alpha, N, M, sample_sizes, integrals, get_normal, "NORMAL");
}


// cauchy with h2=0
void ex3(){
    std::vector<double>  h1_vals{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    double h2 = 0.0;
    double alpha = 0.05;
    std::vector<int> sample_sizes{100, 400, 900};
    int N = 5000;
    int M = 5000;
    std::vector<double> integrals = read_integrals(DistributionType::CAUCHY);
    ex_tmp(h1_vals, h2, alpha, N, M, sample_sizes, integrals, get_cauchy, "CAUCHY");
}


// cauchy with h1=0
void ex4(){
    std::vector<double>  h2_vals{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    double h1 = 0.0;
    double alpha = 0.05;
    std::vector<int> sample_sizes{100, 400, 900};
    int N = 5000;
    int M = 5000;
    std::vector<double> integrals = read_integrals(DistributionType::CAUCHY);
    ex_tmp(h1, h2_vals, alpha, N, M, sample_sizes, integrals, get_cauchy, "CAUCHY");
}


// student with h2=0
void ex5(){
    std::vector<double>  h1_vals{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    double h2 = 0.0;
    double alpha = 0.05;
    std::vector<int> sample_sizes{100, 400, 900};
    int N = 5000;
    int M = 5000;
    int df = 10;
    std::vector<double> integrals = read_integrals(DistributionType::CAUCHY);
    ex_tmp_st(h1_vals, h2, alpha, N, M, df, sample_sizes, integrals);
}

int main() {
    ex1();
    // ex2();
    // ex3();
    // ex4();
}