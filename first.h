#include "boost/random/normal_distribution.hpp"
#include <boost/random.hpp>
#include <iostream>
#include <random>
#include <omp.h>
#include <fstream>


#include <chrono>


class Timer {
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    
public:
    Timer() : start(std::chrono::high_resolution_clock::now()) {}
    
    ~Timer() {
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << "Время выполнения: " << duration.count()/ 1000000.0 << " сек" << std::endl;
    }
};


template <typename T>
void sample(T dist, int sample_size, double *X, boost::random::mt19937& gen){
    for(int i=0; i<sample_size; ++i){
        X[i] = dist(gen);
    }
}


void print_vector(std::vector<double> V);


//
void print_sample(double *sample, int sample_size);

double compute_edf(std::vector<double>& ordered_sample, double t);


double compute_cm(double *X, double *Y, int sample_size);

double compute_etest(double (*g)(double), double *X, double *Y, int sample_size);
double compute_etest2(double (*g)(double), double *X, double *Y, int sample_size);


double compute_wmw(double *X, double *Y, int sample_size);

double compute_ad(double *X, double *Y, int sample_size);

double compute_ks(double *X, double *Y, int sample_size);
//

double compute_asymptotic_power(double alpha, double b, double a);


template <typename T>
double compute_empirical_power(int n, int N, double crit_val, T d1, T d2,
                               std::function<double(double*, double*, int)> compute_test,
                               boost::random::mt19937& gen)
{
    double cnt_reject = 0;

    double *X = new double[n * N];
    double *Y = new double[n * N];

    for(int i=0; i<N; ++i){
        sample(d1, n, X + (i * n), gen);
        sample(d2, n, Y + (i * n), gen);
    }

    #pragma omp parallel for reduction(+:cnt_reject)
    for(int i=0; i<N; ++i){
        cnt_reject += ( compute_test(X + (i * n), Y + (i * n), n) >= crit_val );
    }

    delete[] X;
    delete[] Y;
    
    return cnt_reject / N;
}


//
template<typename T>
double quantile(const std::vector<T>& data, double probability) {
    if (data.empty()) return 0.0;
    if (probability < 0.0 || probability > 1.0) {
        throw std::invalid_argument("Probability must be in [0, 1]");
    }
    
    std::vector<T> sorted_data = data;
    std::sort(sorted_data.begin(), sorted_data.end());
    
    if (probability == 0.0) return sorted_data.front();
    if (probability == 1.0) return sorted_data.back();
    
    double index = probability * (sorted_data.size() - 1);
    size_t lower_index = static_cast<size_t>(std::floor(index));
    size_t upper_index = static_cast<size_t>(std::ceil(index));
    
    if (lower_index == upper_index) {
        return sorted_data[lower_index];
    }
    
    // Линейная интерполяция
    double weight = index - lower_index;
    return (1 - weight) * sorted_data[lower_index] + weight * sorted_data[upper_index];
}


//
template<typename T>
void random_split_direct(std::vector<T>& Z, size_t n, double *X, double *Y, boost::random::mt19937& gen) {
    std::shuffle(Z.begin(), Z.end(), gen);

    for(int i=0; i<n; ++i){
        X[i] = Z[i];
        Y[i] = Z[i+n];
    }

}


//
template <typename T>
double compute_crit_val(int n, int M, double alpha, T d1,
                        std::function<double(double*, double*, int)> compute_test,
                        boost::random::mt19937& gen)
{
    double *Z_ = new double[2*n];
    sample(d1, 2*n, Z_, gen);
    std::vector<double> Z;
    for(int i=0; i<2*n; ++i){
        Z.push_back(Z_[i]);
    }
    delete[] Z_;
    
    std::vector<double> test_vals(M);

    double *X = new double[n * M];
    double *Y = new double[n * M];
    
    for(int i=0; i<M; ++i){
        random_split_direct(Z, n, X + (i * n), Y + (i * n), gen);
    }

    #pragma omp parallel for
    for(int i=0; i<M; ++i){
        double test_val = compute_test(X + (i * n), Y + (i * n), n);
        #pragma omp critical
        {
            test_vals[i] = test_val;
        }
    }

    delete[] X;
    delete[] Y;
    
    return quantile(test_vals, 1 - alpha);
}


template <typename T>
double compute_crit_val_v2(int n, int M, double alpha, T d1,
                        std::function<double(double*, double*, int)> compute_test,
                        boost::random::mt19937& gen)
{
    std::vector<double> test_vals(M);

    double *X = new double[n * M];
    double *Y = new double[n * M];
    
    for(int i=0; i<M; ++i){
        sample(d1, n, X + (i * n), gen);
        sample(d1, n, Y + (i * n), gen);
    }

    #pragma omp parallel for
    for(int i=0; i<M; ++i){
        double test_val = compute_test(X + (i * n), Y + (i * n), n);
        #pragma omp critical
        {
            test_vals[i] = test_val;
        }
    }

    delete[] X;
    delete[] Y;
    
    return quantile(test_vals, 1 - alpha);
}




enum class DistributionType {
    NORMAL = 0,
    CAUCHY = 1
};

std::vector<double> read_numbers_from_file(const std::string& filename);


std::vector<double> read_integrals(DistributionType d_type);


template <typename T>
std::pair<std::vector<double>, double> experiment_step(
                    double h1, double h2, double alpha,
                    int N, int M, std::vector<int> sample_sizes,
                    std::vector<double> integrals,
                    std::vector<double> crit_vals,
                    T dist_template, boost::random::mt19937& gen,
                    std::function<double(double*, double*, int)> compute_test)
{

    auto d1 = dist_template(1,0,0);

    double J1 = integrals[0];
    double J2 = integrals[1];
    double J3 = integrals[2];
    double J1_star = integrals[3] * h1 * h1;
    double J2_star = integrals[4] * h2 * h2;

    double b1 = sqrt(abs(J1_star));
    double b2 = sqrt(abs(J2_star));
    
    double a = pow( J2 + J1 * J1 - 2 * J3, 0.25);
    // double a = pow(fabs(J1), 0.5);
    double b = sqrt( b1 * b1 + b2 * b2);
    
    printf("a = %f  ||  b = %f\n", a, b);
    printf("b1 = %f  ||  b2 = %f\n", b1, b2);
    

    std::vector<double> emp_powers;

    for (int i=0; i<sample_sizes.size(); ++i){
        Timer t1;
        int n = sample_sizes[i];
        std::cout << "n = " << n << "  ||  ";
        auto d2_n = dist_template(n, h1, h2);
        double e_pow = compute_empirical_power(n, N, crit_vals[i], d1, d2_n, compute_test, gen);
        emp_powers.push_back(e_pow);
    }

    double a_pow = compute_asymptotic_power(alpha, b, a);

    return {emp_powers, a_pow};

}



boost::random::normal_distribution<double> get_normal(int n, double h1, double h2);


boost::random::cauchy_distribution<double> get_cauchy(int n, double h1, double h2);


template <typename T>
void run_experiment(double h1, std::vector<double> h2_vals,
                    double alpha, int N, int M,
                    std::vector<int> sample_sizes,
                    T dist_template, std::vector<double> integrals,
                    std::function<double(double*, double*, int)> compute_test,
                    bool compute_AP, boost::random::mt19937& gen)
{   
    std::cout << "N = " << N << "\n"; 
    std::cout << "M = " << M << "\n"; 
    std::cout << "alpha = " << alpha << "\n"; 
    std::cout << "h1 = " << h1 << "\n\n"; 
    auto d1 = dist_template(1,0,0);

    std::vector<double> crit_vals;
    std::cout << "computing crit_vals...\n";
    for(int n : sample_sizes){
        Timer t1;
        double cv = compute_crit_val_v2(n, M, alpha, d1, compute_test, gen);
        std::cout << "n = " << n << "  ||  cv = " <<  cv << " || ";
        crit_vals.push_back(cv);
    }
    std::cout << "\n";
    
    for(double h2 : h2_vals){
        {
        Timer t1;
        std::cout << "h2 = " << h2 << "\n";
        auto [e_pow, a_pow] = experiment_step(h1, h2,
            alpha, N, M, sample_sizes, integrals, crit_vals, dist_template, gen, compute_test);
        std::cout << "emp_powers = ";
        print_vector(e_pow);
        if(compute_AP){
            std::cout << "asp_power = " << a_pow << "\n";
        }
        }
        std::cout << "\n\n";  
    }
}


template <typename T>
void run_experiment(std::vector<double> h1_vals, double h2, 
                    double alpha, int N, int M,
                    std::vector<int> sample_sizes,
                    T dist_template, std::vector<double> integrals,
                    std::function<double(double*, double*, int)> compute_test,
                    bool compute_AP, boost::random::mt19937& gen)
{   

    std::cout << "N = " << N << "\n"; 
    std::cout << "M = " << M << "\n"; 
    std::cout << "alpha = " << alpha << "\n"; 
    std::cout << "h2 = " << h2 << "\n\n"; 
    auto d1 = dist_template(1,0,0);

    std::vector<double> crit_vals;
    std::cout << "computing crit_vals...\n";
    for(int n : sample_sizes){
        Timer t1;
        double cv = compute_crit_val_v2(n, M, alpha, d1, compute_test, gen);
        std::cout << "n = " << n << "  ||  cv = " <<  cv << " || ";
        crit_vals.push_back(cv);
    }
    std::cout << "\n";
    
    for(double h1 : h1_vals){
        {
        Timer t1;
        std::cout << "h1 = " << h1 << "\n";
        auto [e_pow, a_pow] = experiment_step(h1, h2,
            alpha, N, M, sample_sizes, integrals, crit_vals, dist_template, gen, compute_test);
        std::cout << "emp_powers = ";
        print_vector(e_pow);
        if (compute_AP){
            std::cout << "asp_power = " << a_pow << "\n";
        }
        }
        std::cout << "\n\n";  
    }
}




template<typename T = double>
class scaled_student_t_distribution : public boost::random::student_t_distribution<T> {
private:
    T h1, h2;
    int n;
    T denom;
    T shift;
public:
    scaled_student_t_distribution(T df, T h1_, T h2_, int n_)
        : boost::random::student_t_distribution<T>(df), h1(h1_), h2(h2_), n(n_) {
            denom = T(1) + h2 / std::sqrt(static_cast<T>(n));
            shift = h1 / std::sqrt(static_cast<T>(n));
        }

    template<typename Engine>
    T operator()(Engine& eng) {
        T t = boost::random::student_t_distribution<T>::operator()(eng);
        return (t - shift) / denom;
    }
};


template<typename T=double>
std::pair<std::vector<double>, double> experiment_step_st(
                    double h1, double h2, int df, double alpha,
                    int N, int M, std::vector<int> sample_sizes,
                    std::vector<double> integrals,
                    std::vector<double> crit_vals,
                    boost::random::mt19937& gen,
                    std::function<double(double*, double*, int)> compute_test)
{

    auto d1 = scaled_student_t_distribution<double>(df, 0, 0, 1);

    double J1 = integrals[0];
    double J2 = integrals[1];
    double J3 = integrals[2];
    double J1_star = integrals[3] * h1 * h1;
    double J2_star = integrals[4] * h2 * h2;

    double b1 = sqrt(abs(J1_star));
    double b2 = sqrt(abs(J2_star));
    
    double a = pow( J2 + J1 * J1 - 2 * J3, 0.25);
    // double a = pow(fabs(J1), 0.5);
    double b = sqrt( b1 * b1 + b2 * b2);
    
    std::vector<double> emp_powers;

    for (int i=0; i<sample_sizes.size(); ++i){
        Timer t1;
        int n = sample_sizes[i];
        std::cout << "n = " << n << "  ||  ";
        auto d2_n = scaled_student_t_distribution<double>(df, h1, h2, n);
        double e_pow = compute_empirical_power(n, N, crit_vals[i], d1, d2_n, compute_test, gen);
        emp_powers.push_back(e_pow);
    }

    double a_pow = compute_asymptotic_power(alpha, b, a);

    return {emp_powers, a_pow};

}

template<typename T=double>
void run_experiment_st(std::vector<double> h1_vals, double h2, 
                    double alpha, int N, int M, int df,
                    std::vector<int> sample_sizes,
                    std::vector<double> integrals,
                    std::function<double(double*, double*, int)> compute_test,
                    bool compute_AP, boost::random::mt19937& gen)
{   

    std::cout << "N = " << N << "\n"; 
    std::cout << "M = " << M << "\n"; 
    std::cout << "alpha = " << alpha << "\n"; 
    std::cout << "h2 = " << h2 << "\n\n"; 
    auto d1 = boost::random::student_t_distribution<double>(df);

    std::vector<double> crit_vals;
    std::cout << "computing crit_vals...\n";
    for(int n : sample_sizes){
        Timer t1;
        double cv = compute_crit_val_v2(n, M, alpha, d1, compute_test, gen);
        std::cout << "n = " << n << "  ||  cv = " <<  cv << " || ";
        crit_vals.push_back(cv);
    }
    std::cout << "\n";
    
    for(double h1 : h1_vals){
        {
        Timer t1;
        std::cout << "h1 = " << h1 << "\n";
        auto [e_pow, a_pow] = experiment_step_st(h1, h2, df,
            alpha, N, M, sample_sizes, integrals, crit_vals, gen, compute_test);
        std::cout << "emp_powers = ";
        print_vector(e_pow);
        if (compute_AP){
            std::cout << "asp_power = " << a_pow << "\n";
        }
        }
        std::cout << "\n\n";  
    }
}

template<typename T=double>
void run_experiment_st(double h1, std::vector<double> h2_vals,
                    double alpha, int N, int M, int df,
                    std::vector<int> sample_sizes,
                    std::vector<double> integrals,
                    std::function<double(double*, double*, int)> compute_test,
                    bool compute_AP, boost::random::mt19937& gen)
{   
    std::cout << "N = " << N << "\n"; 
    std::cout << "M = " << M << "\n"; 
    std::cout << "alpha = " << alpha << "\n"; 
    std::cout << "h1 = " << h1 << "\n\n"; 
    auto d1 = boost::random::student_t_distribution<double>(df);

    std::vector<double> crit_vals;
    std::cout << "computing crit_vals...\n";
    for(int n : sample_sizes){
        Timer t1;
        double cv = compute_crit_val_v2(n, M, alpha, d1, compute_test, gen);
        std::cout << "n = " << n << "  ||  cv = " <<  cv << " || ";
        crit_vals.push_back(cv);
    }
    std::cout << "\n";
    
    for(double h2 : h2_vals){
        {
        Timer t1;
        std::cout << "h2 = " << h2 << "\n";
        auto [e_pow, a_pow] = experiment_step_st(h1, h2, df,
            alpha, N, M, sample_sizes, integrals, crit_vals, gen, compute_test);
        std::cout << "emp_powers = ";
        print_vector(e_pow);
        if(compute_AP){
            std::cout << "asp_power = " << a_pow << "\n";
        }
        }
        std::cout << "\n\n";  
    }
}



