#include <random>
#include <iostream>



double g(double x){
    return log(1 + x * x);
}

double d2_g(double x){
    return 2 * (1 - x * x) / ((1 + x * x) * (1 + x * x));
}


int main(){
    // Жестко бахаем Монте-Карло с -O3 -march-native

    int N = 1e+8;
    std::mt19937 gen(1339);

    double J1, J2, J3, J1_star, J2_star;
    
    // std::normal_distribution<double> dist(0.0, 1.0);
    // std::cauchy_distribution<double> dist(0.0, 1.0);
    std::student_t_distribution<double> dist(69);
    
    // compute J1
    double curr_mean = 0;
    for(int i=1; i<=N; ++i){
        curr_mean = (curr_mean * (i - 1) + g(dist(gen) - dist(gen)) ) / i;
    }
    J1 = curr_mean;
    printf("J1 = %f\n", J1);


    // compute J2
    curr_mean = 0;
    double v;
    for(int i=1; i<=N; ++i){
        v = g(dist(gen) - dist(gen));
        curr_mean = (curr_mean * (i - 1) + v * v ) / i;
    }
    J2 = curr_mean;
    printf("J2 = %f\n", J2);


    // compute J3
    curr_mean = 0;
    double x,y,z;
    for(int i=1; i<=N; ++i){
        x = dist(gen);
        y = dist(gen);
        z = dist(gen);
        curr_mean = (curr_mean * (i - 1) + g(x - y) * g(x - z) ) / i;
    }
    J3 = curr_mean;
    printf("J3 = %f\n", J3);


    // compute J1_star
    curr_mean = 0;
    for(int i=1; i<=N; ++i){
        curr_mean = (curr_mean * (i - 1) + d2_g(dist(gen) - dist(gen)) ) / i;
    }
    curr_mean /= 2;
    J1_star = curr_mean;
    printf("J1_star = %f\n", J1_star);


    // compute J2_star
    curr_mean = 0;
    for(int i=1; i<=N; ++i){
        x = dist(gen);
        y = dist(gen);
        curr_mean = (curr_mean * (i - 1) + 
                ((y * y - 0.5 * (x - y) * (x - y)) * d2_g(x - y)) ) / i;
    }
    curr_mean /= 2;
    J2_star = curr_mean;
    printf("J2_star = %f\n", J2_star);


    // FILE *f = fopen("./precomputed_integrals_cpp/student_30.txt", "w");
    // fprintf(f, "%f %f %f %f %f", J1, J2, J3, J1_star, J2_star);
    // fclose(f);

}
