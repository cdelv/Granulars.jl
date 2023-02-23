#include <cmath>
#include <iostream>

// Function declaration
double Epsilon(double v_star);
double Target_Epsilon(double v_star, double target_e);
double dEpsilon_dv(double v_star);
double Newton_Raphson(double x0, double target, double (*f)(double, double), double (*fderiv)(double), double eps, int max_iter);
double get_A(double v, double v_star, double E1, double E2, double R1, double R2, double nu1, double nu2, double m1, double m2);
double get_v_star(double v, double A, double E1, double E2, double R1, double R2, double nu1, double nu2, double m1, double m2);

int main(int argc, char *argv[])
{
    // Restitution coefficient at given velocity
    double e = 0.00015;
    double v = 1.5;

    // Density
    double D1 = 2650.0;
    double D2 = 2650.0;

    // Young Modulus
    double E1 = 3.0e8;
    double E2 = 3.0e8;

    // Radius
    double R1 = 0.00000005;
    double R2 = 0.00000005;

    // Poisson Ratio
    double nu1 = 0.3;
    double nu2 = 0.3;

    // Volume
    double V1 = 4*M_PI*std::pow(R1,3)/3;
    double V2 = 4*M_PI*std::pow(R2,3)/3;

    // Mass
    double m1 = D1*V1;
    double m2 = D2*V2;

    double target_v_star = Newton_Raphson(1.0, e, Target_Epsilon, dEpsilon_dv, 1e-12, 200);

    double A = get_A(v, target_v_star, E1, E2, R1, R2, nu1, nu2, m1, m2);
    std::cout << "A = " << A << std::endl;

    std::cout << "Target e = " << e << std::endl;
    std::cout << "Achived e = " << Epsilon(get_v_star(v, A, E1, E2, R1, R2, nu1, nu2, m1, m2)) << std::endl;

    return 0;
}

double Epsilon(double v_star){
    // 1/4 coefficients
    //const double a_i[] = {1.0, 0.501086};
    //const double b_i[] = {1.0, 0.501086, 1.15345, 0.577977, 0.532178};

    // 3/6 coefficients
    const double a_i[] = {1.0, 1.07232, 0.574198, 0.141552};
    const double b_i[] = {1.0, 1.07232, 1.72765, 1.37842, 1.19449, 0.467273, 0.235585};

    // Initialize sum to 0
    double A=0.0, B=0.0, n=0.0;

    for (auto i : a_i){
        A += i*std::pow(v_star, n);
        n += 1.0;
    }

    n = 0.0;
    for (auto i : b_i){
        B += i*std::pow(v_star, n);
        n += 1.0;
    }

    return A/B;
}

double Target_Epsilon(double v_star, double target_e){
    return Epsilon(v_star) - target_e;
}

double dEpsilon_dv(double v_star){
    double dt = 1.0e-8;
    // Compute the derivative with central differences (for simplicity)
    return (Epsilon(v_star + dt) - Epsilon(v_star - dt))/(2*dt);
}

double Newton_Raphson(double x0, double target, double (*f)(double, double), double (*fderiv)(double), double eps, int max_iter)
{   
    double xr = x0;
    int i;

    for (i = 0; i < max_iter; ++i)
    {
        // Check if root was found
        if (std::fabs(f(xr, target)) <= eps)
            break;

        // New root according to the method
        xr = xr - f(xr, target) / fderiv(xr);
    }

    return xr;
}

double get_A(double v, double v_star, double E1, double E2, double R1, double R2, double nu1, double nu2, double m1, double m2){
    double Reff = R1*R2/(R1 + R2);
    double meff = m1*m2/(m1+m2);
    double Eeff = 1.0/( (1.0 - nu1 * nu1)/E1 + (1.0 - nu2 * nu2)/E2 );
    double rho = 4.0*Eeff*std::sqrt(Reff)/3.0;
    double beta = v_star*v_star*std::pow(v, -0.2);  // ^-1/5 = ^-2/10
    double A = 2.0*beta*std::pow(rho/meff, -0.4)/3; // ^-2/5
    return A;
}

double get_v_star(double v, double A, double E1, double E2, double R1, double R2, double nu1, double nu2, double m1, double m2){
    double Reff = R1*R2/(R1 + R2);
    double meff = m1*m2/(m1+m2);
    double Eeff = 1.0/( (1.0 - nu1 * nu1)/E1 + (1.0 - nu2 * nu2)/E2 );
    double rho = 4.0*Eeff*std::sqrt(Reff)/3.0;
    double beta = 3.0*A*std::pow(rho/meff, 0.4)/2; // ^2/5
    return std::sqrt(beta)*std::pow(v, 0.1); // ^1/10
}
