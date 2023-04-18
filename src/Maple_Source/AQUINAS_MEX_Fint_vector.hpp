// Last generated on Thu Mar 09 11:45:00 2023 
#pragma once 
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace bnu = boost::numeric::ublas;

bnu::vector<double> AQUINAS_Fint_vector(double eta, double L, double r, double dphids, double c, double s, double psi, const bnu::vector<double> elDOFs, const bnu::vector<double> Sigma, unsigned int include_nonlinear_G)
{
    double eta2 = eta * eta;
    double eta3 = eta2 * eta;

    bnu::matrix<double> B0a(6,12);
    bnu::matrix<double> Ga(6,12);
    bnu::matrix<double> Ba(6,12);
    bnu::matrix<double> Omega = bnu::zero_matrix<double>(6,6);

    double B0a_1_1 = 0.3e1 / 0.4e1 * c * (-0.1e1 + eta2) / L;
    double B0a_1_2 = c * (-0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
    double B0a_1_3 = 0.0e0;
    double B0a_1_4 = 0.0e0;
    double B0a_1_5 = -0.3e1 / 0.4e1 * s * (-0.1e1 + eta2) / L;
    double B0a_1_6 = -s * (-0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
    double B0a_1_7 = -0.3e1 / 0.4e1 * c * (-0.1e1 + eta2) / L;
    double B0a_1_8 = c * (0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
    double B0a_1_9 = 0.0e0;
    double B0a_1_10 = 0.0e0;
    double B0a_1_11 = 0.3e1 / 0.4e1 * s * (-0.1e1 + eta2) / L;
    double B0a_1_12 = -s * (0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;

    double B0a_2_1 = (eta + 0.2e1) * pow(eta - 0.1e1, 0.2e1) * (c * c + s * s) / r / 0.4e1;
    double B0a_2_2 = L * (eta + 0.1e1) * pow(eta - 0.1e1, 0.2e1) * (c * c + s * s) / r / 0.4e1;
    double B0a_2_3 = 0.0e0;
    double B0a_2_4 = 0.0e0;
    double B0a_2_5 = 0.0e0;
    double B0a_2_6 = 0.0e0;
    double B0a_2_7 = -(eta - 0.2e1) * pow(eta + 0.1e1, 0.2e1) * (c * c + s * s) / r / 0.4e1;
    double B0a_2_8 = L * (eta - 0.1e1) * pow(eta + 0.1e1, 0.2e1) * (c * c + s * s) / r / 0.4e1;
    double B0a_2_9 = 0.0e0;
    double B0a_2_10 = 0.0e0;
    double B0a_2_11 = 0.0e0;
    double B0a_2_12 = 0.0e0;

    double B0a_3_1 = 0.0e0;
    double B0a_3_2 = 0.0e0;
    double B0a_3_3 = 0.0e0;
    double B0a_3_4 = 0.0e0;
    double B0a_3_5 = 0.0e0;
    double B0a_3_6 = 0.0e0;
    double B0a_3_7 = 0.0e0;
    double B0a_3_8 = 0.0e0;
    double B0a_3_9 = 0.0e0;
    double B0a_3_10 = 0.0e0;
    double B0a_3_11 = 0.0e0;
    double B0a_3_12 = 0.0e0;

    double B0a_4_1 = ((-0.3e1 * c * eta2 + 0.3e1 * c) * dphids * L - 0.6e1 * eta * s) * pow(L, -0.2e1) / 0.4e1;
    double B0a_4_2 = (-0.3e1 * (eta + 0.1e1 / 0.3e1) * (eta - 0.1e1) * c * dphids * L - 0.6e1 * eta * s + 0.2e1 * s) / L / 0.4e1;
    double B0a_4_3 = 0.0e0;
    double B0a_4_4 = 0.0e0;
    double B0a_4_5 = ((0.3e1 * eta2 * s - 0.3e1 * s) * dphids * L - 0.6e1 * c * eta) * pow(L, -0.2e1) / 0.4e1;
    double B0a_4_6 = (0.3e1 * (eta + 0.1e1 / 0.3e1) * (eta - 0.1e1) * s * dphids * L - 0.6e1 * c * eta + 0.2e1 * c) / L / 0.4e1;
    double B0a_4_7 = ((0.3e1 * c * eta2 - 0.3e1 * c) * dphids * L + 0.6e1 * eta * s) * pow(L, -0.2e1) / 0.4e1;
    double B0a_4_8 = (-0.3e1 * (eta + 0.1e1) * c * (eta - 0.1e1 / 0.3e1) * dphids * L - 0.6e1 * eta * s - 0.2e1 * s) / L / 0.4e1;
    double B0a_4_9 = 0.0e0;
    double B0a_4_10 = 0.0e0;
    double B0a_4_11 = ((-0.3e1 * eta2 * s + 0.3e1 * s) * dphids * L + 0.6e1 * c * eta) * pow(L, -0.2e1) / 0.4e1;
    double B0a_4_12 = (0.3e1 * (eta + 0.1e1) * s * (eta - 0.1e1 / 0.3e1) * dphids * L - 0.6e1 * c * eta - 0.2e1 * c) / L / 0.4e1;

    double B0a_5_1 = -0.3e1 / 0.4e1 / r * c * s / L * (-0.1e1 + eta2);
    double B0a_5_2 = -0.3e1 / 0.4e1 * (eta + 0.1e1 / 0.3e1) * s * c * (eta - 0.1e1) / r;
    double B0a_5_3 = 0.0e0;
    double B0a_5_4 = 0.0e0;
    double B0a_5_5 = -0.3e1 / 0.4e1 / r * c * c / L * (-0.1e1 + eta2);
    double B0a_5_6 = -c * c * (eta - 0.1e1) / r * (0.3e1 * eta + 0.1e1) / 0.4e1;
    double B0a_5_7 = 0.3e1 / 0.4e1 / r * c * s / L * (-0.1e1 + eta2);
    double B0a_5_8 = -0.3e1 / 0.4e1 * (eta + 0.1e1) * (eta - 0.1e1 / 0.3e1) * s * c / r;
    double B0a_5_9 = 0.0e0;
    double B0a_5_10 = 0.0e0;
    double B0a_5_11 = 0.3e1 / 0.4e1 / r * c * c / L * (-0.1e1 + eta2);
    double B0a_5_12 = -c * c * (eta + 0.1e1) / r * (0.3e1 * eta - 0.1e1) / 0.4e1;

    double B0a_6_1 = 0.0e0;
    double B0a_6_2 = 0.0e0;
    double B0a_6_3 = 0.0e0;
    double B0a_6_4 = 0.0e0;
    double B0a_6_5 = 0.0e0;
    double B0a_6_6 = 0.0e0;
    double B0a_6_7 = 0.0e0;
    double B0a_6_8 = 0.0e0;
    double B0a_6_9 = 0.0e0;
    double B0a_6_10 = 0.0e0;
    double B0a_6_11 = 0.0e0;
    double B0a_6_12 = 0.0e0;


    B0a(0,0) = B0a_1_1;
    B0a(0,1) = B0a_1_2;
    B0a(0,2) = B0a_1_3;
    B0a(0,3) = B0a_1_4;
    B0a(0,4) = B0a_1_5;
    B0a(0,5) = B0a_1_6;
    B0a(0,6) = B0a_1_7;
    B0a(0,7) = B0a_1_8;
    B0a(0,8) = B0a_1_9;
    B0a(0,9) = B0a_1_10;
    B0a(0,10) = B0a_1_11;
    B0a(0,11) = B0a_1_12;

    B0a(1,0) = B0a_2_1;
    B0a(1,1) = B0a_2_2;
    B0a(1,2) = B0a_2_3;
    B0a(1,3) = B0a_2_4;
    B0a(1,4) = B0a_2_5;
    B0a(1,5) = B0a_2_6;
    B0a(1,6) = B0a_2_7;
    B0a(1,7) = B0a_2_8;
    B0a(1,8) = B0a_2_9;
    B0a(1,9) = B0a_2_10;
    B0a(1,10) = B0a_2_11;
    B0a(1,11) = B0a_2_12;

    B0a(2,0) = B0a_3_1;
    B0a(2,1) = B0a_3_2;
    B0a(2,2) = B0a_3_3;
    B0a(2,3) = B0a_3_4;
    B0a(2,4) = B0a_3_5;
    B0a(2,5) = B0a_3_6;
    B0a(2,6) = B0a_3_7;
    B0a(2,7) = B0a_3_8;
    B0a(2,8) = B0a_3_9;
    B0a(2,9) = B0a_3_10;
    B0a(2,10) = B0a_3_11;
    B0a(2,11) = B0a_3_12;

    B0a(3,0) = B0a_4_1;
    B0a(3,1) = B0a_4_2;
    B0a(3,2) = B0a_4_3;
    B0a(3,3) = B0a_4_4;
    B0a(3,4) = B0a_4_5;
    B0a(3,5) = B0a_4_6;
    B0a(3,6) = B0a_4_7;
    B0a(3,7) = B0a_4_8;
    B0a(3,8) = B0a_4_9;
    B0a(3,9) = B0a_4_10;
    B0a(3,10) = B0a_4_11;
    B0a(3,11) = B0a_4_12;

    B0a(4,0) = B0a_5_1;
    B0a(4,1) = B0a_5_2;
    B0a(4,2) = B0a_5_3;
    B0a(4,3) = B0a_5_4;
    B0a(4,4) = B0a_5_5;
    B0a(4,5) = B0a_5_6;
    B0a(4,6) = B0a_5_7;
    B0a(4,7) = B0a_5_8;
    B0a(4,8) = B0a_5_9;
    B0a(4,9) = B0a_5_10;
    B0a(4,10) = B0a_5_11;
    B0a(4,11) = B0a_5_12;

    B0a(5,0) = B0a_6_1;
    B0a(5,1) = B0a_6_2;
    B0a(5,2) = B0a_6_3;
    B0a(5,3) = B0a_6_4;
    B0a(5,4) = B0a_6_5;
    B0a(5,5) = B0a_6_6;
    B0a(5,6) = B0a_6_7;
    B0a(5,7) = B0a_6_8;
    B0a(5,8) = B0a_6_9;
    B0a(5,9) = B0a_6_10;
    B0a(5,10) = B0a_6_11;
    B0a(5,11) = B0a_6_12;


    if (include_nonlinear_G==0)
    {

        Ba = B0a;

    }
    else
    {

        double Ga_1_1 = 0.3e1 / 0.4e1 * s * (-0.1e1 + eta2) / L;
        double Ga_1_2 = s * (-0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
        double Ga_1_3 = 0.0e0;
        double Ga_1_4 = 0.0e0;
        double Ga_1_5 = 0.3e1 / 0.4e1 * c * (-0.1e1 + eta2) / L;
        double Ga_1_6 = c * (-0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
        double Ga_1_7 = -0.3e1 / 0.4e1 * s * (-0.1e1 + eta2) / L;
        double Ga_1_8 = s * (0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
        double Ga_1_9 = 0.0e0;
        double Ga_1_10 = 0.0e0;
        double Ga_1_11 = -0.3e1 / 0.4e1 * c * (-0.1e1 + eta2) / L;
        double Ga_1_12 = c * (0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;

        double Ga_2_1 = 0.0e0;
        double Ga_2_2 = 0.0e0;
        double Ga_2_3 = 0.0e0;
        double Ga_2_4 = 0.0e0;
        double Ga_2_5 = 0.0e0;
        double Ga_2_6 = 0.0e0;
        double Ga_2_7 = 0.0e0;
        double Ga_2_8 = 0.0e0;
        double Ga_2_9 = 0.0e0;
        double Ga_2_10 = 0.0e0;
        double Ga_2_11 = 0.0e0;
        double Ga_2_12 = 0.0e0;

        double Ga_3_1 = 0.3e1 / 0.4e1 * c * (-0.1e1 + eta2) / L;
        double Ga_3_2 = c * (-0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
        double Ga_3_3 = 0.0e0;
        double Ga_3_4 = 0.0e0;
        double Ga_3_5 = -0.3e1 / 0.4e1 * s * (-0.1e1 + eta2) / L;
        double Ga_3_6 = -s * (-0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
        double Ga_3_7 = -0.3e1 / 0.4e1 * c * (-0.1e1 + eta2) / L;
        double Ga_3_8 = c * (0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
        double Ga_3_9 = 0.0e0;
        double Ga_3_10 = 0.0e0;
        double Ga_3_11 = 0.3e1 / 0.4e1 * s * (-0.1e1 + eta2) / L;
        double Ga_3_12 = -s * (0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;

        double Ga_4_1 = 0.0e0;
        double Ga_4_2 = 0.0e0;
        double Ga_4_3 = 0.0e0;
        double Ga_4_4 = 0.0e0;
        double Ga_4_5 = 0.0e0;
        double Ga_4_6 = 0.0e0;
        double Ga_4_7 = 0.0e0;
        double Ga_4_8 = 0.0e0;
        double Ga_4_9 = 0.0e0;
        double Ga_4_10 = 0.0e0;
        double Ga_4_11 = 0.0e0;
        double Ga_4_12 = 0.0e0;

        double Ga_5_1 = (eta + 0.2e1) * pow(eta - 0.1e1, 0.2e1) * (c * c + s * s) / r / 0.4e1;
        double Ga_5_2 = L * (eta + 0.1e1) * pow(eta - 0.1e1, 0.2e1) * (c * c + s * s) / r / 0.4e1;
        double Ga_5_3 = 0.0e0;
        double Ga_5_4 = 0.0e0;
        double Ga_5_5 = 0.0e0;
        double Ga_5_6 = 0.0e0;
        double Ga_5_7 = -(eta - 0.2e1) * pow(eta + 0.1e1, 0.2e1) * (c * c + s * s) / r / 0.4e1;
        double Ga_5_8 = L * (eta - 0.1e1) * pow(eta + 0.1e1, 0.2e1) * (c * c + s * s) / r / 0.4e1;
        double Ga_5_9 = 0.0e0;
        double Ga_5_10 = 0.0e0;
        double Ga_5_11 = 0.0e0;
        double Ga_5_12 = 0.0e0;

        double Ga_6_1 = 0.0e0;
        double Ga_6_2 = 0.0e0;
        double Ga_6_3 = 0.0e0;
        double Ga_6_4 = 0.0e0;
        double Ga_6_5 = 0.0e0;
        double Ga_6_6 = 0.0e0;
        double Ga_6_7 = 0.0e0;
        double Ga_6_8 = 0.0e0;
        double Ga_6_9 = 0.0e0;
        double Ga_6_10 = 0.0e0;
        double Ga_6_11 = 0.0e0;
        double Ga_6_12 = 0.0e0;


        Ga(0,0) = Ga_1_1;
        Ga(0,1) = Ga_1_2;
        Ga(0,2) = Ga_1_3;
        Ga(0,3) = Ga_1_4;
        Ga(0,4) = Ga_1_5;
        Ga(0,5) = Ga_1_6;
        Ga(0,6) = Ga_1_7;
        Ga(0,7) = Ga_1_8;
        Ga(0,8) = Ga_1_9;
        Ga(0,9) = Ga_1_10;
        Ga(0,10) = Ga_1_11;
        Ga(0,11) = Ga_1_12;

        Ga(1,0) = Ga_2_1;
        Ga(1,1) = Ga_2_2;
        Ga(1,2) = Ga_2_3;
        Ga(1,3) = Ga_2_4;
        Ga(1,4) = Ga_2_5;
        Ga(1,5) = Ga_2_6;
        Ga(1,6) = Ga_2_7;
        Ga(1,7) = Ga_2_8;
        Ga(1,8) = Ga_2_9;
        Ga(1,9) = Ga_2_10;
        Ga(1,10) = Ga_2_11;
        Ga(1,11) = Ga_2_12;

        Ga(2,0) = Ga_3_1;
        Ga(2,1) = Ga_3_2;
        Ga(2,2) = Ga_3_3;
        Ga(2,3) = Ga_3_4;
        Ga(2,4) = Ga_3_5;
        Ga(2,5) = Ga_3_6;
        Ga(2,6) = Ga_3_7;
        Ga(2,7) = Ga_3_8;
        Ga(2,8) = Ga_3_9;
        Ga(2,9) = Ga_3_10;
        Ga(2,10) = Ga_3_11;
        Ga(2,11) = Ga_3_12;

        Ga(3,0) = Ga_4_1;
        Ga(3,1) = Ga_4_2;
        Ga(3,2) = Ga_4_3;
        Ga(3,3) = Ga_4_4;
        Ga(3,4) = Ga_4_5;
        Ga(3,5) = Ga_4_6;
        Ga(3,6) = Ga_4_7;
        Ga(3,7) = Ga_4_8;
        Ga(3,8) = Ga_4_9;
        Ga(3,9) = Ga_4_10;
        Ga(3,10) = Ga_4_11;
        Ga(3,11) = Ga_4_12;

        Ga(4,0) = Ga_5_1;
        Ga(4,1) = Ga_5_2;
        Ga(4,2) = Ga_5_3;
        Ga(4,3) = Ga_5_4;
        Ga(4,4) = Ga_5_5;
        Ga(4,5) = Ga_5_6;
        Ga(4,6) = Ga_5_7;
        Ga(4,7) = Ga_5_8;
        Ga(4,8) = Ga_5_9;
        Ga(4,9) = Ga_5_10;
        Ga(4,10) = Ga_5_11;
        Ga(4,11) = Ga_5_12;

        Ga(5,0) = Ga_6_1;
        Ga(5,1) = Ga_6_2;
        Ga(5,2) = Ga_6_3;
        Ga(5,3) = Ga_6_4;
        Ga(5,4) = Ga_6_5;
        Ga(5,5) = Ga_6_6;
        Ga(5,6) = Ga_6_7;
        Ga(5,7) = Ga_6_8;
        Ga(5,8) = Ga_6_9;
        Ga(5,9) = Ga_6_10;
        Ga(5,10) = Ga_6_11;
        Ga(5,11) = Ga_6_12;



        bnu::vector<double> betas = prod(Ga,elDOFs);

        Omega(0,0) = betas(0);
        Omega(0,1) = betas(1);
        Omega(0,2) = betas(2);
        Omega(1,3) = betas(3);
        Omega(1,4) = betas(4);
        Omega(1,5) = betas(5);
        Omega(2,0) = betas(3);
        Omega(2,1) = betas(4);
        Omega(2,2) = betas(5);
        Omega(2,3) = betas(0);
        Omega(2,4) = betas(1);
        Omega(2,5) = betas(2);

        bnu::matrix<double> BLa = prod(Omega,Ga);

        Ba = B0a + BLa;

    }

    bnu::vector<double> Fint = prod(trans(Ba),Sigma);
    return Fint;
}
