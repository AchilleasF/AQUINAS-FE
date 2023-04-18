// Last generated on Thu Mar 09 11:44:18 2023 
#pragma once 
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace bnu = boost::numeric::ublas;

bnu::matrix<double> AQUINAS_GN0G(double eta, double theta, double L, double r, double dphids, double c, double s, double n, double psi, double Nphi, double Ntheta, double Nphitheta)
{
    double eta2 = eta * eta;
    double eta3 = eta2 * eta;
    double cn = cos(n*theta);
    double sn = sin(n*theta);
    
    bnu::matrix<double> G(6,12);
    
    bnu::matrix<double> N0 = bnu::zero_matrix<double>(6,6);
    N0(0,0) = Nphi;
    N0(0,3) = Nphitheta;
    N0(1,1) = Nphi;
    N0(1,4) = Nphitheta;
    N0(2,2) = Nphi;
    N0(2,5) = Nphitheta;
    N0(3,3) = Ntheta;
    N0(3,0) = Nphitheta;
    N0(4,4) = Ntheta;
    N0(4,1) = Nphitheta;
    N0(5,5) = Ntheta;
    N0(5,2) = Nphitheta;

    double G_1_1 = 0.3e1 / 0.4e1 * s * cn * (-0.1e1 + eta2) / L;
    double G_1_2 = s * cn * (-0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
    double G_1_3 = 0.0e0;
    double G_1_4 = 0.0e0;
    double G_1_5 = 0.3e1 / 0.4e1 * c * cn * (-0.1e1 + eta2) / L;
    double G_1_6 = c * cn * (-0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
    double G_1_7 = -0.3e1 / 0.4e1 * s * cn * (-0.1e1 + eta2) / L;
    double G_1_8 = s * cn * (0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
    double G_1_9 = 0.0e0;
    double G_1_10 = 0.0e0;
    double G_1_11 = -0.3e1 / 0.4e1 * c * cn * (-0.1e1 + eta2) / L;
    double G_1_12 = c * cn * (0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;

    double G_2_1 = 0.0e0;
    double G_2_2 = 0.0e0;
    double G_2_3 = 0.3e1 / 0.4e1 * sn * (-0.1e1 + eta2) / L;
    double G_2_4 = sn * (-0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
    double G_2_5 = 0.0e0;
    double G_2_6 = 0.0e0;
    double G_2_7 = 0.0e0;
    double G_2_8 = 0.0e0;
    double G_2_9 = -0.3e1 / 0.4e1 * sn * (-0.1e1 + eta2) / L;
    double G_2_10 = sn * (0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
    double G_2_11 = 0.0e0;
    double G_2_12 = 0.0e0;

    double G_3_1 = 0.3e1 / 0.4e1 * c * cn * (-0.1e1 + eta2) / L;
    double G_3_2 = c * cn * (-0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
    double G_3_3 = 0.0e0;
    double G_3_4 = 0.0e0;
    double G_3_5 = -0.3e1 / 0.4e1 * s * cn * (-0.1e1 + eta2) / L;
    double G_3_6 = -s * cn * (-0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
    double G_3_7 = -0.3e1 / 0.4e1 * c * cn * (-0.1e1 + eta2) / L;
    double G_3_8 = c * cn * (0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
    double G_3_9 = 0.0e0;
    double G_3_10 = 0.0e0;
    double G_3_11 = 0.3e1 / 0.4e1 * s * cn * (-0.1e1 + eta2) / L;
    double G_3_12 = -s * cn * (0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;

    double G_4_1 = -s * n * sn * (eta + 0.2e1) * pow(eta - 0.1e1, 0.2e1) / r / 0.4e1;
    double G_4_2 = -s * n * sn * L * (eta + 0.1e1) * pow(eta - 0.1e1, 0.2e1) / r / 0.4e1;
    double G_4_3 = -sn * (eta + 0.2e1) * pow(eta - 0.1e1, 0.2e1) * s / r / 0.4e1;
    double G_4_4 = -sn * L * (eta + 0.1e1) * pow(eta - 0.1e1, 0.2e1) * s / r / 0.4e1;
    double G_4_5 = -c * n * sn * (eta + 0.2e1) * pow(eta - 0.1e1, 0.2e1) / r / 0.4e1;
    double G_4_6 = -c * n * sn * L * (eta + 0.1e1) * pow(eta - 0.1e1, 0.2e1) / r / 0.4e1;
    double G_4_7 = s * n * sn * (eta - 0.2e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;
    double G_4_8 = -s * n * sn * L * (eta - 0.1e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;
    double G_4_9 = sn * (-0.3e1 * eta - 0.2e1 + eta3) * s / r / 0.4e1;
    double G_4_10 = -sn * L * (eta - 0.1e1) * pow(eta + 0.1e1, 0.2e1) * s / r / 0.4e1;
    double G_4_11 = c * n * sn * (eta - 0.2e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;
    double G_4_12 = -c * n * sn * L * (eta - 0.1e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;

    double G_5_1 = cn * (eta + 0.2e1) * pow(eta - 0.1e1, 0.2e1) * (c * c + s * s) / r / 0.4e1;
    double G_5_2 = cn * L * (eta + 0.1e1) * pow(eta - 0.1e1, 0.2e1) * (c * c + s * s) / r / 0.4e1;
    double G_5_3 = n * cn * (-0.3e1 * eta + 0.2e1 + eta3) / r / 0.4e1;
    double G_5_4 = n * cn * L * (eta + 0.1e1) * pow(eta - 0.1e1, 0.2e1) / r / 0.4e1;
    double G_5_5 = 0.0e0;
    double G_5_6 = 0.0e0;
    double G_5_7 = -cn * (eta - 0.2e1) * pow(eta + 0.1e1, 0.2e1) * (c * c + s * s) / r / 0.4e1;
    double G_5_8 = cn * L * (eta - 0.1e1) * pow(eta + 0.1e1, 0.2e1) * (c * c + s * s) / r / 0.4e1;
    double G_5_9 = -n * cn * (eta - 0.2e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;
    double G_5_10 = n * cn * L * (eta - 0.1e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;
    double G_5_11 = 0.0e0;
    double G_5_12 = 0.0e0;

    double G_6_1 = -c * n * sn * (eta + 0.2e1) * pow(eta - 0.1e1, 0.2e1) / r / 0.4e1;
    double G_6_2 = -c * n * sn * L * (eta + 0.1e1) * pow(eta - 0.1e1, 0.2e1) / r / 0.4e1;
    double G_6_3 = -sn * (eta + 0.2e1) * pow(eta - 0.1e1, 0.2e1) * c / r / 0.4e1;
    double G_6_4 = -sn * L * (eta + 0.1e1) * pow(eta - 0.1e1, 0.2e1) * c / r / 0.4e1;
    double G_6_5 = s * n * sn * (eta + 0.2e1) * pow(eta - 0.1e1, 0.2e1) / r / 0.4e1;
    double G_6_6 = s * n * sn * L * (eta + 0.1e1) * pow(eta - 0.1e1, 0.2e1) / r / 0.4e1;
    double G_6_7 = c * n * sn * (eta - 0.2e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;
    double G_6_8 = -c * n * sn * L * (eta - 0.1e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;
    double G_6_9 = sn * (-0.3e1 * eta - 0.2e1 + eta3) * c / r / 0.4e1;
    double G_6_10 = -sn * L * (eta - 0.1e1) * pow(eta + 0.1e1, 0.2e1) * c / r / 0.4e1;
    double G_6_11 = -s * n * sn * (eta - 0.2e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;
    double G_6_12 = s * n * sn * L * (eta - 0.1e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;


    G(0,0) = G_1_1;
    G(0,1) = G_1_2;
    G(0,2) = G_1_3;
    G(0,3) = G_1_4;
    G(0,4) = G_1_5;
    G(0,5) = G_1_6;
    G(0,6) = G_1_7;
    G(0,7) = G_1_8;
    G(0,8) = G_1_9;
    G(0,9) = G_1_10;
    G(0,10) = G_1_11;
    G(0,11) = G_1_12;

    G(1,0) = G_2_1;
    G(1,1) = G_2_2;
    G(1,2) = G_2_3;
    G(1,3) = G_2_4;
    G(1,4) = G_2_5;
    G(1,5) = G_2_6;
    G(1,6) = G_2_7;
    G(1,7) = G_2_8;
    G(1,8) = G_2_9;
    G(1,9) = G_2_10;
    G(1,10) = G_2_11;
    G(1,11) = G_2_12;

    G(2,0) = G_3_1;
    G(2,1) = G_3_2;
    G(2,2) = G_3_3;
    G(2,3) = G_3_4;
    G(2,4) = G_3_5;
    G(2,5) = G_3_6;
    G(2,6) = G_3_7;
    G(2,7) = G_3_8;
    G(2,8) = G_3_9;
    G(2,9) = G_3_10;
    G(2,10) = G_3_11;
    G(2,11) = G_3_12;

    G(3,0) = G_4_1;
    G(3,1) = G_4_2;
    G(3,2) = G_4_3;
    G(3,3) = G_4_4;
    G(3,4) = G_4_5;
    G(3,5) = G_4_6;
    G(3,6) = G_4_7;
    G(3,7) = G_4_8;
    G(3,8) = G_4_9;
    G(3,9) = G_4_10;
    G(3,10) = G_4_11;
    G(3,11) = G_4_12;

    G(4,0) = G_5_1;
    G(4,1) = G_5_2;
    G(4,2) = G_5_3;
    G(4,3) = G_5_4;
    G(4,4) = G_5_5;
    G(4,5) = G_5_6;
    G(4,6) = G_5_7;
    G(4,7) = G_5_8;
    G(4,8) = G_5_9;
    G(4,9) = G_5_10;
    G(4,10) = G_5_11;
    G(4,11) = G_5_12;

    G(5,0) = G_6_1;
    G(5,1) = G_6_2;
    G(5,2) = G_6_3;
    G(5,3) = G_6_4;
    G(5,4) = G_6_5;
    G(5,5) = G_6_6;
    G(5,6) = G_6_7;
    G(5,7) = G_6_8;
    G(5,8) = G_6_9;
    G(5,9) = G_6_10;
    G(5,10) = G_6_11;
    G(5,11) = G_6_12;



    bnu::matrix<double> N0G = prod(N0,G);
    bnu::matrix<double> GN0G = prod(trans(G),N0G);
    return GN0G;

}