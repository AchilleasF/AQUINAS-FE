// Last generated on Thu Mar 09 11:44:28 2023 
#pragma once 
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace bnu = boost::numeric::ublas;

bnu::matrix<double> AQUINAS_BDTB(double eta, double theta, double L, double r, double dphids, double c, double s, double n, double psi, const bnu::vector<double>& elDofs, const bnu::matrix<double> DT) 
{
    double eta2 = eta * eta;
    double eta3 = eta2 * eta;
    double cn = cos(n*theta);
    double sn = sin(n*theta);
    double U1 = elDofs(0);
    double dUdS_1 = elDofs(1);
    double W1 = elDofs(4);
    double dWdS_1 = elDofs(5);
    double U2 = elDofs(6);
    double dUdS_2 = elDofs(7);
    double W2 = elDofs(10);
    double dWdS_2 = elDofs(11);
    
    bnu::matrix<double> B0(6,12);
    bnu::matrix<double> G(6,12);
    bnu::matrix<double> Omega = bnu::zero_matrix<double>(6,6);
    
    double B0_1_1 = 0.3e1 / 0.4e1 * c * cn * (-0.1e1 + eta2) / L;
    double B0_1_2 = c * cn * (-0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
    double B0_1_3 = 0.0e0;
    double B0_1_4 = 0.0e0;
    double B0_1_5 = -0.3e1 / 0.4e1 * s * cn * (-0.1e1 + eta2) / L;
    double B0_1_6 = -s * cn * (-0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
    double B0_1_7 = -0.3e1 / 0.4e1 * c * cn * (-0.1e1 + eta2) / L;
    double B0_1_8 = c * cn * (0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
    double B0_1_9 = 0.0e0;
    double B0_1_10 = 0.0e0;
    double B0_1_11 = 0.3e1 / 0.4e1 * s * cn * (-0.1e1 + eta2) / L;
    double B0_1_12 = -s * cn * (0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;

    double B0_2_1 = cn * (eta + 0.2e1) * pow(eta - 0.1e1, 0.2e1) * (c * c + s * s) / r / 0.4e1;
    double B0_2_2 = cn * L * (eta + 0.1e1) * pow(eta - 0.1e1, 0.2e1) * (c * c + s * s) / r / 0.4e1;
    double B0_2_3 = n * cn * (-0.3e1 * eta + 0.2e1 + eta3) / r / 0.4e1;
    double B0_2_4 = n * cn * L * (eta + 0.1e1) * pow(eta - 0.1e1, 0.2e1) / r / 0.4e1;
    double B0_2_5 = 0.0e0;
    double B0_2_6 = 0.0e0;
    double B0_2_7 = -cn * (eta - 0.2e1) * pow(eta + 0.1e1, 0.2e1) * (c * c + s * s) / r / 0.4e1;
    double B0_2_8 = cn * L * (eta - 0.1e1) * pow(eta + 0.1e1, 0.2e1) * (c * c + s * s) / r / 0.4e1;
    double B0_2_9 = -n * cn * (eta - 0.2e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;
    double B0_2_10 = n * cn * L * (eta - 0.1e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;
    double B0_2_11 = 0.0e0;
    double B0_2_12 = 0.0e0;

    double B0_3_1 = -c * n * sn * (eta + 0.2e1) * pow(eta - 0.1e1, 0.2e1) / r / 0.4e1;
    double B0_3_2 = -c * n * sn * L * (eta + 0.1e1) * pow(eta - 0.1e1, 0.2e1) / r / 0.4e1;
    double B0_3_3 = -sn * (L * c * eta2 + (c * L - 0.3e1 * r) * eta - 0.2e1 * c * L - 0.3e1 * r) * (eta - 0.1e1) / r / L / 0.4e1;
    double B0_3_4 = -sn * ((-0.3e1 * eta - 0.1e1) * r + L * c * (eta - 0.1e1) * (eta + 0.1e1)) * (eta - 0.1e1) / r / 0.4e1;
    double B0_3_5 = s * n * sn * (eta + 0.2e1) * pow(eta - 0.1e1, 0.2e1) / r / 0.4e1;
    double B0_3_6 = s * n * sn * L * (eta + 0.1e1) * pow(eta - 0.1e1, 0.2e1) / r / 0.4e1;
    double B0_3_7 = c * n * sn * (eta - 0.2e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;
    double B0_3_8 = -c * n * sn * L * (eta - 0.1e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;
    double B0_3_9 = sn * (L * c * eta2 + (-c * L - 0.3e1 * r) * eta - 0.2e1 * c * L + 0.3e1 * r) * (eta + 0.1e1) / r / L / 0.4e1;
    double B0_3_10 = -((-0.3e1 * eta + 0.1e1) * r + L * c * (eta - 0.1e1) * (eta + 0.1e1)) * sn * (eta + 0.1e1) / r / 0.4e1;
    double B0_3_11 = -s * n * sn * (eta - 0.2e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;
    double B0_3_12 = s * n * sn * L * (eta - 0.1e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;

    double B0_4_1 = -0.3e1 / 0.4e1 * (c * dphids * (eta - 0.1e1) * (eta + 0.1e1) * L + 0.2e1 * eta * s) * cn * pow(L, -0.2e1);
    double B0_4_2 = -0.3e1 / 0.4e1 * ((eta + 0.1e1 / 0.3e1) * c * dphids * (eta - 0.1e1) * L + 0.2e1 * s * (eta - 0.1e1 / 0.3e1)) * cn / L;
    double B0_4_3 = 0.0e0;
    double B0_4_4 = 0.0e0;
    double B0_4_5 = 0.3e1 / 0.4e1 * (s * dphids * (eta - 0.1e1) * (eta + 0.1e1) * L - 0.2e1 * c * eta) * cn * pow(L, -0.2e1);
    double B0_4_6 = 0.3e1 / 0.4e1 * (s * (eta + 0.1e1 / 0.3e1) * dphids * (eta - 0.1e1) * L - 0.2e1 * (eta - 0.1e1 / 0.3e1) * c) * cn / L;
    double B0_4_7 = 0.3e1 / 0.4e1 * (c * dphids * (eta - 0.1e1) * (eta + 0.1e1) * L + 0.2e1 * eta * s) * cn * pow(L, -0.2e1);
    double B0_4_8 = -0.3e1 / 0.4e1 * cn * ((eta + 0.1e1) * (eta - 0.1e1 / 0.3e1) * c * dphids * L + 0.2e1 * s * (eta + 0.1e1 / 0.3e1)) / L;
    double B0_4_9 = 0.0e0;
    double B0_4_10 = 0.0e0;
    double B0_4_11 = -0.3e1 / 0.4e1 * (s * dphids * (eta - 0.1e1) * (eta + 0.1e1) * L - 0.2e1 * c * eta) * cn * pow(L, -0.2e1);
    double B0_4_12 = 0.3e1 / 0.4e1 * (s * (eta + 0.1e1) * (eta - 0.1e1 / 0.3e1) * dphids * L - 0.2e1 * (eta + 0.1e1 / 0.3e1) * c) * cn / L;

    double B0_5_1 = s * (L * n * n * eta2 + (n * n * L - 0.3e1 * c * r) * eta - 0.2e1 * n * n * L - 0.3e1 * c * r) * (eta - 0.1e1) * cn * pow(r, -0.2e1) / L / 0.4e1;
    double B0_5_2 = s * (eta - 0.1e1) * cn * ((-0.3e1 * c * eta - c) * r + L * n * n * (eta - 0.1e1) * (eta + 0.1e1)) * pow(r, -0.2e1) / 0.4e1;
    double B0_5_3 = s * n * cn * (eta + 0.2e1) * pow(eta - 0.1e1, 0.2e1) * pow(r, -0.2e1) / 0.4e1;
    double B0_5_4 = s * n * cn * L * (eta + 0.1e1) * pow(eta - 0.1e1, 0.2e1) * pow(r, -0.2e1) / 0.4e1;
    double B0_5_5 = cn * (L * n * n * eta2 + (n * n * L - 0.3e1 * c * r) * eta - 0.2e1 * n * n * L - 0.3e1 * c * r) * (eta - 0.1e1) * c * pow(r, -0.2e1) / L / 0.4e1;
    double B0_5_6 = cn * ((-0.3e1 * eta - 0.1e1) * c * r + L * n * n * (eta - 0.1e1) * (eta + 0.1e1)) * (eta - 0.1e1) * c * pow(r, -0.2e1) / 0.4e1;
    double B0_5_7 = -s * (eta + 0.1e1) * (-L * eta * n * n + L * n * n * eta2 - 0.2e1 * n * n * L - 0.3e1 * c * r * eta + 0.3e1 * c * r) * cn * pow(r, -0.2e1) / L / 0.4e1;
    double B0_5_8 = s * ((-0.3e1 * c * eta + c) * r + L * n * n * (eta - 0.1e1) * (eta + 0.1e1)) * (eta + 0.1e1) * cn * pow(r, -0.2e1) / 0.4e1;
    double B0_5_9 = -s * n * cn * (eta - 0.2e1) * pow(eta + 0.1e1, 0.2e1) * pow(r, -0.2e1) / 0.4e1;
    double B0_5_10 = s * n * cn * L * (eta - 0.1e1) * pow(eta + 0.1e1, 0.2e1) * pow(r, -0.2e1) / 0.4e1;
    double B0_5_11 = -(n * n * (eta + 0.1e1) * (eta - 0.2e1) * L - 0.3e1 * r * c * (eta - 0.1e1)) * (eta + 0.1e1) * c * cn * pow(r, -0.2e1) / L / 0.4e1;
    double B0_5_12 = cn * (eta + 0.1e1) * ((-0.3e1 * eta + 0.1e1) * c * r + L * n * n * (eta - 0.1e1) * (eta + 0.1e1)) * c * pow(r, -0.2e1) / 0.4e1;

    double B0_6_1 = n * sn * (c * (eta + 0.2e1) * (eta - 0.1e1) * (dphids * r - 0.2e1 * s) * L + 0.6e1 * r * s * (eta + 0.1e1)) * (eta - 0.1e1) * pow(r, -0.2e1) / L / 0.4e1;
    double B0_6_2 = n * ((L * c * dphids * eta2 - c * dphids * L + 0.6e1 * eta * s + 0.2e1 * s) * r - 0.2e1 * L * c * s * (eta - 0.1e1) * (eta + 0.1e1)) * (eta - 0.1e1) * sn * pow(r, -0.2e1) / 0.4e1;
    double B0_6_3 = sn * (c * (eta + 0.2e1) * (eta - 0.1e1) * (dphids * r - 0.2e1 * s) * L + 0.3e1 * r * s * (eta + 0.1e1)) * (eta - 0.1e1) * pow(r, -0.2e1) / L / 0.4e1;
    double B0_6_4 = sn * ((L * c * dphids * eta2 - c * dphids * L + 0.3e1 * eta * s + s) * r - 0.2e1 * L * c * s * (eta - 0.1e1) * (eta + 0.1e1)) * (eta - 0.1e1) * pow(r, -0.2e1) / 0.4e1;
    double B0_6_5 = -n * sn * (eta - 0.1e1) * ((eta + 0.2e1) * (eta - 0.1e1) * (s * dphids * r / 0.2e1 + c * c) * L - 0.3e1 * r * c * (eta + 0.1e1)) * pow(r, -0.2e1) / L / 0.2e1;
    double B0_6_6 = -n * sn * ((L * dphids * eta2 * s / 0.2e1 - s * dphids * L / 0.2e1 - 0.3e1 * c * eta - c) * r + L * c * c * (eta - 0.1e1) * (eta + 0.1e1)) * (eta - 0.1e1) * pow(r, -0.2e1) / 0.2e1;
    double B0_6_7 = -n * (eta + 0.1e1) * sn * (c * (eta + 0.1e1) * (eta - 0.2e1) * (dphids * r - 0.2e1 * s) * L + 0.6e1 * r * s * (eta - 0.1e1)) * pow(r, -0.2e1) / L / 0.4e1;
    double B0_6_8 = n * ((L * c * dphids * eta2 - c * dphids * L + 0.6e1 * eta * s - 0.2e1 * s) * r - 0.2e1 * L * c * s * (eta - 0.1e1) * (eta + 0.1e1)) * (eta + 0.1e1) * sn * pow(r, -0.2e1) / 0.4e1;
    double B0_6_9 = -(c * (eta + 0.1e1) * (eta - 0.2e1) * (dphids * r - 0.2e1 * s) * L + 0.3e1 * r * s * (eta - 0.1e1)) * (eta + 0.1e1) * sn * pow(r, -0.2e1) / L / 0.4e1;
    double B0_6_10 = sn * ((L * c * dphids * eta2 - c * dphids * L + 0.3e1 * eta * s - s) * r - 0.2e1 * L * c * s * (eta - 0.1e1) * (eta + 0.1e1)) * (eta + 0.1e1) * pow(r, -0.2e1) / 0.4e1;
    double B0_6_11 = n * (eta + 0.1e1) * sn * ((eta - 0.2e1) * (eta + 0.1e1) * (s * dphids * r / 0.2e1 + c * c) * L - 0.3e1 * r * c * (eta - 0.1e1)) * pow(r, -0.2e1) / L / 0.2e1;
    double B0_6_12 = -n * ((L * dphids * eta2 * s / 0.2e1 - s * dphids * L / 0.2e1 - 0.3e1 * c * eta + c) * r + L * c * c * (eta - 0.1e1) * (eta + 0.1e1)) * (eta + 0.1e1) * sn * pow(r, -0.2e1) / 0.2e1;


    B0(0,0) = B0_1_1;
    B0(0,1) = B0_1_2;
    B0(0,2) = B0_1_3;
    B0(0,3) = B0_1_4;
    B0(0,4) = B0_1_5;
    B0(0,5) = B0_1_6;
    B0(0,6) = B0_1_7;
    B0(0,7) = B0_1_8;
    B0(0,8) = B0_1_9;
    B0(0,9) = B0_1_10;
    B0(0,10) = B0_1_11;
    B0(0,11) = B0_1_12;

    B0(1,0) = B0_2_1;
    B0(1,1) = B0_2_2;
    B0(1,2) = B0_2_3;
    B0(1,3) = B0_2_4;
    B0(1,4) = B0_2_5;
    B0(1,5) = B0_2_6;
    B0(1,6) = B0_2_7;
    B0(1,7) = B0_2_8;
    B0(1,8) = B0_2_9;
    B0(1,9) = B0_2_10;
    B0(1,10) = B0_2_11;
    B0(1,11) = B0_2_12;

    B0(2,0) = B0_3_1;
    B0(2,1) = B0_3_2;
    B0(2,2) = B0_3_3;
    B0(2,3) = B0_3_4;
    B0(2,4) = B0_3_5;
    B0(2,5) = B0_3_6;
    B0(2,6) = B0_3_7;
    B0(2,7) = B0_3_8;
    B0(2,8) = B0_3_9;
    B0(2,9) = B0_3_10;
    B0(2,10) = B0_3_11;
    B0(2,11) = B0_3_12;

    B0(3,0) = B0_4_1;
    B0(3,1) = B0_4_2;
    B0(3,2) = B0_4_3;
    B0(3,3) = B0_4_4;
    B0(3,4) = B0_4_5;
    B0(3,5) = B0_4_6;
    B0(3,6) = B0_4_7;
    B0(3,7) = B0_4_8;
    B0(3,8) = B0_4_9;
    B0(3,9) = B0_4_10;
    B0(3,10) = B0_4_11;
    B0(3,11) = B0_4_12;

    B0(4,0) = B0_5_1;
    B0(4,1) = B0_5_2;
    B0(4,2) = B0_5_3;
    B0(4,3) = B0_5_4;
    B0(4,4) = B0_5_5;
    B0(4,5) = B0_5_6;
    B0(4,6) = B0_5_7;
    B0(4,7) = B0_5_8;
    B0(4,8) = B0_5_9;
    B0(4,9) = B0_5_10;
    B0(4,10) = B0_5_11;
    B0(4,11) = B0_5_12;

    B0(5,0) = B0_6_1;
    B0(5,1) = B0_6_2;
    B0(5,2) = B0_6_3;
    B0(5,3) = B0_6_4;
    B0(5,4) = B0_6_5;
    B0(5,5) = B0_6_6;
    B0(5,6) = B0_6_7;
    B0(5,7) = B0_6_8;
    B0(5,8) = B0_6_9;
    B0(5,9) = B0_6_10;
    B0(5,10) = B0_6_11;
    B0(5,11) = B0_6_12;



    double G_1_1 = 0.3e1 / 0.4e1 * s * cos(n * theta) * (-0.1e1 + eta2) / L;
    double G_1_2 = s * cos(n * theta) * (-0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
    double G_1_3 = 0.0e0;
    double G_1_4 = 0.0e0;
    double G_1_5 = 0.3e1 / 0.4e1 * c * cos(n * theta) * (-0.1e1 + eta2) / L;
    double G_1_6 = c * cos(n * theta) * (-0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
    double G_1_7 = -0.3e1 / 0.4e1 * s * cos(n * theta) * (-0.1e1 + eta2) / L;
    double G_1_8 = s * cos(n * theta) * (0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
    double G_1_9 = 0.0e0;
    double G_1_10 = 0.0e0;
    double G_1_11 = -0.3e1 / 0.4e1 * c * cos(n * theta) * (-0.1e1 + eta2) / L;
    double G_1_12 = c * cos(n * theta) * (0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;

    double G_2_1 = 0.0e0;
    double G_2_2 = 0.0e0;
    double G_2_3 = 0.3e1 / 0.4e1 * sin(n * theta) * (-0.1e1 + eta2) / L;
    double G_2_4 = sin(n * theta) * (-0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
    double G_2_5 = 0.0e0;
    double G_2_6 = 0.0e0;
    double G_2_7 = 0.0e0;
    double G_2_8 = 0.0e0;
    double G_2_9 = -0.3e1 / 0.4e1 * sin(n * theta) * (-0.1e1 + eta2) / L;
    double G_2_10 = sin(n * theta) * (0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
    double G_2_11 = 0.0e0;
    double G_2_12 = 0.0e0;

    double G_3_1 = 0.3e1 / 0.4e1 * c * cos(n * theta) * (-0.1e1 + eta2) / L;
    double G_3_2 = c * cos(n * theta) * (-0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
    double G_3_3 = 0.0e0;
    double G_3_4 = 0.0e0;
    double G_3_5 = -0.3e1 / 0.4e1 * s * cos(n * theta) * (-0.1e1 + eta2) / L;
    double G_3_6 = -s * cos(n * theta) * (-0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
    double G_3_7 = -0.3e1 / 0.4e1 * c * cos(n * theta) * (-0.1e1 + eta2) / L;
    double G_3_8 = c * cos(n * theta) * (0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;
    double G_3_9 = 0.0e0;
    double G_3_10 = 0.0e0;
    double G_3_11 = 0.3e1 / 0.4e1 * s * cos(n * theta) * (-0.1e1 + eta2) / L;
    double G_3_12 = -s * cos(n * theta) * (0.2e1 * eta - 0.1e1 + 0.3e1 * eta2) / 0.4e1;

    double G_4_1 = -s * n * sin(n * theta) * (eta + 0.2e1) * pow(eta - 0.1e1, 0.2e1) / r / 0.4e1;
    double G_4_2 = -s * n * sin(n * theta) * L * (eta + 0.1e1) * pow(eta - 0.1e1, 0.2e1) / r / 0.4e1;
    double G_4_3 = -sin(n * theta) * (eta + 0.2e1) * pow(eta - 0.1e1, 0.2e1) * s / r / 0.4e1;
    double G_4_4 = -sin(n * theta) * L * (eta + 0.1e1) * pow(eta - 0.1e1, 0.2e1) * s / r / 0.4e1;
    double G_4_5 = -c * n * sin(n * theta) * (eta + 0.2e1) * pow(eta - 0.1e1, 0.2e1) / r / 0.4e1;
    double G_4_6 = -c * n * sin(n * theta) * L * (eta + 0.1e1) * pow(eta - 0.1e1, 0.2e1) / r / 0.4e1;
    double G_4_7 = s * n * sin(n * theta) * (eta - 0.2e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;
    double G_4_8 = -s * n * sin(n * theta) * L * (eta - 0.1e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;
    double G_4_9 = sin(n * theta) * (-0.3e1 * eta - 0.2e1 + eta3) * s / r / 0.4e1;
    double G_4_10 = -sin(n * theta) * L * (eta - 0.1e1) * pow(eta + 0.1e1, 0.2e1) * s / r / 0.4e1;
    double G_4_11 = c * n * sin(n * theta) * (eta - 0.2e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;
    double G_4_12 = -c * n * sin(n * theta) * L * (eta - 0.1e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;

    double G_5_1 = cos(n * theta) * (eta + 0.2e1) * pow(eta - 0.1e1, 0.2e1) * (c * c + s * s) / r / 0.4e1;
    double G_5_2 = cos(n * theta) * L * (eta + 0.1e1) * pow(eta - 0.1e1, 0.2e1) * (c * c + s * s) / r / 0.4e1;
    double G_5_3 = n * cos(n * theta) * (-0.3e1 * eta + 0.2e1 + eta3) / r / 0.4e1;
    double G_5_4 = n * cos(n * theta) * L * (eta + 0.1e1) * pow(eta - 0.1e1, 0.2e1) / r / 0.4e1;
    double G_5_5 = 0.0e0;
    double G_5_6 = 0.0e0;
    double G_5_7 = -cos(n * theta) * (eta - 0.2e1) * pow(eta + 0.1e1, 0.2e1) * (c * c + s * s) / r / 0.4e1;
    double G_5_8 = cos(n * theta) * L * (eta - 0.1e1) * pow(eta + 0.1e1, 0.2e1) * (c * c + s * s) / r / 0.4e1;
    double G_5_9 = -n * cos(n * theta) * (eta - 0.2e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;
    double G_5_10 = n * cos(n * theta) * L * (eta - 0.1e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;
    double G_5_11 = 0.0e0;
    double G_5_12 = 0.0e0;

    double G_6_1 = -c * n * sin(n * theta) * (eta + 0.2e1) * pow(eta - 0.1e1, 0.2e1) / r / 0.4e1;
    double G_6_2 = -c * n * sin(n * theta) * L * (eta + 0.1e1) * pow(eta - 0.1e1, 0.2e1) / r / 0.4e1;
    double G_6_3 = -sin(n * theta) * (eta + 0.2e1) * pow(eta - 0.1e1, 0.2e1) * c / r / 0.4e1;
    double G_6_4 = -sin(n * theta) * L * (eta + 0.1e1) * pow(eta - 0.1e1, 0.2e1) * c / r / 0.4e1;
    double G_6_5 = s * n * sin(n * theta) * (eta + 0.2e1) * pow(eta - 0.1e1, 0.2e1) / r / 0.4e1;
    double G_6_6 = s * n * sin(n * theta) * L * (eta + 0.1e1) * pow(eta - 0.1e1, 0.2e1) / r / 0.4e1;
    double G_6_7 = c * n * sin(n * theta) * (eta - 0.2e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;
    double G_6_8 = -c * n * sin(n * theta) * L * (eta - 0.1e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;
    double G_6_9 = sin(n * theta) * (-0.3e1 * eta - 0.2e1 + eta3) * c / r / 0.4e1;
    double G_6_10 = -sin(n * theta) * L * (eta - 0.1e1) * pow(eta + 0.1e1, 0.2e1) * c / r / 0.4e1;
    double G_6_11 = -s * n * sin(n * theta) * (eta - 0.2e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;
    double G_6_12 = s * n * sin(n * theta) * L * (eta - 0.1e1) * pow(eta + 0.1e1, 0.2e1) / r / 0.4e1;


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


    double beta1 = ((((0.3e1 * dWdS_1 + 0.3e1 * dWdS_2) * eta2 + (-0.2e1 * dWdS_1 + 0.2e1 * dWdS_2) * eta - dWdS_1 - dWdS_2) * c - 0.2e1 * s * ((-0.3e1 / 0.2e1 * dUdS_1 - 0.3e1 / 0.2e1 * dUdS_2) * eta2 + (dUdS_1 - dUdS_2) * eta + dUdS_1 / 0.2e1 + dUdS_2 / 0.2e1)) * L + 0.3e1 * (-0.1e1 + eta2) * ((W1 - W2) * c + s * (U1 - U2))) / L / 0.4e1;
    double beta2 = 0.0e0;
    double beta3 = ((((0.3e1 * dUdS_1 + 0.3e1 * dUdS_2) * eta2 + (-0.2e1 * dUdS_1 + 0.2e1 * dUdS_2) * eta - dUdS_1 - dUdS_2) * c + 0.2e1 * s * ((-0.3e1 / 0.2e1 * dWdS_1 - 0.3e1 / 0.2e1 * dWdS_2) * eta2 + (dWdS_1 - dWdS_2) * eta + dWdS_1 / 0.2e1 + dWdS_2 / 0.2e1)) * L + 0.3e1 * ((U1 - U2) * c - s * (W1 - W2)) * (-0.1e1 + eta2)) / L / 0.4e1;
    double beta4 = 0.0e0;
    double beta5 = (((dUdS_1 + dUdS_2) * L + U1 - U2) * pow(eta, 0.3e1) - L * (dUdS_1 - dUdS_2) * eta * eta + ((-dUdS_1 - dUdS_2) * L - 0.3e1 * U1 + 0.3e1 * U2) * eta + L * (dUdS_1 - dUdS_2) + 0.2e1 * U1 + 0.2e1 * U2) * (c * c + s * s) / r / 0.4e1;
    double beta6 = 0.0e0;

    Omega = bnu::zero_matrix<double>(6,6);
    Omega(0,0) = beta1;
    Omega(0,1) = beta2;
    Omega(0,2) = beta3;
    Omega(1,3) = beta4;
    Omega(1,4) = beta5;
    Omega(1,5) = beta6;
    Omega(2,0) = beta4;
    Omega(2,1) = beta5;
    Omega(2,2) = beta6;
    Omega(2,3) = beta1;
    Omega(2,4) = beta2;
    Omega(2,5) = beta3;

    bnu::matrix<double> BL = prod(Omega,G);
    bnu::matrix<double> B = B0 + BL;
    bnu::matrix<double> DTB = prod(DT,B);
    bnu::matrix<double> BDTB = prod(trans(B),DTB);
    return BDTB;

}