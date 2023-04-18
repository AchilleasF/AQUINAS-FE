// N.B. compile with MinGW64 on Windows as: mex CXXFLAGS="-O3 -fopenmp -std=c++17" LDOPTIMFLAGS="$LDOPTIMFLAGS -fopenmp" -IC:\boost_1_75_0 AQUINAS_Cpp_interface_LBA.cpp
// N.B. compile with gcc on Linux as: mex CXXFLAGS="-O3 -fopenmp -fPIC -std=c++17" LDOPTIMFLAGS="$LDOPTIMFLAGS -fopenmp" -lboost_system AQUINAS_Cpp_interface_LBA.cpp
// AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.
#define _USE_MATH_DEFINES
#include <iostream>
#include <unordered_map>
#include <omp.h>
#include <math.h>
#include "mex.hpp"
#include "mexAdapter.hpp"
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "AQUINAS_Cpp_Classes.hpp"
#include "AQUINAS_MEX_B0DTB0_matrix.hpp"
#include "AQUINAS_MEX_BDTB_matrix.hpp"
#include "AQUINAS_MEX_GN0G_matrix.hpp"

namespace bnu = boost::numeric::ublas;
matlab::data::ArrayFactory m_Factory;
#pragma omp declare reduction(+ : matlab::data::TypedArray<double> : \
	std::transform(omp_in.begin(), omp_in.end(), \
				   omp_out.begin(), omp_out.begin(), \
				   std::plus<double>() ) ) \
	initializer (omp_priv(omp_orig))

//////////////////
// MEX Function //
//////////////////
class MexFunction : public matlab::mex::Function
{
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs)
    {
        int nElements = std::move(inputs[0][0]); // No. elements in the segment
        matlab::data::TypedArray<double> daSegNode_r = std::move(inputs[1]); // Array of nodal r coordinates
        matlab::data::TypedArray<double> daSegNode_z = std::move(inputs[2]); // Array of nodal z coordinates
        double dThickness = std::move(inputs[3][0]); // Segment thickness
        unsigned int unC_mode = std::move(inputs[4][0]); // Circumferential harmonic
        unsigned int unNoGaussStations = std::move(inputs[5][0]); // No. of Gauss quadrature integration stations
        unsigned int unElCounter = std::move(inputs[6][0]); // Element counter at start of Segment Object
        matlab::data::TypedArray<double> daElDOFs = std::move(inputs[7]); // Incoming element dof data
        matlab::data::TypedArray<double> daSIGMAS = std::move(inputs[8]); // Incoming stress resultant data
        matlab::data::TypedArray<double> daDTs = std::move(inputs[9]); // Incoming element tangent modulus matrices
        matlab::data::TypedArray<unsigned int> daOffset_vals = std::move(inputs[10]); // Offset values
        unsigned int include_nonlinear_G = std::move(inputs[11][0]); // Option to include or not the nonlinear strains in formulating the material stiffness matrix
        unsigned int include_nonlinear_M = std::move(inputs[12][0]); // Option to include or not the nonlinear tangent modulus matrices in formulating the material stiffness matrix
        matlab::data::TypedArray<double> daMTV_M = std::move(inputs[13]); // Material stiffness matrix coefficients in vector form
        matlab::data::TypedArray<double> daMTV_G = std::move(inputs[14]); // Geometric stiffness matrix coefficients in vector form
        matlab::data::TypedArray<double> daOStorage = std::move(inputs[15]); // Storage matrix for pre-computed variables per element
        unsigned int unNoThreads = std::move(inputs[16][0]); // No. of OpenMP threads

        std::unordered_map<std::string, unsigned int> umOffsets; // Hash table for segment offsets
        umOffsets["D0"] = daOffset_vals[0];
        umOffsets["U1"] = daOffset_vals[1];
        umOffsets["L1"] = daOffset_vals[2];
        umOffsets["U2"] = daOffset_vals[3];
        umOffsets["L2"] = daOffset_vals[4];
        umOffsets["U3"] = daOffset_vals[5];
        umOffsets["L3"] = daOffset_vals[6];
        umOffsets["U4"] = daOffset_vals[7];
        umOffsets["L4"] = daOffset_vals[8];
        umOffsets["U5"] = daOffset_vals[9];
        umOffsets["L5"] = daOffset_vals[10];
        umOffsets["U6"] = daOffset_vals[11];
        umOffsets["L6"] = daOffset_vals[12];
        umOffsets["U7"] = daOffset_vals[13];
        umOffsets["L7"] = daOffset_vals[14];
        umOffsets["TT"] = daOffset_vals[15];

        std::vector<double> dvGauss_Nodes, dvGauss_Weights; // Vectors of nodes and weights for Gaussian quadrature

        // Generate Gauss weights and create storage containers for uncondensed K_M and K_G
        this->obtain_Gauss_Legendre_nodes_weights(unNoGaussStations, dvGauss_Nodes, dvGauss_Weights);

        matlab::data::TypedArray<double> daKMStorage = m_Factory.createArray<double>({12,12,(unsigned int)nElements},{0.0});
        matlab::data::TypedArray<double> daKGStorage = m_Factory.createArray<double>({12,12,(unsigned int)nElements},{0.0});


        // Main element loop
        #pragma omp parallel shared(daKMStorage,daKGStorage,daSegNode_r,daSegNode_z,nElements,dThickness,unNoGaussStations,dvGauss_Nodes,dvGauss_Weights,\
            include_nonlinear_G,include_nonlinear_M,unElCounter,unC_mode,umOffsets,daElDOFs,daSIGMAS,daDTs,daOStorage)\
            num_threads(unNoThreads) reduction(+:daMTV_M) reduction(+:daMTV_G)
        {
            // All variables declared within the parallel region are private by default
            double r1, r2, z1, z2;
            double drds1, dzds1, drds2, dzds2;
            double phi1, phi2, dphids1, dphids2, d2phids2_1, d2phids2_2;
            double L, psi, eta, eta2, eta3;
            double Nphi, Ntheta, Nphitheta;
            double N01, N11, N02, N12, r, z, phi, dphids, c, s;
            unsigned int el_off, el_off1, el_off2, el_off3;
            unsigned int no_cond_DOFs = (unC_mode == 0) ? 2 : 4;

            bnu::matrix<double> DT(6,6);
            bnu::matrix<double> dm_KBARE_M(12,12);
            bnu::matrix<double> dm_KBARE_G(12,12);
            double c1, s1, c2, s2;
            double dV1, dV5, dV7, dV11;
            bnu::vector<double> dc_KBARE_R1, dc_KBARE_R5, dc_KBARE_R7, dc_KBARE_R11;
            bnu::vector<double> dr_KBARE_R1, dr_KBARE_R5, dr_KBARE_R7, dr_KBARE_R11;
            bnu::matrix<double> dm_KRR_M(8,8), dm_KRC_M(8,no_cond_DOFs), dm_KCC_M(no_cond_DOFs,no_cond_DOFs);
            bnu::matrix<double> dm_KRR_G(8,8), dm_KRC_G(8,no_cond_DOFs), dm_KCC_G(no_cond_DOFs,no_cond_DOFs);
            bnu::matrix<double> dm_K_M(8,8), dm_K_G(8,8), dm_Kpr1(8,8), dm_Kpr2(8,8);
            bnu::vector<double> elDOFs(12);

            if (include_nonlinear_M==0)
            {
                for (unsigned int I = 0; I < 6; ++I)
                {
                    for (unsigned int J = 0; J < 6; ++J)
                    {
                        DT(I,J) = daDTs[I][J];
                    }
                }
            }

            #pragma omp for schedule(dynamic,1)
            for (int E = 0; E < nElements; ++E)
            {
                // Recall element geometry
                r1 = daSegNode_r[E]; z1 = daSegNode_z[E]; // Element top node r-z coordinates
                r2 = daSegNode_r[E+1]; z2 = daSegNode_z[E+1]; // Element bottom node r-z coordinates
                drds1 = daOStorage[1][E]; drds2 = daOStorage[0][E];
                dzds1 = daOStorage[3][E]; dzds2 = daOStorage[2][E];
                phi1 = daOStorage[5][E]; phi2 = daOStorage[4][E];
                dphids1 = daOStorage[7][E]; dphids2 = daOStorage[6][E];
                d2phids2_1 = daOStorage[9][E]; d2phids2_2 = daOStorage[8][E];
                L = daOStorage[10][E]; psi = daOStorage[11][E];

                // Nodal dof values from axisymmetric prebuckling state
                for (unsigned int I = 0; I < 12; ++I)
                {
                    elDOFs(I) = daElDOFs[I][E+unElCounter];
                }

                // Integrate element stiffness matrix and nodal load vector
                for (unsigned int I = 0; I < 12; ++I)
                {
                    for (unsigned int J = 0; J < 12; ++J)
                    {
                        dm_KBARE_M(I,J) = 0.0;
                        dm_KBARE_G(I,J) = 0.0;
                    }
                }
                for (unsigned int P = 0; P < unNoGaussStations; ++P)
                {
                    if (include_nonlinear_M == 1)
                    {
                        for (unsigned int I = 0; I < 6; ++I)
                        {
                            for (unsigned int J = 0; J < 6; ++J)
                            {
                                DT(I,J) = daDTs[E+unElCounter][P][I][J];
                            }
                        }
                    }
                    eta = dvGauss_Nodes[P];
                    eta2 = eta * eta;
                    eta3 = eta2 * eta;
                    N01 = (2. - 3. * eta + eta3) * 0.25;
                    N11 = L * (1. - eta - eta2 + eta3) * 0.25;
                    N02 = (2. + 3. * eta - eta3) * 0.25;
                    N12 = L * (-1. - eta + eta2 + eta3) * 0.25;
                    r = N01 * r1 + N02 * r2 + N11 * drds1 + N12 * drds2;
                    z = N01 * z1 + N02 * z2 + N11 * dzds1 + N12 * dzds2;
                    phi = N01 * phi1 + N02 * phi2 + N11 * dphids1 + N12 * dphids2;
                    dphids = N01 * dphids1 + N02 * dphids2 + N11 * d2phids2_1 + N12 * d2phids2_2;
                    Nphi = daSIGMAS[0][E+unElCounter][P]; Ntheta = daSIGMAS[1][E+unElCounter][P]; Nphitheta = daSIGMAS[2][E+unElCounter][P];
                    c = cos(phi);
                    s = sin(phi);
                    if (unC_mode == 0)
                    {
                        if (include_nonlinear_G == 0)
                        {
                            dm_KBARE_M += dvGauss_Weights[P] * L * r * AQUINAS_B0DTB0(eta, 0.0, L, r, dphids, c, s, (double)unC_mode, psi, DT);
                        }
                        else if (include_nonlinear_G == 1)
                        {
                            dm_KBARE_M += dvGauss_Weights[P] * L * r * AQUINAS_BDTB(eta, 0.0, L, r, dphids, c, s, (double)unC_mode, psi, elDOFs, DT);
                        }
                        dm_KBARE_G += dvGauss_Weights[P] * L * r * AQUINAS_GN0G(eta, 0.0, L, r, dphids, c, s, (double)unC_mode, psi, Nphi, Ntheta, Nphitheta);
                    }
                    else
                    {
                        if (include_nonlinear_G == 0)
                        {
                            dm_KBARE_M += dvGauss_Weights[P] * L * r * (AQUINAS_B0DTB0(eta, 0.0, L, r, dphids, c, s, (double)unC_mode, psi, DT) \
                                + AQUINAS_B0DTB0(eta, M_PI_2/(double)unC_mode, L, r, dphids, c, s, (double)unC_mode, psi, DT));
                        }
                        else if (include_nonlinear_G == 1)
                        {
                            dm_KBARE_M += dvGauss_Weights[P] * L * r * (AQUINAS_BDTB(eta, 0.0, L, r, dphids, c, s, (double)unC_mode, psi, elDOFs, DT) \
                                + AQUINAS_BDTB(eta, M_PI_2/(double)unC_mode, L, r, dphids, c, s, (double)unC_mode, psi, elDOFs, DT));
                        }
                        dm_KBARE_G += dvGauss_Weights[P] * L * r * (AQUINAS_GN0G(eta, 0.0, L, r, dphids, c, s, (double)unC_mode, psi, Nphi, Ntheta, Nphitheta) \
                            + AQUINAS_GN0G(eta, M_PI_2/(double)unC_mode, L, r, dphids, c, s, (double)unC_mode, psi, Nphi, Ntheta, Nphitheta));
                    }
                }
                for (unsigned int I = 0; I < 12; ++I)
                {
                    for (unsigned int J = 0; J < 12; ++J)
                    {
                        daKMStorage[I][J][E] = dm_KBARE_M(I,J);
                        daKGStorage[I][J][E] = dm_KBARE_G(I,J);
                    }
                }

                // Transform element matrices
                c1 = cos(phi1); s1 = sin(phi1); c2 = cos(phi2); s2 = sin(phi2);
                dc_KBARE_R1 = column(dm_KBARE_M, 1);
                dc_KBARE_R5 = column(dm_KBARE_M, 5);
                dc_KBARE_R7 = column(dm_KBARE_M, 7);
                dc_KBARE_R11 = column(dm_KBARE_M, 11);
                column(dm_KBARE_M, 1) = -s1 * dc_KBARE_R1 - c1 * dc_KBARE_R5;
                column(dm_KBARE_M, 5) = c1 * dc_KBARE_R1 - s1 * dc_KBARE_R5;
                column(dm_KBARE_M, 7) = -s2 * dc_KBARE_R7 - c2 * dc_KBARE_R11;
                column(dm_KBARE_M,11) = c2 * dc_KBARE_R7 - s2 * dc_KBARE_R11;
                dc_KBARE_R1 = column(dm_KBARE_G, 1);
                dc_KBARE_R5 = column(dm_KBARE_G, 5);
                dc_KBARE_R7 = column(dm_KBARE_G, 7);
                dc_KBARE_R11 = column(dm_KBARE_G, 11);
                column(dm_KBARE_G, 1) = -s1 * dc_KBARE_R1 - c1 * dc_KBARE_R5;
                column(dm_KBARE_G, 5) = c1 * dc_KBARE_R1 - s1 * dc_KBARE_R5;
                column(dm_KBARE_G, 7) = -s2 * dc_KBARE_R7 - c2 * dc_KBARE_R11;
                column(dm_KBARE_G,11) = c2 * dc_KBARE_R7 - s2 * dc_KBARE_R11;

                dr_KBARE_R1 = row(dm_KBARE_M, 1);
                dr_KBARE_R5 = row(dm_KBARE_M, 5);
                dr_KBARE_R7 = row(dm_KBARE_M, 7);
                dr_KBARE_R11 = row(dm_KBARE_M, 11);
                row(dm_KBARE_M, 1) = -s1 * dr_KBARE_R1 - c1 * dr_KBARE_R5;
                row(dm_KBARE_M, 5) = c1 * dr_KBARE_R1 - s1 * dr_KBARE_R5;
                row(dm_KBARE_M, 7) = -s2 * dr_KBARE_R7 - c2 * dr_KBARE_R11;
                row(dm_KBARE_M,11) = c2 * dr_KBARE_R7 - s2 * dr_KBARE_R11;
                dr_KBARE_R1 = row(dm_KBARE_G, 1);
                dr_KBARE_R5 = row(dm_KBARE_G, 5);
                dr_KBARE_R7 = row(dm_KBARE_G, 7);
                dr_KBARE_R11 = row(dm_KBARE_G, 11);
                row(dm_KBARE_G, 1) = -s1 * dr_KBARE_R1 - c1 * dr_KBARE_R5;
                row(dm_KBARE_G, 5) = c1 * dr_KBARE_R1 - s1 * dr_KBARE_R5;
                row(dm_KBARE_G, 7) = -s2 * dr_KBARE_R7 - c2 * dr_KBARE_R11;
                row(dm_KBARE_G,11) = c2 * dr_KBARE_R7 - s2 * dr_KBARE_R11;

                // Partition element matrix
                if (unC_mode == 0)
                {
                    std::vector<unsigned int> vR_IDs{0,2,4,1,6,8,10,7}, vC_IDs{5,11};
                    for (unsigned int I = 0; I < 8; ++I)
                    {
                        for (unsigned int J = 0; J < 8; ++J)
                        {
                            dm_KRR_M(I,J) = dm_KBARE_M(vR_IDs[I],vR_IDs[J]);
                            dm_KRR_G(I,J) = dm_KBARE_G(vR_IDs[I],vR_IDs[J]);
                        }
                        for (unsigned int J = 0; J < no_cond_DOFs; ++J)
                        {
                            dm_KRC_M(I,J) = dm_KBARE_M(vR_IDs[I],vC_IDs[J]);
                            dm_KRC_G(I,J) = dm_KBARE_G(vR_IDs[I],vC_IDs[J]);
                            if (I < no_cond_DOFs)
                            {
                                dm_KCC_M(I,J) = dm_KBARE_M(vC_IDs[I],vC_IDs[J]);
                                dm_KCC_G(I,J) = dm_KBARE_G(vC_IDs[I],vC_IDs[J]);
                            }
                        }
                    }
                }
                else
                {
                    std::vector<unsigned int> vR_IDs{0,2,4,1,6,8,10,7}, vC_IDs{3,5,9,11};
                    for (unsigned int I = 0; I < 8; ++I)
                    {
                        for (unsigned int J = 0; J < 8; ++J)
                        {
                            dm_KRR_M(I,J) = dm_KBARE_M(vR_IDs[I],vR_IDs[J]);
                            dm_KRR_G(I,J) = dm_KBARE_G(vR_IDs[I],vR_IDs[J]);
                        }
                        for (unsigned int J = 0; J < no_cond_DOFs; ++J)
                        {
                            dm_KRC_M(I,J) = dm_KBARE_M(vR_IDs[I],vC_IDs[J]);
                            dm_KRC_G(I,J) = dm_KBARE_G(vR_IDs[I],vC_IDs[J]);
                            if (I < no_cond_DOFs)
                            {
                                dm_KCC_M(I,J) = dm_KBARE_M(vC_IDs[I],vC_IDs[J]);
                                dm_KCC_G(I,J) = dm_KBARE_G(vC_IDs[I],vC_IDs[J]);
                            }
                        }
                    }
                }

                bnu::matrix<double> dm_KCC_M_inv = bnu::identity_matrix<float>(dm_KCC_M.size1());
                bnu::permutation_matrix<size_t> pm_M(dm_KCC_M.size1());
                bnu::lu_factorize(dm_KCC_M, pm_M);
                bnu::lu_substitute(dm_KCC_M, pm_M, dm_KCC_M_inv);

                // Condense element matrix
                dm_Kpr1 = prod(dm_KRC_M, dm_KCC_M_inv);
                dm_K_M = dm_KRR_M - prod(dm_Kpr1, trans(dm_KRC_M));
                dm_Kpr1 = prod(dm_KCC_M_inv, trans(dm_KRC_M));
                dm_Kpr1 = prod(dm_KRC_G, dm_Kpr1);
                dm_Kpr2 = prod(dm_KCC_M_inv,trans(dm_KRC_M));
                dm_Kpr2 = prod(dm_KCC_G, dm_Kpr2);
                dm_Kpr2 = prod(dm_KCC_M_inv, dm_Kpr2);
                dm_Kpr2 = prod(dm_KRC_M, dm_Kpr2);
                dm_K_G = dm_KRR_G - dm_Kpr1 - trans(dm_Kpr1) + dm_Kpr2;

                // Insert condensed matrix / vector coefficients into vector form containers
                el_off = (E + 1 - 1) * 4; el_off1 = (E + 1 - 1) * 3; el_off2 = (E + 1 - 1) * 2; el_off3 = (E + 1 - 1);
                for (unsigned int I = 0; I < 8; ++I)
                {
                    daMTV_M[I + el_off + umOffsets["D0"]] += dm_K_M(I,I); // Main diagonal
                    daMTV_G[I + el_off + umOffsets["D0"]] += dm_K_G(I,I); // Main diagonal
                }
                for (unsigned int I = 0; I < 7; ++I)
                {
                    daMTV_M[I + el_off + umOffsets["U1"]] += dm_K_M(I,I+1); // First upper off-diagonal
                    daMTV_M[I + el_off + umOffsets["L1"]] += dm_K_M(I+1,I); // First lower off-diagonal
                    daMTV_G[I + el_off + umOffsets["U1"]] += dm_K_G(I,I+1); // First upper off-diagonal
                    daMTV_G[I + el_off + umOffsets["L1"]] += dm_K_G(I+1,I); // First lower off-diagonal
                }
                for (unsigned int I = 0; I < 6; ++I)
                {
                    daMTV_M[I + el_off + umOffsets["U2"]] += dm_K_M(I,I+2); // Second upper off-diagonal
                    daMTV_M[I + el_off + umOffsets["L2"]] += dm_K_M(I+2,I); // Second lower off-diagonal
                    daMTV_G[I + el_off + umOffsets["U2"]] += dm_K_G(I,I+2); // Second upper off-diagonal
                    daMTV_G[I + el_off + umOffsets["L2"]] += dm_K_G(I+2,I); // Second lower off-diagonal
                }
                for (unsigned int I = 0; I < 5; ++I)
                {
                    daMTV_M[I + el_off + umOffsets["U3"]] += dm_K_M(I,I+3); // Third upper off-diagonal
                    daMTV_M[I + el_off + umOffsets["L3"]] += dm_K_M(I+3,I); // Third lower off-diagonal
                    daMTV_G[I + el_off + umOffsets["U3"]] += dm_K_G(I,I+3); // Third upper off-diagonal
                    daMTV_G[I + el_off + umOffsets["L3"]] += dm_K_G(I+3,I); // Third lower off-diagonal
                }
                for (unsigned int I = 0; I < 4; ++I)
                {
                    daMTV_M[I + el_off + umOffsets["U4"]] += dm_K_M(I,I+4); // Fourth upper off-diagonal
                    daMTV_M[I + el_off + umOffsets["L4"]] += dm_K_M(I+4,I); // Fourth lower off-diagonal
                    daMTV_G[I + el_off + umOffsets["U4"]] += dm_K_G(I,I+4); // Fourth upper off-diagonal
                    daMTV_G[I + el_off + umOffsets["L4"]] += dm_K_G(I+4,I); // Fourth lower off-diagonal
                }
                for (unsigned int I = 0; I < 3; ++I)
                {
                    daMTV_M[I + el_off1 + umOffsets["U5"]] += dm_K_M(I,I+5); // Fifth upper off-diagonal
                    daMTV_M[I + el_off1 + umOffsets["L5"]] += dm_K_M(I+5,I); // Fifth lower off-diagonal
                    daMTV_G[I + el_off1 + umOffsets["U5"]] += dm_K_G(I,I+5); // Fifth upper off-diagonal
                    daMTV_G[I + el_off1 + umOffsets["L5"]] += dm_K_G(I+5,I); // Fifth lower off-diagonal
                }
                for (unsigned int I = 0; I < 2; ++I)
                {
                    daMTV_M[I + el_off2 + umOffsets["U6"]] += dm_K_M(I,I+6); // Sixth upper off-diagonal
                    daMTV_M[I + el_off2 + umOffsets["L6"]] += dm_K_M(I+6,I); // Sixth lower off-diagonal
                    daMTV_G[I + el_off2 + umOffsets["U6"]] += dm_K_G(I,I+6); // Sixth upper off-diagonal
                    daMTV_G[I + el_off2 + umOffsets["L6"]] += dm_K_G(I+6,I); // Sixth lower off-diagonal
                }
                daMTV_M[0 + el_off3 + umOffsets["U7"]] += dm_K_M(0,0+7); // Seventh upper off-diagonal
                daMTV_M[0 + el_off3 + umOffsets["L7"]] += dm_K_M(0+7,0); // Seventh lower off-diagonal
                daMTV_G[0 + el_off3 + umOffsets["U7"]] += dm_K_G(0,0+7); // Seventh upper off-diagonal
                daMTV_G[0 + el_off3 + umOffsets["L7"]] += dm_K_G(0+7,0); // Seventh lower off-diagonal
            } // End of parallel loop
        } // End of parallel region

        outputs[0] = daMTV_M;
        outputs[1] = daMTV_G;
        outputs[2] = daKMStorage;
        outputs[3] = daKGStorage;
    }

    void obtain_Gauss_Legendre_nodes_weights(unsigned int unNoGaussStations, std::vector<double>& dvGauss_Nodes, std::vector<double>& dvGauss_Weights)
    {
        switch (unNoGaussStations)
        {
            case 1:
                dvGauss_Nodes.push_back( 0.0 );
                dvGauss_Weights.push_back( 2.0 );
                break;
            case 2:
                dvGauss_Nodes.push_back(-0.577350269189626 );
                dvGauss_Nodes.push_back( 0.577350269189626 );
                dvGauss_Weights.push_back( 1.0 );
                dvGauss_Weights.push_back( 1.0 );
                break;
            case 3:
                dvGauss_Nodes.push_back(-0.774596669241483 );
                dvGauss_Nodes.push_back( 0.0 );
                dvGauss_Nodes.push_back( 0.774596669241483 );
                dvGauss_Weights.push_back( 0.555555555555556 );
                dvGauss_Weights.push_back( 0.888888888888889 );
                dvGauss_Weights.push_back( 0.555555555555556 );
                break;
            case 4:
                dvGauss_Nodes.push_back(-0.861136311594053 );
                dvGauss_Nodes.push_back(-0.339981043584856 );
                dvGauss_Nodes.push_back( 0.339981043584856 );
                dvGauss_Nodes.push_back( 0.861136311594053 );
                dvGauss_Weights.push_back( 0.347854845137454 );
                dvGauss_Weights.push_back( 0.652145154862546 );
                dvGauss_Weights.push_back( 0.652145154862546 );
                dvGauss_Weights.push_back( 0.347854845137454 );
                break;
            case 5:
                dvGauss_Nodes.push_back(-0.906179845938664 );
                dvGauss_Nodes.push_back(-0.538469310105683 );
                dvGauss_Nodes.push_back( 0.0 );
                dvGauss_Nodes.push_back( 0.538469310105683 );
                dvGauss_Nodes.push_back( 0.906179845938664 );
                dvGauss_Weights.push_back( 0.236926885056189 );
                dvGauss_Weights.push_back( 0.478628670499366 );
                dvGauss_Weights.push_back( 0.568888888888889 );
                dvGauss_Weights.push_back( 0.478628670499366 );
                dvGauss_Weights.push_back( 0.236926885056189 );
                break;
            case 6:
                dvGauss_Nodes.push_back(-0.9324695142031521 );
                dvGauss_Nodes.push_back(-0.6612093864662645 );
                dvGauss_Nodes.push_back(-0.2386191860831969 );
                dvGauss_Nodes.push_back( 0.2386191860831969 );
                dvGauss_Nodes.push_back( 0.6612093864662645 );
                dvGauss_Nodes.push_back( 0.9324695142031521 );
                dvGauss_Weights.push_back( 0.1713244923791704 );
                dvGauss_Weights.push_back( 0.3607615730481386 );
                dvGauss_Weights.push_back( 0.4679139345726910 );
                dvGauss_Weights.push_back( 0.4679139345726910 );
                dvGauss_Weights.push_back( 0.3607615730481386 );
                dvGauss_Weights.push_back( 0.1713244923791704 );
                break;
            case 7:
                dvGauss_Nodes.push_back(-0.9491079123427585 );
                dvGauss_Nodes.push_back(-0.7415311855993945 );
                dvGauss_Nodes.push_back(-0.4058451513773972 );
                dvGauss_Nodes.push_back( 0.0 );
                dvGauss_Nodes.push_back( 0.4058451513773972 );
                dvGauss_Nodes.push_back( 0.7415311855993945 );
                dvGauss_Nodes.push_back( 0.9491079123427585 );
                dvGauss_Weights.push_back( 0.1294849661688697 );
                dvGauss_Weights.push_back( 0.2797053914892766 );
                dvGauss_Weights.push_back( 0.3818300505051189 );
                dvGauss_Weights.push_back( 0.4179591836734694 );
                dvGauss_Weights.push_back( 0.3818300505051189 );
                dvGauss_Weights.push_back( 0.2797053914892766 );
                dvGauss_Weights.push_back( 0.1294849661688697 );
                break;
        }
    }
};


// BSD 3-Clause License
//  
// Copyright (c) 2023, Mr Achilleas Filippidis and Dr Adam Jan Sadowski of 
// Imperial College London. All rights reserved.
//  
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//  
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//  
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//  
// 3. Neither the name of the copyright holder, nor of Imperial College, nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//  
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//  
// THE USE OF THE MATLAB LANGUAGE DOES NOT IMPLY ENDORSEMENT BY MATHWORKS.
