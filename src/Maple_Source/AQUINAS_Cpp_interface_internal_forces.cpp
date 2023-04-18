// N.B. compile with MinGW64 on Windows as: mex CXXFLAGS="-O3 -fopenmp -std=c++17" LDOPTIMFLAGS="$LDOPTIMFLAGS -fopenmp" -IC:\boost_1_75_0 AQUINAS_Cpp_interface_internal_forces.cpp
// N.B. compile with gcc on Linux as: mex CXXFLAGS="-O3 -fopenmp -fPIC -std=c++17" LDOPTIMFLAGS="$LDOPTIMFLAGS -fopenmp" -lboost_system AQUINAS_Cpp_interface_internal_forces.cpp
// AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
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
#include "AQUINAS_MEX_Fint_vector.hpp"
#include "AQUINAS_MEX_strains_vector.hpp"

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
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
    void displayOnMATLAB(std::ostringstream& stream) {
        // Pass stream content to MATLAB fprintf function
        matlabPtr->feval(u"fprintf", 0,
            std::vector<matlab::data::Array>({ m_Factory.createScalar(stream.str()) }));
        // Clear stream buffer
        stream.str("");
    }


    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs)
    {
        int nElements = std::move(inputs[0][0]); // No. elements in the segment
        matlab::data::TypedArray<double> daSegNode_r = std::move(inputs[1]); // Array of nodal r coordinates
        matlab::data::TypedArray<double> daSegNode_z = std::move(inputs[2]); // Array of nodal z coordinates
        double dThickness = std::move(inputs[3][0]); // Segment thickness
        double dnu = std::move(inputs[4][0]); // Poisson Ratio
        double dE = std::move(inputs[4][1]); // Elastic stiffness modulus
        double dG = std::move(inputs[4][2]); // Shear stiffness modulus
        unsigned int unNoMatCurvePoints = std::move(inputs[5][0]); // No. of Gauss quadrature integration stations
        matlab::data::TypedArray<double> daeps = std::move(inputs[6]); // Array holding the plastic strains for the ture stress - strain curve of the material
        matlab::data::TypedArray<double> dasys = std::move(inputs[7]); // Array holding the yield stresses for the ture stress - strain curve of the material
        unsigned int unNoGaussStations = std::move(inputs[8][0]); // No. of Gauss quadrature integration stations
        unsigned int unNoSimpsonStations = std::move(inputs[9][0]); // No. of stations for through thickness integration according to Simpson's 1/3 rule
        unsigned int unElCounter = std::move(inputs[10][0]); // Element counter at start of Segment Object
        matlab::data::TypedArray<double> daElDOFs = std::move(inputs[11]); // Incoming element dof data
        matlab::data::TypedArray<unsigned int> daAlphasIn = std::move(inputs[12]); // Incoming stress resultant data
        matlab::data::TypedArray<double> daSigmasIn = std::move(inputs[13]); // Incoming stresses per element, per Gauss point and per Simpson station
        matlab::data::TypedArray<double> daSigma_ysIn = std::move(inputs[14]); // Incoming yield stresses, for comparison with the von Mises criterion, per element, per Gauss point and per Simpson station
        matlab::data::TypedArray<double> daEpsilonsIn = std::move(inputs[15]); // Incoming strains per element, per Gauss point and per Simpson station
        matlab::data::TypedArray<double> daDepnsIn = std::move(inputs[16]); // Incoming incremental equivalent plastic strains per element, per Gauss point and per Simpson station
        matlab::data::TypedArray<double> daEpnsIn = std::move(inputs[17]); // Incoming cummulative equivalent plastic strains per element, per Gauss point and per Simpson station
        double epsilon_s = std::move(inputs[18][0]); // Reference strain to determine number of sub increments for materially nonlinear analyses
        double max_epsilon_bar = std::move(inputs[18][1]); // Maximum effective strain increment that may be attempted to be solved for with the sub-increment method (corresponding to a number of sub-increments Nsb = max_epsilon_bar/epsilon_s sub-increments)
        unsigned int include_nonlinear_G = std::move(inputs[19][0]); // Option to include or not the nonlinear strains in formulating the material stiffness matrix
        unsigned int include_nonlinear_M = std::move(inputs[20][0]); // Option to include or not the nonlinear tangent modulus matrices in formulating the material stiffness matrix
        matlab::data::TypedArray<double> daVCV = std::move(inputs[21]); // Nodal load vector coefficients
        matlab::data::TypedArray<double> daMTV = std::move(inputs[22]); // Matrix stiffness coefficients in vector form
        matlab::data::TypedArray<double> daFDistr = std::move(inputs[23]); // Storage matrix for element internal force vectors
        matlab::data::TypedArray<double> daKStorage = std::move(inputs[24]); // Storage matrix for pre-computed stiffness matrices per element
        matlab::data::TypedArray<double> daOStorage = std::move(inputs[25]); // Storage matrix for pre-computed variables per element
        matlab::data::TypedArray<unsigned int> daOffset_vals = std::move(inputs[26]); // Offset values
        unsigned int unNoThreads = std::move(inputs[27][0]); // No. of OpenMP threads

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

        // Create Material Object, Gauss weights and create storage containers for uncondensed K and F
        cMaterial* p_cMaterial = new cMaterial(0, dE, dnu, dG, unNoMatCurvePoints, dasys, daeps, epsilon_s, max_epsilon_bar); // Pointer to a cMaterial object
        this->obtain_Gauss_Legendre_nodes_weights(unNoGaussStations, dvGauss_Nodes, dvGauss_Weights);

        matlab::data::TypedArray<int> daSubFlags = m_Factory.createArray<int>({(unsigned int)nElements,unNoGaussStations,unNoSimpsonStations},{0});
        matlab::data::TypedArray<double> daFINTER = m_Factory.createArray<double>({12,(unsigned int)nElements},{0.0});
        matlab::data::TypedArray<double> daSIGMAS = m_Factory.createArray<double>({6,(unsigned int)nElements,(unsigned int)unNoGaussStations},{0.0});
        matlab::data::TypedArray<double> daEpsilons = m_Factory.createArray<double>({6,(unsigned int)nElements,(unsigned int)unNoGaussStations},{0.0});
        matlab::data::TypedArray<double> daSigmas = m_Factory.createArray<double>({3,(unsigned int)nElements,(unsigned int)unNoGaussStations,(unsigned int)unNoSimpsonStations},{0.0});
        matlab::data::TypedArray<unsigned int> daAlphas = m_Factory.createArray<unsigned int>({(unsigned int)nElements,(unsigned int)unNoGaussStations,(unsigned int)unNoSimpsonStations},{0});
        matlab::data::TypedArray<double> daSigma_ys = m_Factory.createArray<double>({(unsigned int)nElements,(unsigned int)unNoGaussStations,(unsigned int)unNoSimpsonStations},{0});
        matlab::data::TypedArray<double> daEpns = m_Factory.createArray<double>({(unsigned int)nElements,(unsigned int)unNoGaussStations,(unsigned int)unNoSimpsonStations},{0});
        matlab::data::TypedArray<double> daDepns = m_Factory.createArray<double>({(unsigned int)nElements,(unsigned int)unNoGaussStations,(unsigned int)unNoSimpsonStations},{0});
        for (int E = 0; E < nElements; ++E)
        {
            for (int I = 0; I < unNoGaussStations; ++I)
            {
                for (int J = 0; J < unNoSimpsonStations; ++J)
                {
                    daAlphas[E][I][J] = daAlphasIn[E][I][J];
                    daEpns[E][I][J] = daEpnsIn[E][I][J];
                    daDepns[E][I][J] = daDepnsIn[E][I][J];
                    daSigma_ys[E][I][J] = daSigma_ysIn[E][I][J];
                    for (int i = 0; i < 3; ++i)
                        daSigmas[i][E][I][J] = daSigmasIn[i][E][I][J];
                }
                for (int i = 0; i < 6; ++i)
                    daEpsilons[i][E][I] = daEpsilonsIn[i][E][I];
            }
        }
        // Main element loop
        #pragma omp parallel shared(daFINTER,daSegNode_z,daSegNode_r,nElements,dThickness,unNoGaussStations,unNoSimpsonStations,p_cMaterial,unElCounter,umOffsets,\
            daElDOFs,daAlphas,daSIGMAS,daSigmas,daEpsilons,daDepns,daEpns,daSigma_ys,dvGauss_Nodes,dvGauss_Weights,daKStorage,daFDistr,daOStorage,include_nonlinear_G,include_nonlinear_M,daSubFlags)\
            num_threads(unNoThreads) reduction(+:daVCV) reduction(+:daMTV)
        {
            // All variables declared within the parallel region are private by default
            double z1, z2, r1, r2;
            double drds1, dzds1, drds2, dzds2;
            double phi1, phi2, dphids1, dphids2, d2phids2_1, d2phids2_2;
            double L, psi, eta, eta2, eta3;
            double N01, N11, N02, N12, r, z, phi, dphids, c, s;
            unsigned int el_off, el_off1, el_off2, el_off3;

            bnu::matrix<double> DT(6, 6);
            bnu::vector<double> epsilons(6), epsilons_prev(6), Delta_epsilon(6), epsilon0Z(3), DeltaepsilonZ(3), SIGMAintegral(6);
            bnu::matrix<double> Z_J, Z_a, Z_m, Z_b;
            bnu::matrix<double> dm_KBARE(12,12);
            bnu::vector<double> dv_FDISTR(12);
            bnu::vector<double> dv_FINTER(12);
            bnu::vector<double> dv_Fel(12);
            double c1, s1, c2, s2;
            double dV1, dV5, dV7, dV11;
            bnu::vector<double> dc_KBARE_R1, dc_KBARE_R5, dc_KBARE_R7, dc_KBARE_R11;
            bnu::vector<double> dr_KBARE_R1, dr_KBARE_R5, dr_KBARE_R7, dr_KBARE_R11;
            bnu::matrix<double> dm_KRR(8,8), dm_KRC(8,2), dm_KCC(2,2);
            bnu::vector<double> dv_FR(8), dv_FC(2);
            bnu::matrix<double> dm_K(8,8), dm_Kpr1(8,8);
            bnu::vector<double> dv_F(8);
            bnu::vector<double> elDOFs(12);

            if (include_nonlinear_M == 0) { DT = p_cMaterial->DTe(dThickness); }

            #pragma	omp for schedule(static,1)
            for (int E = 0; E < nElements; ++E)
            {
                // Recall element geometry
                r1 = daSegNode_r[E]; z1 = daSegNode_z[E]; // Element node 1 r-z coordinates
                r2 = daSegNode_r[E+1]; z2 = daSegNode_z[E+1]; // Element node 2 r-z coordinates
                drds1 = daOStorage[1][E]; drds2 = daOStorage[0][E];
                dzds1 = daOStorage[3][E]; dzds2 = daOStorage[2][E];
                phi2 = daOStorage[4][E]; phi1 = daOStorage[5][E];
                dphids2 = daOStorage[6][E]; dphids1 = daOStorage[7][E];
                d2phids2_2 = daOStorage[8][E]; d2phids2_1 = daOStorage[9][E];
                L = daOStorage[10][E]; psi = daOStorage[11][E];

                // Nodal dof values at element ends for axisymmetric pre-buckling state
                for (unsigned int I = 0; I < 12; ++I)
                    elDOFs(I) = daElDOFs[I][unElCounter + E];

                // Initialize internal element forces vector to zeros
                dv_Fel = bnu::zero_vector<double>(12);
                dv_FINTER = bnu::zero_vector<double>(12);
                // Integrate element stiffness matrix and nodal load vector
                for (unsigned int I = 0; I < unNoGaussStations; ++I)
                {
                    eta = dvGauss_Nodes[I];
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
                    c = cos(phi);
                    s = sin(phi);

                    if (include_nonlinear_M == 1)
                    {
                        for (int i = 0; i < 6; ++i) { epsilons_prev(i) = daEpsilons[i][E][I]; } // Store previous strains
                        epsilons = AQUINAS_strains_vector(eta, L, r, dphids, c, s, psi, elDOFs, include_nonlinear_G);
                        Delta_epsilon = epsilons - epsilons_prev;
                        for (int J = 0; J < unNoSimpsonStations; ++J)
                        {
                            double zJ = -dThickness/2 + J*dThickness/(unNoSimpsonStations-1);
                            Z_J = p_cMaterial->Z(zJ);
                            epsilon0Z = prod(Z_J,epsilons_prev);
                            DeltaepsilonZ = prod(Z_J,Delta_epsilon);
                            p_cMaterial->sub_incremental_computations(daSubFlags, epsilon0Z, DeltaepsilonZ, daSigmas, daSigma_ys, daEpns, daDepns,  daAlphas, E, I, J);
                        }
                        for (int i = 0; i < 6; ++i) { daEpsilons[i][E][I] = epsilons(i); } // Update strains
                        SIGMAintegral = bnu::zero_vector<double>(6);
                        int noLayers = unNoSimpsonStations - 1;
                        // Through thickness integration of stresses to compute stress resultants
                        for (int J = 0; J < (int)(noLayers / 2); ++J)
                        {
                            int a = 2 * J; int m = 2 * J + 1; int b = 2 * J + 2;
                            double z_a = -dThickness / 2 + (double)a*dThickness/(double)noLayers; double z_m = -dThickness / 2 + (double)m * dThickness / (double)noLayers; double z_b = -dThickness / 2 + (double)b * dThickness / (double)noLayers;
                            Z_a = p_cMaterial->Z(z_a);
                            Z_m = p_cMaterial->Z(z_m);
                            Z_b = p_cMaterial->Z(z_b);
                            bnu::vector<double> sigma_a = bnu::zero_vector<double>(3);
                            bnu::vector<double> sigma_m = bnu::zero_vector<double>(3);
                            bnu::vector<double> sigma_b = bnu::zero_vector<double>(3);
                            for (int i = 0; i < 3; ++i)
                            {
                                sigma_a(i) = daSigmas[i][E][I][a];
                                sigma_m(i) = daSigmas[i][E][I][m];
                                sigma_b(i) = daSigmas[i][E][I][b];
                            }
                            bnu::vector<double> SIGMA_a = prod(trans(sigma_a),Z_a);
                            bnu::vector<double> SIGMA_m = prod(trans(sigma_m),Z_m);
                            bnu::vector<double> SIGMA_b = prod(trans(sigma_b),Z_b);
                            SIGMAintegral += (z_b-z_a)*trans(SIGMA_a + 4*SIGMA_m + SIGMA_b)/6;
                        }
                    }
                    else
                    {
                        epsilons = AQUINAS_strains_vector(eta, L, r, dphids, c, s, psi, elDOFs, include_nonlinear_G);
                        SIGMAintegral = prod(DT, epsilons);
                    }
                    for (int i = 0; i < 6; ++i) { daSIGMAS[i][E][I] = SIGMAintegral(i); } // Update stress resultants
                    dv_FINTER += dvGauss_Weights[I] * L * r * AQUINAS_Fint_vector(eta, L, r, dphids, c, s, psi, elDOFs, SIGMAintegral, include_nonlinear_G);
                }
                // Store element internal force vector
                for (unsigned int I = 0; I < 12; ++I)
                {
                    daFINTER[I][E] = dv_FINTER(I);
                }
                // Copy 12 by 12 stiffness matrix, before static condensation
                for (unsigned int I = 0; I < 12; ++I)
                {
                    dv_FDISTR(I) = daFDistr[I][E];
                    for (unsigned int J = 0; J < 12; ++J)
                    {
                        dm_KBARE(I,J) = daKStorage[I][J][E];
                    }
                }
                // Transform element matrix & nodal load vector
                c1 = cos(phi1); s1 = sin(phi1); c2 = cos(phi2); s2 = sin(phi2);
                dc_KBARE_R1 = column(dm_KBARE, 1);
                dc_KBARE_R5 = column(dm_KBARE, 5);
                dc_KBARE_R7 = column(dm_KBARE, 7);
                dc_KBARE_R11 = column(dm_KBARE, 11);
                column(dm_KBARE, 1) = -s1 * dc_KBARE_R1 - c1 * dc_KBARE_R5;
                column(dm_KBARE, 5) = c1 * dc_KBARE_R1 - s1 * dc_KBARE_R5;
                column(dm_KBARE, 7) = -s2 * dc_KBARE_R7 - c2 * dc_KBARE_R11;
                column(dm_KBARE,11) = c2 * dc_KBARE_R7 - s2 * dc_KBARE_R11;
                dr_KBARE_R1 = row(dm_KBARE, 1);
                dr_KBARE_R5 = row(dm_KBARE, 5);
                dr_KBARE_R7 = row(dm_KBARE, 7);
                dr_KBARE_R11 = row(dm_KBARE, 11);
                row(dm_KBARE, 1) = -s1 * dr_KBARE_R1 - c1 * dr_KBARE_R5;
                row(dm_KBARE, 5) = c1 * dr_KBARE_R1 - s1 * dr_KBARE_R5;
                row(dm_KBARE, 7) = -s2 * dr_KBARE_R7 - c2 * dr_KBARE_R11;
                row(dm_KBARE,11) = c2 * dr_KBARE_R7 - s2 * dr_KBARE_R11;
                dV1 = dv_FDISTR(1) - dv_FINTER(1);
                dV5 = dv_FDISTR(5) - dv_FINTER(5);
                dV7 = dv_FDISTR(7) - dv_FINTER(7);
                dV11 = dv_FDISTR(11) - dv_FINTER(11);
                dv_Fel(0) = dv_FDISTR(0) - dv_FINTER(0);
                dv_Fel(1) = -s1 * dV1 - c1 * dV5;
                dv_Fel(4) = dv_FDISTR(4) - dv_FINTER(4);
                dv_Fel(5) = c1 * dV1 - s1 * dV5;
                dv_Fel(6) = dv_FDISTR(6) - dv_FINTER(6);
                dv_Fel(7) = -s2 * dV7 - c2 * dV11;
                dv_Fel(10) = dv_FDISTR(10) - dv_FINTER(10);
                dv_Fel(11) = c2 * dV7 - s2 * dV11;
                // Partition element matrix
                std::vector<unsigned int> vR_IDs{0,2,4,1,6,8,10,7}, vC_IDs{5,11};
                for (unsigned int I = 0; I < 8; ++I)
                {
                    for (unsigned int J = 0; J < 8; ++J)
                        dm_KRR(I,J) = dm_KBARE(vR_IDs[I],vR_IDs[J]);
                    for (unsigned int J = 0; J < 2; ++J)
                    {
                        dm_KRC(I,J) = dm_KBARE(vR_IDs[I],vC_IDs[J]);
                        if (I < 2) dm_KCC(I,J) = dm_KBARE(vC_IDs[I],vC_IDs[J]);
                    }
                    dv_FR(I) = dv_Fel(vR_IDs[I]);
                    if (I< 2) dv_FC(I) = dv_Fel(vC_IDs[I]);
                }
                bnu::matrix<double> dm_KCC_inv = bnu::identity_matrix<float>(dm_KCC.size1());
                bnu::permutation_matrix<size_t> pm(dm_KCC.size1());
                bnu::lu_factorize(dm_KCC, pm);
                bnu::lu_substitute(dm_KCC, pm, dm_KCC_inv);
                // Condense element matrix
                dm_Kpr1 = prod(dm_KRC, dm_KCC_inv);
                dm_K = dm_KRR - prod(dm_Kpr1, trans(dm_KRC));
                dv_F = dv_FR - prod(dm_Kpr1, trans(dv_FC));
                // Insert condensed matrix / vector coefficients into vector form containers
                el_off = (E + 1 - 1) * 4; el_off1 = (E + 1 - 1) * 3; el_off2 = (E + 1 - 1) * 2; el_off3 = (E + 1 - 1);
                for (unsigned int I = 0; I < 8; ++I)
                {
                    daMTV[I + el_off + umOffsets["D0"]] += dm_K(I, I); // Main diagonal
                    daVCV[I + el_off + umOffsets["D0"]] += dv_F(I); // Nodal load vector
                }
                for (unsigned int I = 0; I < 7; ++I)
                {
                    daMTV[I + el_off + umOffsets["U1"]] += dm_K(I, I + 1); // First upper off-diagonal
                    daMTV[I + el_off + umOffsets["L1"]] += dm_K(I + 1, I); // First lower off-diagonal
                }
                for (unsigned int I = 0; I < 6; ++I)
                {
                    daMTV[I + el_off + umOffsets["U2"]] += dm_K(I, I + 2); // Second upper off-diagonal
                    daMTV[I + el_off + umOffsets["L2"]] += dm_K(I + 2, I); // Second lower off-diagonal
                }
                for (unsigned int I = 0; I < 5; ++I)
                {
                    daMTV[I + el_off + umOffsets["U3"]] += dm_K(I, I + 3); // Third upper off-diagonal
                    daMTV[I + el_off + umOffsets["L3"]] += dm_K(I + 3, I); // Third lower off-diagonal
                }
                for (unsigned int I = 0; I < 4; ++I)
                {
                    daMTV[I + el_off + umOffsets["U4"]] += dm_K(I, I + 4); // Fourth upper off-diagonal
                    daMTV[I + el_off + umOffsets["L4"]] += dm_K(I + 4, I); // Fourth lower off-diagonal
                }
                for (unsigned int I = 0; I < 3; ++I)
                {
                    daMTV[I + el_off1 + umOffsets["U5"]] += dm_K(I, I + 5); // Fifth upper off-diagonal
                    daMTV[I + el_off1 + umOffsets["L5"]] += dm_K(I + 5, I); // Fifth lower off-diagonal
                }
                for (unsigned int I = 0; I < 2; ++I)
                {
                    daMTV[I + el_off2 + umOffsets["U6"]] += dm_K(I, I + 6); // Sixth upper off-diagonal
                    daMTV[I + el_off2 + umOffsets["L6"]] += dm_K(I + 6, I); // Sixth lower off-diagonal
                }
                daMTV[0 + el_off3 + umOffsets["U7"]] += dm_K(0, 0 + 7); // Seventh upper off-diagonal
                daMTV[0 + el_off3 + umOffsets["L7"]] += dm_K(0 + 7, 0); // Seventh lower off-diagonal
            } // End of parallel loop
        } // End of parallel region

        outputs[0] = daVCV;
        outputs[1] = daMTV;
        outputs[2] = daFINTER;
        outputs[3] = daAlphas;
        outputs[4] = daSIGMAS;
        outputs[5] = daSigmas;
        outputs[6] = daSigma_ys;
        outputs[7] = daEpsilons;
        outputs[8] = daDepns;
        outputs[9] = daSubFlags;
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
