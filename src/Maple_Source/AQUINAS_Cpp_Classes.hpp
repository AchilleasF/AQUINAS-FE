// AQUINAS is distributed under a BSD 3-Clause License, included at the end of this file.
#pragma once
namespace bnu = boost::numeric::ublas;


/////////////////////////////////
// Abstract Segment Base Class //
/////////////////////////////////
class cSegment
{
public:
    double mer_curv_sign; // Sign of meridional curvature of segment with respect to the axis of revolution. Only relevant for elliptical segments, zero otherwise.
    double rbot; // Global cylindrical radial reference surface r coordinate of base node (LHS node for ring)
    double zbot; // Global cylindrical vertical reference surface z coordinate of base node (LHS node for ring)
    double rtop; // Global cylindrical radial reference surface r coordinate of top node (RHS node for ring)
    double ztop; // Global cylindrical vertical reference surface z coordinate of top node (RHS node for ring)

public:
    void set_bot_top(matlab::data::TypedArray<double> daBotTop)
    {
        this->rbot = daBotTop[0];
        this->zbot = daBotTop[1];
        this->rtop = daBotTop[2];
        this->ztop = daBotTop[3];
    }

public:
    virtual double r(double z) = 0; // Global r coordinate of reference meridian i.t.o. global z coordinate
    virtual double drdz(double z) = 0; // dr/dz of reference meridian i.t.o. global z coordinate
    virtual double d2rdz2(double z) = 0; // d2r/dz2 of reference meridian i.t.o. global z coordinate
    virtual double d3rdz3(double z) = 0; // d3r/dz3 of reference meridian i.t.o. global z coordinate

public:
    virtual double S(double B, double T) // Arc length coordinate of refernece meridian from lower to upper global z coordinate positions
    {
        // For cCone and cEllipse, B and T are zB and zT
        auto arc = [this](const double &z) -> double { return sqrt(1.0 + this->drdz(z) * this->drdz(z)); };
        double error;
        return boost::math::quadrature::gauss_kronrod<double, 61>::integrate(arc, B, T, 20, 1e-9, &error);
    }

    virtual double PHI(double z) // Meridional angle of normal to reference meridian i.t.o. r(z) i.t.o. global z coordinate
    {
        return M_PI_2 - atan2(-this->drdz(z), 1.0);
    }

    double dRdS(double z) // dr/ds i.t.o. global z coordinate
    {
        if (isinf(this->drdz(z)))
            if (this->rbot > this->rtop)
                return 1.0;
            else
                return -1.0;
        else
            return (-this->drdz(z) / sqrt(1.0 + this->drdz(z) * this->drdz(z)));
    }

    double dZdS(double z) // dz/ds i.t.o. global z coordinate
    {
        if (isinf(this->drdz(z)))
            return 0.0;
        else
            return (-1.0 / sqrt(1.0 + this->drdz(z) * this->drdz(z)));
    }

    double dPHIdS(double z) // dphi/ds curvature of reference meridian i.t.o. r(z) i.t.o. global z coordinate
    {
        if (isinf(this->drdz(z)))
            return 0.0;
        else
            return (this->dZdS(z) * this->d2rdz2(z) / (1.0 + this->drdz(z) * this->drdz(z)));
    }

    double d2PHIdS2(double z) // d2phi/ds2 aberrancy of reference meridian i.t.o. r(z) i.t.o. global z coordinate
    {
        if (isinf(this->drdz(z)))
            return 0.0;
        else
            return -this->dZdS(z) * (this->d3rdz3(z) / pow(1.0 + this->drdz(z) * this->drdz(z), 1.5) - 3.0*(this->drdz(z) * this->d2rdz2(z) * this->d2rdz2(z) / pow(1.0 + this->drdz(z) * this->drdz(z), 2.5)));
    }
};

/////////////////////////////////
// Plate Segment Derived Class //
/////////////////////////////////
class cPlate : public cSegment
{
public:
    cPlate(matlab::data::TypedArray<double> daBotTop)
    {
        this->set_bot_top(daBotTop);
        this->mer_curv_sign = 0;
    }

    ~cPlate() {}

private:
    double S(double B, double T) override
    {
        // For cPlate, B and T are rB and rT
        return abs(T-B);
    }

    double PHI(double z) override
    {
        return 0;
    }

    double r(double z) override
    {
        return ((this->rbot - this->rtop)*z - (this->ztop * this->rbot - this->zbot * this->rtop)) / (this->zbot - this->ztop);
    }

    double drdz(double z) override
    {
        return INFINITY;
    }

    double d2rdz2(double z) override
    {
        return 0.0;
    }

    double d3rdz3(double z) override
    {
        return 0.0;
    }
};

///////////////////////////////////
// Conical Segment Derived Class //
///////////////////////////////////
class cCone : public cSegment
{
public:
    cCone(matlab::data::TypedArray<double> daBotTop)
    {
        this->set_bot_top(daBotTop);
        this->mer_curv_sign = 0;
    }

    ~cCone() {}

private:
    double r(double z) override
    {
        return ((this->rbot - this->rtop)*z - (this->ztop * this->rbot - this->zbot * this->rtop)) / (this->zbot - this->ztop);
    }

    double drdz(double z) override
    {
        return (this->rbot - this->rtop) / (this->zbot - this->ztop);
    }

    double d2rdz2(double z) override
    {
        return 0.0;
    }

    double d3rdz3(double z) override
    {
        return 0.0;
    }
};

//////////////////////////////////////
// Elliptical Segment Derived Class //
//////////////////////////////////////
class cEllipse : public cSegment
{
public:
    double cenr;
    double cenz;
    double rhor;
    double rhoz;

public:
    cEllipse(matlab::data::TypedArray<double> daBotTop, matlab::data::TypedArray<double> daGeom, double merCurvSign)
    {
        this->set_bot_top(daBotTop);
        this->cenr = daGeom[0];
        this->cenz = daGeom[1];
        this->rhor = daGeom[2];
        this->rhoz = daGeom[3];
        this->mer_curv_sign = merCurvSign;
    }

    ~cEllipse() {}

private:
    double r(double z) override
    {
        return this->cenr + this->mer_curv_sign * (this->rhor / this->rhoz) * pow((this->rhoz + this->cenz - z) * (this->rhoz - this->cenz + z), 0.5);
    }

    double drdz(double z) override
    {
        return this->mer_curv_sign * (this->rhor / this->rhoz) * (this->cenz - z) * pow((this->rhoz + this->cenz - z) * (this->rhoz - this->cenz + z), -0.5);
    }

    double d2rdz2(double z) override
    {
        return -this->mer_curv_sign * this->rhor * this->rhoz * pow((this->rhoz + this->cenz - z) * (this->rhoz - this->cenz + z), -1.5);
    }

    double d3rdz3(double z) override
    {
        return 3.0 * this->mer_curv_sign * this->rhor * this->rhoz * (this->cenz - z) * pow((this->rhoz + this->cenz - z) * (this->rhoz - this->cenz + z), -2.5);
    }
};

////////////////////////////////
// Distributed Pressure Class //
////////////////////////////////
class cDistributedPressure
{
public :
    int type; // Distribute Pressure type (1: normal / 2: tangent)

public :
    cDistributedPressure(int type)
    {
        this->type = type;
    };

    ~cDistributedPressure() {}

    double eta_integral_equivalent_nodal_load_vector(unsigned int I, double eta, double p1, double p2, double N01, double N11, double N02, double N12,double L, double c, double s)
    {
        // Returns the vector of equivalent nodal loads due to Distributed Pressure along the meridian of the element, term by term
        // This vector is evaluated through numerical integration of the Distributed Pressure multiplied by the the shape functions
        double out;
        double p = p1 * (1 - eta) / 2 + p2 * (1 + eta) / 2;
        switch (this->type)
        {
        case 1: // Normal distributed load (pn)
            // F = Integral of : magnitude * r * [N01*s N11*s 0 0 N01*c N11*c N02*s N12*s 0 0 N02*c N12*c] for eta = -1 -> 1
            switch (I)
            {
            case 0:
                out = p * N01 * s;
                break;
            case 1:
                out = p * N11 * s;
                break;
            case 4:
                out = p * N01 * c;
                break;
            case 5:
                out = p * N11 * c;
                break;
            case 6:
                out = p * N02 * s;
                break;
            case 7:
                out = p * N12 * s;
                break;
            case 10:
                out = p * N02 * c;
                break;
            case 11:
                out = p * N12 * c;
                break;
            default:
                out = 0.0;
            }
            break;
        case 2: // Tangential distributed load (pt)
            // F = Integral of : magnitude * r * [N01*c N11*c 0 0 N01*s N11*s N02*c N12*c 0 0 N02*s N12*s] for eta = -1 -> 1
            switch (I)
            {
            case 0:
                out = p * N01 * c;
                break;
            case 1:
                out = p * N11 * c;
                break;
            case 4:
                out = - p * N01 * s;
                break;
            case 5:
                out = - p * N11 * s;
                break;
            case 6:
                out = p * N02 * c;
                break;
            case 7:
                out = p * N12 * c;
                break;
            case 10:
                out = - p * N02 * s;
                break;
            case 11:
                out = - p * N12 * s;
                break;
            default:
                out = 0.0;
            }
            break;
        default:
            out = 0.0;
            break;
        }
        return out;
    }

};

////////////////////
// Material Class //
////////////////////
class cMaterial
{

public:
    int type; // type of material, 0 for isotropic or 1 for orthotropic (orthotropic material not currently fully supported)
    double E; // Isotropic Young's modulus
    double nu; // Isotropic Poisson's ratio
    double G; // Shear modulus
    double yield_tol;
    double epsilon_s;
    double max_epsilon_bar;
    bnu::matrix<double> De;
    std::vector<double> sy;
    std::vector<double> ep;

public:
    cMaterial(int type, double E, double nu, double G, int noMatCurvePoints, matlab::data::TypedArray<double>& sigy, matlab::data::TypedArray<double>& epsp, double epsilon_s, double max_epsilon_bar)
    {
        this->type = type;
        this->E = E;
        this->nu = nu;
        this->G = G;
        this->yield_tol = 1e-8;
        this->epsilon_s = epsilon_s;
        this->max_epsilon_bar = max_epsilon_bar;
        this->De = bnu::zero_matrix<double>(3, 3);
        De(0, 0) = this->E / (1 - this->nu * this->nu);
        De(0, 1) = this->nu * this->E / (1 - this->nu * this->nu);
        De(1, 1) = this->E / (1 - this->nu * this->nu);
        De(1, 0) = this->nu * this->E / (1 - this->nu * this->nu);
        De(2, 2) = this->G;
        for (int i = 0; i < noMatCurvePoints; i++)
        {
            sy.push_back(sigy[i]);
            ep.push_back(epsp[i]);
        };
    };

    // Auxiliary method for computing the Z matrix, accoring to eq. 12 of Teng and Rotter (1989a)
    bnu::matrix<double> Z(double z)
    {
        bnu::matrix<double> Z = bnu::zero_matrix<double>(3,6);
        for (int i = 0; i < 3; ++i)
        {
            Z(i,i) = 1.0; Z(i,i+3) = z;
        }
        return Z;
    }

    // Compute strain hardening parameter H, according to eq. 27 of Teng and Rotter(1989a)
    double H(double epn)
    {
        if ( (this->ep.empty()) || (this->ep.size() == 1) ) { return 0.0; }
        // Find the range that of the stres-strain curve that the current equivalnent plastic strain belongs to
        int I;
        for (int i = 0; i < this->ep.size(); ++i)
        {
            if ( i == (this->ep.size()-1)) { I = i; break; }
            if ((this->ep[i] < epn) & (this->ep[i+1] > epn)) { I = i; break; }
        }
        if (I == (this->ep.size()-1)) { I = I - 1; }
        double ep_p = this->ep[I]; double ep_n = this->ep[I+1];
        double sy_p = this->sy[I]; double sy_n = this->sy[I+1];
        // Compute isotropic strain hardening parameter
        double Hi = (sy_n - sy_p)/(ep_n - ep_p);
        return Hi;
    }

    // Compute yield stress of material point, according to par. 10-(8)-(v)-(b) of Teng and Rotter (1989a)
    double sigmay(double sigy,double epn,double depn)
    {
        if ((epn + depn) < 1e-16) { return sigy; } // if there is no equivalnet plastic strain developing, do not alter the value of the yield stress
        return (sigy + this->H(epn)*depn);
    }

    // Compute(1 - ksi) portion of plastic strain out of the current strain increment, according to eq. 85 of Teng and Rotter(1989a)
    double elastic_strain_portion(bnu::vector<double> sigma, bnu::vector<double> Deltasigma, double sigmay)
    {
        // E - element counter | I - Gauss point | J - Simpson point
        double ksi;
        double A = Deltasigma(0) * Deltasigma(0) + Deltasigma(1) * Deltasigma(1) - Deltasigma(0) * Deltasigma(1) + 3 * Deltasigma(2) * Deltasigma(2);
        double B = Deltasigma(0) * (2 * sigma(0) - sigma(1)) + Deltasigma(1) * (2 * sigma(1) - sigma(0)) + 6 * sigma(2) * Deltasigma(2);
        double C = sigma(0) * sigma(0) + sigma(1) * sigma(1) - sigma(0) * sigma(1) + 3 * sigma(2) * sigma(2) - sigmay * sigmay;
        double D = B*B - 4 * A * C;
        if (D < 0) {return 0;}
        ksi = (-B + sqrt(B*B - 4 * A * C)) / 2 / A;
        if (ksi < 0) { ksi = 0; }
        return ksi;
    };

    // Compute Dep matrix at a through - thickness station of coordinate z, according to eq. 34 of Teng and Rotter(1989a)
    // Method overloaded for matlab::data::TypedArray<double> sigmas
    bnu::matrix<double> Dep(matlab::data::TypedArray<double>& sigmas, double sigmay, matlab::data::TypedArray<double>& epns, int E, int I, int J)
    {
        // E - element counter | I - Gauss point | J - Simpson point
        // von Mises equivalent stress, according to eq. 24 of Teng and Rotter(1989a)
        double sigmabar = sqrt(sigmas[0][E][I][J] * sigmas[0][E][I][J] + sigmas[1][E][I][J] * sigmas[1][E][I][J] - sigmas[0][E][I][J] * sigmas[1][E][I][J] + 3 * sigmas[2][E][I][J] * sigmas[2][E][I][J]);
        if (sigmabar > sigmay * (1 - this->yield_tol))
        {
            // Deviatoric sigmas, according to eqs. 30 of Teng and Rotter(1989a)
            double sphi = (2 * sigmas[0][E][I][J] - sigmas[1][E][I][J]) / 3; double stheta = (2 * sigmas[1][E][I][J] - sigmas[0][E][I][J]) / 3;
            // S1, S2, S4and S5 parameters, according to eqs. 35 of Teng and Rotter(1989a)
            double S1 = sphi + this->nu * stheta; double S2 = stheta + this->nu * sphi; double S4 = sphi * sphi + stheta * stheta + 2 * this->nu * sphi * stheta + 2 * (1 - this->nu) * sigmas[2][E][I][J] * sigmas[2][E][I][J]; double S5 = 2 * (sigmabar * sigmabar) * (1 - this->nu) * this->H(epns[E][I][J]) / (9 * this->G) + S4;
            bnu::matrix<double> Dp(3,3);
            Dp(0, 0) = this->E / ((1 - this->nu * this->nu) * S5) * (S1 * S1);
            Dp(0, 1) = this->E / ((1 - this->nu * this->nu) * S5) * (S1 * S2);
            Dp(0, 2) = this->E / ((1 - this->nu * this->nu) * S5) * (S1 * (1 - this->nu) * sigmas[2][E][I][J]);
            Dp(1, 0) = this->E / ((1 - this->nu * this->nu) * S5) * (S1 * S2);
            Dp(1, 1) = this->E / ((1 - this->nu * this->nu) * S5) * (S2 * S2);
            Dp(1, 2) = this->E / ((1 - this->nu * this->nu) * S5) * (S2 * (1 - this->nu) * sigmas[2][E][I][J]);
            Dp(2, 0) = this->E / ((1 - this->nu * this->nu) * S5) * (S1 * (1 - this->nu) * sigmas[2][E][I][J]);
            Dp(2, 1) = this->E / ((1 - this->nu * this->nu) * S5) * (S2 * (1 - this->nu) * sigmas[2][E][I][J]);
            Dp(2, 2) = this->E / ((1 - this->nu * this->nu) * S5) * (sigmas[2][E][I][J] * sigmas[2][E][I][J]) * ((1 - this->nu) * (1 - this->nu));
            return this->De - Dp;
        }
        else
        {
            return this->De;
        }
    };

    // Compute Dep matrix at a through - thickness station of coordinate z, according to eq. 34 of Teng and Rotter(1989a)
    // Method overloaded for boost::numeric::ublas::vector<double> sigma
    bnu::matrix<double> Dep(bnu::vector<double> sigma, double sigmay, double epn)
    {
        // E - element counter | I - Gauss point | J - Simpson point
        // von Mises equivalent stress, according to eq. 24 of Teng and Rotter(1989a)
        double sigmabar = sqrt(sigma(0) * sigma(0) + sigma(1) * sigma(1) - sigma(0) * sigma(1) + 3 * sigma(2) * sigma(2));
        if (sigmabar > sigmay * (1 - this->yield_tol))
        {
            // Deviatoric sigmas, according to eqs. 30 of Teng and Rotter(1989a)
            double sphi = (2 * sigma(0) - sigma(1)) / 3; double stheta = (2 * sigma(1) - sigma(0)) / 3;
            // S1, S2, S4and S5 parameters, according to eqs. 35 of Teng and Rotter(1989a)
            double S1 = sphi + this->nu * stheta; double S2 = stheta + this->nu * sphi; double S4 = sphi * sphi + stheta * stheta + 2 * this->nu * sphi * stheta + 2 * (1 - this->nu) * sigma(2) * sigma(2); double S5 = 2 * (sigmabar * sigmabar) * (1 - this->nu) * this->H(epn) / (9 * this->G) + S4;
            bnu::matrix<double> Dp(3,3);
            Dp(0, 0) = this->E / ((1 - this->nu * this->nu) * S5) * (S1 * S1);
            Dp(0, 1) = this->E / ((1 - this->nu * this->nu) * S5) * (S1 * S2);
            Dp(0, 2) = this->E / ((1 - this->nu * this->nu) * S5) * (S1 * (1 - this->nu) * sigma(2));
            Dp(1, 0) = this->E / ((1 - this->nu * this->nu) * S5) * (S1 * S2);
            Dp(1, 1) = this->E / ((1 - this->nu * this->nu) * S5) * (S2 * S2);
            Dp(1, 2) = this->E / ((1 - this->nu * this->nu) * S5) * (S2 * (1 - this->nu) * sigma(2));
            Dp(2, 0) = this->E / ((1 - this->nu * this->nu) * S5) * (S1 * (1 - this->nu) * sigma(2));
            Dp(2, 1) = this->E / ((1 - this->nu * this->nu) * S5) * (S2 * (1 - this->nu) * sigma(2));
            Dp(2, 2) = this->E / ((1 - this->nu * this->nu) * S5) * (sigma(2) * sigma(2)) * ((1 - this->nu) * (1 - this->nu));
            return this->De - Dp;
        }
        else
        {
            return this->De;
        }
    };

    // Compute elastic modulus matrix DT, according to eq. 40 of Teng and Rotter(1989a), without the inclusion of the plastic part Dp
    bnu::matrix<double> DTe(double t)
    {
        double zmin = -t / 2; double zmax = t / 2;
        double zm = zmax - zmin; double zm2 = (zmax * zmax - zmin * zmin) * 0.5; double zm3 = (zmax * zmax * zmax - zmin * zmin * zmin) / 3;
        bnu::matrix<double> DTe(6, 6);
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                DTe(i, j) = zm * this->De(i,j);
                DTe(i+3, j) = zm2 * this->De(i, j);
                DTe(i, j+3) = zm2 * this->De(i, j);
                DTe(i+3, j+3) = zm3 * this->De(i, j);
            }
        }
        return DTe;
    }

    // Compute tangent modulus matrix DT according to eq. 40 of Teng and Rotter(1989a)
    bnu::matrix<double> DTep(double t, matlab::data::TypedArray<double>& sigmas, matlab::data::TypedArray<double>& sigma_ys, matlab::data::TypedArray<double>& epns, int E, int I, int noSimpson)
    {
        // E - element counter | I - Gauss point | J - Simpson point
        double zmin = -t / 2; double zmax = t / 2;
        int noLayers = noSimpson - 1;
        bnu::matrix<double> DTelastoplastic = bnu::zero_matrix<double>(6,6);
        for (int J = 0; J < (int)(noLayers / 2); ++J)
        {
            int a = 2 * J; int m = 2 * J + 1; int b = 2 * J + 2;
            double z_a = -t / 2 + (double)a*t/(double)noLayers; double z_m = -t / 2 + (double)m * t / (double)noLayers; double z_b = -t / 2 + (double)b * t / (double)noLayers;
            bnu::matrix<double> Z_a = this->Z(z_a);
            bnu::matrix<double> Z_m = this->Z(z_m);
            bnu::matrix<double> Z_b = this->Z(z_b);
            bnu::matrix<double> Dep_a = this->Dep(sigmas,sigma_ys[E][I][a],epns,E,I,a);
            bnu::matrix<double> Dep_m = this->Dep(sigmas,sigma_ys[E][I][m],epns,E,I,m);
            bnu::matrix<double> Dep_b = this->Dep(sigmas,sigma_ys[E][I][b],epns,E,I,b);
            bnu::matrix<double> DepZ_a = prod(Dep_a, Z_a);
            bnu::matrix<double> DepZ_m = prod(Dep_m, Z_m);
            bnu::matrix<double> DepZ_b = prod(Dep_b, Z_b);
            bnu::matrix<double> DTep_a = prod(trans(Z_a), DepZ_a);
            bnu::matrix<double> DTep_m = prod(trans(Z_m), DepZ_m);
            bnu::matrix<double> DTep_b = prod(trans(Z_b), DepZ_b);
            for (int i = 0; i < 6; ++i)
            {
                for (int j = 0; j < 6; ++j)
                {
                    DTelastoplastic(i,j) += (z_b - z_a) * (DTep_a(i, j) + 4 * DTep_m(i, j) + DTep_b(i, j)) / 6;
                }
            }
        }
        return DTelastoplastic;
    };

    // Compute plastic strains - sigmas using the sub incremental technique
    void sub_incremental_computations(matlab::data::TypedArray<int>& daSubFlags, bnu::vector<double> epsilon0, bnu::vector<double> DeltaepsilonZ, matlab::data::TypedArray<double>& sigmas, matlab::data::TypedArray<double>& sigma_ys, matlab::data::TypedArray<double>& epns, matlab::data::TypedArray<double>& Depns, matlab::data::TypedArray<unsigned int>& alphas, int E, int I, int J)
    {
        bnu::vector<double> sigman(3), deltasigma;
        sigman(0) = sigmas[0][E][I][J]; sigman(1) = sigmas[1][E][I][J]; sigman(2) = sigmas[2][E][I][J];
        double ksi = 0; double sigmabar = 0;
        double Depn = Depns[E][I][J]; double sigmayn = sigma_ys[E][I][J];
        if ((alphas[E][I][J] == 0) || (Depn < 0))
        {
            bnu::vector<double> deltasigma = prod(this->De, DeltaepsilonZ);
            sigmabar = sqrt((sigman(0) + deltasigma(0)) * (sigman(0) + deltasigma(0)) + (sigman(1) + deltasigma(1)) * (sigman(1) + deltasigma(1)) - (sigman(0) + deltasigma(0)) * (sigman(1) + deltasigma(1)) + 3 * (sigman(2) + deltasigma(2)) * (sigman(2) + deltasigma(2))); // von Mises stress
            if (sigmabar > sigma_ys[E][I][J] * (1 - this->yield_tol)) // The point has yielded
            {
                // (1 - ksi) is the portion of the current strain increment that causes plasticity to occur for the material
                ksi = this->elastic_strain_portion(sigman, deltasigma, sigma_ys[E][I][J]);
                alphas[E][I][J] = 1; // This point is now plastic hardening
            }
            else
            {
                sigman = sigman + deltasigma;
                sigmas[0][E][I][J] = sigman(0); sigmas[1][E][I][J] = sigman(1); sigmas[2][E][I][J] = sigman(2);
                return;
            }
        }
        sigman = sigman + ksi * prod(this->De, DeltaepsilonZ);
        bnu::vector<double> epsilonn = epsilon0 + ksi * DeltaepsilonZ;
        // Effective strain increment, according to par. 10 - (8) - (iv)of Teng and Rotter(1989a)
        double Delta_epsilon_bar = (2 * (1 - ksi) / sqrt(3)) * sqrt(DeltaepsilonZ(0) * DeltaepsilonZ(0) + DeltaepsilonZ(1) * DeltaepsilonZ(1) + DeltaepsilonZ(0) * DeltaepsilonZ(1) + (DeltaepsilonZ(2) * DeltaepsilonZ(2)) / 4);
        // An effective strain increment that is higher than the max_epsilon_bar will not be attempted
        if (Delta_epsilon_bar > this->max_epsilon_bar) { daSubFlags[E][I][J] = 1; return; }
        // Number of required sub increments, according to par. 10 - (8) - (iv)of Teng and Rotter(1989a)
        int Nsb = (int)(Delta_epsilon_bar / this->epsilon_s) + 1;
        // Strain increment to be applied per sub - increment, according to par. 10 - (8) - (iv)of Teng and Rotter(1989a)
        bnu::vector<double> depsilon_Z = ((1 - ksi) / (double)Nsb) * DeltaepsilonZ;
        // Initializations to avoid constant re-assignment of memory space
        bnu::matrix<double> Dep(3,3), Dep_halfstep(3,3), Dep_fullstep(3,3);
        bnu::vector<double> dsigma_n(3), deenZ(3), depnZ(3), deviatoric_stresses(3), sigma_n_halfstep(3);
        for (int n = 0; n < Nsb; n++)
        {
            if (alphas[E][I][J] == 0) // The point is behaving elastically(unloaded during the n - 1 step)
            {
                sigmas[0][E][I][J] = sigman(0); sigmas[1][E][I][J] = sigman(1); sigmas[2][E][I][J] = sigman(2); Depns[E][I][J] = Depn;
                this->sub_incremental_computations(daSubFlags, epsilonn, depsilon_Z, sigmas, sigma_ys, epns, Depns, alphas, E, I, J); // recursive call to the function it self, in order to return to step(ii) of methodology of paragraph 10 - (8) of Teng and Rotter 1989a
                if (daSubFlags[E][I][J] == 1) { return; }
                epsilonn += depsilon_Z;
            }
            else
            {
                // von Mises stress, according to eq. 24 of Teng and Rotter(1989a)
                sigmabar = sqrt(sigman(0) * sigman(0) + sigman(1) * sigman(1) - sigman(0) * sigman(1) + 3 * sigman(2) * sigman(2));
                // Deviatoric sigmas, according to eqs. 30 of Teng and Rotter(1989a)
                double sphi = (2 * sigman(0) - sigman(1)) / 3; double stheta = (2 * sigman(1) - sigman(0)) / 3; double sphitheta = 2 * sigman(2);
                // S1, S2, S3, S4and S5 parameters, according to eqs. 35 of Teng and Rotter(1989a)
                double S1 = sphi + this->nu * stheta; double S2 = stheta + this->nu * sphi; double S4 = sphi * sphi + stheta * stheta + 2 * this->nu * sphi * stheta + 2 * (1 - this->nu) * sigman(2) * sigman(2); double S5 = 2 * (sigmabar * sigmabar) * (1 - this->nu) * this->H(epns[E][I][J]) / (9 * this->G) + S4;
                double S3 = S1 * depsilon_Z[0] + S2 * depsilon_Z[1] + (1 - this->nu) * sigman(2) * depsilon_Z[2];
                // Sub - incremental equivalent plastic strain, according to eq. 36 of of Teng and Rotter(1989a)
                double depn = (2 * sigmabar / 3) * S3 / S5;
                // Accumulated sub - incremental equivalent plastic strain, according to par. 10 - (8) - (v)-(b)of Teng and Rotter(1989a)
                Depn += depn;

                if (Depn >= 0)
                {
                    // Update yield stress, according to par. 10 - (8) - (v)-(b)of Teng and Rotter(1989a)
                    sigmayn = this->sigmay(sigma_ys[E][I][J], epns[E][I][J], Depns[E][I][J]);
                    // Stresses at half the step, according to eq. 72 of Teng and Rotter(1989a)
                    Dep_halfstep = this->Dep(sigman, sigmayn, epns[E][I][J]);
                    sigma_n_halfstep = sigman + prod(Dep_halfstep, depsilon_Z);
                    // Current stress increment, according to eq. 73 of Teng and Rotter(1989a)
                    Dep_fullstep = this->Dep(sigma_n_halfstep, sigmayn, epns[E][I][J]);
                    Dep = 0.5 * (Dep_halfstep + Dep_fullstep);
                    dsigma_n = prod(Dep, depsilon_Z);
                }
                else if (Depn < 0)
                {
                    alphas[E][I][J] = 0;
                    deviatoric_stresses(0) = sphi; deviatoric_stresses(1) = stheta; deviatoric_stresses(2) = sphitheta;
                    depnZ = Depn * (3 / 2 / sigmabar) * deviatoric_stresses;
                    deenZ = (depsilon_Z - depnZ);
                    dsigma_n = prod(this->De, deenZ);
                }
                sigman += dsigma_n;
                epsilonn += depsilon_Z;
            }
        }
        // von Mises stress, according to eq. 24 of Teng and Rotter(1989a)
        sigmabar = sqrt(sigman(0) * sigman(0) + sigman(1) * sigman(1) - sigman(0) * sigman(1) + 3 * sigman(2) * sigman(2));
        // Yield surface function value
        if (Depn < 0)
        {
            sigmayn = this->sigmay(sigma_ys[E][I][J], epns[E][I][J], 0);
        }
        else
        {
            sigmayn = this->sigmay(sigma_ys[E][I][J], epns[E][I][J], Depns[E][I][J]);
        }
        double F1 = sigmabar - sigmayn; int iters = 0;
        if (F1 > -this->yield_tol * sigmayn) { alphas[E][I][J] = 1; }
        while (F1 > this->yield_tol * sigmayn)
        {
            // Deviatoric sigmas, according to eqs. 30 of Teng and Rotter(1989a)
            double sphi = (2 * sigman(0) - sigman(1)) / 3; double stheta = (2 * sigman(1) - sigman(0)) / 3; double sphitheta = 2 * sigman(2);
            // Compute stress corrections, according to eq. 79 of Teng and Rotter(1989a)
            deviatoric_stresses(0) = sphi; deviatoric_stresses(1) = stheta; deviatoric_stresses(2) = sphitheta;
            deltasigma = -2 * sigmabar * F1 / (3 * (sphi * sphi + stheta * stheta + 4 * sigman(2) * sigman(2))) * deviatoric_stresses;
            // Update sigmas of current sub increment n with the stress corrections
            sigman += deltasigma;
            // von Mises stress, according to eq. 24 of Teng and Rotter(1989a)
            sigmabar = sqrt(sigman(0) * sigman(0) + sigman(1) * sigman(1) - sigman(0) * sigman(1) + 3 * sigman(2) * sigman(2));
            // Yield surface function value
            F1 = sigmabar - sigmayn;
            iters = iters + 1;
            if (iters > 1000) { daSubFlags[E][I][J] = 1; return; }
        }
        sigmas[0][E][I][J] = sigman(0); sigmas[1][E][I][J] = sigman(1); sigmas[2][E][I][J] = sigman(2); Depns[E][I][J] = Depn;
        return;
    };

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