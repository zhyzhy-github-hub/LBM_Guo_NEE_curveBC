#include "CDE2019.h"
#include "NSE.h"

using std::cout;
using std::endl;
CDE2019::CDE2019(int NX_, int NY_)
{

    NX = NX_;
    NY = NY_;
    const unsigned int NN = NX * NY;
    C = new double[NN];
    C_old = new double[NN];
    F = new double[NN];

    f1 = new double[NN * Q9];
    f0 = new double[NN * Q9];

    s = new double[Q9];

    B_ade = new BADE[Q9];
}

CDE2019::~CDE2019()
{
    delete[] C;
    delete[] C_old;
    delete[] F;

    delete[] f0;
    delete[] f1;

    delete[] s;
    delete[] B_ade;
}

size_t CDE2019::scalar_index(size_t i, size_t j) const
{
    return (NX * j + i);
}
size_t CDE2019::field_index(size_t x, size_t y, size_t d) const
{
    // return (NX * (NY * (d) + y) + x);
    return (Q9 * (NX * y + x) + d);
}

double CDE2019::get_C(const size_t i, const size_t j) const
{
    return C[scalar_index(i, j)];
}

CDE2019::CDE2019(const int NX_, const int NY_, const double mu_, const double D_, const double c_, const double s_add)
    : NX(NX_), NY(NY_), mu(mu_), D(D_), c(c_)
{
    const unsigned int NN = NX * NY;
    C = new double[NN];
    C_old = new double[NN];
    F = new double[NN];

    f1 = new double[NN * Q9];
    f0 = new double[NN * Q9];

    s = new double[Q9];
    B_ade = new BADE[NN];

    //-----------------------------------

    dx = 1.0 / (NY - 1);
    dt = mu * dx * dx;

    Cs = c / sqrt(3);
    Cs2 = pow(Cs, 2);

    double s3 = 2.0 / (6 * mu * D + 1.0);
    tau = 1.0 / s3;

    for (int i = 0; i < Q9; ++i)
    {
        if (i == 3 || i == 5)
        {
            s[i] = s3;
        }
        else
        {
            // s[i] = 1.0 + s_add;
            s[i] = s3; //1.0 + s_add;
        }
    }
    cout << "--------------ADE paramater----------------" << endl;
    cout << "dx = " << dx << endl;
    cout << "dt = " << dt << endl;
    cout << "mu = " << mu << endl;
    cout << "tau = " << tau << endl;
    cout << " c = " << c << endl;
    for (int i = 0; i < Q9; ++i)
    {
        cout << "s[" << i << "] " << s[i] << " | ";
    }
    cout << endl;
}

#ifdef DIFF_SCALING
void CDE2019::initCDE2019_Diff(double mu_, double D_, double s_add)
{
    mu = mu_;
    dx = 1.0 / (NY - 1);
    dt = mu * dx * dx;
    // dt = mu * dx ;
    c = dx / dt;
    Cs = c / sqrt(3);

    double s3 = 2.0 / (6 * mu * D_ + 1);
    for (int i = 0; i < Q9; ++i)
    {
        if (i == 3 || i == 5)
        {
            s[i] = s3;
        }
        else
        {
            s[i] = s3; //1.0 + s_add;
        }
    }

    using std::cout;
    using std::endl;
    cout << "dx = " << dx << endl;
    cout << "dt = " << dt << endl;
    cout << "mu = " << mu << endl;
    cout << " c = " << c << endl;
    for (int i = 0; i < Q9; ++i)
    {
        cout << "s[" << i << "] " << s[i] << " | ";
    }
    cout << endl;
}
#elif
void CDE2019::initCDE2019(double c_, double D_, double s_add)
{
    c = c_;
    dt = dx / c;
    Cs = c / sqrt(3);
    Cs2 = c * c / 3;

    D = D_;
    tau = 0.5 + D / (Cs2 * dt);
    double tau_inv = 1.0 / tau;

    for (int i = 0; i < Q9; ++i)
    {
        if (i == 3 || i == 5)
        {
            s[i] = tau_inv;
        }
        else
        {
            s[i] = 1.0 + s_add;
        }
    }
}
#endif

void CDE2019::initScalar(size_t n, double *x, double *y)
{
    double t = dt * n;
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t index = scalar_index(i, j);
            double y_ = y[index];
            double x_ = x[index];

            // C[index] = C0 * exp(-(p2(x_ - 0.5) + p2(y_ - 0.5)) / (2 * p2(delta0)));
            F[index] = 0;
            C[index] = (t + 1) * sin(2 * PI * x_) * cos(2 * PI * y_);
        }
    }
}

void CDE2019::getCOld()
{
    size_t NN = NX * NY;
    for (int i = 0; i < NN; ++i)
    {
        C_old[i] = C[i];
    }
}

void CDE2019::getFi(double *Fi, double F)
{
    for (int i = 0; i < Q9; ++i)
    {
        if (i == 0)
        {
            Fi[i] = dt * w0 * F;
        }
        // else if (i < 5 && i > 0)
        else if (i == 1 || i == 2 || i == 3 || i == 4)
        {
            Fi[i] = dt * wc * F;
        }
        else
        {
            Fi[i] = dt * wd * F;
        }
    }
}

double CDE2019::meqADE2019(int k, double phi, double D, BADE B)
{
    switch (k)
    {
    case 0:
        return phi;
    case 1:
        return 2 * D - 4 * phi;
    case 2:
        return 3 * phi - 2 * D;
    case 3:
        return mu * dx * B.Bx;
    case 4:
        return -mu * dx * B.Bx;
    case 5:
        return mu * dx * B.By;
    case 6:
        return -mu * dx * B.By;
    case 7:
        return 0;
    case 8:
        return 0;
    }
}

// double CDE2019::meqADE2019(int k, double phi, double D, BADE B)
// {
//     double feq[9];
//     double meq[9];

//     feq[0] = w0 * (2 * phi - D);
//     feq[1] = wc * (2 * phi - D + 3 * B.Bx / c + (3 * (D - phi) * 1) / 2);
//     feq[2] = wc * (2 * phi - D + 3 * B.By / c + (3 * (D - phi) * 1) / 2);
//     feq[3] = wc * (2 * phi - D - 3 * B.Bx / c + (3 * (D - phi) * 1) / 2);
//     feq[4] = wc * (2 * phi - D - 3 * B.By / c + (3 * (D - phi) * 1) / 2);
//     feq[5] = wd * (2 * phi - D + (3 * (B.Bx + B.By)) / c + (3 * (D - phi) * 2) / 2);
//     feq[6] = wd * (2 * phi - D + (3 * (-B.Bx + B.By)) / c + (3 * (D - phi) * 2) / 2);
//     feq[7] = wd * (2 * phi - D - (3 * (B.Bx + B.By )) / c + (3 * (D - phi) * 2) / 2);
//     feq[8] = wd * (2 * phi - D - (3 * (-B.Bx + B.By)) / c + (3 * (D - phi) * 2) / 2);

//     getMR(meq, feq);
//     return meq[k];
// }

double CDE2019::feqADE2019(int k, double phi, double D, BADE B)
{
    switch (k)
    {
    case 0:
        return (w0 * (-D + 2 * phi));
    case 1:
        return (wc * (6 * B.Bx + c * (D + phi)) / (2 * c));
    case 2:
        return (wc * (6 * B.By + c * (D + phi)) / (2 * c));
    case 3:
        return (wc * (-6 * B.Bx + c * (D + phi)) / (2 * c));
    case 4:
        return (wc * (-6 * B.By + c * (D + phi)) / (2 * c));
    case 5:
        return (wd * (3 * B.Bx + 3 * B.By + c * (2 * D - phi)) / c);
    case 6:
        return (wd * (-3 * B.Bx + 3 * B.By + c * (2 * D - phi)) / c);
    case 7:
        return (-wd * (3 * B.Bx + 3 * B.By - c * (2 * D - phi)) / c);
    case 8:
        return (wd * (3 * B.Bx - 3 * B.By + c * (2 * D - phi)) / c);
        // case 1:
        //     return wc * (2 * phi - D + 3 * B.Bx / c + (3 * (D - phi) * 1) / 2);
        // case 2:
        //     return wc * (2 * phi - D + 3 * B.By / c + (3 * (D - phi) * 1) / 2);
        // case 3:
        //     return wc * (2 * phi - D - 3 * B.Bx / c + (3 * (D - phi) * 1) / 2);
        // case 4:
        //     return wc * (2 * phi - D - 3 * B.By / c + (3 * (D - phi) * 1) / 2);
        // case 5:
        //     return wd * (2 * phi - D + (3 * (B.Bx + B.By)) / c + (3 * (D - phi) * 2) / 2);
        // case 6:
        //     return wd * (2 * phi - D + (3 * (-B.Bx + B.By)) / c + (3 * (D - phi) * 2) / 2);
        // case 7:
        //     return wd * (2 * phi - D - (3 * (B.Bx + B.By)) / c + (3 * (D - phi) * 2) / 2);
        // case 8:
        //     return wd * (2 * phi - D - (3 * (-B.Bx + B.By)) / c + (3 * (D - phi) * 2) / 2);
    }
}

void CDE2019::get_each_phi_D_B(const double C_, const BADE B_, double &phi, double &D, BADE &B)
{
    phi = C_;
    D = C_;
    B = B_;
}

void CDE2019::init_scalar(const size_t n, const double *x, const double *y, const int *flag)
{
    double t = dt * n;
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t index = scalar_index(i, j);
            double y_ = y[index];
            double x_ = x[index];

            if (flag[index] == COMP)
            {
                C[index] = 0.0;
            }
            else
            {
                C[index] = 0.0;
                // C[index] = 1.0;
            }

            F[index] = 0;
        }
    }
}

void CDE2019::AnsScalar(size_t n, double *x, double *y)
{
    // initScalar(t, x, y);
    double t = dt * n;
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t index = scalar_index(i, j);
            double y_ = y[index];
            double x_ = x[index];
            C_old[index] = (t + 1) * sin(2 * PI * x_) * cos(2 * PI * y_);

            // double deltaD = sqrt(0.1 * t * 2);
            // C_old[index] = p2(delta0) / (p2(delta0) + p2(deltaD)) * C0 * exp(-(p2(x_ - 0.5) + p2(y_ - 0.5)) / (2 * (p2(delta0) + p2(deltaD))));
        }
    }
}
void CDE2019::initFeqAndMeq(const int *flag)
{
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t s_idx = scalar_index(i, j);
            // double phi = C[s_idx];
            // double D_ = sin(phi);
            // BADE B_ = {phi, phi};

            // double C_{C[s_idx]};
            // double phi{C_};
            // double D_{C_};
            // BADE B_{B_ade[s_idx].Bx, B_ade[s_idx].By};

            double phi, D_;
            BADE B_;

            get_each_phi_D_B(C[s_idx], B_ade[s_idx], phi, D_, B_);

            for (int k = 0; k < Q9; ++k)
            {
                size_t f_idx = field_index(i, j, k);
                f0[f_idx] = feqADE2019(k, phi, D_, B_);
                f1[f_idx] = meqADE2019(k, phi, D_, B_);
                // f1[f_idx] = feqADE2019(k, phi, D_, B_);
                if (flag[s_idx] == OBSTACLE)
                {
                    f0[f_idx] = 0;
                }
            }
        }
    }
}
void CDE2019::initFeqAndMeq()
{
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t s_idx = scalar_index(i, j);
            // double phi = C[s_idx];
            // double D_ = sin(phi);
            // BADE B_ = {phi, phi};

            // double C_{C[s_idx]};
            // double phi{C_};
            // double D_{C_};
            // BADE B_{B_ade[s_idx].Bx, B_ade[s_idx].By};

            double phi, D_;
            BADE B_;

            get_each_phi_D_B(C[s_idx], B_ade[s_idx], phi, D_, B_);

            for (int k = 0; k < Q9; ++k)
            {
                size_t f_idx = field_index(i, j, k);
                f0[f_idx] = feqADE2019(k, phi, D_, B_);
                f1[f_idx] = meqADE2019(k, phi, D_, B_);
                // f1[f_idx] = feqADE2019(k, phi, D_, B_);
            }
        }
    }
}

void CDE2019::getF(size_t t_, double *x, double *y)
{
    double t = t_ * dt;
    // std::cout << t << std::endl;
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t index = scalar_index(i, j);
            double y_ = y[index];
            double x_ = x[index];
            F[index] = sin(2 * PI * x_) * cos(2 * PI * y_)                                                                                           //
                       + 2 * PI * (t + 1) * cos(2 * PI * x_ + 2 * PI * y_)                                                                           //
                       + 0.4 * p2(PI) * p2(t + 1) * sin((t + 1) * sin(2 * PI * x_) * cos(2 * PI * y_)) * p2(cos(2 * PI * x_)) * p2(cos(2 * PI * y_)) //
                       + 0.4 * p2(PI) * p2(t + 1) * sin((t + 1) * sin(2 * PI * x_) * cos(2 * PI * y_)) * p2(sin(2 * PI * x_)) * p2(sin(2 * PI * y_)) //
                       + 0.8 * p2(PI) * (t + 1) * cos((t + 1) * sin(2 * PI * x_) * cos(2 * PI * y_)) * (sin(2 * PI * x_)) * (cos(2 * PI * y_));      //
        }
    }
}

void CDE2019::getMR(double *m, double *R)
{
    const double f0 = R[0];
    const double f1 = R[1];
    const double f2 = R[2];
    const double f3 = R[3];
    const double f4 = R[4];
    const double f5 = R[5];
    const double f6 = R[6];
    const double f7 = R[7];
    const double f8 = R[8];
    m[0] = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
    m[1] = (-4 * f0 - f1 - f2 - f3 - f4 + 2 * f5 + 2 * f6 + 2 * f7 + 2 * f8);
    m[2] = (4 * f0 - 2 * f1 - 2 * f2 - 2 * f3 - 2 * f4 + f5 + f6 + f7 + f8);
    m[3] = (f1 - f3 + f5 - f6 - f7 + f8);
    m[4] = (-2 * f1 + 2 * f3 + f5 - f6 - f7 + f8);
    m[5] = (f2 - f4 + f5 + f6 - f7 - f8);
    m[6] = (-2 * f2 + 2 * f4 + f5 + f6 - f7 - f8);
    m[7] = (f1 - f2 + f3 - f4);
    m[8] = (f5 - f6 + f7 - f8);
}

void CDE2019::m2f()
{
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            double F_ = F[scalar_index(i, j)];

            double m0 = f1[field_index(i, j, 0)];
            double m1 = f1[field_index(i, j, 1)];
            double m2 = f1[field_index(i, j, 2)];
            double m3 = f1[field_index(i, j, 3)];
            double m4 = f1[field_index(i, j, 4)];
            double m5 = f1[field_index(i, j, 5)];
            double m6 = f1[field_index(i, j, 6)];
            double m7 = f1[field_index(i, j, 7)];
            double m8 = f1[field_index(i, j, 8)];

            f1[field_index(i, j, 0)] = m0 / 9 - m1 / 9 + m2 / 9;
            f1[field_index(i, j, 1)] = m0 / 9 - m1 / 36 - m2 / 18 + m3 / 6 - m4 / 6 + m7 / 4;
            f1[field_index(i, j, 2)] = m0 / 9 - m1 / 36 - m2 / 18 + m5 / 6 - m6 / 6 - m7 / 4;
            f1[field_index(i, j, 3)] = m0 / 9 - m1 / 36 - m2 / 18 - m3 / 6 + m4 / 6 + m7 / 4;
            f1[field_index(i, j, 4)] = m0 / 9 - m1 / 36 - m2 / 18 - m5 / 6 + m6 / 6 - m7 / 4;
            f1[field_index(i, j, 5)] = m0 / 9 + m1 / 18 + m2 / 36 + m3 / 6 + m4 / 12 + m5 / 6 + m6 / 12 + m8 / 4;
            f1[field_index(i, j, 6)] = m0 / 9 + m1 / 18 + m2 / 36 - m3 / 6 - m4 / 12 + m5 / 6 + m6 / 12 - m8 / 4;
            f1[field_index(i, j, 7)] = m0 / 9 + m1 / 18 + m2 / 36 - m3 / 6 - m4 / 12 - m5 / 6 - m6 / 12 + m8 / 4;
            f1[field_index(i, j, 8)] = m0 / 9 + m1 / 18 + m2 / 36 + m3 / 6 + m4 / 12 - m5 / 6 - m6 / 12 - m8 / 4;
            // f1[field_index(i, j, 0)] += +dt * w0 * F_;
            // f1[field_index(i, j, 1)] += +dt * wc * F_;
            // f1[field_index(i, j, 2)] += +dt * wc * F_;
            // f1[field_index(i, j, 3)] += +dt * wc * F_;
            // f1[field_index(i, j, 4)] += +dt * wc * F_;
            // f1[field_index(i, j, 5)] += +dt * wd * F_;
            // f1[field_index(i, j, 6)] += +dt * wd * F_;
            // f1[field_index(i, j, 7)] += +dt * wd * F_;
            // f1[field_index(i, j, 8)] += +dt * wd * F_;
        }
    }
}

void CDE2019::streaming_cruve(const int *flag)
{
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t s_idx = scalar_index(i, j);
            if (flag[s_idx] == COMP)
            {
                unsigned int ip1 = (i + 1) % (NX);
                unsigned int jp1 = (j + 1) % (NY);
                unsigned int im1 = (NX + i - 1) % (NX);
                unsigned int jm1 = (NY + j - 1) % (NY);

                double ft0 = f1[field_index(i, j, 0)];
                double ft1 = f1[field_index(im1, j, 1)];
                double ft2 = f1[field_index(i, jm1, 2)];
                double ft3 = f1[field_index(ip1, j, 3)];
                double ft4 = f1[field_index(i, jp1, 4)];
                double ft5 = f1[field_index(im1, jm1, 5)];
                double ft6 = f1[field_index(ip1, jm1, 6)];
                double ft7 = f1[field_index(ip1, jp1, 7)];
                double ft8 = f1[field_index(im1, jp1, 8)];

                const size_t f_idx = field_index(i, j, 0);
                f0[f_idx + 0] = ft0;
                f0[f_idx + 1] = ft1;
                f0[f_idx + 2] = ft2;
                f0[f_idx + 3] = ft3;
                f0[f_idx + 4] = ft4;
                f0[f_idx + 5] = ft5;
                f0[f_idx + 6] = ft6;
                f0[f_idx + 7] = ft7;
                f0[f_idx + 8] = ft8;
                // f0[field_index(i, j, 0)] = ft0;
                // f0[field_index(i, j, 1)] = ft1;
                // f0[field_index(i, j, 2)] = ft2;
                // f0[field_index(i, j, 3)] = ft3;
                // f0[field_index(i, j, 4)] = ft4;
                // f0[field_index(i, j, 5)] = ft5;
                // f0[field_index(i, j, 6)] = ft6;
                // f0[field_index(i, j, 7)] = ft7;
                // f0[field_index(i, j, 8)] = ft8;

                //--------------------------------------------------------------------
                C[s_idx] = ft0                     //
                           + ft1 + ft2 + ft3 + ft4 //
                           + ft5 + ft6 + ft7 + ft8;
            }
            else
            {
                C[s_idx] = 1.0; //
            }
        }
    }
}

void CDE2019::BC_curve_NEE(const double *delta, const int *flag,
                           const double *x, const double *y)
{

    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t s_idx = scalar_index(i, j);
            double C_ = C[s_idx];
        }
    }
}

void CDE2019::streaming()
{
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            unsigned int ip1 = (i + 1) % (NX);
            unsigned int jp1 = (j + 1) % (NY);
            unsigned int im1 = (NX - 1 + i) % (NX);
            unsigned int jm1 = (NY - 1 + j) % (NY);
            // unsigned int ip1 = (i + 1) % (NX - 1);
            // unsigned int jp1 = (j + 1) % (NY - 1);
            // unsigned int im1 = (NX - 1 + i - 1) % (NX - 1);
            // unsigned int jm1 = (NY - 1 + j - 1) % (NY - 1);

            double ft0 = f1[field_index(i, j, 0)];
            double ft1 = f1[field_index(im1, j, 1)];
            double ft2 = f1[field_index(i, jm1, 2)];
            double ft3 = f1[field_index(ip1, j, 3)];
            double ft4 = f1[field_index(i, jp1, 4)];
            double ft5 = f1[field_index(im1, jm1, 5)];
            double ft6 = f1[field_index(ip1, jm1, 6)];
            double ft7 = f1[field_index(ip1, jp1, 7)];
            double ft8 = f1[field_index(im1, jp1, 8)];

            C[scalar_index(i, j)] = ft0                     //
                                    + ft1 + ft2 + ft3 + ft4 //
                                    + ft5 + ft6 + ft7 + ft8;

            f0[field_index(i, j, 0)] = ft0;
            f0[field_index(i, j, 1)] = ft1;
            f0[field_index(i, j, 2)] = ft2;
            f0[field_index(i, j, 3)] = ft3;
            f0[field_index(i, j, 4)] = ft4;
            f0[field_index(i, j, 5)] = ft5;
            f0[field_index(i, j, 6)] = ft6;
            f0[field_index(i, j, 7)] = ft7;
            f0[field_index(i, j, 8)] = ft8;
        }
    }
}

void CDE2019::collision()
{
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t index = scalar_index(i, j);

            double Fi[9];
            getFi(Fi, F[index]);

            double MF[9];
            getMR(MF, Fi);

            // // double phi = C_old[index];
            // double phi = C[index];
            // // double D_ = sin(phi);
            // double D_ = sin(phi);
            // BADE B_; // = {phi, phi};
            // B_.By = phi;
            // B_.Bx = phi;

            // double C_ = C[index];
            // double phi{C_};
            // double D_{C_};
            // BADE B_{B_ade[index].Bx, B_ade[index].By};
            // // BADE B_{B_ade[index].Bx * C_, B_ade[index].By * C_};

            double phi, D_;
            BADE B_;

            get_each_phi_D_B(C[index], B_ade[index], phi, D_, B_);

            double f0_[9];
            double m0[9];
            for (int k = 0; k < Q9; ++k)
            {
                f0_[k] = f0[field_index(i, j, k)];
            }
            getMR(m0, f0_);
            size_t fidx = field_index(i, j, 0);
            f1[fidx + 0] = m0[0] - s[0] * (m0[0] - meqADE2019(0, phi, D_, B_)) + MF[0];
            f1[fidx + 1] = m0[1] - s[1] * (m0[1] - meqADE2019(1, phi, D_, B_)) + MF[1];
            f1[fidx + 2] = m0[2] - s[2] * (m0[2] - meqADE2019(2, phi, D_, B_)) + MF[2];
            f1[fidx + 3] = m0[3] - s[3] * (m0[3] - meqADE2019(3, phi, D_, B_)) + MF[3];
            f1[fidx + 4] = m0[4] - s[4] * (m0[4] - meqADE2019(4, phi, D_, B_)) + MF[4];
            f1[fidx + 5] = m0[5] - s[5] * (m0[5] - meqADE2019(5, phi, D_, B_)) + MF[5];
            f1[fidx + 6] = m0[6] - s[6] * (m0[6] - meqADE2019(6, phi, D_, B_)) + MF[6];
            f1[fidx + 7] = m0[7] - s[7] * (m0[7] - meqADE2019(7, phi, D_, B_)) + MF[7];
            f1[fidx + 8] = m0[8] - s[8] * (m0[8] - meqADE2019(8, phi, D_, B_)) + MF[8];

            // for (int k = 0; k < Q9; ++k)
            // {
            //     size_t fidx = field_index(i, j, k);
            //     f1[fidx] = m0[k] - s[k] * (m0[k] - meqADE2019(k, phi, D_, B_)) + MF[k];
            // }
        }
    }
}

void CDE2019::collisionBGK()
{
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t index = scalar_index(i, j);
            double phi = C[index];
            // double phi = C_old[index];
            double D_ = sin(phi);
            // double D_ = sin(phi);
            BADE B_ = {phi, phi};

            for (int k = 0; k < Q9; ++k)
            {
                size_t f_idx = field_index(i, j, k);
                // f1[f_idx] = f0[f_idx] - s[3] * (f0[f_idx] - feqADE2019(k, phi, D_, B_)) + w[k] * dt * F[index];
                f1[f_idx] = (1 - s[3]) * f0[f_idx] + s[3] * feqADE2019(k, phi, D_, B_) + w[k] * dt * F[index];
                // f1[f_idx] = feqADE2019(k, phi, D_, B_)  + w[k] * dt * F[index];
            }
        }
    }
}

void CDE2019::streamingBGK()
{
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {

            // double C_ = 0;
            // for (int k = 0; k < Q9; ++k)
            // {
            //     unsigned int xmd = (NX + i - dirx[k]) % NX;
            //     unsigned int ymd = (NY + j - diry[k]) % NY;

            //     f0[field_index(i, j, k)] = f1[field_index(xmd, ymd, k)];

            //     C_ += f0[field_index(i, j, k)];
            // }
            // C[scalar_index(i, j)] = C_;

            unsigned int ip1 = (i + 1) % (NX - 1);
            unsigned int jp1 = (j + 1) % (NY - 1);
            unsigned int im1 = (NX - 1 + i - 1) % (NX - 1);
            unsigned int jm1 = (NY - 1 + j - 1) % (NY - 1);

            double ft0 = f1[field_index(i, j, 0)];
            double ft1 = f1[field_index(im1, j, 1)];
            double ft2 = f1[field_index(i, jm1, 2)];
            double ft3 = f1[field_index(ip1, j, 3)];
            double ft4 = f1[field_index(i, jp1, 4)];
            double ft5 = f1[field_index(im1, jm1, 5)];
            double ft6 = f1[field_index(ip1, jm1, 6)];
            double ft7 = f1[field_index(ip1, jp1, 7)];
            double ft8 = f1[field_index(im1, jp1, 8)];

            C[scalar_index(i, j)] = ft0 + ft1 + ft2 + ft3 + ft4 + ft5 + ft6 + ft7 + ft8;
            // C[scalar_index(i, j)] = C_old[scalar_index(i, j)];

            f0[field_index(i, j, 0)] = ft0;
            f0[field_index(i, j, 1)] = ft1;
            f0[field_index(i, j, 2)] = ft2;
            f0[field_index(i, j, 3)] = ft3;
            f0[field_index(i, j, 4)] = ft4;
            f0[field_index(i, j, 5)] = ft5;
            f0[field_index(i, j, 6)] = ft6;
            f0[field_index(i, j, 7)] = ft7;
            f0[field_index(i, j, 8)] = ft8;
        }
    }
}

void get_B(const NSE &nse, CDE2019 &cde)
{
    const size_t NX_ = cde.NX;
    const size_t NY_ = cde.NY;
    for (size_t j = 0; j < NY_; ++j)
    {
        for (size_t i = 0; i < NX_; ++i)
        {
            size_t s_idx = cde.scalar_index(i, j);
            cde.B_ade[s_idx].Bx = nse.get_ux(i, j) * cde.C[s_idx];
            cde.B_ade[s_idx].By = nse.get_uy(i, j) * cde.C[s_idx];

            // cout << nse.get_ux(i, j) << " , " << nse.get_uy(i, j) << endl;
            // cout << cde.B_ade[s_idx].Bx << " , " << cde.B_ade[s_idx].By << endl;
        }
    }
}

double CDE2019::errC()
{
    double temp1 = 0;
    double temp2 = 0;
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t idx = scalar_index(i, j);
            temp1 += p2(C[idx] - C_old[idx]);
            temp2 += p2(C_old[idx]);
        }
    }
    double err = sqrt(temp1) / sqrt(temp2);
    return err;
}

// void CDE2019::streaming_cruve(const int *flag)
// {
//     for (int j = 0; j < NY; ++j)
//     {
//         for (int i = 0; i < NX; ++i)
//         {
//             size_t s_idx = scalar_index(i, j);
//             if (flag[s_idx] == COMP)
//             {
//                 unsigned int ip1 = (i + 1) % (NX);
//                 unsigned int jp1 = (j + 1) % (NY);
//                 unsigned int im1 = (NX + i - 1) % (NX);
//                 unsigned int jm1 = (NY + j - 1) % (NY);

//                 double ft0 = f1[field_index(i, j, 0)];
//                 double ft1 = f1[field_index(im1, j, 1)];
//                 double ft2 = f1[field_index(i, jm1, 2)];
//                 double ft3 = f1[field_index(ip1, j, 3)];
//                 double ft4 = f1[field_index(i, jp1, 4)];
//                 double ft5 = f1[field_index(im1, jm1, 5)];
//                 double ft6 = f1[field_index(ip1, jm1, 6)];
//                 double ft7 = f1[field_index(ip1, jp1, 7)];
//                 double ft8 = f1[field_index(im1, jp1, 8)];

//                 f0[field_index(i, j, 0)] = ft0;
//                 f0[field_index(i, j, 1)] = ft1;
//                 f0[field_index(i, j, 2)] = ft2;
//                 f0[field_index(i, j, 3)] = ft3;
//                 f0[field_index(i, j, 4)] = ft4;
//                 f0[field_index(i, j, 5)] = ft5;
//                 f0[field_index(i, j, 6)] = ft6;
//                 f0[field_index(i, j, 7)] = ft7;
//                 f0[field_index(i, j, 8)] = ft8;

//                 //----------------------- computing macro valiables-----------------
//                 double C_ = ft0                     //
//                               + ft1 + ft2 + ft3 + ft4 //
//                               + ft5 + ft6 + ft7 + ft8;

//                 C[s_idx] = C_;
//             }
//             else
//             {
//                 C[s_idx] = 1;
//             }
//         }
//     }
// }

void CDE2019::BC_curve_Guo(const double *delta, const int *flag,
                           const double *x, const double *y)
{
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            const size_t s_idx = scalar_index(i, j);

            const double C_f = C[s_idx];
            // const double phi_f{C_f};
            // const double D_f{C_f};
            // const BADE B_f{B_ade[s_idx]};

            double phi_f, D_f;
            BADE B_f;

            get_each_phi_D_B(C_f, B_ade[s_idx], phi_f, D_f, B_f);

            if (flag[s_idx] == COMP)
            {
                for (int k = 0; k < Q9; k++)
                {
                    const int ip = i - dirx[k];
                    const int jp = j - diry[k];
                    const size_t f_idx = field_index(i, j, k);

                    const size_t s_p = scalar_index(ip, jp);

                    if (flag[s_p] == OBSTACLE)
                    {
                        const int xff = i + dirx[k];
                        const int yff = j + diry[k];
                        const size_t sff_idx = scalar_index(xff, yff);
                        const size_t f_ff = field_index(xff, yff, k);

                        //***********fluid
                        double C_ghost, phi_ghost, D_ghost;
                        BADE B_ghost;

                        double fneq_ghost;

                        const double C_ff{C[sff_idx]};
                        // const double phi_ff{C_ff}, D_ff{C_ff};

                        // const BADE B_ff{B_ade[sff_idx]};

                        double phi_ff, D_ff;
                        BADE B_ff;

                        get_each_phi_D_B(C_ff, B_ade[sff_idx], phi_ff, D_ff, B_ff);

                        const double Th = 1.0;
                        const double q_delta = delta[f_idx];
                        // cout << q_delta << " " << phi_ff << " " << D_ff << " " << B_ff.Bx << " " << B_ff.By << endl;

                        const double Tb1 = (Th + (q_delta - 1.0) * C_f) / q_delta;
                        const double Tb2 = (2.0 * Th + (q_delta - 1.0) * C_ff) / (1.0 + q_delta);

                        // cout << Tb1 << " " << Tb2 << " " << C_f << " " << C_ff << endl;

                        const BADE B_b{0.0, 0.0};

                        double ux_b1 = (B_b.Bx + (q_delta - 1.0) * B_f.Bx) / q_delta;
                        double uy_b1 = (B_b.By + (q_delta - 1.0) * B_f.By) / q_delta;
                        double ux_b2 = (2.0 * B_b.Bx + (q_delta - 1.0) * B_ff.Bx) / (1.0 + q_delta);
                        double uy_b2 = (2.0 * B_b.By + (q_delta - 1.0) * B_ff.By) / (1.0 + q_delta);

                        const BADE B_b1{ux_b1, uy_b1};
                        const BADE B_b2{ux_b2, uy_b2};

                        double Tneq1 = f1[f_idx] - feqADE2019(k, phi_f, D_f, B_f);
                        double Tneq2 = f1[f_ff] - feqADE2019(k, phi_ff, D_ff, B_ff);

                        // if (q_delta > 0)
                        // cout << q_delta << " " << Tb1 << " " << Tb2 << " " << ux_b1 << " " << uy_b1 << endl;

                        // C_ghost = Tb1;
                        // B_ghost.Bx = 0;
                        // B_ghost.By = 0;
                        // fneq_ghost = Tneq1;
                        if (q_delta >= 0.75)
                        {
                            C_ghost = Tb1;
                            B_ghost.Bx = B_b1.Bx;
                            B_ghost.By = B_b1.By;
                            fneq_ghost = Tneq1;
                        }
                        else
                        {
                            C_ghost = q_delta * Tb1 + (1.0 - q_delta) * Tb2;
                            B_ghost.Bx = q_delta * B_b1.Bx + (1.0 - q_delta) * B_b2.Bx;
                            B_ghost.By = q_delta * B_b1.By + (1.0 - q_delta) * B_b2.By;
                            fneq_ghost = q_delta * Tneq1 + (1.0 - q_delta) * Tneq2;
                        }
                        phi_ghost = D_ghost = C_ghost;

                        // cout << phi_ghost << " " << D_ghost << " " << B_ghost.Bx << " " << B_ghost.By << endl;
                        //
                        f1[field_index(ip, jp, k)] = feqADE2019(k, phi_ghost, D_ghost, B_ghost) + (1.0 - 1.0 / tau) * fneq_ghost;
                        // f1[field_index(ip, jp, k)] = feqADE2019(k, phi_ghost, D_ghost, B_ghost);
                    }
                }
            }
        }
    }
}

void CDE2019::BC_NEE_Dirichlet()
{
    //------------------------------top and bottom ------------------------------------
    for (int i = 0; i < NX; ++i)
    {
        size_t bottom = scalar_index(i, 0);
        size_t bottom_neighborhood = scalar_index(i, 0 + 1);

        size_t top = scalar_index(i, NY - 1);
        size_t top_neighborhood = scalar_index(i, NY - 2);

        C[top] = 0;
        C[bottom] = 0.0;
        // C[top] = C[top_neighborhood];
        // C[bottom] = C[bottom_neighborhood];

        size_t flag_this;
        size_t flag_neigh;
        double C_this;
        double C_neigh;
        double phi_this, D_this;
        double phi_neigh, D_neigh;
        BADE B_this;
        BADE B_neigh;

        // -------------top--------
        flag_this = NY - 1;
        flag_neigh = NY - 2;

        C_this = C[top];


        get_each_phi_D_B(C_this, B_ade[top], phi_this, D_this, B_this);

        // B_this.Bx = B_ade[top].Bx;
        // B_this.By = B_ade[top].By;
        // phi_this = D_this = C_this;

        C_neigh = C[top_neighborhood];

        get_each_phi_D_B(C_neigh, B_ade[top_neighborhood], phi_neigh, D_neigh, B_neigh);


        // B_neigh.Bx = B_ade[top_neighborhood].Bx;
        // B_neigh.By = B_ade[top_neighborhood].By;
        // phi_neigh = D_neigh = C_neigh;
        for (int k = 0; k < Q9; k += 3)
        {
            f0[field_index(i, flag_this, k + 0)] = feqADE2019(k + 0, phi_this, D_this, B_this) + f0[field_index(i, flag_neigh, k + 0)] - feqADE2019(k + 0, phi_neigh, D_neigh, B_neigh);
            f0[field_index(i, flag_this, k + 1)] = feqADE2019(k + 1, phi_this, D_this, B_this) + f0[field_index(i, flag_neigh, k + 1)] - feqADE2019(k + 1, phi_neigh, D_neigh, B_neigh);
            f0[field_index(i, flag_this, k + 2)] = feqADE2019(k + 2, phi_this, D_this, B_this) + f0[field_index(i, flag_neigh, k + 2)] - feqADE2019(k + 2, phi_neigh, D_neigh, B_neigh);
        }

        // ----------------------------------------------bottom------------------------------------------------
        flag_this = 0;
        flag_neigh = 1;
        C_this = C[bottom];

        get_each_phi_D_B(C_this, B_ade[bottom], phi_this, D_this, B_this);


        // B_this.Bx = B_ade[bottom].Bx;
        // B_this.By = B_ade[bottom].By;
        // phi_this = D_this = C_this;

        C_neigh = C[bottom_neighborhood];

        get_each_phi_D_B(C_neigh, B_ade[bottom_neighborhood], phi_neigh, D_neigh, B_neigh);

        // B_neigh.Bx = B_ade[bottom_neighborhood].Bx;
        // B_neigh.By = B_ade[bottom_neighborhood].By;
        // phi_neigh = D_neigh = C_neigh;

        for (int k = 0; k < Q9; k += 3)
        {
            f0[field_index(i, flag_this, k + 0)] = feqADE2019(k + 0, phi_this, D_this, B_this) + f0[field_index(i, flag_neigh, k + 0)] - feqADE2019(k + 0, phi_neigh, D_neigh, B_neigh);
            f0[field_index(i, flag_this, k + 1)] = feqADE2019(k + 1, phi_this, D_this, B_this) + f0[field_index(i, flag_neigh, k + 1)] - feqADE2019(k + 1, phi_neigh, D_neigh, B_neigh);
            f0[field_index(i, flag_this, k + 2)] = feqADE2019(k + 2, phi_this, D_this, B_this) + f0[field_index(i, flag_neigh, k + 2)] - feqADE2019(k + 2, phi_neigh, D_neigh, B_neigh);
        }
    }

    //-----------------------------------------------------------------left and right -------------------------------------------------------------
    for (size_t j = 0; j < NY; ++j)
    {
        size_t left = scalar_index(0, j);
        size_t left_neighborhood = scalar_index(0 + 1, j);

        size_t right = scalar_index(NX - 1, j);
        size_t right_neighborhood = scalar_index(NX - 2, j);

        C[right] = 0;
        C[left] = 0.0;
        // C[right] = C[right_neighborhood];
        // C[left] = C[left_neighborhood];

        size_t flag_this;
        size_t flag_neigh;
        double C_this;
        double C_neigh;
        double phi_this, D_this;
        double phi_neigh, D_neigh;
        BADE B_this;
        BADE B_neigh;

        // -----------------------------------------------------------------------left-----------------------------------------------
        flag_this = 0;
        flag_neigh = 1;

        C_this = C[left];

        get_each_phi_D_B(C_this, B_ade[left], phi_this, D_this, B_this);


        // B_this.Bx = B_ade[left].Bx;
        // B_this.By = B_ade[left].By;
        // phi_this = D_this = C_this;

        C_neigh = C[left_neighborhood];

        get_each_phi_D_B(C_neigh, B_ade[left_neighborhood], phi_neigh, D_neigh, B_neigh);

        // B_neigh.Bx = B_ade[left_neighborhood].Bx;
        // B_neigh.By = B_ade[left_neighborhood].By;
        // phi_neigh = D_neigh = C_neigh;

        for (int k = 0; k < 9; k += 3)
        {
            f0[field_index(flag_this, j, k + 0)] = feqADE2019(k + 0, phi_this, D_this, B_this) + f0[field_index(flag_neigh, j, k + 0)] - feqADE2019(k + 0, phi_neigh, D_neigh, B_neigh);
            f0[field_index(flag_this, j, k + 1)] = feqADE2019(k + 1, phi_this, D_this, B_this) + f0[field_index(flag_neigh, j, k + 1)] - feqADE2019(k + 1, phi_neigh, D_neigh, B_neigh);
            f0[field_index(flag_this, j, k + 2)] = feqADE2019(k + 2, phi_this, D_this, B_this) + f0[field_index(flag_neigh, j, k + 2)] - feqADE2019(k + 2, phi_neigh, D_neigh, B_neigh);
        }

        // ---------------------------------------------------------------- right---------------------------------------------------------------------
        flag_this = NX - 1;
        flag_neigh = NX - 2;
        C_this = C[right];

        get_each_phi_D_B(C_this, B_ade[right], phi_this, D_this, B_this);


        // B_this.Bx = B_ade[right].Bx;
        // B_this.By = B_ade[right].By;
        // phi_this = D_this = C_this;

        C_neigh = C[right_neighborhood];
        get_each_phi_D_B(C_neigh, B_ade[right_neighborhood], phi_neigh, D_neigh, B_neigh);

        // B_neigh.Bx = B_ade[right_neighborhood].Bx;
        // B_neigh.By = B_ade[right_neighborhood].By;
        // phi_neigh = D_neigh = C_neigh;

        for (int k = 0; k < 9; k += 3)
        {
            f0[field_index(flag_this, j, k + 0)] = feqADE2019(k + 0, phi_this, D_this, B_this) + f0[field_index(flag_neigh, j, k + 0)] - feqADE2019(k + 0, phi_neigh, D_neigh, B_neigh);
            f0[field_index(flag_this, j, k + 1)] = feqADE2019(k + 1, phi_this, D_this, B_this) + f0[field_index(flag_neigh, j, k + 1)] - feqADE2019(k + 1, phi_neigh, D_neigh, B_neigh);
            f0[field_index(flag_this, j, k + 2)] = feqADE2019(k + 2, phi_this, D_this, B_this) + f0[field_index(flag_neigh, j, k + 2)] - feqADE2019(k + 2, phi_neigh, D_neigh, B_neigh);
        }
    }
}

void CDE2019::outputTec(double t, double *x, double *y)
{
    std::ostringstream name;
    name << "CDE_t_" << t << "_iter"
         << ".dat";
    std::ofstream out(name.str().c_str());
    out << "Title= \"Lid Driven Flow\"\n";
    out << "VARIABLES = \"X\", \"Y\", \"C\",  \"Co\", \"Bx\", \"By\" \n";
    out << "ZONE T= \"BOX\",I=" << NX << ",J=" << NY << ",F=	POINT" << std::endl;

    for (size_t j = 0; j < NY; ++j)
    {
        for (size_t i = 0; i < NX; ++i)
        {
            const size_t index = scalar_index(i, j);
            const size_t findex = field_index(i, j, 0);
            out << std::fixed << std::setprecision(10)
                << " " << x[index] << " " << y[index]
                // << " " << f0[findex]
                // << " " << f1[findex]
                << " " << C[index]
                << " " << C_old[index]
                << " " << B_ade[index].Bx
                << " " << B_ade[index].By
                << std::endl;
        }
    }
    out.close();
}
