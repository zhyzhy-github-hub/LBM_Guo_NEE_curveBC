#include "NSE.h"
#include "CDE2019.h"

using std::cout;
using std::endl;

// I added this
// 123456
NSE::NSE(const int NX_, const int NY_, const double nu_, const double c_)
    : NX(NX_), NY(NY_), nu(nu_), c(c_)
{
    const unsigned int NN = NX * NY;
    s = new double[Q9];
    ux = new double[NN];
    uy = new double[NN];
    rho = new double[NN];
    ux_old = new double[NN];
    uy_old = new double[NN];
    force_x = new double[NN];
    force_y = new double[NN];
    f_col = new double[NN * Q9];
    f_str = new double[NN * Q9];

    dx = 1.0 / (NY - 1);
    // dx = 1.0 ;
    dt = dx / c;
    Cs = c / sqrt(3);
    Cs2 = pow(Cs, 2);

    tau = 0.5 + nu / (Cs2 * dt);
    // tau = 0.5 + nu / (Cs2);

    cout << "tau = " << tau << endl;
    cout << "nu = " << nu << endl;

    const double s_nu = 1.0 / tau;
    const double s_e = s_nu;
    const double s_q = 8 * (2 - s_nu) / (8 - s_nu);
    const double s_eps = s_nu;

    s[0] = 0.0;
    s[1] = s_e;
    s[2] = s_eps;
    s[3] = 0.0;
    s[4] = s_q;
    s[5] = 0.0;
    s[6] = s_q;
    s[7] = s_nu;
    s[8] = s_nu;

    cout << "--------------NSE paramater----------------" << endl;
    cout << "dx = " << dx << endl;
    cout << "dt = " << dt << endl;
    cout << "nu = " << nu << endl;
    cout << "tau = " << tau << endl;
    cout << " c = " << c << endl;
    for (int i = 0; i < Q9; ++i)
    {
        cout << "s[" << i << "]=" << s[i] << " | ";
    }
    cout << endl;
}

NSE::~NSE()
{
    delete[] s;
    delete[] ux;
    delete[] uy;
    delete[] rho;
    delete[] ux_old;
    delete[] uy_old;
    delete[] force_x;
    delete[] force_y;
    delete[] f_col;
    delete[] f_str;
}

size_t NSE::scalar_index(const size_t i, const size_t j) const
{
    return (NX * j + i);
}

size_t NSE::field_index(const size_t x, const size_t y, const size_t d) const
{
    // return (NX * (NY * (d) + y) + x);
    // return (9 * (NX * y + x) + d);
    return (Q9 * (NX * y + x) + d);
}
double NSE::meqNSE2019(int k, double rho, double ux, double uy)
{
    switch (k)
    {
    case 0:
        return (rho);
        break;
    case 1:
        return (rho / (c * c) * (3 * ux * ux + 3 * uy * uy - 2 * c * c));
        break;
    case 2:
        return (rho / (c * c) * (-3 * ux * ux - 3 * uy * uy + 1 * c * c));
        break;
    case 3:
        return (rho * ux) / c;
        break;
    case 4:
        return ((-rho * ux) / c);
        break;
    case 5:
        return (rho * uy) / c;
        break;
    case 6:
        return ((-rho * uy) / c);
        break;
    case 7:
        return (rho / (c * c) * (ux * ux - uy * uy));
        break;
    case 8:
        return (rho / (c * c) * ux * uy);
        break;
    default:
        break;
    }
}

double NSE::feqNSE2019(int k, double rho, double ux, double uy)
{
    double feq;
    //double CU;
    double CU_Cs2inv;
    const double Cs2_inv = 1.0 / Cs2;
    const double Omusq = 1.0 - (ux * ux + uy * uy) * 0.5 * Cs2_inv;
    switch (k)
    {
    case 0:
        feq = w0 * rho * Omusq;
        break;
    case 1:
        CU_Cs2inv = ux * Cs2_inv * c;
        feq = wc * rho * (CU_Cs2inv + 0.5 * (CU_Cs2inv) * (CU_Cs2inv) + Omusq);
        break;
    case 2:
        CU_Cs2inv = uy * Cs2_inv * c;
        feq = wc * rho * (CU_Cs2inv + 0.5 * (CU_Cs2inv) * (CU_Cs2inv) + Omusq);
        break;
    case 3:
        CU_Cs2inv = -ux * Cs2_inv * c;
        feq = wc * rho * (CU_Cs2inv + 0.5 * (CU_Cs2inv) * (CU_Cs2inv) + Omusq);
        break;
    case 4:
        CU_Cs2inv = -uy * Cs2_inv * c;
        feq = wc * rho * (CU_Cs2inv + 0.5 * (CU_Cs2inv) * (CU_Cs2inv) + Omusq);
        break;
    case 5:
        CU_Cs2inv = (ux + uy) * Cs2_inv * c;
        feq = wd * rho * (CU_Cs2inv + 0.5 * (CU_Cs2inv) * (CU_Cs2inv) + Omusq);
        break;
    case 6:
        CU_Cs2inv = (-ux + uy) * Cs2_inv * c;
        feq = wd * rho * (CU_Cs2inv + 0.5 * (CU_Cs2inv) * (CU_Cs2inv) + Omusq);
        break;
    case 7:
        CU_Cs2inv = (-ux - uy) * Cs2_inv * c;
        feq = wd * rho * (CU_Cs2inv + 0.5 * (CU_Cs2inv) * (CU_Cs2inv) + Omusq);
        break;
    case 8:
        CU_Cs2inv = (ux - uy) * Cs2_inv * c;
        feq = wd * rho * (CU_Cs2inv + 0.5 * (CU_Cs2inv) * (CU_Cs2inv) + Omusq);
        break;
    }
    return feq;
}

void NSE::initScalar(const size_t n, const double *x, const double *y)
{
    double t = dt * n;
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t index = scalar_index(i, j);
            double y_ = y[index];
            double x_ = x[index];

            rho[index] = rho0;
            ux[index] = 0.0;
            uy[index] = 0.0;
            ux_old[index] = 0.0;
            uy_old[index] = 0.0;
            force_x[index] = 0.0;
            force_y[index] = 0.0;
        }
    }
}

void NSE::BC_NEE_Guo()
{
    //------------------------------top and bottom ------------------------------------
    for (int i = 0; i < NX; ++i)
    {
        size_t bottom = scalar_index(i, 0);
        size_t bottom_neighborhood = scalar_index(i, 0 + 1);

        size_t top = scalar_index(i, NY - 1);
        size_t top_neighborhood = scalar_index(i, NY - 2);

        rho[top] = rho[top_neighborhood];
        ux[top] = u_max * 1;
        uy[top] = 0.0;
        rho[bottom] = rho[bottom_neighborhood];
        ux[bottom] = 0.0;
        uy[bottom] = 0.0;

        size_t flag_this;
        size_t flag_neigh;
        double rho_this, ux_this, uy_this;
        double rho_neigh, ux_neigh, uy_neigh;

        // -------------top--------
        flag_this = NY - 1;
        flag_neigh = NY - 2;

        rho_this = rho[top];
        ux_this = ux[top];
        uy_this = uy[top];
        rho_neigh = rho[top_neighborhood];
        ux_neigh = ux[top_neighborhood];
        uy_neigh = uy[top_neighborhood];
        for (int k = 0; k < 9; k += 3)
        {
            f_str[field_index(i, flag_this, k + 0)] = feqNSE2019(k + 0, rho_this, ux_this, uy_this) + f_str[field_index(i, flag_neigh, k + 0)] - feqNSE2019(k + 0, rho_neigh, ux_neigh, uy_neigh);
            f_str[field_index(i, flag_this, k + 1)] = feqNSE2019(k + 1, rho_this, ux_this, uy_this) + f_str[field_index(i, flag_neigh, k + 1)] - feqNSE2019(k + 1, rho_neigh, ux_neigh, uy_neigh);
            f_str[field_index(i, flag_this, k + 2)] = feqNSE2019(k + 2, rho_this, ux_this, uy_this) + f_str[field_index(i, flag_neigh, k + 2)] - feqNSE2019(k + 2, rho_neigh, ux_neigh, uy_neigh);
        }
        // for (int k = 0; k < 3; ++k)
        // {
        //     f_str[field_index(i, flag_this, 3 * k + 0)] = feqNSE2019(3 * k + 0, rho_this, ux_this, uy_this) + f_str[field_index(i, flag_neigh, 3 * k + 0)] - feqNSE2019(3 * k + 0, rho_neigh, ux_neigh, uy_neigh);
        //     f_str[field_index(i, flag_this, 3 * k + 1)] = feqNSE2019(3 * k + 1, rho_this, ux_this, uy_this) + f_str[field_index(i, flag_neigh, 3 * k + 1)] - feqNSE2019(3 * k + 1, rho_neigh, ux_neigh, uy_neigh);
        //     f_str[field_index(i, flag_this, 3 * k + 2)] = feqNSE2019(3 * k + 2, rho_this, ux_this, uy_this) + f_str[field_index(i, flag_neigh, 3 * k + 2)] - feqNSE2019(3 * k + 2, rho_neigh, ux_neigh, uy_neigh);
        // }
        // for (int k = 0; k < Q9; ++k)
        // {
        //     f_str[field_index(i, flag_this, k)] = feqNSE2019(k, rho_this, ux_this, uy_this) + f_str[field_index(i, flag_neigh, k)] - feqNSE2019(k, rho_neigh, ux_neigh, uy_neigh);
        // }

        // ----------------------------------------------bottom------------------------------------------------
        flag_this = 0;
        flag_neigh = 1;
        rho_this = rho[bottom];
        ux_this = ux[bottom];
        uy_this = uy[bottom];
        rho_neigh = rho[bottom_neighborhood];
        ux_neigh = ux[bottom_neighborhood];
        uy_neigh = uy[bottom_neighborhood];
        // for (int k = 0; k < Q9; ++k)
        // {
        //     f_str[field_index(i, flag_this, k)] = feqNSE2019(k, rho_this, ux_this, uy_this) + f_str[field_index(i, flag_neigh, k)] - feqNSE2019(k, rho_neigh, ux_neigh, uy_neigh);
        // }
        for (int k = 0; k < 9; k += 3)
        {
            f_str[field_index(i, flag_this, k + 0)] = feqNSE2019(k + 0, rho_this, ux_this, uy_this) + f_str[field_index(i, flag_neigh, k + 0)] - feqNSE2019(k + 0, rho_neigh, ux_neigh, uy_neigh);
            f_str[field_index(i, flag_this, k + 1)] = feqNSE2019(k + 1, rho_this, ux_this, uy_this) + f_str[field_index(i, flag_neigh, k + 1)] - feqNSE2019(k + 1, rho_neigh, ux_neigh, uy_neigh);
            f_str[field_index(i, flag_this, k + 2)] = feqNSE2019(k + 2, rho_this, ux_this, uy_this) + f_str[field_index(i, flag_neigh, k + 2)] - feqNSE2019(k + 2, rho_neigh, ux_neigh, uy_neigh);
        }
        // for (int k = 0; k < 3; ++k)
        // {
        //     f_str[field_index(i, flag_this, 3 * k + 0)] = feqNSE2019(3 * k + 0, rho_this, ux_this, uy_this) + f_str[field_index(i, flag_neigh, 3 * k + 0)] - feqNSE2019(3 * k + 0, rho_neigh, ux_neigh, uy_neigh);
        //     f_str[field_index(i, flag_this, 3 * k + 1)] = feqNSE2019(3 * k + 1, rho_this, ux_this, uy_this) + f_str[field_index(i, flag_neigh, 3 * k + 1)] - feqNSE2019(3 * k + 1, rho_neigh, ux_neigh, uy_neigh);
        //     f_str[field_index(i, flag_this, 3 * k + 2)] = feqNSE2019(3 * k + 2, rho_this, ux_this, uy_this) + f_str[field_index(i, flag_neigh, 3 * k + 2)] - feqNSE2019(3 * k + 2, rho_neigh, ux_neigh, uy_neigh);
        // }
    }

    //-----------------------------------------------------------------left and right -------------------------------------------------------------
    for (size_t j = 0; j < NY; ++j)
    {
        size_t left = scalar_index(0, j);
        size_t left_neighborhood = scalar_index(0 + 1, j);

        size_t right = scalar_index(NX - 1, j);
        size_t right_neighborhood = scalar_index(NX - 2, j);

        rho[right] = rho[right_neighborhood];
        ux[right] = 0.0;
        uy[right] = 0.0;
        rho[left] = rho[left_neighborhood];
        ux[left] = 0.0;
        uy[left] = 0.0;

        size_t flag_this;
        size_t flag_neigh;
        double rho_this, ux_this, uy_this;
        double rho_neigh, ux_neigh, uy_neigh;

        // -----------------------------------------------------------------------left-----------------------------------------------
        flag_this = 0;
        flag_neigh = 1;

        rho_this = rho[left];
        ux_this = ux[left];
        uy_this = uy[left];
        rho_neigh = rho[left_neighborhood];
        ux_neigh = ux[left_neighborhood];
        uy_neigh = uy[left_neighborhood];
        for (int k = 0; k < 9; k += 3)
        {
            f_str[field_index(flag_this, j, k + 0)] = feqNSE2019(k + 0, rho_this, ux_this, uy_this) + f_str[field_index(flag_neigh, j, k + 0)] - feqNSE2019(k + 0, rho_neigh, ux_neigh, uy_neigh);
            f_str[field_index(flag_this, j, k + 1)] = feqNSE2019(k + 1, rho_this, ux_this, uy_this) + f_str[field_index(flag_neigh, j, k + 1)] - feqNSE2019(k + 1, rho_neigh, ux_neigh, uy_neigh);
            f_str[field_index(flag_this, j, k + 2)] = feqNSE2019(k + 2, rho_this, ux_this, uy_this) + f_str[field_index(flag_neigh, j, k + 2)] - feqNSE2019(k + 2, rho_neigh, ux_neigh, uy_neigh);
        }
        // for (int k = 0; k < 3; ++k)
        // {
        //     f_str[field_index(flag_this, j, 3 * k + 0)] = feqNSE2019(3 * k + 0, rho_this, ux_this, uy_this) + f_str[field_index(flag_neigh, j, 3 * k + 0)] - feqNSE2019(3 * k + 0, rho_neigh, ux_neigh, uy_neigh);
        //     f_str[field_index(flag_this, j, 3 * k + 1)] = feqNSE2019(3 * k + 1, rho_this, ux_this, uy_this) + f_str[field_index(flag_neigh, j, 3 * k + 1)] - feqNSE2019(3 * k + 1, rho_neigh, ux_neigh, uy_neigh);
        //     f_str[field_index(flag_this, j, 3 * k + 2)] = feqNSE2019(3 * k + 2, rho_this, ux_this, uy_this) + f_str[field_index(flag_neigh, j, 3 * k + 2)] - feqNSE2019(3 * k + 2, rho_neigh, ux_neigh, uy_neigh);
        // }
        // for (int k = 0; k < Q9; ++k)
        // {
        //     f_str[field_index(flag_this, j, k)] = feqNSE2019(k, rho_this, ux_this, uy_this) + f_str[field_index(flag_neigh, j, k)] - feqNSE2019(k, rho_neigh, ux_neigh, uy_neigh);
        // }
        // ---------------------------------------------------------------- right---------------------------------------------------------------------
        flag_this = NX - 1;
        flag_neigh = NX - 2;

        rho_this = rho[right];
        ux_this = ux[right];
        uy_this = uy[right];
        rho_neigh = rho[right_neighborhood];
        ux_neigh = ux[right_neighborhood];
        uy_neigh = uy[right_neighborhood];
        // for (int k = 0; k < Q9; ++k)
        // {
        //     f_str[field_index(flag_this, j, k)] = feqNSE2019(k, rho_this, ux_this, uy_this) + f_str[field_index(flag_neigh, j, k)] - feqNSE2019(k, rho_neigh, ux_neigh, uy_neigh);
        // }
        for (int k = 0; k < 9; k += 3)
        {
            f_str[field_index(flag_this, j, k + 0)] = feqNSE2019(k + 0, rho_this, ux_this, uy_this) + f_str[field_index(flag_neigh, j, k + 0)] - feqNSE2019(k + 0, rho_neigh, ux_neigh, uy_neigh);
            f_str[field_index(flag_this, j, k + 1)] = feqNSE2019(k + 1, rho_this, ux_this, uy_this) + f_str[field_index(flag_neigh, j, k + 1)] - feqNSE2019(k + 1, rho_neigh, ux_neigh, uy_neigh);
            f_str[field_index(flag_this, j, k + 2)] = feqNSE2019(k + 2, rho_this, ux_this, uy_this) + f_str[field_index(flag_neigh, j, k + 2)] - feqNSE2019(k + 2, rho_neigh, ux_neigh, uy_neigh);
        }
        // for (int k = 0; k < 3; ++k)
        // {
        //     f_str[field_index(flag_this, j, 3 * k + 0)] = feqNSE2019(3 * k + 0, rho_this, ux_this, uy_this) + f_str[field_index(flag_neigh, j, 3 * k + 0)] - feqNSE2019(3 * k + 0, rho_neigh, ux_neigh, uy_neigh);
        //     f_str[field_index(flag_this, j, 3 * k + 1)] = feqNSE2019(3 * k + 1, rho_this, ux_this, uy_this) + f_str[field_index(flag_neigh, j, 3 * k + 1)] - feqNSE2019(3 * k + 1, rho_neigh, ux_neigh, uy_neigh);
        //     f_str[field_index(flag_this, j, 3 * k + 2)] = feqNSE2019(3 * k + 2, rho_this, ux_this, uy_this) + f_str[field_index(flag_neigh, j, 3 * k + 2)] - feqNSE2019(3 * k + 2, rho_neigh, ux_neigh, uy_neigh);
        // }
    }
}

void NSE::streaming_curve(const int *flag)
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

                double ft0 = f_col[field_index(i, j, 0)];
                double ft1 = f_col[field_index(im1, j, 1)];
                double ft2 = f_col[field_index(i, jm1, 2)];
                double ft3 = f_col[field_index(ip1, j, 3)];
                double ft4 = f_col[field_index(i, jp1, 4)];
                double ft5 = f_col[field_index(im1, jm1, 5)];
                double ft6 = f_col[field_index(ip1, jm1, 6)];
                double ft7 = f_col[field_index(ip1, jp1, 7)];
                double ft8 = f_col[field_index(im1, jp1, 8)];

                const size_t f_idx = field_index(i, j, 0);
                f_str[f_idx + 0] = ft0;
                f_str[f_idx + 1] = ft1;
                f_str[f_idx + 2] = ft2;
                f_str[f_idx + 3] = ft3;
                f_str[f_idx + 4] = ft4;
                f_str[f_idx + 5] = ft5;
                f_str[f_idx + 6] = ft6;
                f_str[f_idx + 7] = ft7;
                f_str[f_idx + 8] = ft8;
                // f_str[field_index(i, j, 0)] = ft0;
                // f_str[field_index(i, j, 1)] = ft1;
                // f_str[field_index(i, j, 2)] = ft2;
                // f_str[field_index(i, j, 3)] = ft3;
                // f_str[field_index(i, j, 4)] = ft4;
                // f_str[field_index(i, j, 5)] = ft5;
                // f_str[field_index(i, j, 6)] = ft6;
                // f_str[field_index(i, j, 7)] = ft7;
                // f_str[field_index(i, j, 8)] = ft8;

                //----------------------- computing macro valiables-----------------
                double rho_ = ft0                     //
                              + ft1 + ft2 + ft3 + ft4 //
                              + ft5 + ft6 + ft7 + ft8;

                double rhoinv = 1.0 / rho_;
                double ux_ = rhoinv * (ft1 + ft5 + ft8 - (ft3 + ft6 + ft7)) * c;
                double uy_ = rhoinv * (ft2 + ft5 + ft6 - (ft4 + ft7 + ft8)) * c;

                rho[s_idx] = rho_;

                double force_x_ = force_x[s_idx];
                double force_y_ = force_y[s_idx];

                // ux[s_idx] = ux_ + force_x_ * 0.5 * rhoinv;
                // uy[s_idx] = uy_ + force_y_ * 0.5 * rhoinv;
                ux[s_idx] = ux_ + force_x_ * dt * 0.5 * rhoinv;
                uy[s_idx] = uy_ + force_y_ * dt * 0.5 * rhoinv;
            }
            else
            {
                rho[s_idx] = rho0;
                ux[s_idx] = 0.0;
                uy[s_idx] = 0.0;
            }
        }
    }
}

void NSE::BC_curve_Guo(const double *delta, const int *flag,
                       const double *x, const double *y)
{
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t s_idx = scalar_index(i, j);
            double rho_ = rho[s_idx];
            double ux_ = ux[s_idx];
            double uy_ = uy[s_idx];

            if (flag[s_idx] == COMP)
            {
                for (int k = 0; k < Q9; k++)
                {
                    int ip = i - dirx[k];
                    int jp = j - diry[k];
                    size_t f_idx = field_index(i, j, k);

                    if (flag[scalar_index(ip, jp)] == OBSTACLE)
                    {
                        int xff = i + dirx[k];
                        int yff = j + diry[k];
                        //***********fluid
                        double rho_b = rho[s_idx];
                        size_t sff_idx = scalar_index(xff, yff);

                        double ux_ghost, uy_ghost, fneq_ghost;

                        const double ux_b = 0.0;
                        const double uy_b = 0.0;
                        const double rho_ff = rho[sff_idx];
                        const double ux_ff = ux[sff_idx];
                        const double uy_ff = uy[sff_idx];
                        const double q_delta = delta[f_idx];
                        double ux_b1 = (ux_b + (q_delta - 1.0) * ux_) / q_delta;
                        double uy_b1 = (uy_b + (q_delta - 1.0) * uy_) / q_delta;
                        double ux_b2 = (2.0 * ux_b + (q_delta - 1.0) * ux_ff) / (1.0 + q_delta);
                        double uy_b2 = (2.0 * uy_b + (q_delta - 1.0) * uy_ff) / (1.0 + q_delta);

                        double fneq1 = f_col[f_idx] - feqNSE2019(k, rho_, ux_, uy_);
                        double fneq2 = f_col[field_index(xff, yff, k)] - feqNSE2019(k, rho_ff, ux_ff, uy_ff);

                        if (q_delta >= 0.75)
                        {
                            ux_ghost = ux_b1;
                            uy_ghost = uy_b1;
                            fneq_ghost = fneq1;
                        }
                        else
                        {
                            ux_ghost = q_delta * ux_b1 + (1.0 - q_delta) * ux_b2;
                            uy_ghost = q_delta * uy_b1 + (1.0 - q_delta) * uy_b2;
                            // ub[1] = Delta[x][y][k] * ub1[1] + (1.0 - Delta[x][y][k]) * ub2[1];
                            fneq_ghost = q_delta * fneq1 + (1.0 - q_delta) * fneq2;
                        }

                        f_col[field_index(ip, jp, k)] = feqNSE2019(k, rho_b, ux_ghost, uy_ghost) + (1.0 - 1.0 / tau) * fneq_ghost;
                        // f_col[field_index(ip, jp, k)] = feq(k, rhob, ub) + (1.0 - 1.0 / TAUF[x][y]) * fneq;
                        // ub1[0] = (0.0 + (Delta[x][y][k] - 1.0) * u[x][y][0]) / Delta[x][y][k];
                        // ub1[1] = (0.0 + (Delta[x][y][k] - 1.0) * u[x][y][1]) / Delta[x][y][k];
                        // ub2[0] = (2.0 * 0.0 + (Delta[x][y][k] - 1.0) * u[xff][yff][0]) / (1.0 + Delta[x][y][k]);
                        // ub2[1] = (2.0 * 0.0 + (Delta[x][y][k] - 1.0) * u[xff][yff][1]) / (1.0 + Delta[x][y][k]);
                        // fneq1 = f[x][y][k] - feq(k, rho[x][y], u[x][y]);
                        // fneq2 = f[xff][yff][k] - feq(k, rho[xff][yff], u[xff][yff]);
                    }
                }
            }
        }
    }
}

void NSE::initFeqAndMeq()
{
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t s_idx = scalar_index(i, j);
            double rho_ = rho[s_idx];
            double ux_ = ux[s_idx];
            double uy_ = uy[s_idx];

            for (int k = 0; k < Q9; ++k)
            {
                size_t f_idx = field_index(i, j, k);
                f_str[f_idx] = feqNSE2019(k, rho_, ux_, uy_);
                f_col[f_idx] = meqNSE2019(k, rho_, ux_, uy_);
            }
        }
    }
}

void NSE::streaming()
{
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            unsigned int ip1 = (i + 1) % (NX);
            unsigned int jp1 = (j + 1) % (NY);
            unsigned int im1 = (NX + i - 1) % (NX);
            unsigned int jm1 = (NY + j - 1) % (NY);
            // unsigned int ip1 = (i + 1) % (NX - 1);
            // unsigned int jp1 = (j + 1) % (NY - 1);
            // unsigned int im1 = (NX - 1 + i - 1) % (NX - 1);
            // unsigned int jm1 = (NY - 1 + j - 1) % (NY - 1);

            double ft0 = f_col[field_index(i, j, 0)];
            double ft1 = f_col[field_index(im1, j, 1)];
            double ft2 = f_col[field_index(i, jm1, 2)];
            double ft3 = f_col[field_index(ip1, j, 3)];
            double ft4 = f_col[field_index(i, jp1, 4)];
            double ft5 = f_col[field_index(im1, jm1, 5)];
            double ft6 = f_col[field_index(ip1, jm1, 6)];
            double ft7 = f_col[field_index(ip1, jp1, 7)];
            double ft8 = f_col[field_index(im1, jp1, 8)];

            f_str[field_index(i, j, 0)] = ft0;
            f_str[field_index(i, j, 1)] = ft1;
            f_str[field_index(i, j, 2)] = ft2;
            f_str[field_index(i, j, 3)] = ft3;
            f_str[field_index(i, j, 4)] = ft4;
            f_str[field_index(i, j, 5)] = ft5;
            f_str[field_index(i, j, 6)] = ft6;
            f_str[field_index(i, j, 7)] = ft7;
            f_str[field_index(i, j, 8)] = ft8;

            //-------------------------------- computing macro valiables
            double rho_ = ft0                     //
                          + ft1 + ft2 + ft3 + ft4 //
                          + ft5 + ft6 + ft7 + ft8;

            double rhoinv = 1.0 / rho_;
            double ux_ = rhoinv * (ft1 + ft5 + ft8 - (ft3 + ft6 + ft7)) * c;
            double uy_ = rhoinv * (ft2 + ft5 + ft6 - (ft4 + ft7 + ft8)) * c;

            size_t s_idx = scalar_index(i, j);
            rho[s_idx] = rho_;

            double force_x_ = force_x[s_idx];
            double force_y_ = force_y[s_idx];

            ux[s_idx] = ux_ + force_x_ * 0.5 * rhoinv;
            uy[s_idx] = uy_ + force_y_ * 0.5 * rhoinv;
            // ux[s_idx] = ux_ + force_x_ * dt * 0.5 * rhoinv;
            // uy[s_idx] = uy_ + force_y_ * dt * 0.5 * rhoinv;
        }
    }
}

void NSE::getMR(double *m, double *R)
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

void NSE::get_MF_Guo(double *MF, const double ux, const double uy,
                     const double forceX, const double forceY)
{
     MF[0] = 0;
     MF[1] = 6.0 * (forceX * ux + forceY * uy);
     MF[2] = -6.0 * (forceX * ux + forceY * uy);
     MF[3] = 1.0 * forceX;
     MF[4] = -1.0 * forceX;
     MF[5] = 1.0 * forceY;
     MF[6] = -1.0 * forceY;
     MF[7] = 2.0 * forceX * ux - 2.0 * forceY * uy;
     MF[8] = 1.0 * forceX * uy + 1.0 * forceY * ux;
}

void NSE::collisionMRT()
{
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t index = scalar_index(i, j);

            double rho_ = rho[index];
            double ux_ = ux[index];
            double uy_ = uy[index];
            double force_x_ = force_x[index];
            double force_y_ = force_y[index];

            double f_str_[9];
            double m_str_[9];
            for (int k = 0; k < Q9; ++k)
            {
                f_str_[k] = f_str[field_index(i, j, k)];
            }
            getMR(m_str_, f_str_);

            double MF_force[9];
            get_MF_Guo(MF_force, ux_, uy_, force_x_, force_y_);

            size_t fidx = field_index(i, j, 0);
            // f_col[fidx + 0] = m_str_[0] - s[0] * (m_str_[0] - meqNSE2019(0, rho_, ux_, uy_));
            // f_col[fidx + 1] = m_str_[1] - s[1] * (m_str_[1] - meqNSE2019(1, rho_, ux_, uy_));
            // f_col[fidx + 2] = m_str_[2] - s[2] * (m_str_[2] - meqNSE2019(2, rho_, ux_, uy_));
            // f_col[fidx + 3] = m_str_[3] - s[3] * (m_str_[3] - meqNSE2019(3, rho_, ux_, uy_)) + force_x_ * dt / c;
            // f_col[fidx + 4] = m_str_[4] - s[4] * (m_str_[4] - meqNSE2019(4, rho_, ux_, uy_));
            // f_col[fidx + 5] = m_str_[5] - s[5] * (m_str_[5] - meqNSE2019(5, rho_, ux_, uy_)) + force_y_ * dt / c;
            // f_col[fidx + 6] = m_str_[6] - s[6] * (m_str_[6] - meqNSE2019(6, rho_, ux_, uy_));
            // f_col[fidx + 7] = m_str_[7] - s[7] * (m_str_[7] - meqNSE2019(7, rho_, ux_, uy_));
            // f_col[fidx + 8] = m_str_[8] - s[8] * (m_str_[8] - meqNSE2019(8, rho_, ux_, uy_));

               //-----------------------------------------------------------------------------------------
                for (int k = 0; k < Q9; ++k)
                {
                    size_t fidx = field_index(i, j, k);
                    f_col[fidx] = m_str_[k] - s[k] * (m_str_[k] - meqNSE2019(k, rho_, ux_, uy_)) //
                                //   + MF_force[k] * (2.0 - s[k]) * 0.5;
                      + dt * MF_force[k] * (2.0 - s[k]) * 0.5;
                }

            // double MF0 = 0;
            // double MF1 = 6.0 * (force_x_ * ux_ + force_y_ * uy_);
            // double MF2 = -6.0 * (force_x_ * ux_ + force_y_ * uy_);
            // double MF3 = 1.0 * force_x_;
            // double MF4 = -1.0 * force_x_;
            // double MF5 = 1.0 * force_y_;
            // double MF6 = -1.0 * force_y_;
            // double MF7 = 2.0 * force_x_ * ux_ - 2.0 * force_y_ * uy_;
            // double MF8 = 1.0 * force_x_ * uy_ + 1.0 * force_y_ * ux_;

            // f_col[fidx + 0] = m_str_[0] - s[0] * (m_str_[0] - meqNSE2019(0, rho_, ux_, uy_)) + dt * MF0 * (2.0 - s[0]) * 0.5;
            // f_col[fidx + 1] = m_str_[1] - s[1] * (m_str_[1] - meqNSE2019(1, rho_, ux_, uy_)) + dt * MF1 * (2.0 - s[1]) * 0.5;
            // f_col[fidx + 2] = m_str_[2] - s[2] * (m_str_[2] - meqNSE2019(2, rho_, ux_, uy_)) + dt * MF2 * (2.0 - s[2]) * 0.5;
            // f_col[fidx + 3] = m_str_[3] - s[3] * (m_str_[3] - meqNSE2019(3, rho_, ux_, uy_)) + dt * MF3 * (2.0 - s[3]) * 0.5;
            // f_col[fidx + 4] = m_str_[4] - s[4] * (m_str_[4] - meqNSE2019(4, rho_, ux_, uy_)) + dt * MF4 * (2.0 - s[4]) * 0.5;
            // f_col[fidx + 5] = m_str_[5] - s[5] * (m_str_[5] - meqNSE2019(5, rho_, ux_, uy_)) + dt * MF5 * (2.0 - s[5]) * 0.5;
            // f_col[fidx + 6] = m_str_[6] - s[6] * (m_str_[6] - meqNSE2019(6, rho_, ux_, uy_)) + dt * MF6 * (2.0 - s[6]) * 0.5;
            // f_col[fidx + 7] = m_str_[7] - s[7] * (m_str_[7] - meqNSE2019(7, rho_, ux_, uy_)) + dt * MF7 * (2.0 - s[7]) * 0.5;
            // f_col[fidx + 8] = m_str_[8] - s[8] * (m_str_[8] - meqNSE2019(8, rho_, ux_, uy_)) + dt * MF8 * (2.0 - s[8]) * 0.5;
        }
    }
}

void NSE::m2f()
{
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            double m0 = f_col[field_index(i, j, 0)];
            double m1 = f_col[field_index(i, j, 1)];
            double m2 = f_col[field_index(i, j, 2)];
            double m3 = f_col[field_index(i, j, 3)];
            double m4 = f_col[field_index(i, j, 4)];
            double m5 = f_col[field_index(i, j, 5)];
            double m6 = f_col[field_index(i, j, 6)];
            double m7 = f_col[field_index(i, j, 7)];
            double m8 = f_col[field_index(i, j, 8)];

            f_col[field_index(i, j, 0)] = m0 / 9 - m1 / 9 + m2 / 9;
            f_col[field_index(i, j, 1)] = m0 / 9 - m1 / 36 - m2 / 18 + m3 / 6 - m4 / 6 + m7 / 4;
            f_col[field_index(i, j, 2)] = m0 / 9 - m1 / 36 - m2 / 18 + m5 / 6 - m6 / 6 - m7 / 4;
            f_col[field_index(i, j, 3)] = m0 / 9 - m1 / 36 - m2 / 18 - m3 / 6 + m4 / 6 + m7 / 4;
            f_col[field_index(i, j, 4)] = m0 / 9 - m1 / 36 - m2 / 18 - m5 / 6 + m6 / 6 - m7 / 4;
            f_col[field_index(i, j, 5)] = m0 / 9 + m1 / 18 + m2 / 36 + m3 / 6 + m4 / 12 + m5 / 6 + m6 / 12 + m8 / 4;
            f_col[field_index(i, j, 6)] = m0 / 9 + m1 / 18 + m2 / 36 - m3 / 6 - m4 / 12 + m5 / 6 + m6 / 12 - m8 / 4;
            f_col[field_index(i, j, 7)] = m0 / 9 + m1 / 18 + m2 / 36 - m3 / 6 - m4 / 12 - m5 / 6 - m6 / 12 + m8 / 4;
            f_col[field_index(i, j, 8)] = m0 / 9 + m1 / 18 + m2 / 36 + m3 / 6 + m4 / 12 - m5 / 6 - m6 / 12 - m8 / 4;
        }
    }
}

// BADE get_u(const size_t i, const size_t j)
// {
//     size_t s_idx = scalar_
//     BADE u{ux,uy}
// }

double NSE::get_ux(const size_t i, const size_t j) const
{
    return ux[scalar_index(i, j)];
};
double NSE::get_uy(const size_t i, const size_t j) const
{
    return uy[scalar_index(i, j)];
};
double NSE::get_rho(const size_t i, const size_t j) const
{
    return rho[scalar_index(i, j)];
};

void get_NSE_force(const CDE2019 &cde, NSE &nse, const double alpha_g)
{
    const size_t NX_ = nse.NX;
    const size_t NY_ = nse.NY;
    for (size_t j = 0; j < NY_; ++j)
    {
        for (size_t i = 0; i < NX_; ++i)
        {
            size_t s_idx = nse.scalar_index(i, j);
            // nse.force_x[s_idx] = alpha_g;
            // nse.force_y[s_idx] = 0;
            nse.force_x[s_idx] = 0.0;
            nse.force_y[s_idx] = rho0 * alpha_g * (cde.get_C(i, j) - 0.5);
        }
    }
}

void NSE::outputTec(double t, double *x, double *y)
{
    std::ostringstream name;
    name << "t_" << t << "_"
         << ".dat";
    std::ofstream out(name.str().c_str());
    out << "Title= \"Lid Driven Flow\"\n";
    out << "VARIABLES = \"X\", \"Y\", \"rho\", \"ux\", \"uy\", \"fx\", \"fy\" \n";
    out << "ZONE T= \"BOX\",I=" << NX << ",J=" << NY << ",F=	POINT" << std::endl;

    for (size_t j = 0; j < NY; ++j)
    {
        for (size_t i = 0; i < NX; ++i)
        {
            const size_t index = scalar_index(i, j);
            const size_t findex = field_index(i, j, 0);
            out << std::fixed << std::setprecision(10)
                << " " << x[index] << " " << y[index]
                << " " << rho[index]

                << " " << ux[index]
                << " " << uy[index]
                << " " << force_x[index]
                << " " << force_y[index]
                << std::endl;
        }
    }
    out.close();
}