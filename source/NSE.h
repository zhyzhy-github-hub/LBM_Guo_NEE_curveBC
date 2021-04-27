#pragma once
#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "Mesh.h"

const double u_max = 0.0;

class CDE2019;

class NSE
{

public:
    NSE(const int NX_, const int NY_, const double nu_, const double c_);
    ~NSE();

    void initScalar(const size_t n, const double *x, const double *y);
    void initFeqAndMeq();
    void streaming();
    void streaming_curve(const int *flag);

    void m2f();

    void BC_NEE_Guo();

    void BC_curve_Guo(const double *delta, const int *flag,
                      const double *x, const double *y);

    void collisionMRT();

    void outputTec(double t, double *x, double *y);

    void get_force();

    double get_ux(const size_t i, const size_t j) const;
    double get_uy(const size_t i, const size_t j) const;
    double get_rho(const size_t i, const size_t j) const;

    friend void get_NSE_force(const CDE2019 &cde, NSE &nse, const double alpha_g);

private:
    // const int Q9 = 9;

    int NX;
    int NY;

    double dx;

    double dt, c, Cs, Cs2;
    double nu;
    double tau;

    double *s;

    double *ux, *uy, *rho;
    double *ux_old, *uy_old;
    double *force_x, *force_y;

    double *f_col; // 1 is post collision
    double *f_str; // 0 is post streaming

    size_t scalar_index(const size_t i, const size_t j) const;
    size_t field_index(const size_t x, const size_t y, const size_t d) const;

    double meqNSE2019(int k, double rho, double ux, double uy);
    double feqNSE2019(int k, double rho, double ux, double uy);

    void getMR(double *m, double *R);

    void get_MF_Guo(double *MF, const double ux, const double uy,
                    const double forceX, const double forceY);
};