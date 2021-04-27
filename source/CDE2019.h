#pragma once

#include <cmath>
#include <iostream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "Mesh.h"

#define DIFF_SCALING

class NSE;

class CDE2019
{

public:
    CDE2019(int NX_, int NY_);
    CDE2019(const int NX_, const int NY_, const double mu_, const double D_, const double c_, const double s_add);

    ~CDE2019();

    double get_time_CDE(const size_t n) const {return n * dt * D / pow(1, 2);};

    void get_each_phi_D_B(const double C_, const BADE B_, double& phi, double& D, BADE & B);

    friend void get_B(const NSE &nse, CDE2019 &cde);

#ifdef DIFF_SCALING
    void initCDE2019_Diff(double mu, double D_, double s_add);
#elif
    void initCDE2019(double c_, double D_, double s_add);
#endif

    void getCOld();

    void initScalar(size_t n, double *x, double *y);

    void init_scalar(const size_t n, const double *x, const double *y, const int *flag);

    void initFeqAndMeq();
    void initFeqAndMeq(const int *flag);

    void AnsScalar(size_t n, double *x, double *y);

    void outputTec(double t, double *x, double *y);

    void getF(size_t t_, double *x, double *y);

    double get_C(const size_t i, const size_t j) const;

    void collision();
    void streaming();
    // void streaming_cruve(const int *flag);

    void BC_NEE_Dirichlet();

    void streaming_cruve(const int *flag);

    void BC_curve_NEE(const double *delta, const int *flag,
                      const double *x, const double *y);

    void BC_curve_Guo(const double *delta, const int *flag,
                      const double *x, const double *y);

    void streamingBGK();
    void collisionBGK();

    void m2f();
    double errC();

    double time(size_t n) { return n * dt; }

private:
    int NX;
    int NY;

    double dx;

#ifdef DIFF_SCALING
    double dt;
    double D, tau;
    double mu;
    double c;
    double Cs, Cs2;
#elif
    double dt, c, Cs, Cs2;
    double D;
    double tau;

#endif

    double *s;

    double *C;
    double *C_old;
    double *F;
    double *f1; // 1 is post collision
    double *f0; // 0 is post streaming

    BADE *B_ade;

    size_t scalar_index(size_t i, size_t j) const;
    size_t field_index(size_t x, size_t y, size_t d) const;

    void getMR(double *m, double *R);
    void getFi(double *Fi, double F);
    double meqADE2019(int k, double phi, double D, BADE B);
    double feqADE2019(int k, double phi, double D, BADE B);
};
