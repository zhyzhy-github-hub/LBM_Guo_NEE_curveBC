#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#define PI 3.14159265359

#define w0 4.0 / 9
#define wc 1.0 / 9
#define wd 1.0 / 36

#define delta0 0.0001
#define C0 0.01

#define rho0 1.0

const int Q9 = 9;

const double w[9] = {w0, wc, wc, wc, wc, wd, wd, wd, wd};

const int dirx[] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
const int diry[] = {0, 0, 1, 0, -1, 1, 1, -1, -1};

const double xCen = 0.5;
const double yCen = 0.5;
const double rIn = 0.2;

// const double xCen1 = 0.2;
// const double yCen1 = 0.2;
// const double rIn_1 = 0.2;

// #define cirIn(x, y) ((x - xCen) * (x - xCen) + (y - yCen) * (y - yCen) - (rIn) * (rIn))

enum
{
    OBSTACLE = -2,
    SILID = -1,
    COMP,
    RIGHT,
    TOP,
    LEFT,
    BOTTOM,
    RIGHT_TOP,
    LEFT_TOP,
    LEFT_BOTTOM,
    RIGHT_BOTTOM
};

struct BADE
{
    double Bx;
    double By;
};

class Mesh
{
public:
    double *xNorm, *yNorm;
    int *solid_flag;
    double *delta_Q9;

    Mesh();
    Mesh(int NX_, int NY_);
    ~Mesh();
    void InitMesh();
    void mesh_flag();
    void mesh_Tec();
    void mesh_delta_Tec();
    void init_delta_Q9();

    size_t get_scalar_index(const size_t i, const size_t j) const { return scalar_index(i, j); };

    size_t get_NX() const { return NX; };
    size_t get_NY() const { return NY; };

private:
    size_t NX;
    size_t NY;

    double *x, *y;

    size_t scalar_index(int i, int j) const { return (NX * j + i); }
    size_t field_index(const size_t x, const size_t y, const size_t d) const { return (Q9 * (NX * y + x) + d); }

    double getDelta(const double tr, const double x0, const double y0,
                    const double x1, const double y1, const double x2, const double y2);
};

template <typename T>
inline T p2(const T x)
{
    return (x * x);
}

inline int circle_obstruct(const double R, const double x, const double y, const double xCen, const double yCen)
{
    //  double len = pow((x - xCen), 2) + pow((y - yCen), 2);
    double len = p2((x - xCen)) + p2((y - yCen));
    int result = (len > R * R ? 1 : 0);
    return result;
    // return (len > R * R);
}