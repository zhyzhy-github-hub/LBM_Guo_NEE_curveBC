#include "Mesh.h"
Mesh::~Mesh()
{
    delete[] x;
    delete[] y;
    delete[] solid_flag;
    delete[] xNorm;
    delete[] yNorm;
    delete[] delta_Q9;
}

Mesh::Mesh(int NX_, int NY_)
{
    NX = NX_;
    NY = NY_;
    const unsigned int NN = NX * NY;
    x = new double[NN];
    y = new double[NN];
    xNorm = new double[NN];
    yNorm = new double[NN];
    solid_flag = new int[NN];
    delta_Q9 = new double[NN * Q9];
}

void Mesh::InitMesh()
{
    const double LC = NY - 1;
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t index = scalar_index(i, j);
            x[index] = i;
            y[index] = j;
            solid_flag[index] = COMP;
            xNorm[index] = double(i) / LC;
            yNorm[index] = double(j) / LC;
        }
    }
}

// size_t Mesh::field_index(const size_t x, const size_t y, const size_t d)
// {
//     // return (NX * (NY * (d) + y) + x);
//     // return (9 * (NX * y + x) + d);
//     return (Q9 * (NX * y + x) + d);
// }

void Mesh::mesh_flag()
{
    for (size_t y = 0; y < NY; ++y)
    {
        for (size_t x = 0; x < NX; ++x)
        {
            size_t index = scalar_index(x, y);
            if (circle_obstruct(rIn, xNorm[index], yNorm[index], xCen, yCen) == 0)
            {
                solid_flag[index] = OBSTACLE;
            }
            // else if (circle_obstruct(rIn_1, xNorm[index], yNorm[index], xCen1, yCen1) == 0)
            // {
            //     solid_flag[index] = OBSTACLE;
            // }
            else if (y == 0 || y == NY - 1 || x == 0 || x == NX - 1)
            {
                // size_t x_in, y_in;

                if (y == 0 && x != 0 && x != NX - 1)
                {
                    // --------------------------------------lower plate
                    solid_flag[index] = BOTTOM;
                }
                else if (x == 0 && y == 0)
                {
                    // ***********************************left bottom
                    solid_flag[index] = LEFT_BOTTOM;
                }
                else if (x == NX - 1 && y == 0)
                {
                    // ******************************right bottom
                    solid_flag[index] = RIGHT_BOTTOM;
                }
                else if (x == 0 && y == NY - 1)
                {
                    //********************************** left top
                    solid_flag[index] = LEFT_TOP;
                }
                else if (x == NX - 1 && y == NY - 1)
                {
                    // *********************************right top
                    solid_flag[index] = RIGHT_TOP;
                }
                else if (x == 0 && y != 0 && y != NY - 1)
                {
                    // -----------------------------------left plate
                    solid_flag[index] = LEFT;
                }
                else if (x == NX - 1 && y != 0 && y != NY - 1)
                {
                    // ---------------------------------right plate
                    solid_flag[index] = RIGHT;
                }
                else if (y == NY - 1 && x != 0 && x != NX - 1)
                {
                    // ---------------------------------upper plate
                    solid_flag[index] = TOP;
                }
            }
            else
            {
                solid_flag[index] = COMP;
            }
        }
    }
}
//************************************************
double Mesh::getDelta(const double tr, const double x0, const double y0,
                      const double x1, const double y1, const double x2, const double y2)
{
    double a = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
    double b = 2.0 * (x1 - x2) * (x1 - x0) + 2.0 * (y1 - y2) * (y1 - y0);
    double c = (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) - tr * tr;

    double t1 = (-b + sqrt(b * b - 4.0 * a * c)) / (2.0 * a);
    double t2 = (-b - sqrt(b * b - 4.0 * a * c)) / (2.0 * a);

    double xx1 = x1 + (x1 - x2) * t1;
    double yy1 = y1 + (y1 - y2) * t1;

    double xx2 = x1 + (x1 - x2) * t2;
    double yy2 = y1 + (y1 - y2) * t2;

    double q1 = sqrt((x1 - xx1) * (x1 - xx1) + (y1 - yy1) * (y1 - yy1)) / sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
    double q2 = sqrt((x1 - xx2) * (x1 - xx2) + (y1 - yy2) * (y1 - yy2)) / sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
    if (q1 <= q2)
    {
        return q1;
    }
    else
    {
        return q2;
    }
}

void Mesh::init_delta_Q9()
{
    int L = NY - 1;
    for (int x = 0; x < NX; ++x)
    {
        for (int y = 0; y < NY; ++y)
        {
            size_t index = scalar_index(x, y);

            if (solid_flag[index] == COMP)
            {
                for (int k = 0; k < Q9; ++k)
                {
                    int xp = x - dirx[k];
                    int yp = y - diry[k];
                    size_t index_p = scalar_index(xp, yp);

                    if (solid_flag[index_p] == OBSTACLE)
                    {
                        delta_Q9[field_index(x, y, k)] = getDelta(rIn * L, xCen * L, yCen * L, x, y, xp, yp);
                    }
                }
            }
        }
    }
}

void Mesh::mesh_delta_Tec()
{
    std::ofstream out("Mesh_MSG_delta.dat");
    out << "Title= \"delta\"\n";
    out << "VARIABLES = \"X\", \"Y\", \"D0\",\"D1\",\"D2\",\"D3\",\"D4\",\"D5\",\"D6\",\"D7\",\"D8\" \n";
    out << "ZONE T= \"BOX\",I=" << NX << ",J=" << NY << ",F=	POINT" << std::endl;

    for (size_t j = 0; j < NY; ++j)
    {
        for (size_t i = 0; i < NX; ++i)
        {
            const size_t index = scalar_index(i, j);
            out << std::fixed << std::setprecision(10)
                << " " << xNorm[index]
                << " " << yNorm[index]
                << " " << delta_Q9[field_index(i, j, 0)]
                << " " << delta_Q9[field_index(i, j, 1)]
                << " " << delta_Q9[field_index(i, j, 2)]
                << " " << delta_Q9[field_index(i, j, 3)]
                << " " << delta_Q9[field_index(i, j, 4)]
                << " " << delta_Q9[field_index(i, j, 5)]
                << " " << delta_Q9[field_index(i, j, 6)]
                << " " << delta_Q9[field_index(i, j, 7)]
                << " " << delta_Q9[field_index(i, j, 8)]
                << std::endl;
        }
    }
}

void Mesh::mesh_Tec()
{
    std::ofstream out("Mesh_MSG.dat");
    out << "Title= \"Lid Driven Flow\"\n";
    out << "VARIABLES = \"X\", \"Y\", \"xNorm\",  \"yNorm\", \"Solid_Flag\" \n";
    out << "ZONE T= \"BOX\",I=" << NX << ",J=" << NY << ",F=	POINT" << std::endl;

    for (size_t j = 0; j < NY; ++j)
    {
        for (size_t i = 0; i < NX; ++i)
        {
            const size_t index = scalar_index(i, j);
            out << std::fixed << std::setprecision(10)
                << " " << x[index] << " " << y[index]
                // << " " << f0[findex]
                // << " " << f1[findex]
                << " " << xNorm[index]
                << " " << yNorm[index]
                << " " << solid_flag[index]
                << std::endl;
        }
    }
}