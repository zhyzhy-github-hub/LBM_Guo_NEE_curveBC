#include "source/CDE2019.h"
#include "source/Mesh.h"
#include "source/NSE.h"
#include "seconds.h"

using std::cout;
using std::endl;

void output_Tec_result(const double t, const Mesh &mesh0,
                       NSE &nse0, CDE2019 &cde);

int main()
{
    const int NX = 201;
    const int NY = 201;
    // const int NX = 161;
    // const int NY = 161;

    // double Re = 1000;
    // double c_nu = 1.0;

    // // double nu = u_max * (NY-1) / Re;
    // double nu = u_max * 1 / Re;

    double c_nu = 1.0;

    double s_add = 0.4;
    double mu = NY - 1;
    double c_kappa = 1.0; // c_D * 1.0 / dx_Drm

    double L = 1.0;
    double Ra = 1e5;
    double Pr = 0.71;
    // double Pr = 1;
    double Ma = 0.08;
    double nu = Ma * L / sqrt(3) * sqrt(Pr / Ra);

    // double nu = 0.001;
    // double dp = 0.001;
    // double u_poi = dp * L * L / 8 / nu;
    // cout << "u_poi = " << u_poi << endl;

    double kappa = nu / Pr;

    // double alpha_g = kappa * kappa * Pr * Ra / 1.0 / (L * L * L);
    double alpha_g = kappa * kappa * Pr * Ra / 1.0 / (L * L * L);

    cout << " nu = " << nu << ",  Kappa = " << kappa << ", alpla_g = " << alpha_g << endl;

    Mesh mesh0(NX, NY);
    mesh0.InitMesh();
    mesh0.mesh_flag();
    mesh0.mesh_Tec();
    mesh0.init_delta_Q9();
    mesh0.mesh_delta_Tec();

    NSE NSE0(NX, NY, nu, c_nu);

    NSE0.initScalar(0, mesh0.xNorm, mesh0.yNorm);
    NSE0.initFeqAndMeq();
    NSE0.outputTec(0, mesh0.xNorm, mesh0.yNorm);

    CDE2019 T0(NX, NY, mu, kappa, c_kappa, s_add);
    T0.init_scalar(0, mesh0.xNorm, mesh0.yNorm, mesh0.solid_flag);
    T0.initFeqAndMeq(mesh0.solid_flag);
    // T0.initFeqAndMeq();
    get_B(NSE0, T0);
    T0.outputTec(0, mesh0.xNorm, mesh0.yNorm);

    double start_time = seconds();

    size_t top_step;
    for (int n = 1; n < 100000; ++n)
    {
        //--------------NSE--------------------
        // get_NSE_force(T0, NSE0, dp);
        get_NSE_force(T0, NSE0, alpha_g);

        NSE0.m2f();
        NSE0.BC_curve_Guo(mesh0.delta_Q9, mesh0.solid_flag, mesh0.xNorm, mesh0.yNorm);
        NSE0.streaming_curve(mesh0.solid_flag);
        // const double *temp_ux = get_ux(NSE0);
        // cout << "-------------\n";
        // cout << temp_ux << endl;
        // NSE0.streaming();

        NSE0.BC_NEE_Guo();
        NSE0.collisionMRT();

        // //--------------CDE--------------------
        get_B(NSE0, T0);
        T0.m2f();
        T0.BC_curve_Guo(mesh0.delta_Q9, mesh0.solid_flag, mesh0.xNorm, mesh0.yNorm);
        T0.streaming_cruve(mesh0.solid_flag);

        // T0.streaming();
        T0.BC_NEE_Dirichlet();
        T0.collision();

        top_step = n;

        if (n % 100 == 0)
        {
            cout << "Step = " << n << ", time = " << T0.get_time_CDE(n) << endl;
        }
        if (n % 2000 == 0 || n == 500)
        {
            cout << "Out Tec file-------------" << endl;
            // NSE0.outputTec(top_step, mesh0.xNorm, mesh0.yNorm);

            // T0.outputTec(n, mesh0.xNorm, mesh0.yNorm);

            output_Tec_result(n, mesh0, NSE0, T0);
            // break;
        }
    }
    double end_time = seconds();
    cout << "The all time is : ----------------" << end_time - start_time << endl;
    NSE0.outputTec(top_step, mesh0.xNorm, mesh0.yNorm);
}

void output_Tec_result(const double t, const Mesh &mesh0,
                       NSE &nse0, CDE2019 &cde)
{
    const size_t NX = mesh0.get_NX();
    const size_t NY = mesh0.get_NY();

    std::ostringstream name;
    name << "Natural_convection_" << t << ".dat";
    std::ofstream out(name.str().c_str());
    out << "Title= \"NC cylinder Flow\"\n";
    out << "VARIABLES = \"X\", \"Y\", \"rho\", \"ux\", \"uy\", \"C\"  \n";
    out << "ZONE T= \"BOX\",I=" << NX << ",J=" << NY << ",F=	POINT" << std::endl;

    for (size_t j = 0; j < NY; ++j)
    {
        for (size_t i = 0; i < NX; ++i)
        {
            const size_t index = mesh0.get_scalar_index(i, j);
            out << std::fixed << std::setprecision(16)
                << " " << mesh0.xNorm[index] << " " << mesh0.yNorm[index]
                << " " << nse0.get_rho(i, j)
                << " " << nse0.get_ux(i, j)
                << " " << nse0.get_uy(i, j)
                << " " << cde.get_C(i, j)
                << std::endl;
        }
    }
    out.close();
}
