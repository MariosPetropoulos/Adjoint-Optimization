#include <cmath>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm>


static const double U   = 16.0/3.0;
double Q = -0.05 * U;
static const double mi  = 1.75e-5;
static const double rho = 1.25;
static const double nu  = mi / rho;
static const double eps = 1e-6;
static const double A0 = M_PI*nu/2;
static const double C = M_PI*mi*U/2;


static inline double v0(double x, double a, double b, double c)
{
    if (x <= 0.2 || x >= 0.8) return 0.0;
    double s = std::sin(M_PI * x / 6.0);
    double g = a*x*x + b*x + c*s;
    return -(g*g);
}


static inline double ddelta_dx(double x, double delta,
                               double a, double b, double c)
{
    const double K = (2.0 / M_PI - 0.5);

    delta = std::max(delta, 1e-14);
    double v = v0(x, a, b, c);

    return ( (0.5 * M_PI) * nu - v * delta ) / (K * U * delta);
}





double computeF(double a, double b, double c,
                int N,
                bool writeFiles = true,
                const std::string& deltaFile = "delta.txt",
                const std::string& fricFile  = "friction.txt")
{
    const double x0 = 0.0, x1 = 1.0;
    const double dx = (x1 - x0) / (N - 1);

    std::vector<double> x(N), delta(N), tauw(N);
    for (int i = 0; i < N; ++i) x[i] = x0 + i * dx;

    delta[0] = eps;

    for (int i = 0; i < N - 1; ++i)
    {
        double xi = x[i];

        double k1 = ddelta_dx(xi,          delta[i],              a, b, c);
        double k2 = ddelta_dx(xi + dx/2.0, delta[i] + dx*k1/2.0,  a, b, c);

        delta[i+1] = delta[i] + dx * k2;
        if (delta[i+1] < eps) delta[i+1] = eps;
    }

    for (int i = 0; i < N; ++i) {
        delta[i] = std::max(delta[i], 1e-14);
        tauw[i]  = (M_PI/2.0) * mi * U / delta[i];
    }

    double F = 0.0;
    for (int i = 0; i < N - 1; ++i) {
        F += 0.5 * (tauw[i] + tauw[i+1]) * dx;
    }

    if (writeFiles) {
        std::ofstream fdelta(deltaFile);
        fdelta << std::scientific << std::setprecision(17);
        for (int i = 0; i < N; ++i) fdelta << x[i] << " " << delta[i] << "\n";

        std::ofstream ffric(fricFile);
        ffric << std::scientific << std::setprecision(17);
        for (int i = 0; i < N; ++i) ffric << x[i] << " " << tauw[i] << "\n";
    }

    return F;
}




static inline double dv0_da(double x, double a, double b, double c)
{
    if (x <= 0.2 || x >= 0.8) return 0.0;

    double s = std::sin(M_PI * x / 6.0);
    double g = a*x*x + b*x + c*s;

    return -2.0 * g * (x*x);
}

static inline double dv0_db(double x, double a, double b, double c)
{
    if (x <= 0.2 || x >= 0.8) return 0.0;

    double s = std::sin(M_PI * x / 6.0);
    double g = a*x*x + b*x + c*s;

    return -2.0 * g * x;
}

static inline double dv0_dc(double x, double a, double b, double c)
{
    if (x <= 0.2 || x >= 0.8) return 0.0;

    double s = std::sin(M_PI * x / 6.0);
    double g = a*x*x + b*x + c*s;

    return -2.0 * g * s;
}





static inline double rhs(double delta_theta, double delta, double v0_theta,
                         double A0, double K, double U)
{
    const double KU = K * U;
    const double eps = 1e-14;
    if (std::fabs(delta) < eps) delta = (delta >= 0.0 ? eps : -eps);

    return -(A0 / KU) * (delta_theta / (delta * delta)) - (1.0 / KU) * v0_theta;
}



void solve_sens_RK2(
    int N, double x0, double h,
    const double* x,           
    const double* delta,       
    double a, double b, double c,
    double A0, double K, double U,
    double da0, double db0, double dc0,  
    double* da, double* db, double* dc   
){
    da[0] = da0;
    db[0] = db0;
    dc[0] = dc0;

    for (int i = 0; i < N-1; ++i)
    {
        double xi   = x ? x[i] : (x0 + i*h);
        double xmid = xi + 0.5*h;

        double di   = delta[i];
        double dmid = 0.5 * (delta[i] + delta[i+1]);   

        // --- a ---
        double k1a = rhs(da[i], di,   dv0_da(xi,   a,b,c), A0,K,U);
        double k2a = rhs(da[i] + 0.5*h*k1a, dmid, dv0_da(xmid, a,b,c), A0,K,U);
        da[i+1] = da[i] + h*k2a;

        // --- b ---
        double k1b = rhs(db[i], di,   dv0_db(xi,   a,b,c), A0,K,U);
        double k2b = rhs(db[i] + 0.5*h*k1b, dmid, dv0_db(xmid, a,b,c), A0,K,U);
        db[i+1] = db[i] + h*k2b;

        // --- c ---
        double k1c = rhs(dc[i], di,   dv0_dc(xi,   a,b,c), A0,K,U);
        double k2c = rhs(dc[i] + 0.5*h*k1c, dmid, dv0_dc(xmid, a,b,c), A0,K,U);
        dc[i+1] = dc[i] + h*k2c;
    }
}


static void solve_delta_RK2(int N, double a, double b, double c,
                            std::vector<double>& x, std::vector<double>& delta)
{
    const double x0 = 0.0, x1 = 1.0;
    const double dx = (x1 - x0) / (N - 1);

    x.resize(N);
    delta.resize(N);

    for (int i = 0; i < N; ++i) x[i] = x0 + i * dx;

    delta[0] = eps;

    for (int i = 0; i < N - 1; ++i)
    {
        double xi = x[i];

        double k1 = ddelta_dx(xi,          delta[i],              a, b, c);
        double k2 = ddelta_dx(xi + dx/2.0, delta[i] + dx*k1/2.0,  a, b, c);

        delta[i+1] = delta[i] + dx * k2;
        if (delta[i+1] < eps) delta[i+1] = eps;
    }
}

void compute_df_ddelta(const std::vector<double>& delta,
                       double dx,
                       std::vector<double>& df_ddelta)
{
    int N = (int)delta.size();
    df_ddelta.resize(N);

    const double C = (M_PI/2.0) * mi * U;

    for (int i = 0; i < N; ++i)
    {
        double di = std::max(delta[i], 1e-14);

        double w = (i == 0 || i == N-1) ? 0.5 : 1.0;
        df_ddelta[i] = -C * dx * w / (di * di);
    }
}

double dot(const std::vector<double>& v1,
           const std::vector<double>& v2)
{
    double s = 0.0;
    int N = (int)v1.size();
    for (int i = 0; i < N; ++i)
        s += v1[i] * v2[i];
    return s;
}


// ADJOINT

static inline double adj_rhs(double x, double lam,
                             double delta, double A, double B, double C, double a, double b, double c)
{
    return ((v0(x,a,b,c))/(A*delta))*lam - (C/A)*(1.0/(delta*delta*delta));
}



void solve_adjoint_RK2_backward(
    int N,
    const double* x,           
    const double* delta,      
    double A, double B, double C,
    double a, double b, double c,
    double lamL,               // boundary condition at x=L
    double* lam               
)
{
    
    lam[N-1] = lamL;

    
    for (int i = N - 1; i > 0; --i)
    {
        double xi   = x[i];
        double him1 = x[i-1];
        double h    = him1 - xi;        

        double xmid = xi + 0.5 * h;

        
        double d_i   = delta[i];
        double d_mid = 0.5 * (delta[i] + delta[i-1]);

        
        double k1 = adj_rhs(xi,          lam[i],
                            d_i, A, B, C, a, b, c);

        double k2 = adj_rhs(xmid,        lam[i] + 0.5*h*k1,
                            d_mid, A, B, C, a, b, c);

        
        lam[i-1] = lam[i] + h * k2;
    }
}


double trapz_lam_delta_q(const std::vector<double>& x,
                         const std::vector<double>& lam,
                         const std::vector<double>& delta,
                         const std::vector<double>& q)
{
    int N = (int)x.size();
    double I = 0.0;

    for (int i = 0; i < N - 1; ++i)
    {
        double dx = x[i+1] - x[i];

        double di   = std::max(delta[i],   1e-14);
        double dip1 = std::max(delta[i+1], 1e-14);

        double gi   = lam[i]   * di   * q[i];
        double gip1 = lam[i+1] * dip1 * q[i+1];

        I += 0.5 * (gi + gip1) * dx;
    }

    return I;
}


// ALM

double compute_F_and_grad_CA(double a, double b, double c, int N,
                             double& dFda, double& dFdb, double& dFdc)
{
    
    std::vector<double> x, delta;
    solve_delta_RK2(N, a, b, c, x, delta);

    
    std::vector<double> lam(N);
    double lamL = 0.0;

    double A = (2/M_PI - 0.5)*U;
    double B = (M_PI*nu/2);

    solve_adjoint_RK2_backward(
        N, x.data(), delta.data(),
        A, B, C, a, b, c,
        lamL, lam.data()
    );

    
    std::vector<double> q_a(N), q_b(N), q_c(N);
    for (int i = 0; i < N; ++i) {
        q_a[i] = dv0_da(x[i], a, b, c);
        q_b[i] = dv0_db(x[i], a, b, c);
        q_c[i] = dv0_dc(x[i], a, b, c);
    }

    double Ia = trapz_lam_delta_q(x, lam, delta, q_a);
    double Ib = trapz_lam_delta_q(x, lam, delta, q_b);
    double Ic = trapz_lam_delta_q(x, lam, delta, q_c);

    
    dFda = Ia;
    dFdb = Ib;
    dFdc = Ic;

    
    return computeF(a, b, c, N, false);
}


double compute_c_and_grad(double a, double b, double c,
                          double& dcda, double& dcdb, double& dcdc)
{
    // Constraint:
    // 0.0444450027370472 c^2
    // + (0.104737806763403 a + 0.172815991859899 b) c
    // + (0.065472 a^2 + 0.204 a b + 0.168 b^2 - 0.05 U) = 0

    constexpr double K2   = 0.0444450027370472;   
    constexpr double Ka_c = 0.104737806763403;    
    constexpr double Kb_c = 0.172815991859899;    

    constexpr double Kaa  = 0.065472;             
    constexpr double Kab  = 0.204;               
    constexpr double Kbb  = 0.168;               

    // residual cval 
    const double cval =
        K2 * c * c
        + (Ka_c * a + Kb_c * b) * c
        + (Kaa * a * a + Kab * a * b + Kbb * b * b - 0.05 * U);

    
    dcda = Ka_c * c + (2.0 * Kaa * a + Kab * b);
    dcdb = Kb_c * c + (Kab * a + 2.0 * Kbb * b);
    dcdc = 2.0 * K2 * c + (Ka_c * a + Kb_c * b);

    return cval;
}








int main()
{
    int N =80001;
    double h = 1e-6;

    double a = 11.0;
    double b = -53.0;
    double c = 89.86;

    double F0 = computeF(a,b,c,N,true);


    // ΠΕΠΕΡΑΣΜΕΝΕΣ ΔΙΑΦΟΡΕΣ

    double Fa_p = computeF(a+h,b,c,N,false);
    double Fa_m = computeF(a-h,b,c,N,false);
    double Fb_p = computeF(a,b+h,c,N,false);
    double Fb_m = computeF(a,b-h,c,N,false);
    double Fc_p = computeF(a,b,c+h,N,false);
    double Fc_m = computeF(a,b,c-h,N,false);

    double dF_da = (Fa_p - Fa_m) / (2*h);
    double dF_db = (Fb_p - Fb_m) / (2*h);
    double dF_dc = (Fc_p - Fc_m) / (2*h);

    std::cout << std::scientific << std::setprecision(12);
    std::cout << "F = " << F0 << "\n";
    std::cout << "dF/da (FD) = " << dF_da << "\n";
    std::cout << "dF/db (FD) = " << dF_db << "\n";
    std::cout << "dF/dc (FD) = " << dF_dc << "\n";
    std::cout << "N  = " << N << "\n";
    std::cout << "dx = " << 1.0/(N-1) << "\n";


    // DD
    std::vector<double> x, delta;
    solve_delta_RK2(N, a, b, c, x, delta);

    std::vector<double> da(N), db(N), dc(N);

    
    double da0 = 0.0, db0 = 0.0, dc0 = 0.0;

    const double dx = 1.0 / (N - 1);
   
    const double Kc = (2.0 / M_PI - 0.5); 

    solve_sens_RK2(
        N, 0.0, dx,
        x.data(),
        delta.data(),
        a, b, c,
        A0, Kc, U,
        da0, db0, dc0,
        da.data(), db.data(), dc.data()
    );

    // π.χ. έλεγχος ότι όντως βγήκαν:
    // std::cout << "delta_a(1) = " << da.back() << "\n";
    // std::cout << "delta_b(1) = " << db.back() << "\n";
    // std::cout << "delta_c(1) = " << dc.back() << "\n";


    std::vector<double> df_ddelta;
    compute_df_ddelta(delta, dx, df_ddelta);

    double dF_da_DD = dot(df_ddelta, da);
    double dF_db_DD = dot(df_ddelta, db);
    double dF_dc_DD = dot(df_ddelta, dc);

    std::cout << "---------------------------- " <<  "\n";
    std::cout << "dF/da (DD) = " << dF_da_DD << "\n";
    std::cout << "dF/db (DD) = " << dF_db_DD << "\n";
    std::cout << "dF/dc (DD) = " << dF_dc_DD << "\n";


    //ADJOINT

    std::vector<double> lam(N);
    double lamL = 0;
    double A = (2/M_PI - 0.5)*U;
    double B = (M_PI*nu/2);

    solve_adjoint_RK2_backward(
        N,
        x.data(),
        delta.data(),
        A, B, C,
        a, b, c,
        lamL,
        lam.data()
    );

    std::vector<double> q_a(N), q_b(N), q_c(N);

    for (int i = 0; i < N; ++i)
    {
        q_a[i] = dv0_da(x[i], a, b, c);
        q_b[i] = dv0_db(x[i], a, b, c);
        q_c[i] = dv0_dc(x[i], a, b, c);
    }


    double Ia = trapz_lam_delta_q(x, lam, delta, q_a);
    double Ib = trapz_lam_delta_q(x, lam, delta, q_b);
    double Ic = trapz_lam_delta_q(x, lam, delta, q_c);

    std::cout << "---------------------------- " <<  "\n";
    std::cout << "Integral_a = " << Ia << "\n";
    std::cout << "Integral_b = " << Ib << "\n";
    std::cout << "Integral_c = " << Ic << "\n";

    std::ofstream flam("adjoint_lambda.txt");
    flam << std::scientific << std::setprecision(17);

    for (int i = 0; i < N; ++i)
    {
        flam << x[i] << " " << lam[i] << "\n";
    }

    flam.close();


    // DISCRETE ADJOINT
    std::vector<double> diag(N);   
    std::vector<double> psi(N);

    diag[0]=1;

    for (int i = 1; i < N; ++i)
    {
        double di = delta[i];
        double xi = x[i];

        diag[i] = -1.0 - 2*B*dx *((-2)*A*di*di + B*dx)/((2*A*di*di - v0(xi,a,b,c)*dx*di + B*dx)*(2*A*di*di - v0(xi,a,b,c)*dx*di + B*dx)) ;
    }

    
    psi[N-1] = df_ddelta[N-1];

   
    for (int i = N-2; i >= 1; --i)
    {
        psi[i] = (df_ddelta[i] - psi[i+1])/(diag[i]);
        }

    
    

        
        psi[0] = df_ddelta[0]; 

        std::ofstream psi_file("adjoint_psi.txt");

    if (!psi_file.is_open())
    {
        std::cerr << "Error opening adjoint_psi.txt" << std::endl;
    }
    else
    {
        for (int i = 0; i < N; ++i)
        {
            psi_file << i << " " << psi[i] << "\n";
        }
        psi_file.close();
    }
   

    
    std::vector<double> dR_da(N);
    std::vector<double> dR_db(N);
    std::vector<double> dR_dc(N);

    dR_da[0] = 0;
    dR_db[0] = 0;
    dR_dc[0] = 0;

    for (int i = 1; i < N; ++i)
    {
        
        double xi   = x[i];
        double di   = delta[i];
        double xmid = x[i] + dx/2;
        

        
        dR_da[i] = -2*B*dx*di*(1.0 / ((2*A*di*di + B*dx - dx*di*v0(xi,a,b,c))*(2*A*di*di + B*dx - dx*di*v0(xi,a,b,c))))*dv0_da(xi,a,b,c)*dx*di - (dx/A)*dv0_da(xmid,a,b,c);
        dR_db[i] = -2*B*dx*di*(1.0 / ((2*A*di*di + B*dx - dx*di*v0(xi,a,b,c))*(2*A*di*di + B*dx - dx*di*v0(xi,a,b,c))))*dv0_db(xi,a,b,c)*dx*di - (dx/A)*dv0_db(xmid,a,b,c);
        dR_dc[i] = -2*B*dx*di*(1.0 / ((2*A*di*di + B*dx - dx*di*v0(xi,a,b,c))*(2*A*di*di + B*dx - dx*di*v0(xi,a,b,c))))*dv0_dc(xi,a,b,c)*dx*di - (dx/A)*dv0_dc(xmid,a,b,c);
    }


    double psiT_dR_da = 0.0;
    double psiT_dR_db = 0.0;
    double psiT_dR_dc = 0.0;

    for (int i = 0; i < N; ++i)
    {
        psiT_dR_da += -psi[i] * dR_da[i];
        psiT_dR_db += -psi[i] * dR_db[i];
        psiT_dR_dc += -psi[i] * dR_dc[i];
    }


    std::cout << std::scientific << std::setprecision(12);

    std::cout << "----------------------------\n";
    std::cout << "Discrete Adjoint Gradients\n";
    std::cout << "dF/da (Adjoint) = " << psiT_dR_da << "\n";
    std::cout << "dF/db (Adjoint) = " << psiT_dR_db<< "\n";
    std::cout << "dF/dc (Adjoint) = " << psiT_dR_dc << "\n";



    //ALM

 

    // ALM params
    int outer_max = 20;
    int inner_max = 40;

    double lambdaALM = -0.01;  
    double omega = 3.0;     
    double omega_mult = 1.1;

    double eps_c = 1e-8;
    double eps_g = 1e-8;

    
    double alpha = 3;

    for (int kout = 0; kout < outer_max; ++kout)
    {
        // ===== INNER=====
        for (int kin = 0; kin < inner_max; ++kin)
        {
            
            double dFda, dFdb, dFdc;
            double F = compute_F_and_grad_CA(a, b, c, N, dFda, dFdb, dFdc);

          
            double dcda, dcdb, dcdc;
            double ccon = compute_c_and_grad(a, b, c, dcda, dcdb, dcdc);

            
            double g_a = dFda - lambdaALM*dcda + 2.0*omega*ccon*dcda;
            double g_b = dFdb - lambdaALM*dcdb + 2.0*omega*ccon*dcdb;
            double g_c = dFdc - lambdaALM*dcdc + 2.0*omega*ccon*dcdc;

            double gnorm = std::sqrt(g_a*g_a + g_b*g_b + g_c*g_c);
            if (gnorm < eps_g) break;

            // gradient step
            a -= alpha * g_a;
            b -= alpha * g_b;
            c -= alpha * g_c;
        }

        // ===== OUTER: update lambda
        double dcda, dcdb, dcdc;
        double ccon = compute_c_and_grad(a, b, c, dcda, dcdb, dcdc);

        // ALM update
        lambdaALM = lambdaALM - 2.0 * omega * ccon;

        
        double dFda_tmp, dFdb_tmp, dFdc_tmp;
        double F_now = compute_F_and_grad_CA(a, b, c, N, dFda_tmp, dFdb_tmp, dFdc_tmp);

        std::cout << "OUTER " << kout
            << "  F=" << std::scientific << F_now
            << "  c=" << ccon
            << "  lambda=" << lambdaALM
            << "  omega=" << omega
            << "  a=" << a << " b=" << b << " c=" << c
            << "\n";

        if (std::abs(ccon) < eps_c) break;

        omega *= omega_mult; 
    }

    // final 
    double Ffinal = computeF(a, b, c, N, true);
    
    double dcda_f, dcdb_f, dcdc_f;
    double c_res = compute_c_and_grad(a, b, c, dcda_f, dcdb_f, dcdc_f);

    std::cout << "FINAL: F=" << Ffinal << "  a=" << a << " b=" << b << " c=" << c
              << "  lambda=" << lambdaALM << "\n";

    
    std::cout << "FINAL CONSTRAINT residual c_con=" << std::scientific << c_res
          << "  (target: 0)\n";

    return 0;
}
