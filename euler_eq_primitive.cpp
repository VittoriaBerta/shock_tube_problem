#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

//------------------------------------------------------------------------------------------------------------------------------------
// Macros:
#define NX 4096                    // # Active cells
#define NGHOST 3                   // # Ghost zones per bound
#define NTOT ((NX) + (2 * NGHOST)) // Column array dimension
#define NVAR 3                     // Line array dimension
#define GAMMA (5.0/3.0)            // Heat capacity ratio for ideal gas
#define RHO 0                      // Fluid density: rho
#define MX1 1                      // Fluid momentum along x: mx = rho*vx
#define ENG 2                      // Fluid energy density: E = p/(gamma - 1)
#define VX1 MX1                    // Fluid velocity: vx
#define PRS ENG                    // Fluid pressure: p

// Essential macros:
#define MAX(a,b) ((a) > (b) ? (a):(b))
#define MIN(a,b) ((a) < (b) ? (a):(b))
#define SIGN(a)  ((a) > 0.0 ? 1.0:-1.0)
// Second-order limiters for PLM reconstruction:
#define MINMOD_LIMITER(a,b)  ((a)*(b) > 0.0 ? (fabs(a) < fabs(b) ? (a):(b)):0.0)
#define VANLEER_LIMITER(a,b) ((a)*(b) > 0.0 ? 2.0*(a)*(b)/((a)+(b)):0.0)
#define PPM_LIMITER(a,b)     ((a)*(b) > 0.0 ? MIN(2.0*MIN(fabs(a),fabs(b)),0.5*fabs(a+b))*SIGN(0.5*(a+b)):0.0)
// Spatial reconstruction order:
#define LINEAR 1
#define PARABOLIC 2
#define RECONSTRUCTION PARABOLIC
// Boundary conditions:
#define PERIODIC 0
#define OUTFLOW  1
#define BOUNDARY OUTFLOW

// -----------------------------------------------------------------------------------------------------------------------------------

void Init(double *, double **);
void Init_Array(int, int, double **);
void Boundary(double **, int, int);
void Reconstruct(double **, double **, double **, int, int);
void Flux(double **, double **, double **, double *, int, int);
void RK2Trapezoidal(double **, double **, double *, double, double, int, int);
void Prim2Cons(double **, double **, int, int);
void Cons2Prim(double **, double **, int, int);
void Write(double *, double **, int, int);

// -----------------------------------------------------------------------------------------------------------------------------------

int main()
{
  cout << setiosflags(ios::scientific);

  double    xL = -1.0,    xR = 1.0; // Space interval
  double t_min =  0.0, t_max = 0.5; // Time interval
  double     c =   0.8;             // CFL number
  double    dt = 1e-04;             // Initial timestep
  double     t = t_min;
  double         c_max;             // Maximal eigenvlue from spectral decomposition (overall)
  double dx = fabs((xR - xL))/(NX); // Spatial step

  int      ibeg =          NGHOST;  // Initial active cells domain
  int      iend =   ibeg + NX - 1;  // Final   active cells domain
  int last_step =               0;  // Ensure evolution until exactly t_max
  int   counter =               0;  // Controls data outflow on disk

  // Defining solution and grid arrays:
  // GRID array:
  double   *x = new double[NTOT]; 
  // Initial condition PRIMITIVE:
  double **V0 = new double *[NVAR]; 
  Init_Array(NVAR,NTOT,V0);
  // Initial condition CONSERVATIVE:
  double **Q0 = new double *[NVAR];
  Init_Array(NVAR,NTOT,Q0);
  // 2nd order solution:
  double **V2 = new double *[NVAR];
  Init_Array(NVAR,NTOT,V2);

  // Initializing grid array:
  for (int i = 0; i < NTOT; i++)
    x[i] = xL + (i - ibeg + 0.5) * dx;
  // Initial Conditions:
  Init(x, V0);
  Write(x, V0, ibeg, iend);

  // ---------------------------------------------------------------------------------------------------------------------------------

  while (!last_step)
  { // Evaluating dt to end up at exactly t_max:
    if (t + dt >= t_max)
    {
      last_step++;
      dt = t_max - t;
    }
    // Evolving solution over time:
    RK2Trapezoidal(V0, V2, &c_max, dt, dx, ibeg, iend);
    t += dt;
    // Redefining the timestep:
    dt = c*(dx/c_max);
    // Recording data on disk at t = t_max
    if (counter % 100 == 0)
    {
      cout << "Step = " << counter << "; t = " << t << "; dt = " << dt <<  "; cmax = " << c_max << endl;
      Write(x, V2, ibeg, iend);
    }    
    counter++;
  }

  delete x;
  delete V0[0];
  delete V0;
  delete Q0[0];
  delete Q0;
  delete V2[0];
  delete V2;

  return 0;
}

// -----------------------------------------------------------------------------------------------------------------------------------

void Init(double *x, double **A0)
{
  for (int i = 0; i < NTOT; i++)
  {
    A0[RHO][i] = x[i] < 0.0 ? 1.0:0.1;
    A0[VX1][i] = 0.0;
    A0[PRS][i] = x[i] < 0.0 ? 1.0:0.1;       
  }
}

void Boundary(double **A0, int beg, int end)
{
#if BOUNDARY == PERIODIC
  for(int m = 0; m < NVAR; m++)
  {
    for (int i = 0; i <= beg - 1; i++)
      A0[m][i] = A0[m][i + NX];

    for (int i = end; i <= end + NGHOST; i++)
      A0[m][i] = A0[m][i - NX];
  }
#endif

#if BOUNDARY == OUTFLOW 
  for(int m = 0; m < NVAR; m++)
  {
    for (int i = 0; i <= beg - 1; i++)
      A0[m][i] = A0[m][beg];

    for (int i = end; i <= end + NGHOST; i++)
      A0[m][i] = A0[m][end];
  }
#endif
}

void Reconstruct(double **A0, double **AL, double **AR, int beg, int end)
{
  // Second order PLM spatial reconstruction:
#if RECONSTRUCTION == LINEAR
  double dA_fwd, dA_bck, dA;
  for(int m = 0; m < NVAR; m ++)
  {
    for(int i = beg - 1; i <= end + 1; i++)
    {
      dA_fwd = A0[m][i+1] - A0[m][i];
      dA_bck = A0[m][i] - A0[m][i-1];
      dA = VANLEER_LIMITER(dA_fwd, dA_bck);
      AL[m][i] = A0[m][i] + 0.5*dA;
      AR[m][i - 1] = A0[m][i] - 0.5*dA;
    }
  }
#endif
  // Third order PPM spatial reconstruction:
#if RECONSTRUCTION == PARABOLIC
  double dAp, dAm, dApp, dmA, dmA_nxt;
  double A_dif, A_med, A_sqr;
  for(int m = 0; m < NVAR; m ++)
  {
    for(int i = beg - 1; i <= end + 1; i++)
    {
      dApp = (A0[m][i+2] - A0[m][i+1]);
      dAp  = (A0[m][i+1] - A0[m][i]);
      dAm  = (A0[m][i] - A0[m][i-1]);
      dmA     = PPM_LIMITER(dAp, dAm);
      dmA_nxt = PPM_LIMITER(dApp, dAp);
      // Leftmost interface value:
      AL[m][i]   = 0.5*(A0[m][i+1] + A0[m][i]) + (1.0/6.0)*(dmA - dmA_nxt);
      // Rightmost interface value:
      AR[m][i] = AL[m][i];
      // monotonicity check:
      if((AL[m][i]-A0[m][i])*(A0[m][i]-AR[m][i-1]) <= 0.0)
      {
        AR[m][i-1] = A0[m][i];
        AL[m][i]   = AL[m][i];
      }
      A_dif = (AL[m][i] - AR[m][i-1]);
      A_med = 0.5*(AL[m][i] + AR[m][i-1]);
      A_sqr = (1.0/6.0)*(A_dif)*(A_dif);
      if(A_dif*(A0[m][i] - A_med) > A_sqr)
        AR[m][i-1] = 3*A0[m][i] - 2*AL[m][i];
      else if(A_sqr < -(A_dif*(A0[m][i] - A_med)))
        AL[m][i] = 3*A0[m][i] - 2*AR[m][i-1];
    }  
  }
#endif
}

void Flux(double **AL, double **AR, double **F, double *c_max, int beg, int end)
{
  double           cs; // Adiabatic speed of sound
  double   rhoL, rhoR; // Fluid density
  double vx, vxL, vxR; // Fluid speed along x
  double       pL, pR; // Fluid pressure
  double       EL, ER; // Fluid internal energy
  double   lambda_max; // Maximal eigenvlue from spectral decomposition for each active cell

  double FR[NVAR];     // Rightmost flux
  double FL[NVAR];     // Leftmost  flux

  *c_max = 0;

  for(int i = beg - 1; i <= end; i++)
  {
    rhoL = AL[RHO][i];
    rhoR = AR[RHO][i];
    vxL  = AL[VX1][i];
    vxR  = AR[VX1][i];
    vx   = 0.5*(vxL + vxR);
    pL   = AL[PRS][i];
    pR   = AR[PRS][i];
    ER   =   pR/(GAMMA - 1) + 0.5*rhoR*vxR*vxR;
    EL   =   pL/(GAMMA - 1) + 0.5*rhoL*vxL*vxL;
    cs   = sqrt(GAMMA*(pL + pR)/(rhoL + rhoR));
    lambda_max = fabs(vx) + cs;
    FR[RHO] =         rhoR*vxR;
    FR[MX1] = FR[RHO]*vxR + pR;
    FR[ENG] =    (ER + pR)*vxR;
    FL[RHO] =         rhoL*vxL;
    FL[MX1] = FL[RHO]*vxL + pL;
    FL[ENG] =    (EL + pL)*vxL;
    // Lax-Friedrichs Riemann solver:       
    F[RHO][i] = 0.5*((FR[RHO] + FL[RHO]) - lambda_max*(rhoR - rhoL));
    F[MX1][i] = 0.5*((FR[MX1] + FL[MX1]) - lambda_max*(rhoR*vxR - rhoL*vxL));
    F[ENG][i] = 0.5*((FR[ENG] + FL[ENG]) - lambda_max*(ER - EL));
    *c_max = MAX(*c_max,lambda_max);
  }
}

void RK2Trapezoidal(double **A0, double **A2, double *c_max, double dt, double dx, int beg, int end)
{
  // Nonlinear SSP 3rd order Runge-Kutta Method: 
  // RHS: Lu = -(1/dx)*(F[m][i] - F[m][i-1]) (using Gottlieb-Shu_Tadmor's paper notation)
  double Lu;
  static double **B0, **B2, **AL, **AR, **Apred1, **Apred2, **Fpred, **Bpred1, **Bpred2;
  // Initial condition CONSERVATIVE:
  if(B0 == NULL)
  {
    B0 = new double *[NVAR];
    Init_Array(NVAR,NTOT,B0);
  }
  // 2nd order solution CONSERVATIVE:
  if(B2 == NULL)
  {
    B2 = new double *[NVAR];
    Init_Array(NVAR,NTOT,B2);
  }
  // Leftmost spatial reconstruction PRIMITIVE:
  if(AL == NULL)
  {
    AL = new double *[NVAR];
    Init_Array(NVAR,NTOT,AL);
  }
  // Rightmost spatial reconstruction PRIMITIVE:
  if(AR == NULL)
  {
    AR = new double *[NVAR];
    Init_Array(NVAR,NTOT,AR);
  }
  // 1st step solution predictor PRIMITIVE
  if(Apred1 == NULL)
  {
    Apred1 = new double *[NVAR];
    Init_Array(NVAR,NTOT,Apred1);
  }
  // 2nd step solution predictor PRIMITIVE
  if(Apred2 == NULL)
  {
    Apred2 = new double *[NVAR];
    Init_Array(NVAR,NTOT,Apred2);
  }
  // 1st step solution predictor CONSERVATIVE
  if(Bpred1 == NULL)
  {
    Bpred1 = new double *[NVAR];
    Init_Array(NVAR,NTOT,Bpred1);
  }
  // 2nd step solution predictor CONSERVATIVE
  if(Bpred2 == NULL)
  {
    Bpred2 = new double *[NVAR];
    Init_Array(NVAR,NTOT,Bpred2);
  }
  // Flux predictor PRIMITIVE
  if(Fpred == NULL)
  {
    Fpred = new double *[NVAR];
    Init_Array(NVAR,NTOT,Fpred);
  }
  //----------------------------------------------------------------------------------------------------------------------------------
  // Predictor 1:
  Boundary(A0, beg, end);
  Reconstruct(A0, AL, AR, beg, end);
  Flux(AL, AR, Fpred, c_max, beg, end);
  // Converting primitive into conservative variables to evolve the solution in time
  Prim2Cons(A0, B0, beg, end);
  for(int m = 0; m < NVAR; m ++)
  {
    for (int i = beg; i <= end; i++)
    {
      Lu = -(1/dx)*(Fpred[m][i] - Fpred[m][i-1]);
      Bpred1[m][i] = B0[m][i] + (dt) * (Lu);
    }
  }
  // Converting conservative into primitive variables to compute fluxes
  Cons2Prim(Bpred1, Apred1, beg, end);
  // Predictor 2:
  Boundary(Apred1, beg, end);
  Reconstruct(Apred1, AL, AR, beg, end);
  Flux(AL, AR, Fpred, c_max, beg, end);
  Prim2Cons(Apred1, Bpred1, beg, end);
  for(int m = 0; m < NVAR; m ++)
  {
    for (int i = beg; i <= end; i++)
    {
      Lu = -(1/dx)*(Fpred[m][i] - Fpred[m][i-1]);
      Bpred2[m][i] = 0.25*(3*B0[m][i] + Bpred1[m][i] + (dt) * (Lu));
    }
  }
  // Converting conservative into primitive variables to compute fluxes
  Cons2Prim(Bpred2, Apred2, beg, end);
  // Corrector
  Boundary(Apred2, beg, end);
  Reconstruct(Apred2, AL, AR, beg, end);
  Flux(AL, AR, Fpred, c_max, beg, end);
  Prim2Cons(Apred2, Bpred2, beg, end);
  for(int m = 0; m < NVAR; m++)
  {
    for (int i = beg; i <= end; i++)
    {
      Lu = -(1/dx)*(Fpred[m][i] - Fpred[m][i-1]);
      B2[m][i] = (1.0/3.0)*(B0[m][i] + 2*Bpred2[m][i] + 2*(dt) * (Lu));
      // Updating previous condition:
      B0[m][i] = B2[m][i];
    }
  }
  Cons2Prim(B2, A2, beg, end);
  Cons2Prim(B0, A0, beg, end);
}

void Prim2Cons(double **V, double ** Q, int beg, int end)
{
  for(int i = beg; i <= end; i++)
  {
    Q[RHO][i] = V[RHO][i];
    Q[MX1][i] = V[VX1][i]*V[RHO][i];
    Q[ENG][i] = V[PRS][i]/(GAMMA - 1) + 0.5*V[RHO][i]*V[VX1][i]*V[VX1][i];
  }
}

void Cons2Prim(double **Q, double ** V, int beg, int end)
{
  for(int i = beg; i <= end; i++)
  {
    V[RHO][i] = Q[RHO][i];
    V[VX1][i] = Q[MX1][i]/Q[RHO][i];
    V[PRS][i] = (GAMMA - 1)*(Q[ENG][i] - 0.5*V[RHO][i]*V[VX1][i]*V[VX1][i]);
  }
}

void Init_Array(int nrow, int ncol, double **matrix)
{
  matrix[0] = new double[ncol*nrow];
  for (int i = 1; i < nrow; i++)
    matrix[i] = matrix[i - 1] + ncol;
}

void Write(double *x, double **A, int beg, int end)
{
  ofstream fdata;
  char fname[32];
  static int f_counter;
  sprintf(fname, "euler_eq_primitive_%02d.dat", f_counter);

  fdata.open(fname, std::ios_base::out);
  fdata << setiosflags(ios::scientific);

    for (int i = beg; i <= end; i++)
      fdata << x[i] << " " << A[RHO][i] << " " << A[MX1][i] << " " << A[ENG][i] << endl;
    fdata << endl
          << endl;

  fdata.close();
  f_counter++;
}
