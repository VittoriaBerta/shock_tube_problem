# shock_tube_problem
Resolution of a Sod Shock tube problem w/ a third-order numerical scheme
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                EULER EQUATIONS:                                                                  //
// System of Conservation laws for mass, momentum, energy (per unit volume).                                                        //
// Euler Eq. are intrinsically NON-linear hyperbolic PDEs.                                                                          //
//                                                                                                                                  //
// Considering an ideal, adiabatic, force-free gas/fluid w/                                                                         //
// 1) Initial conditions: rho = rho0, mx = rho*vx = rho0*vx0, E = E0;                                                               //
// 2) EoS: E = p/(gamma - 1) w/ gamma = Cp/Cv = 5/3, E = Fluid energy density = Etot - T                                            //
//                                                                                                                                  //
// The system of conservation laws is written as:                                                                                   //
// dt(rho)   + Div(rho*vx)           = 0       (mass)                                                                               //
// dt(rho*vx) + Div[(rho*vx)*vx + I*p] = 0     (momentum)                                                                           //
// dt(E)     + Div[(E + p)*vx]       = 0       (energy)                                                                             //
// Total energy density: Etot = p/(gamma - 1) + 0.5*rho*vx*vx                                                                       //
// Q = (rho, rho*vx, E) array of CONSERVATIVE VARIABLES                                                                             //
//                                                                                                                                  //
// The system can be written in the so-called "primitive form" using primitive variables:                                           //
// A = | vx    rho      0    |                                                                                                      //
//     | 0     vx      1/rho |    - A is a nxn matrix w/ COEFFICIENTS varying in space and time                                     //
//     | 0   gamma*p    vx   |    - A can be rewritten only using rho, v and cs                                                     //
// V = (rho, vx, p) = array of PRIMITIVE VARIABLES;                                                                                 //
//                                                                                                                                  //
// By solving the characteristic polynomial eq for A: Det(A - lambda*I) = 0 the right eigenvector decomposition is obtained giving: //
// lambda1 = (vx - |cs|),  r1 = (1, - |cs|/rho, cs^2)        SHOCK/RAREFACTION WAVE                                                 //
// lambda2 = (vx),         r2 = (1,0,0)                      CONTACT DISCONTINUITY                                                  //
// lambda3 = (vx + |cs|),  r3 = (1, |cs|/rho, cs^2)          SHOCK/RAREFACTION WAVE                                                 //
// |cs| = sqrt(gamma*p/rho) is the ADIABATIC SPEED OF SOUND in the fluid.                                                           //
// The eigenvalues represent the speed at which information (= the solution waves) travel.                                          //
//                                                                                                                                  //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
