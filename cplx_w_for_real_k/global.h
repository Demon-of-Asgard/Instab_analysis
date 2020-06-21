#define NW 1   // number of omega[:=1/(neutrino energy)] grids
#define NU 400 // number of costheta grids u=cos\theta
#define NP 400 // # of phi grids
#define NR 289 //# of x-grids
#define NT 289 //# of z-grids
//
//
#define nv0 6.1573851
//#define XST 1.8 // ignore
//#define ZST 0.6 // ignore
#define IEND 10 // scan kz space
//#define alpha 1.0
#define beta 0.0
//
//
#define GFIHBARC 1.16637e-11 // G_F*(hbar c)^-3 [MeV-2]
#define HBARC 197.3269631    // MeV*fm
#define VC 2.99792458e5      // speed of light in km
#define ERG 6.24150934e5     // in MeV
#define AMU 1.660538921e-24  // gram
//
//
#define PER_CM3_TO_MEV3 7.761e-33     //(cm)^-3 to MeV^3
#define GM_PER_CM3_TO_MEV4 4.3610E-54 //(gm)(cm)^-3 to (MeV)^4
#define AMU_MEV 931.4941              //(MEV)
#define GF_PER_MEV2 1.1664e-11        //(MEV)^-2
//
//
#define NX 1
#define NZ 1
#define LNX 289
#define LNZ 289
#define ND1 78
#define ND2 60
//{
#define ORE 1.0213e+01  //initial guess of real part of Omega
#define OIM 52.6210e-02 // initial guess of imag. part of Omega
#define DX 0.1
double KZ = -5.58; // kz, we assume kx=ky=0
double KZ_MAX = 11.52;
double xx = 0.40;
double zz = 0.25;
double dz = 0.001;
double dx = 0.001;

//}
double wg[NW] = {0.0};         // value of omega grids
double ug[NU] = {0.0};         // ... costheta
double pg[NP] = {0.0};         // ... phi
double gp[NW][NU][NP] = {0.0}; // g(omega,u,phi) for omega > 0
double gn[NW][NU][NP] = {0.0}; // g(omega,u,phi) for omega < 0
double dw = 0.0;               // grid size
double du = 0.0;               // ..
double dp = 0.0;               // ..
double eps = 0.0;              // \int dw du dp (gp-gn)
double epsx = 0.0;             // \int dw du dp (gp-gn) vx
double epsy = 0.0;             //... vy
double epsz = 0.0;             // ...vz
double C_Hnu = 0.0;            //sqrt(2.0) * GFIHBARC * HBARC * HBARC * 1e-26;      // sqrt(2)*G_F*n_nu [cm^-1]
long double C_Hm = 0.0;        //sqrt(2.0) * GFIHBARC * HBARC * HBARC / AMU * 1e-26;
double alpha = 0.0;

char profname[100] = " "; //"nu_anu_profile_35_40.dat";
char out_fname[100] = " ";
char params_fname[100] = " ";

struct rparams
{
  double mur;  // mu
  double lamr; // lambda
};
