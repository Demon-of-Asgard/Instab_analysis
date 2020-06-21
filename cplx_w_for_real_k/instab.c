#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include "global.h"
#include "profile_generator.h"
#include "subroutines.h"
#include "DR_function.h"

int main(void)
{
  int counter = 0;
  double Rd1, Rd2, Hd1, Hd2;
  double umax1, umax2, pmax1, pmax2, umin1, umin2, phase1 = 0.0, phase2 = 0.0;
  long double lambda_xz = 0.0, mu_xz = 0.0;
  long double n_nu = 0.0;
  long double rho_xz = 0.0, ye_xz = 0.0;

  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
  int status, i, j, k, ii;
  size_t iter = 0;
  FILE *fp1, *fp2, *fparams;

  //{
  char out_fname_tmp[100] = "";
  strcpy(out_fname_tmp, out_fname);
  int xx_tmp = round(xx * 100); //(100/115 = 0.87) Renormalizing to actual size from initial normalization (by 100).
  int zz_tmp = round(zz * 100);

  sprintf(out_fname_tmp, "../CplxWForRealK/OutInstab/InstabOut%d%d.dat", xx_tmp, zz_tmp);
  sprintf(params_fname, "../CplxWForRealK/OutInstab//params_%d_%d.dat", xx_tmp, zz_tmp);
  strcpy(out_fname, out_fname_tmp);

  sprintf(profname, "../CplxWForRealK/OutProf/ElnProf%d%d.dat", xx_tmp, zz_tmp);

  printf("saved to %s\n", profname);
  printf("Dispersion data saved to %s\n", out_fname);
  //}

  fp1 = fopen(out_fname, "w");
  if (fp1 == NULL)
  {
    printf("fp1 == NULL \n");
    exit(0);
  }
  double x_init[2] = {ORE, OIM};
  double ORE_count_0 = ORE;
  double OIM_count_0 = OIM;

  const size_t n = 2; // Number of equations

  phase1 = 0.0;
  phase2 = 0.0;
  lambda_xz = 0.0;
  mu_xz = 0.0;
  n_nu = 0.0;
  rho_xz = 0.0;
  ye_xz = 0.0;

  //{2.0, 0.05}; // roots vector
  struct rparams p = {1.0, 0.0}; // {mur,lamr,eps} = {effective neutrino strength, matter strength,lepton asymmetry}
  //
  //
  C_Hnu = sqrt(2.0) * GFIHBARC * HBARC * HBARC * 1e-26; // sqrt(2)*G_F*n_nu [cm^-1]
  C_Hm = sqrt(2.0) * GFIHBARC * HBARC * HBARC / AMU * 1e-26;

  //  Rd1=1.0;Rd2=0.75*Rd1;Hd1=0.25*Rd1;Hd2=0.25*Rd2;
  //  xx=XST*Rd1;zz=Hd1+ZST*Rd1;

  dw = 2.0;                       //2.0 / ((double)NW); //*nv0;//6.1573851;
  du = 2.0 / ((double)NU);        //
  dp = 2.0 * M_PI / ((double)NP); //

  // assign the grid values and the spectrum

  //READING THE PROFILE.

  profile_generator(xx, zz);
  n_nu = read_nu_anu_profile();
  ye_xz = read_ye(xx, zz);
  rho_xz = read_rho(xx, zz);
  printf("%LE\t%LE\t%LE\n", ye_xz, rho_xz, C_Hm);

  //
  //

  lambda_xz = rho_xz * ye_xz * C_Hm;
  mu_xz = n_nu * C_Hnu;

  printf("\n \tx_init[0] = %lf\n", x_init[0]);
  printf(" \tx_init[1] = %lf\n", x_init[1]);
  printf(" \t*********************************\n");
  printf(" \t* lambda_xz = %LE\t*\n \t* mu_xz = %LE\t\t*\n", lambda_xz, mu_xz);
  printf(" \t*********************************\n\n");

  fparams = fopen(params_fname, "w");
  fprintf(fparams, "(xx = %lf\tzz = %lf)\n lambda = %LE\t mu = %LE\t alpha = %lf \t (KZ = %lf, KZ_MAX =%lf \t ORE = %lf\tOIM = %lf)\n", xx * 0.87, zz * 0.87, lambda_xz, mu_xz, alpha, KZ, KZ_MAX, ORE, OIM);

  //
  //

  p.lamr = lambda_xz;
  p.mur = mu_xz;

  eps = 0.0;
  epsx = 0.0;
  epsy = 0.0;
  epsz = 0.0;
  phase1 = 0.0;

  for (i = 0; i < NW; i++)
  {
    for (k = 0; k < NP; k++)
    {
      //fprintf(fp1,"%4.3e %4.3e %4.3e %4.3e %4.3e\n",pg[k],umax1,umax2,umin1,umin2);
      for (j = 0; j < NU; j++)
      {
        eps = eps + dw * du * dp / (4.0 * M_PI) * (gp[i][j][k] + gn[i][j][k]);                                            //[MeV]
        epsx = epsx + dw * du * dp / (4.0 * M_PI) * (gp[i][j][k] + gn[i][j][k]) * cos(pg[k]) * sqrt(1.0 - ug[j] * ug[j]); //[MeV]
        epsy = epsy + dw * du * dp / (4.0 * M_PI) * (gp[i][j][k] + gn[i][j][k]) * sin(pg[k]) * sqrt(1.0 - ug[j] * ug[j]); //[MeV]
        epsz = epsz + dw * du * dp / (4.0 * M_PI) * (gp[i][j][k] + gn[i][j][k]) * ug[j];                                  //[MeV]
        phase1 = phase1 + dw * du * dp / (4.0 * M_PI) * gp[i][j][k];                                                      //[MeV]
        //        phase2=phase2-dw*du*dp/(4.0*M_PI)*gn[i][j][k]/alpha;
      }
    }
  }

  fprintf(fparams, " (eps)\t\t (epsx)\t\t (epsy)\t\t (epsz)\t\t (n_nu)\n");
  fprintf(fparams, "%4.3E\t%4.3E\t%4.3E\t%4.3E\t%LE\n", eps, epsx, epsy, epsz, n_nu);
  //printf("%4.3e %4.3e %4.3e %4.3e %4.3e %4.3e\n", eps, epsx, epsy, epsz, phase1, phase2);
  fclose(fparams);
  // searching for the solution
  while (KZ < KZ_MAX)
  {
    gsl_multiroot_function f = {&rosenbrock_f, n, &p};
    iter = 0;
    gsl_vector *x = gsl_vector_alloc(n);
    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc(T, 2);
    gsl_vector_set(x, 0, x_init[0]);
    gsl_vector_set(x, 1, x_init[1]);
    gsl_multiroot_fsolver_set(s, &f, x);

    do
    { // finding the root from an initial guess
      iter++;
      status = gsl_multiroot_fsolver_iterate(s);
      if (status)
      {
        break;
      }
      status = gsl_multiroot_test_residual(s->f, 1e-5);
    } while (status == GSL_CONTINUE && iter < 1000);
    printf("KZ = %lf ", KZ);
    print_state(iter, s);
    printf("status = %s\n", gsl_strerror(status));
    if (!status)
    {
      x_init[0] = gsl_vector_get(s->x, 0);
      x_init[1] = gsl_vector_get(s->x, 1); //{2.0, 0.05}; // roots vector

      if (counter == 0)
      {
        ORE_count_0 = x_init[0];
        OIM_count_0 = x_init[1];
      }

      //fprintf(fp1, "%4.3e\t%4.3e\t%4.3e\t%4.3e\n",xx, zz , x_init[0], fabs(x_init[1]));
      fprintf(fp1, "%4.3e\t%4.3e\t%4.3e\n", KZ, x_init[0], fabs(x_init[1]));
    }
    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);
    KZ = KZ + DX;
  }
  fclose(fp1);
  return (0);
}