

  /*******************************
   **  haloSampling.c		**
   **  Chieh-An Lin		**
   **  Version 2016.03.03	**
   *******************************/


#include "haloSampling.h"


//----------------------------------------------------------------------
//-- Functions related to halo_t, halo_node, halo_list

void set_halo_t(cosmo_hm *cmhm, halo_t *h, double z, double M, double pos[2], error **err)
{
  testErrorRet(h==NULL, peak_null, "halo_t *h is NULL", *err, __LINE__,);
  
  h->z        = z;
  double a    = 1.0 / (1.0 + z);
  h->a        = a;
  h->w        = w(cmhm->cosmo, a, 0, err);       forwardError(*err, __LINE__,); //-- wOmegar = 0
  h->M        = M;
  double c    = concentration(cmhm, M, a, err);  forwardError(*err, __LINE__,);
  h->c        = c;
  h->c_sq     = SQ(c);
  h->f        = 1.0 / (log(1.0 + c) - c/(1.0 + c));
  if (pos != NULL) {
    h->pos[0] = pos[0];
    h->pos[1] = pos[1];
  }
  
  //-- r_vir is halo's physical size
  //-- M = (4 pi / 3) * rho_vir(z) * r_vir^3
  //--   = (4 pi / 3) * (rho_vir(z) / rho_m(z)) * rho_m(z) * r_vir^3
  //--   = (4 pi / 3) * Delta_vir(z) * rho_m(0) * (1+z)^3 * r_vir^3
  //--   = (4 pi / 3) * Delta_vir(z) * rho_crit(0) * Omega_m * (r_vir / a)^3
  h->r_vir        = a * cbrt(M / (FOUR_PI_OVER_THREE * Delta_vir(cmhm, a) * CRITICAL_DENSITY * cmhm->cosmo->Omega_m));
  h->r_vir_sq     = SQ(h->r_vir);
  //-- theta_vir = r_vir / D_A
  //-- given by angular diameter distance D_A = f_K / (1+z)
  double D_l      = a * f_K(cmhm->cosmo, h->w, err); forwardError(*err, __LINE__,);
  h->theta_vir    = h->r_vir / D_l * RADIAN_TO_ARCMIN;
  h->theta_vir_sq = SQ(h->theta_vir);
  
  //-- factor = FOUR_PI_G_OVER_C2 * (D_l M f c_NFW^2 / 2 pi r_vir^2)
  //-- See rayTracing.c for details
  h->factor = FOUR_PI_G_OVER_C2 * D_l * M * h->f * h->c_sq / (TWO_PI * h->r_vir_sq);
  return;
}

halo_node *initialize_halo_node(error **err)
{
  halo_node *hNode = (halo_node*)malloc_err(sizeof(halo_node), err); forwardError(*err, __LINE__, NULL);
  hNode->h         = (halo_t*)malloc_err(sizeof(halo_t), err);       forwardError(*err, __LINE__, NULL);
  hNode->next      = NULL;
  return hNode;
}

halo_list *initialize_halo_list(error **err)
{
  halo_list *hList = (halo_list*)malloc_err(sizeof(halo_list), err); forwardError(*err, __LINE__, NULL);
  hList->length    = 0;
  hList->size      = 0;
  hList->first     = NULL;
  hList->last      = NULL;
  return hList;
}

void free_halo_list(halo_list *hList)
{
  halo_node *hNode;
  while (hList->first != NULL) {
    hNode        = hList->first;
    hList->first = hNode->next;
    if (hNode->h) {free(hNode->h); hNode->h = NULL;}
    free(hNode); hNode = NULL;
  }
  free(hList); hList = NULL;
  return;
}

void reset_halo_list(halo_list *hList)
{
  hList->size = 0;
  return;
}

void append_halo_list(cosmo_hm *cmhm, halo_list *hList, double z, double M, double pos[2], error **err)
{
  halo_t *h;
  if (hList->length == 0) {
    hList->first = initialize_halo_node(err);      forwardError(*err, __LINE__,);
    hList->last  = hList->first;
    hList->length++;
  }
  else if (hList->length == hList->size) {
    hList->last->next = initialize_halo_node(err); forwardError(*err, __LINE__,);
    hList->last       = hList->last->next;
    hList->length++;
  }
  else if (hList->size == 0) hList->last = hList->first;
  else                       hList->last = hList->last->next;
  
  h = hList->last->h;
  set_halo_t(cmhm, h, z, M, pos, err);             forwardError(*err, __LINE__,);
  hList->size++;
  return;
}

//----------------------------------------------------------------------
//-- Functions related to halo_map

halo_map *initialize_halo_map(int N1, int N2, double theta_pix, error **err)
{
  halo_map *hMap  = (halo_map*)malloc_err(sizeof(halo_map), err);                    forwardError(*err, __LINE__, NULL);
  hMap->N1        = N1;
  hMap->N2        = N2;
  hMap->length    = N1 * N2;
  hMap->total     = 0;
  hMap->theta_pix = theta_pix;
  hMap->theta_pix_inv = 1.0 / theta_pix;
  
  hMap->limits[0] = 0;
  hMap->limits[1] = N1 * theta_pix;
  hMap->limits[2] = 0;
  hMap->limits[3] = N2 * theta_pix;
  hMap->center[0] = 0.5 * hMap->limits[1];
  hMap->center[1] = 0.5 * hMap->limits[3];
  
  hMap->map       = (halo_list**)malloc_err(hMap->length * sizeof(halo_list*), err); forwardError(*err, __LINE__, NULL);
  int i;
  for (i=0; i<hMap->length; i++) hMap->map[i] = initialize_halo_list(err);
  return hMap;
}

void free_halo_map(halo_map *hMap)
{
  int i;
  if (hMap->map) {
    for (i=0; i<hMap->length; i++) {free_halo_list(hMap->map[i]); hMap->map[i] = NULL;}
    free(hMap->map); hMap->map = NULL;
  }
  free(hMap); hMap = NULL;
  return;
}

void reset_halo_map(halo_map *hMap)
{
  hMap->total = 0;
  int i;
  for (i=0; i<hMap->length; i++) reset_halo_list(hMap->map[i]);
  return;
}

void append_halo_map(cosmo_hm *cmhm, halo_map *hMap, double z, double M, double pos[2], error **err)
{
  double theta_pix = hMap->theta_pix;
  int i = (int)((pos[0] - hMap->limits[0]) / theta_pix);
  int j = (int)((pos[1] - hMap->limits[2]) / theta_pix);
  append_halo_list(cmhm, hMap->map[i+j*hMap->N1], z, M, pos, err); forwardError(*err, __LINE__,);
  hMap->total++;
  return;
}

void read_halo_map(char name[], cosmo_hm *cmhm, halo_map *hMap, error **err)
{
  FILE *file = fopen_err(name, "r", err); forwardError(*err, __LINE__,);
  printf("Reading halo map...\r");
  fflush(stdout);
  
  char buffer[STRING_LENGTH_MAX], *buffer1;
  int buffer2, count = 0;
  double pos[2], w, z, M;
  
  int c = fgetc(file);
  while (c != EOF) {
    if (c == (int)'#') buffer1 = fgets(buffer, STRING_LENGTH_MAX, file);
    else {
      ungetc(c, file);
      buffer2 = fscanf(file, "%lf %lf %*f %lf %lf\n", &pos[0], &pos[1], &z, &M);
      
      append_halo_map(cmhm, hMap, z, M, pos, err); forwardError(*err, __LINE__,);
      count++;
    }
    c = fgetc(file);
  }

  printf("Done\n");
  
  fclose(file);
  testErrorRet(count!=hMap->total, peak_match, "Halo number match error", *err, __LINE__,);
  printf("\"%s\" read       \n", name);
  printf("%d halos generated\n", count);
  return;
}

void output_halo_map(FILE *file, peak_param *peak, halo_map *hMap)
{
  fprintf(file, "# Number of halos = %d\n", hMap->total);
  fprintf(file, "#\n");
  
  if (peak->field == aardvark_hPatch04 || peak->field == aardvark_gPatch086) { //-- For aardvark, positions are RA, DEC in [deg]
    fprintf(file, "#       RA        DEC         w        z          M\n");
    fprintf(file, "#     [deg]      [deg]   [Mpc/h]      [-]  [M_sol/h]\n");
  }
  else {
    fprintf(file, "#   theta_x    theta_y        w        z          M\n");
    fprintf(file, "#  [arcmin]   [arcmin]   [Mpc/h]      [-]  [M_sol/h]\n");
  }
  
  halo_list *hList;
  halo_node *hNode;
  halo_t *h;
  int i, j;
  for (i=0; i<hMap->length; i++) {
    hList = hMap->map[i];
    for (j=0, hNode=hList->first; j<hList->size; j++, hNode=hNode->next) {
      h = hNode->h;
      fprintf(file, "  %9.3f  %9.3f  %8.3f  %7.5f  %9.3e \n", h->pos[0], h->pos[1], h->w, h->z, h->M);
    }
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to mass function

double massFct(cosmo_hm_params *cANDp, double mass, error **err)
{
  double dn = LN_10 * dn_dlnM(mass, (void*)cANDp, err); forwardError(*err, __LINE__, -1.0);
  return dn;
}

void fillMassFct(cosmo_hm_params *cANDp, sampler_t *samp, error **err)
{
  //-- WARNING: The pdf of this function is a "normalized" mass function, not the exact one!
  
  double *x    = samp->x;
  double *pdf  = samp->pdf;
  int i;
  
  //-- Fill pdf
  for (i=0; i<samp->length; i++) {
    pdf[i] = massFct(cANDp, x[i], err);
    forwardError(*err, __LINE__,);
  }
  
  //-- Precompute cdf, total pdf, <x>
  int setTotalToOne = 0; //-- This is to keep the total of pdf.
  set_sampler_t(samp, setTotalToOne);
  return;
}

void outputMassFct(char name[], cosmo_hm *cmhm, peak_param *peak, double z, error **err)
{
  cosmo_hm_params *cANDp = (cosmo_hm_params*)malloc_err(sizeof(cosmo_hm_params), err); forwardError(*err, __LINE__,);
  cANDp->model = cmhm;
  cANDp->asymptotic = 0;
  cANDp->a = 1.0/(1.0+z);
  
  FILE *file;
  if (name != NULL) file = fopen(name, "w");
  else              file = stdout;
  
  fprintf(file, "# Halo mass function\n");
  fprintf(file, "# Model = %s, redshift = %g\n", smassfct_t(cmhm->massfct), z);
  outputCosmoParam(file, cmhm, peak);
  fprintf(file, "#\n");
  //fprintf(file, "#          M      dn/dlogM         deltac     sigma*D+    nu    f(nu)   dnu/dlnM\n");
  fprintf(file, "#          M      dn/dlogM\n");
  fprintf(file, "#   [M_sol/h]  [(Mpc/h)^-3]\n");
  
  double M, dn;
  for (M=peak->M_min; M<=peak->M_max; M*=1.01) {
//       double deltac = 0.161000E+01; //delta_c(cmhm->cosmo, cANDp->a, err);
//       double     dp = D_plus(cmhm->cosmo, cANDp->a, 1, err);
//       double     sM = sqrt(sigmasqr_M(cmhm, M, err));
//       double     nu = deltac/(dp*sM); 
//       double    nfn = nufnu_val(nu);
//       double dnudlnM= dnu_dlnM(cmhm, M, cANDp->a, err);
    
    dn = massFct(cANDp, M, err); forwardError(*err, __LINE__,);
    //fprintf(file, " %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n", M, dn, deltac, sM, nu, nfn, dnudlnM/nu);
    fprintf(file, " %12.6e  %12.6e\n", M, dn);
  }
  
  fclose(file);
  free(cANDp);
  if (name != NULL) printf("\"%s\" made\n", name);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to fast simulations

double dVol(cosmo_hm *cmhm, peak_param *peak, double z, double ww, double dz, error **err)
{
  //-- Comoving volume between z and z+dz, peak->area in [arcmin^2]
  double dV, A;
  double fK = f_K(cmhm->cosmo, ww, err);                   forwardError(*err, __LINE__, -1.0);
  dV  = HUBBLE_DISTANCE * SQ(fK);
  dV *= peak->area * ARCMIN_SQ_TO_RADIAN_SQ;
  dV *= dz / sqrt(Esqr(cmhm->cosmo, 1.0/(1.0+z), 0, err)); forwardError(*err, __LINE__, -1.0); //-- wOmegar = 0
  return dV;
}

void randomizePosition(peak_param *peak, double *pos, error **err)
{
  double *Omega      = peak->Omega;
  gsl_rng *generator = peak->generator;
  
  //-- Set random position n in solid angle Omega
  double theta_x, theta_y, theta_sq;
  switch (peak->field) {
    case circle:
      theta_sq = SQ(Omega[0]);
      do {
	theta_x = gsl_ran_flat(generator, -1.0, 1.0) * Omega[0];
	theta_y = gsl_ran_flat(generator, -1.0, 1.0) * Omega[0];
      } while (SUM_SQ_2(theta_x, theta_y) > theta_sq);
      break;
      
    case rectangle:
      theta_x = gsl_ran_flat(generator, 0.0, 1.0) * Omega[0];
      theta_y = gsl_ran_flat(generator, 0.0, 1.0) * Omega[1];
      break;
      
    default:
      theta_x = gsl_ran_flat(generator, -0.5, 0.5) * Omega[0];
      theta_y = gsl_ran_flat(generator, -0.5, 0.5) * Omega[1];
      return;
   }

   pos[0] = theta_x;
   pos[1] = theta_y;
   return;
}

void setMassSamplers(cosmo_hm *cmhm, peak_param *peak, sampler_arr *sampArr, error **err)
{
  if (peak->printMode == 0) {printf("Computing mass function...\r"); fflush(stdout);}
  
  //-- Structure for mass function
  cosmo_hm_params *cANDp = (cosmo_hm_params*)malloc_err(sizeof(cosmo_hm_params), err); forwardError(*err, __LINE__,);
  cANDp->model = cmhm;
  cANDp->asymptotic = 0;
  
  double dz = peak->dz_halo; 
  double r  = pow(10, peak->dlogM); //-- Convert bin width to ratio
  sampler_t *samp;
  double *x, z, ww, dV;
  int i, j;
  
  //-- Loop over redshifts
  for (i=0, z=0.5*dz; i<sampArr->length; i++, z+=dz) {
    samp     = sampArr->array[i];
    samp->dx = peak->dlogM;
    x        = samp->x;
    x[0]     = peak->M_min;
    for (j=1; j<peak->N_M; j++) x[j] = x[j-1] * r; //-- Fill masses
    
    //-- Compute the volume
    cANDp->a = 1.0/(1.0+z);
    ww = w(cmhm->cosmo, cANDp->a, 0, err); forwardError(*err, __LINE__,); //-- wOmegar = 0
    dV = dVol(cmhm, peak, z, ww, dz, err); forwardError(*err, __LINE__,);
    
    //-- Fill mass function and cdf
    fillMassFct(cANDp, samp, err);         forwardError(*err, __LINE__,);
    samp->x_mean *= dV; //-- x_mean becomes total mass in dV
  }
  
  free(cANDp);
  if (peak->printMode < 2) printf("Mass function computation done\n");
  return;
}

void sampleHalos(cosmo_hm *cmhm, peak_param *peak, sampler_t *samp, halo_map *hMap, double z, error **err)
{
  gsl_rng *generator = peak->generator;
  double totalMass   = samp->x_mean;
  double sum = 0.0;
  double p, M, pos[2];
  
  //-- Sample masses
  while (sum < totalMass) {
    p = gsl_ran_flat(generator, 0.0, 1.0);
    M = execute_sampler_t(samp, p);
    
    randomizePosition(peak, pos, err);           forwardError(*err, __LINE__,);
    append_halo_map(cmhm, hMap, z, M, pos, err); forwardError(*err, __LINE__,);
    sum += M;
  }
  return;
}

void makeFastSimul(cosmo_hm *cmhm, peak_param *peak, sampler_arr *sampArr, halo_map *hMap, error **err)
{
  //-- Reset
  reset_halo_map(hMap);
  
  double dz = peak->dz_halo;
  double z;
  int i;
  
  for (i=0, z=0.5*dz; i<peak->N_z_halo; i++, z+=dz) {
    sampleHalos(cmhm, peak, sampArr->array[i], hMap, z, err);
    forwardError(*err, __LINE__,);
  }
  
  if      (peak->printMode < 2)   printf("%d halos generated         \n", hMap->total);
  else if (peak->printMode == 2) {printf("%6d halos, ", hMap->total); fflush(stdout);}
  return;
}

void outputFastSimul(char name[], cosmo_hm *cmhm, peak_param *peak, halo_map *hMap)
{
  FILE *file = fopen(name, "w");
  
  fprintf(file, "# Halo list, fast simulation\n");
  fprintf(file, "# Model = %s, field = %s, Omega = (%g, %g) [arcmin]\n", smassfct_t(cmhm->massfct), STR_FIELD_T(peak->field), peak->Omega[0], peak->Omega[1]);
  fprintf(file, "# z_halo_max = %g, N_z_halo = %d, M_min = %8.2e [M_sol/h], M_max = %8.2e\n", peak->z_halo_max, peak->N_z_halo, peak->M_min, peak->M_max);
  fprintf(file, "#\n");
  outputCosmoParam(file, cmhm, peak);
  fprintf(file, "#\n");
  
  output_halo_map(file, peak, hMap);
  
  fclose(file);
  printf("\"%s\" made\n", name);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to mass sheet

void makeAInterpolator(cosmo_hm *cmhm, interpolator_t *a_inter, double z_max, error **err)
{
  //-- Give a(w)
  
  double a_f = 1.0 / (1.0 + z_max);
  double a_0 = 1.0;
  int N_a = a_inter->length - 1;
  int i;
  
  //-- Fill scale factor interpolator
  a_inter->dx = -(a_0 - a_f) / (double)N_a;
  for (i=0; i<a_inter->length; i++) {
    a_inter->value[i] = a_0 + i*a_inter->dx;
    a_inter->x[i]     = w(cmhm->cosmo, a_inter->value[i], 0, err); forwardError(*err, __LINE__,);
  }
  return;
}

void makeKappa1Interpolator(cosmo_hm *cmhm, interpolator_t *a_inter, interpolator_t *kappa1_inter, double z_max, error **err)
{
  //-- TODO
  //-- Give kappa_1(a)
  return;
}

void massSheet(cosmo_hm *cmhm, peak_param *peak, sampler_t *samp, interpolator_t *a_inter, double sheet[2], error **err)
{
  double r  = pow(10, peak->dlogM); //-- Convert bin width to ratio
  double *x;
  int i, j;
  
  //-- Fill mass
  for (i=0; i<peak->N_z_halo; i++) {
    samp->dx = peak->dlogM;
    x        = samp->x;
    x[0]     = peak->M_min;
    for (j=1; j<peak->N_M; j++) x[j] = x[j-1] * r;
  }
  
  //-- Structure for mass function
  cosmo_hm_params *cANDp = (cosmo_hm_params*)malloc_err(sizeof(cosmo_hm_params), err); forwardError(*err, __LINE__,);
  cANDp->model = cmhm;
  cANDp->asymptotic = 0;
  
  double w_s   = a_inter->x[a_inter->length-1];
  double dw    = w_s / (double)peak->N_z_halo;
  double rho_0 = CRITICAL_DENSITY * cmhm->cosmo->Omega_m;
  double kappa_0 = 0.0;
  double kappa_1 = 0.0;
  double w_curr, a_curr, lambda, factor;
  
  //-- Compute kappa_0 and kappa_1
  for (i=0; i<peak->N_z_halo; i++) {
    //-- Fill mass function and cdf
    w_curr   = (i + 0.5) * dw;
    a_curr   = execute_interpolator_t(a_inter, w_curr);
    cANDp->a = a_curr;
    fillMassFct(cANDp, samp, err);                                                         forwardError(*err, __LINE__,);
    
    lambda = samp->x_mean / rho_0; //-- samp->x_mean is the average mass density.
    factor = f_K(cmhm->cosmo, w_s - w_curr, err) * f_K(cmhm->cosmo, w_curr, err) / a_curr; forwardError(*err, __LINE__,);
    kappa_0 += factor;
    kappa_1 += factor * lambda;
  }
  free(cANDp);
  
  double fK_s  = f_K(cmhm->cosmo, w_s, err); forwardError(*err, __LINE__,);
  factor = FOUR_PI_G_OVER_C2 * rho_0 * dw / fK_s;
  sheet[0] = kappa_0 * factor;
  sheet[1] = kappa_1 * factor;
  return;
}

//----------------------------------------------------------------------
//-- Main functions

void doMassFct(cosmo_hm *cmhm, peak_param *peak, error **err)
{
  char name[STRING_LENGTH_MAX];
  double z;
  for (z=0.05; z<1.0; z+=0.1) {
    sprintf(name, "massFct_%s_z%.3f", smassfct_t(cmhm->massfct), z);
    outputMassFct(name, cmhm, peak, z, err);
  }
  printf("------------------------------------------------------------------------\n");
  return;
}

void doFastSimulation(cosmo_hm *cmhm, peak_param *peak, error **err)
{
  sampler_arr *sampArr = initialize_sampler_arr(peak->N_z_halo, peak->N_M);
  setMassSamplers(cmhm, peak, sampArr, err);                                                        forwardError(*err, __LINE__,);
  halo_map *hMap       = initialize_halo_map(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
  
  makeFastSimul(cmhm, peak, sampArr, hMap, err); forwardError(*err, __LINE__,);
  outputFastSimul("haloCat", cmhm, peak, hMap);
  
  free_sampler_arr(sampArr);
  free_halo_map(hMap);
  printf("------------------------------------------------------------------------\n");
  return;
}

void doMassSheet(cosmo_hm *cmhm, peak_param *peak, double z_halo_max, double M_min, error **err)
{
  peak->z_halo_max = z_halo_max;
  peak->dz_halo    = peak->z_halo_max / (double)(peak->N_z_halo);
  peak->M_min      = M_min;
  peak->N_M        = (int)ceil(log10(peak->M_max / peak->M_min) / peak->dlogM); //-- Number of mass bins
  
  printf("z_s      = %f\n", peak->z_halo_max);
  printf("N_z_halo = %d\n", peak->N_z_halo);
  printf("logMmin  = %f\n", log10(peak->M_min));
  printf("logMmax  = %f\n", log10(peak->M_max));
  
  int N_a = 4000;
  double sheet[2];
  sampler_t *samp = initialize_sampler_t(peak->N_M);
  interpolator_t *a_inter = initialize_interpolator_t(N_a+1);
  makeAInterpolator(cmhm, a_inter, peak->z_halo_max, err); forwardError(*err, __LINE__,);
  massSheet(cmhm, peak, samp, a_inter, sheet, err);        forwardError(*err, __LINE__,);
  
  printf("kappa_0 = %f\n", sheet[0]);
  printf("kappa_1 = %f\n", sheet[1]);
  //printf("sigma_noise = %f\n", peak->sigma_noise[0]);
  
  free_sampler_t(samp);
  free_interpolator_t(a_inter);
  printf("------------------------------------------------------------------------\n");
  return;
}
//----------------------------------------------------------------------

// NEW HOD

void outputFastSimul_HOD(char name_cmhm[],char name[], cosmo_hm *cmhm, peak_param *peak, halo_map *hMap)
{

  FILE *file = fopen(name, "w");


  fprintf(file, "# Halo list, fast simulation\n");
  fprintf(file, "# Model = %s, field = %s, Omega = (%g, %g) [arcmin]\n", smassfct_t(cmhm->massfct), STR_FIELD_T(peak->field), peak->Omega[0], peak->Omega[1]);
  fprintf(file, "# z_halo_max = %g, N_z_halo = %d, M_min = %8.2e [M_sol/h], M_max = %8.2e\n", peak->z_halo_max, peak->N_z_halo, peak->M_min, peak->M_max);
  fprintf(file, "#\n");
  outputCosmoParam(file, cmhm, peak);
  fprintf(file, "#\n");
  
  output_halo_map_HOD(name_cmhm,file,cmhm, peak, hMap);
  
  fclose(file);
  printf("\"%s\" made\n", name);
  return;
}



void output_halo_map_HOD(char name_cmhm[],FILE *file, cosmo_hm *cmhm, peak_param *peak, halo_map *hMap)
{
  double ng,ngc,ngs,Mh,zz,ngal_all ;
  error *myerr = NULL, **err = &myerr;


  fprintf(file, "# Number of halos = %d\n", hMap->total);
  fprintf(file, "#\n");
  
  if (peak->field == aardvark_hPatch04 || peak->field == aardvark_gPatch086) { //-- For aardvark, positions are RA, DEC in [deg]
    fprintf(file,"#  theta_x   theta_y      w          z          M         Ngal_c    Ngal_s      Rv   \n");
    fprintf(file, "# [deg]  [deg]    [Mpc/h]     [-]     [M_sol/h]      [-]       [-]     [arcmin]    \n");
  }
  else {
    fprintf(file, "#  theta_x   theta_y      w          z          M         Ngal_c    Ngal_s      Rv \n");
    fprintf(file,"# [arcmin]  [arcmin]    [Mpc/h]    [-]     [M_sol/h]       [-]       [-]     [arcmin]  \n");
  }

  halo_list *hList;
  halo_node *hNode;
  halo_t *h;
  int i,ii, j,ngtot_cat;
  ngtot_cat=0;
  ii=0;
  for (i=0; i<hMap->length; i++) {
    hList = hMap->map[i];

	
    for (j=0, hNode=hList->first; j<hList->size; j++, hNode=hNode->next) {

      h = hNode->h;
	  Mh = h->M ;
	  zz = h->z ;
  	// read_cosmo_hm(name_cmhm, &cmhm, err);   
  	 // forwardError(*err, __LINE__,);
	  //if(zz<0.4) {name_cmhm=name2_cmhm;};
	  //if((zz>0.4)&&(zz<0.6)) {name_cmhm=name3_cmhm;};
	  //if((zz>0.6)&&(zz<0.8)) {name_cmhm=name4_cmhm;};
	  //if((zz>0.8)) {name_cmhm=name5_cmhm;};
	 // read_cosmo_hm(name_cmhm, &cmhm, err); 
     // quitOnError(*err, __LINE__, stderr);

	ngal_all=12.*180.*180. ;

//	if(zz<0.4){
		cmhm->log10M_min=13.17;
		cmhm->log10M1=14.53;
		cmhm->log10M0=11.09;
		cmhm->sigma_log_M=0.39;
		cmhm->alpha=1.27;
//		}
//	else if(zz<0.6){
//		cmhm->log10M_min=13.18;
	//	cmhm->log10M1=14.47;
	//	cmhm->log10M0=10.93;
	//	cmhm->sigma_log_M=0.3;
	//	cmhm->alpha=1.36;
	//	}
	//else if(zz<0.8){
	//	cmhm->log10M_min=12.96;
	//	cmhm->log10M1=14.1;
	//	cmhm->log10M0=12.47;
	//	cmhm->sigma_log_M=0.38;
	//	cmhm->alpha=1.28;
	//	}
	//else if(zz<1){
	//	cmhm->log10M_min=12.8;
	//	cmhm->log10M1=13.94;
	//	cmhm->log10M0=12.15;
	//	cmhm->sigma_log_M=0.33;
	//	cmhm->alpha=1.52;
	//	}
	//else if(zz>1){
	//	cmhm->log10M_min=12.62;
	//	cmhm->log10M1=13.79;
	//	cmhm->log10M0=8.67;
	//	cmhm->sigma_log_M=0.3;
	//	cmhm->alpha=1.5;
	//	};

	  ngc = Ngal_c(cmhm, Mh, cmhm->log10Mstar_min, cmhm->log10Mstar_max, err);
  	  forwardError(*err, __LINE__,);
	  ngs = Ngal_s(cmhm, Mh, cmhm->log10Mstar_min, cmhm->log10Mstar_max, err);
  	  forwardError(*err, __LINE__,);
	  ng=ngc*(1+ngs);
	  ii=ii+1;
		//printf("ok ii %d %f %f \n",ii,zz,ng);
		//return ;

	  ngtot_cat=ng+ngtot_cat;
	  fprintf(file, "%9.3f  %9.3f    %8.3f  %7.5f  %9.3e   %8.3f  %8.3f   %9.3f  \n", h->pos[0], h->pos[1], h->w, h->z, h->M,ngc,ngs,h->r_vir);
    }
  }
  printf("OK \n");
  printf("Nb galaxies created : %i \n",ngtot_cat);
  return;
}

void set_galaxy(double pos[2], double z, int halo_id, galaxy *gal, error **err) {
  gal->z = z;
  gal->halo_id = halo_id;
  if (pos != NULL) {
    gal->pos[0] = pos[0];
    gal->pos[1] = pos[1];
  }
  return;
}

galaxy_node *initialize_galaxy_node(error **err)
{
  galaxy_node *gNode = (galaxy_node*)malloc_err(sizeof(galaxy_node), err); forwardError(*err, __LINE__, NULL);
  gNode->gal         = (galaxy*)malloc_err(sizeof(galaxy), err);       forwardError(*err, __LINE__, NULL);
  gNode->next      = NULL;
  return gNode;
}

galaxy_list *initialize_galaxy_list(error **err)
{
  galaxy_list *gList = (galaxy_list*)malloc_err(sizeof(galaxy_list), err); forwardError(*err, __LINE__, NULL);
  gList->length    = 0;
  gList->size      = 0;
  gList->first     = NULL;
  gList->last      = NULL;
  return gList;
}

void free_galaxy_list(galaxy_list *gList)
{
  galaxy_node *gNode;
  while (gList->first != NULL) {
    gNode        = gList->first;
    gList->first = gNode->next;
    if (gNode->gal) {free(gNode->gal); gNode->gal = NULL;}
    free(gNode); gNode = NULL;
  }
  free(gList); gList = NULL;
  return;
}

void reset_galaxy_list(galaxy_list *gList)
{
  gList->size = 0;
  return;
}

void append_galaxy_list(double pos[2], double z, int halo_id, galaxy_list *gList, error **err)
{
  galaxy *gal;
  if (gList->length == 0) {
    gList->first = initialize_galaxy_node(err);      forwardError(*err, __LINE__,);
    gList->last  = gList->first;
    gList->length++;
  }
  else if (gList->length == gList->size) {
    gList->last->next = initialize_galaxy_node(err); forwardError(*err, __LINE__,);
    gList->last       = gList->last->next;
    gList->length++;
  }
  else if (gList->size == 0) gList->last = gList->first;
  else                       gList->last = gList->last->next;
  
  gal = gList->last->gal;
  set_galaxy(pos, z, halo_id, gal, err);             forwardError(*err, __LINE__,);
  gList->size++;
  return;
}



void outputFastSimul_galaxies(char name_cmhm[], char name[], char name2[], cosmo_hm *cmhm, peak_param *peak, halo_map *hMap)
{
  error *myerr = NULL, **err = &myerr;
  galaxy_list *gList = initialize_galaxy_list(err); forwardError(*err, __LINE__,);
  FILE *file = fopen(name, "w");
  FILE *file2 = fopen(name2, "w");

  fprintf(file, "# Halo list, fast simulation\n");
  fprintf(file, "# Model = %s, field = %s, Omega = (%g, %g) [arcmin]\n", smassfct_t(cmhm->massfct), STR_FIELD_T(peak->field), peak->Omega[0], peak->Omega[1]);
  fprintf(file, "# z_halo_max = %g, N_z_halo = %d, M_min = %8.2e [M_sol/h], M_max = %8.2e\n", peak->z_halo_max, peak->N_z_halo, peak->M_min, peak->M_max);
  fprintf(file, "#\n");
  outputCosmoParam(file, cmhm, peak);
  fprintf(file, "#\n");

  fprintf(file2, "# Galaxy catalogue\n");
  fprintf(file2,  "# Field = %s, Omega = (%g, %g) [arcmin]\n", STR_FIELD_T(peak->field), peak->Omega[0], peak->Omega[1]);
  fprintf(file2, "# n_gal = %g [arcmin^-2], z_s = %g\n", peak->n_gal, peak->z_s);
  fprintf(file2, "#\n");
  outputCosmoParam(file2, cmhm, peak);
  fprintf(file2, "#\n");

  //printf("test3 \n");
  output_halo_map_galaxies(file,file2,cmhm, peak, hMap, gList);
  //printf("test4 \n");
  free_galaxy_list(gList);
  //printf("Gfreen \n");
  fclose(file);
  fclose(file2);
  //printf("close \n");
  printf("\"%s\" made\n", name);
  printf("\"%s\" made\n", name2);
  return;
}

void output2_halo_map_galaxies(FILE *file,FILE *file2, cosmo_hm *cmhm, peak_param *peak, halo_map *hMap, galaxy_list *gList)
{
  halo_list *hList;
  halo_node *hNode;
  error *myerr = NULL, **err = &myerr;
  int i,j,k;
  srand(time(NULL));
  
  fprintf(file, "# Number of halos = %d\n", hMap->total);
  fprintf(file, "#\n");
  
  if (peak->field == aardvark_hPatch04 || peak->field == aardvark_gPatch086) { //-- For aardvark, positions are RA, DEC in [deg]
    fprintf(file,"#  theta_x   theta_y      w          z          M         Ngal_c    Ngal_s      Rv   \n");
    fprintf(file, "# [deg]  [deg]    [Mpc/h]     [-]     [M_sol/h]      [-]       [-]     [arcmin]    \n");
  }
  else {
    fprintf(file, "#  theta_x   theta_y      w          z          M         Ngal_c    Ngal_s      Rv \n");
    fprintf(file,"# [arcmin]  [arcmin]    [Mpc/h]    [-]     [M_sol/h]       [-]       [-]     [arcmin]  \n");
  }
  fprintf(file2, "#\n");
  fprintf(file2, "#\n");
  
 
  fprintf(file2, "#  theta_x   theta_y      z      halo_id    \n  ");
  fprintf(file2,"# [arcmin]  [arcmin]      [-]       [-]      \n");
  
  //printf("test5 \n");

  halo_t *h;
  //printf("l = %i  \n",hMap->length);
  for (i=0; i<hMap->length; i++) {
    hList = hMap->map[i];
    //printf("s = %i  \n",hList->size);
    //printf("i = %i  \n",i);
    for (j=0, hNode=hList->first; j<hList->size; j++, hNode=hNode->next) {
      h=hNode->h;
      double Mh=h->M;
      
      //printf("j = %i  \n",j);
      double ngc = Ngal_c(cmhm, Mh, cmhm->log10Mstar_min, cmhm->log10Mstar_max, err);
      forwardError(*err, __LINE__,);
      double ngs = Ngal_s(cmhm, Mh, cmhm->log10Mstar_min, cmhm->log10Mstar_max, err);
      forwardError(*err, __LINE__,);
      
      fprintf(file, "%9.3f  %9.3f    %8.3f  %7.5f  %9.3e   %8.3f  %8.3f   %9.3f  \n", h->pos[0], h->pos[1], h->w, h->z, Mh,ngc,ngs,h->r_vir);
     
      
      if( (rand()/(double)RAND_MAX)<ngc) {
	append_galaxy_list(h->pos, h->z, j, gList, err); forwardError(*err, __LINE__,);
	fprintf(file2, "%9.3f  %9.3f   %7.5f  %i  \n", h->pos[0], h->pos[1], h->z, j);
      }
      
  
      //printf("ngs = %8.2f \n",ngs);
      
      for (k = 0;k<ngc*ngs+0.5;k++) {
	//printf("k = %i  \n",k);
	int bool = 1;
	double r;
	while(bool){
	  double rtest = (rand()/(double)RAND_MAX);
	  if( NFW(5*rtest) > (rand()/(double)RAND_MAX)) {
	    bool = 0;
	    r = rtest*h->r_vir;
	  }
	}
	double theta = 2*M_PI*(rand()/(double)RAND_MAX);
	double phi = acos(2*(rand()/(double)RAND_MAX)-1);
	double pos[2];
        pos[0] = cos(theta) * sin(phi) * r + h->pos[0];
        pos[1] = sin(theta) * sin(phi) * r + h->pos[1];
	append_galaxy_list(pos, h->z, j, gList, err); forwardError(*err, __LINE__,);
	fprintf(file2, "%9.3f  %9.3f   %7.5f  %i  \n", pos[0], pos[1], h->z, j);
      }
      
      //printf("ok \n");
    }
  }
  return;
}



void output_halo_map_galaxies(FILE *file,FILE *file2, cosmo_hm *cmhm, peak_param *peak, halo_map *hMap, galaxy_list *gList)
{
  halo_list *hList;
  halo_node *hNode;
  error *myerr = NULL, **err = &myerr;
  int i,j,k;

  srand(time(NULL));
  
  fprintf(file, "# Number of halos = %d\n", hMap->total);
  fprintf(file, "#\n");
  
  if (peak->field == aardvark_hPatch04 || peak->field == aardvark_gPatch086) { //-- For aardvark, positions are RA, DEC in [deg]
    fprintf(file,"#  theta_x   theta_y      w          z          M         Ngal_c    Ngal_s      Rv   \n");
    fprintf(file, "# [deg]  [deg]    [Mpc/h]     [-]     [M_sol/h]      [-]       [-]     [arcmin]    \n");
  }
  else {
    fprintf(file, "#  theta_x   theta_y      w          z          M         Ngal_c    Ngal_s      Rv \n");
    fprintf(file,"# [arcmin]  [arcmin]    [Mpc/h]    [-]     [M_sol/h]       [-]       [-]     [arcmin]  \n");
  }
  fprintf(file2, "#\n");
  fprintf(file2, "#\n");
  
 
  fprintf(file2, "#  theta_x   theta_y      z      halo_id    \n  ");
  fprintf(file2,"# [arcmin]  [arcmin]      [-]       [-]      \n");
  
  //printf("test5 \n");

  halo_t *h;

		cmhm->zmin=0.01 ;
		
        double vol_z=vc(cmhm,cmhm->zmin,cmhm->zmax,err);
		double deg2rad=3.14/180.;
		double arc2dec=1./60.;
		double ng= peak->n_gal/(arc2dec*deg2rad)/(arc2dec*deg2rad);
		double n_gal_obs=ng/vol_z ;

		//printf("n_gal_obs volz %9.3E %9.3E   \n",n_gal_obs,vol_z);
		double Mmin = pow(10.,10.)*pow(n_gal_obs,-0.24);
		double M1 = 17.*Mmin;
		cmhm->log10M1=log10(M1) ;
		cmhm->log10M_min = log10(Mmin) ;


  //printf("l = %i  \n",hMap->length);
  for (i=0; i<hMap->length; i++) {
    hList = hMap->map[i];
    //printf("s = %i  \n",hList->size);
    //printf("i = %i  \n",i);
    for (j=0, hNode=hList->first; j<hList->size; j++, hNode=hNode->next) {
      h=hNode->h;
      double Mh=h->M;

		//printf("Mmin %9.3E    \n",Mmin);
        double ngc = Ngal_c(cmhm, Mh, cmhm->log10Mstar_min, cmhm->log10Mstar_max, err);
        forwardError(*err, __LINE__,);
        double ngs = Ngal_s(cmhm, Mh, cmhm->log10Mstar_min, cmhm->log10Mstar_max, err);
        forwardError(*err, __LINE__,);
      
      fprintf(file, "%9.3f  %9.3f    %8.3f  %7.5f  %9.3e   %8.3f  %8.3f   %9.3f  \n", h->pos[0], h->pos[1], h->w, h->z, Mh,ngc,ngs,h->r_vir);
     
      
      if( (rand()/(double)RAND_MAX)<ngc) {
	append_galaxy_list(h->pos, h->z, j, gList, err); forwardError(*err, __LINE__,);
	fprintf(file2, "%9.3f  %9.3f   %7.5f  %i  \n", h->pos[0], h->pos[1], h->z, j);
      }
      
  
      //printf("ngs = %8.2f \n",ngs);
      
      for (k = 0;k<ngc*ngs+0.5;k++) {
	//printf("k = %i  \n",k);
	int bool = 1;
	double r;
	while(bool){
	  double rtest = (rand()/(double)RAND_MAX);
	  if( NFW(5*rtest) > (rand()/(double)RAND_MAX)) {
	    bool = 0;
	    r = rtest*h->r_vir;
	  }
	}
	double theta = 2*M_PI*(rand()/(double)RAND_MAX);
	double phi = acos(2*(rand()/(double)RAND_MAX)-1);
	double pos[2];
        pos[0] = cos(theta) * sin(phi) * r + h->pos[0];
        pos[1] = sin(theta) * sin(phi) * r + h->pos[1];
	append_galaxy_list(pos, h->z, j, gList, err); forwardError(*err, __LINE__,);
	fprintf(file2, "%9.3f  %9.3f   %7.5f  %i  \n", pos[0], pos[1], h->z, j);
      }
      
      //printf("ok \n");
    }
  }
  return;
}


double NFW(double x)
{
  return 1/(x*(1+pow(x,2)));
}
