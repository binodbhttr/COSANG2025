#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "allvars.h"
#include "core_proto.h"



void init_s(void)
{
  int i;
  Snaplistlen = fileset+1;
  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, 42);	 // start-up seed 
 // printf("Before set units\n");
  set_units_s();
  srand((unsigned) time(NULL));
 // printf("before read output snaps\n");
  read_output_snaps();
 // printf("after read output snaps\n");
  read_snap_list();
 // printf("After read snap list\n");

  for(i = 0; i < Snaplistlen; i++)
  {
   // printf("ZZ = %d, AA = %d \n", ZZ[i], AA[i]);
    ZZ[i] = 1 / AA[i] - 1;
    Age[i] = time_to_present(ZZ[i]);
   // printf("ZZ[%d] = %f, AA[%d] = %f, Age[%d] = %f\n", i, ZZ[i], i, AA[i], i, Age[i]);
   // printf("Snaplistlen = %d\n", Snaplistlen); 
  }

  a0 = 1.0 / (1.0 + Reionization_z0);
  ar = 1.0 / (1.0 + Reionization_zr);
 // printf("Before cooling\n");
 // read_cooling_functions();  //Called in main.c so it is only called once.

}



void set_units_s(void)
{

  /*UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  UnitTime_in_Megayears = UnitTime_in_s / SEC_PER_MEGAYEAR;
  G = GRAVITY / pow(UnitLength_in_cm, 3) * UnitMass_in_g * pow(UnitTime_in_s, 2);
  UnitDensity_in_cgs = UnitMass_in_g / pow(UnitLength_in_cm, 3);
  UnitPressure_in_cgs = UnitMass_in_g / UnitLength_in_cm / pow(UnitTime_in_s, 2);
  UnitCoolingRate_in_cgs = UnitPressure_in_cgs / UnitTime_in_s;
  UnitEnergy_in_cgs = UnitMass_in_g * pow(UnitLength_in_cm, 2) / pow(UnitTime_in_s, 2);
*/
 /* if(ThisTask == 0)
  {
     printf("\nHubble (internal units) = %g\n", All.Hubble);
     printf("G (internal units) = %g\n", All.G);
     printf("UnitMass_in_g = %g \n", All.UnitMass_in_g);
     printf("UnitTime_in_s = %g \n", All.UnitTime_in_s);
     printf("UnitVelocity_in_cm_per_s = %g \n", All.UnitVelocity_in_cm_per_s);
     printf("UnitDensity_in_cgs = %g \n", All.UnitDensity_in_cgs);
     printf("UnitEnergy_in_cgs = %g \n", All.UnitEnergy_in_cgs);
     printf("UnitPressure_in_cgs = %g\n", All.UnitPressure_in_cgs);
     printf("UnitTime_in_Megayears = %g\n", All.UnitTime_in_Megayears);
     printf("UnitCoolingRate_in_cgs = %g\n", All.UnitCoolingRate_in_cgs);
  }*/
  EnergySNcode = EnergySN_s / All.UnitEnergy_in_cgs * Hubble_h_s;
  EtaSNcode = EtaSN * (All.UnitMass_in_g / SOLAR_MASS) / Hubble_h_s;
 // if(ThisTask == 0)
   // printf("\nEnergySNcode*EtaSNcode= %g\n", EnergySNcode * EtaSNcode);

  // convert some physical input parameters to internal units 
  //Hubble = HUBBLE * All.UnitTime_in_s;

  // compute a few quantitites 
  RhoCrit = 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

}



void read_output_snaps(void)
{
   int i;

   for(i = 0; i < NOUT; i++)
   {
      ListOutputSnaps[i] = fileset;
   }


  /*int i;

  char buf[1000];
  FILE *fd;

  sprintf(buf, "%s", FileWithOutputSnaps_s);

  if(!(fd = fopen(buf, "r")))
  {
    printf("file `%s' not found.\n", buf);
    exit(1);
  }

  for(i = 0; i < NOUT; i++)
  {
    if(fscanf(fd, " %d ", &ListOutputSnaps[i]) != 1)
    {
      printf("I/O error in file '%s'\n", buf);
      exit(1);
    }
  }
  fclose(fd);*/
}



void read_snap_list(void)
{
  int i;
  double a;
  a = 0.0;

  a = All.TimeOfFirstSnapshot; //* All.TimeBetSnapshot;
  
  for(i = 0; i<fileset+1; i++)
  {
     AA[i] = a;
     a = a * All.TimeBetSnapshot;
  }

  if(All.Ti_Current >= TIMEBASE)
  {
     AA[fileset] = All.TimeMax; 
  }
  
  


/* FILE *fd;
  char fname[1000];

  sprintf(fname, "%s", FileWithSnapList_s);

  if(!(fd = fopen(fname, "r")))
  {
    printf("can't read output list in file '%s'\n", fname);
    ABORT(1);
  }

  Snaplistlen = 0;
  do
  {
    if(fscanf(fd, " %lg ", &AA[Snaplistlen]) == 1)
      Snaplistlen++;
    else
      break;
  }
  while(Snaplistlen < MAXSNAPS);

  fclose(fd);

  if(ThisTask == 0)
    printf("found %d defined times in snaplist.\n", Snaplistlen);
*/

}



double time_to_present(double z)
{
#define WORKSIZE 1000
  gsl_function F;
  gsl_integration_workspace *workspace;
  double time, result, abserr;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &integrand_time_to_present;

  gsl_integration_qag(&F, 1.0 / (z + 1), 1.0, 1.0 / All.Hubble,
    1.0e-8, WORKSIZE, GSL_INTEG_GAUSS21, workspace, &result, &abserr);

  time = 1 / All.Hubble * result;
 // printf("TIMETOPRESENT: time = %f, result = %f\n", time, result);
  gsl_integration_workspace_free(workspace);
  workspace = NULL;  
  // return time to present as a function of redshift 
  return time;
}



double integrand_time_to_present(double a, void *param)
{
  return 1 / sqrt(Omega / a + (1 - Omega - OmegaLambda_s) + OmegaLambda_s * a * a);
}



