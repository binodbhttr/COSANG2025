/*
** @author Krista McCord, 2016
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>  
#include <mpi.h> 

#include "proto.h"
#include "allvars.h"

void OldGal_create()
{
  malOldgal();

  Second_fill();

  if(ThisTask == 0)
  {
     OldGal_write();
  }

  return;
}


void malOldgal()
{   
    OldNumGalaxies = NumGalaxies;

    //Designate size of file for the array of structs
    if((OldGal = (struct Old_Gal_info*)malloc(NumGalaxies * sizeof(struct Old_Gal_info))) == NULL)
    {
        printf("Failed to allocate for OldGal...");
        return;
    }
  return;
}


void First_fill()
{
  int i;
//  printf("First_fill\n");
  for(i=0; i<NumGalaxies; i++)
  { 
     OldGal[i].oldColdGas = AllGal[i].ColdGas;
     OldGal[i].oldStellarMass = AllGal[i].StellarMass;
     OldGal[i].oldBulgeMass = AllGal[i].BulgeMass;
     OldGal[i].oldBlackHoleMass = AllGal[i].BlackHoleMass;
     OldGal[i].oldDiskScaleRadius = AllGal[i].DiskScaleRadius;
     OldGal[i].oldmostbound = AllGal[i].mostbound;
     OldGal[i].timestep = All.Time;
     OldGal[i].oldJ[0] = AllGal[i].J[0];
     OldGal[i].oldJ[1] = AllGal[i].J[1];
     OldGal[i].oldJ[2] = AllGal[i].J[2];
     OldGal[i].oldGalVel[0] = AllGal[i].Vel[0];
     OldGal[i].oldGalVel[1] = AllGal[i].Vel[1];
     OldGal[i].oldGalVel[2] = AllGal[i].Vel[2];


     AllGal[i].timestep = All.Time;
     AllGal[i].NewColdGas = 0.0;
     AllGal[i].NewStellarMass = 0.0;
     AllGal[i].NewBulgeMass = 0.0;
     AllGal[i].NewBlackHoleMass = 0.0;
     AllGal[i].NewDiskScaleRadius = 0.0;
     AllGal[i].NewGalVel[0] = 0.0;
     AllGal[i].NewGalVel[1] = 0.0;
     AllGal[i].NewGalVel[2] = 0.0;
     
     AllGal[i].NewJ[0] = 0.0;
     AllGal[i].NewJ[1] = 0.0;
     AllGal[i].NewJ[2] = 0.0;

  }

  return;
}

void Second_fill()
{
  int i;
 // printf("second_fill\n");
  for(i=0; i<NumGalaxies; i++)
  {
    
   
     OldGal[i].oldColdGas = AllGal[i].ColdGas;
     OldGal[i].oldStellarMass = AllGal[i].StellarMass;
     OldGal[i].oldBulgeMass = AllGal[i].BulgeMass;
     OldGal[i].oldBlackHoleMass = AllGal[i].BlackHoleMass;
     OldGal[i].oldDiskScaleRadius = AllGal[i].DiskScaleRadius;
     OldGal[i].oldmostbound = AllGal[i].mostbound;
     if(AllGal[i].sub_len >= 1000)
     {
        OldGal[i].oldGalPos[0] = AllGal[i].CM_Pos[0];
        OldGal[i].oldGalPos[1] = AllGal[i].CM_Pos[1];
        OldGal[i].oldGalPos[2] = AllGal[i].CM_Pos[2];
     }else{
        OldGal[i].oldGalPos[0] = AllGal[i].Pos[0];
        OldGal[i].oldGalPos[1] = AllGal[i].Pos[1];
        OldGal[i].oldGalPos[2] = AllGal[i].Pos[2];
     }
     OldGal[i].oldGalVel[0] = AllGal[i].Vel[0];
     OldGal[i].oldGalVel[1] = AllGal[i].Vel[1];
     OldGal[i].oldGalVel[2] = AllGal[i].Vel[2];
     OldGal[i].timestep = AllGal[i].timestep;
     OldGal[i].oldJ[0] = AllGal[i].J[0];
     OldGal[i].oldJ[1] = AllGal[i].J[1];
     OldGal[i].oldJ[2] = AllGal[i].J[2]; 
     OldGal[i].oldsub_len = AllGal[i].sub_len;
     AllGal[i].NewColdGas = 0.0;
     AllGal[i].NewStellarMass = 0.0;
     AllGal[i].NewBulgeMass = 0.0;
     AllGal[i].NewBlackHoleMass = 0.0;
     AllGal[i].NewDiskScaleRadius = 0.0;
     AllGal[i].NewGalVel[0] = 0.0;
     AllGal[i].NewGalVel[1] = 0.0;
     AllGal[i].NewGalVel[2] = 0.0;

   
     AllGal[i].NewJ[0] = 0.0;
     AllGal[i].NewJ[1] = 0.0;
     AllGal[i].NewJ[2] = 0.0;
     
  }


  return;
}

void Gal_search()
{

  int i, j;
  //printf("Gal_search\n");
  for(i=0; i<NumGalaxies; i++)
  {

    AllGal[i].NewColdGas = AllGal[i].ColdGas;
    AllGal[i].NewStellarMass = AllGal[i].StellarMass;
    AllGal[i].NewBulgeMass = AllGal[i].BulgeMass;
    AllGal[i].NewBlackHoleMass = AllGal[i].BlackHoleMass;
    AllGal[i].NewDiskScaleRadius = AllGal[i].DiskScaleRadius;
    AllGal[i].coldslope = 0.0;
    AllGal[i].coldb = AllGal[i].ColdGas;
    AllGal[i].stellarslope = 0.0;
    AllGal[i].stellarb = AllGal[i].StellarMass;
    AllGal[i].bulgeslope = 0.0;
    AllGal[i].bulgeb = AllGal[i].BulgeMass;
    AllGal[i].blackholeslope = 0.0;
    AllGal[i].blackholeb = AllGal[i].BlackHoleMass;
    AllGal[i].scaleslope = 0.0;
    AllGal[i].scaleb = AllGal[i].DiskScaleRadius;
    AllGal[i].NewJ[0] = AllGal[i].J[0];
    AllGal[i].NewJ[1] = AllGal[i].J[1];
    AllGal[i].NewJ[2] = AllGal[i].J[2];
    AllGal[i].Jb[0] = AllGal[i].J[0];
    AllGal[i].Jb[1] = AllGal[i].J[1];
    AllGal[i].Jb[2] = AllGal[i].J[2];
    AllGal[i].Jslope[0] = 0.0;
    AllGal[i].Jslope[1] = 0.0;
    AllGal[i].Jslope[2] = 0.0;
    AllGal[i].NewGalVel[0] = AllGal[i].Vel[0];
    AllGal[i].NewGalVel[1] = AllGal[i].Vel[1];
    AllGal[i].NewGalVel[2] = AllGal[i].Vel[2]; 
    AllGal[i].Velslope[0] = 0.0;
    AllGal[i].Velslope[1] = 0.0;
    AllGal[i].Velslope[2] = 0.0;
    AllGal[i].Velb[0] = AllGal[i].Vel[0];
    AllGal[i].Velb[1] = AllGal[i].Vel[1];
    AllGal[i].Velb[2] = AllGal[i].Vel[2];

    for(j=0; j<NumGalaxies; j++)
    {
        rotation_angles(j);
    }


    for(j=0; j<OldNumGalaxies; j++)
    {
      if(AllGal[i].sub_len >= 1000)
      {
       //if(AllGal[i].CM_Pos[0] < (OldGal[j].oldGalPos[0] + 0.06) && AllGal[i].CM_Pos[1] < (OldGal[j].oldGalPos[1] + 0.06) && AllGal[i].CM_Pos[2] < (OldGal[j].oldGalPos[2] + 0.06) && AllGal[i].CM_Pos[0] > (OldGal[j].oldGalPos[0] - 0.06) && AllGal[i].CM_Pos[1] > (OldGal[j].oldGalPos[1] - 0.06) && AllGal[i].CM_Pos[2] > (OldGal[j].oldGalPos[2] - 0.06))
       if( fabs(AllGal[i].CM_Pos[0] - OldGal[j].oldGalPos[0]) < 0.06 && fabs(AllGal[i].CM_Pos[1] - OldGal[j].oldGalPos[1]) < 0.06 && fabs(AllGal[i].CM_Pos[2] - OldGal[j].oldGalPos[2]) < 0.06)
       {
        // if(ThisTask == 0)
        // {
           //  printf("Extrap:  NewPosx = %f, NewPosy = %f, NewPosz = %f\n", AllGal[i].NewGalPos[0], AllGal[i].NewGalPos[1], AllGal[i].NewGalPos[2]);           
         //}
 
         calc_Mass_slopes(i,j);
         
       }
      }
    }


  }

  return;
}

void calc_Mass_slopes(int newgal, int oldgal)
{
  float slopecg, slopesm, slopebulge, slopebhm;
  float bintcg, bintsm, bintbulge, bintbhm;
  float diskmass, totalgalmass;
  float olddiskmass, oldtotgalmass;
  float oldscale, newscale, ascale;
  float scalediff, scalediff2;
  int i; 
  FILE *data2;
  char filelen2[1000];

  ascale = All.TimeOfFirstSnapshot; //* All.TimeBetSnapshot;
 

  for(i = 0; i<All.SnapshotFileCount; i++)

  {
     ascale = ascale * All.TimeBetSnapshot;

     if(i == (All.SnapshotFileCount - 2))
     {
        oldscale = ascale;
        newscale = ascale * All.TimeBetSnapshot;
        scalediff = fabs(newscale - oldscale);
     }

  }

  scalediff2 = AllGal[newgal].timestep - OldGal[oldgal].timestep;



  olddiskmass = (OldGal[oldgal].oldStellarMass - OldGal[oldgal].oldBulgeMass) + OldGal[oldgal].oldColdGas;
  oldtotgalmass = olddiskmass + OldGal[oldgal].oldBlackHoleMass;

  diskmass = (AllGal[newgal].NewStellarMass - AllGal[newgal].NewBulgeMass) + AllGal[newgal].NewColdGas;
  totalgalmass = diskmass + AllGal[newgal].NewBlackHoleMass;


  /*if(ThisTask == 0)
  {
 

    sprintf(filelen2, "%s/data2_step.txt", All.OutputDir);
    data2 = fopen(filelen2, "a+");

    

    fprintf(data2, "%d %f %f %f \n", All.SnapshotFileCount-1, All.Time, oldtotgalmass, totalgalmass);

    
    fclose(data2);
  
   }*/



  if(totalgalmass <= ((0.20 * oldtotgalmass) + oldtotgalmass))
  {
     
     slopecg = (AllGal[newgal].NewColdGas - OldGal[oldgal].oldColdGas) / scalediff;

     bintcg = AllGal[newgal].NewColdGas - (slopecg * All.Time);

     slopesm = (AllGal[newgal].NewStellarMass - OldGal[oldgal].oldStellarMass) / scalediff;

     bintsm = AllGal[newgal].NewStellarMass - (slopesm * All.Time);

     slopebulge = (AllGal[newgal].NewBulgeMass - OldGal[oldgal].oldBulgeMass) / scalediff;

     bintbulge = AllGal[newgal].NewBulgeMass - (slopebulge * All.Time);

     slopebhm = (AllGal[newgal].NewBlackHoleMass - OldGal[oldgal].oldBlackHoleMass) / scalediff;

     bintbhm = AllGal[newgal].NewBlackHoleMass - (slopebhm * All.Time);

     /*if(ThisTask == 0)
     {
       printf("slopecg = %f, bintcg = %f, slopesm = %f, bintsm = %f, All[Gal].NewColdGas = %g, OldGal[oldgal].oldColdGas = %g, NewStellarMass = %g, oldStellarMass=%g, All.Time = %f, scalediff = %f\n", slopecg, bintcg, slopesm, bintsm, AllGal[newgal].NewColdGas, OldGal[oldgal].oldColdGas, AllGal[newgal].NewStellarMass, OldGal[oldgal].oldStellarMass, All.Time, scalediff2);

     }*/

     AllGal[newgal].coldslope = slopecg;
     AllGal[newgal].coldb = bintcg;
     AllGal[newgal].stellarslope = slopesm;
     AllGal[newgal].stellarb = bintsm;
     AllGal[newgal].bulgeslope = slopebulge;
     AllGal[newgal].bulgeb = bintbulge;
     AllGal[newgal].blackholeslope = slopebhm;
     AllGal[newgal].blackholeb = bintbhm;

     calc_Radius_slopes(newgal, oldgal);

     calc_J_slopes(newgal, oldgal);
  }

  return;
}

void calc_Radius_slopes(int newgal, int oldgal)
{

  float slope, slope2, slope3, slope4;
  float bint, bint2, bint3, bint4;
  int i;
  float oldscale, newscale, ascale;
  float scalediff, scalediff2;

  ascale = All.TimeOfFirstSnapshot; //* All.TimeBetSnapshot;


  for(i = 0; i<All.SnapshotFileCount; i++)

  {
     ascale = ascale * All.TimeBetSnapshot;

     if(i == (All.SnapshotFileCount - 2))
     {
        oldscale = ascale;
        newscale = ascale * All.TimeBetSnapshot;
        scalediff = fabs(newscale - oldscale);
     }

  }

  scalediff2 = AllGal[i].timestep - OldGal[i].timestep;
 /* if(ThisTask == 0)
  {
   printf("DiskScaleRadius = %f, NewDiskRadius = %f\n", AllGal[newgal].DiskScaleRadius, AllGal[newgal].NewDiskScaleRadius);
  }*/
  
  slope = (AllGal[newgal].NewDiskScaleRadius - OldGal[oldgal].oldDiskScaleRadius) / scalediff;

  bint = AllGal[newgal].NewDiskScaleRadius - (slope * All.Time);

  slope2 = (AllGal[newgal].NewGalVel[0] - OldGal[oldgal].oldGalVel[0]) / scalediff;
  bint2 = AllGal[newgal].NewGalVel[0] - (slope2 * All.Time);
  
  slope3 = (AllGal[newgal].NewGalVel[1] - OldGal[oldgal].oldGalVel[1]) / scalediff;
  bint3 = AllGal[newgal].NewGalVel[1] - (slope3 * All.Time);

  slope4 = (AllGal[newgal].NewGalVel[2] - OldGal[oldgal].oldGalVel[2]) / scalediff;
  bint4 = AllGal[newgal].NewGalVel[2] - (slope4 * All.Time);


  AllGal[newgal].scaleslope = slope;
  AllGal[newgal].scaleb = bint;

  AllGal[newgal].Velslope[0] = slope2;
  AllGal[newgal].Velslope[1] = slope3;
  AllGal[newgal].Velslope[2] = slope4;
  AllGal[newgal].Velb[0] = bint2;
  AllGal[newgal].Velb[1] = bint3;
  AllGal[newgal].Velb[2] = bint4;


  return;
}

void calc_J_slopes(int newgal, int oldgal)
{
  float slope1, slope2, slope3;
  float bint1, bint2, bint3;
  int i;
  float oldscale, newscale, ascale;
  float scalediff;

  ascale = All.TimeOfFirstSnapshot; //* All.TimeBetSnapshot;


  for(i = 0; i<All.SnapshotFileCount; i++)

  {
     ascale = ascale * All.TimeBetSnapshot;

     if(i == (All.SnapshotFileCount - 2))
     {
        oldscale = ascale;
        newscale = ascale * All.TimeBetSnapshot;
        scalediff = fabs(newscale - oldscale);
     }

  }


  slope1 = (AllGal[newgal].NewJ[0] - OldGal[oldgal].oldJ[0]) / scalediff;

  bint1 = AllGal[newgal].NewJ[0] - (slope1 * All.Time);

  slope2 = (AllGal[newgal].NewJ[1] - OldGal[oldgal].oldJ[1]) / scalediff;

  bint2 = AllGal[newgal].NewJ[1] - (slope2 * All.Time);

  slope3 = (AllGal[newgal].NewJ[2] - OldGal[oldgal].oldJ[2]) / scalediff;

  bint3 = AllGal[newgal].NewJ[2] - (slope3 * All.Time);


  AllGal[newgal].Jslope[0] = slope1;
  AllGal[newgal].Jb[0] = bint1;
  AllGal[newgal].Jslope[1] = slope2;
  AllGal[newgal].Jb[1] = bint2;
  AllGal[newgal].Jslope[2] = slope3;
  AllGal[newgal].Jb[2] = bint3;

  /*if(ThisTask == 0)
  {
   printf("Jslope1= %f, Jb1 = %f\n", slope1, bint1);
  }
*/


}


void cal_new_galparams()
{
  int i;
  float cg, sm, bm, bhm, dsr, J1, J2, J3, Vel1, Vel2, Vel3;

 // printf("cal_new_galparams\n");
  for(i=0; i<NumGalaxies; i++)
  {
    AllGal[i].timestep = All.Time;      
    cg = (All.Time * AllGal[i].coldslope) + AllGal[i].coldb;
    sm = (All.Time * AllGal[i].stellarslope) + AllGal[i].stellarb;
    bm = (All.Time * AllGal[i].bulgeslope) + AllGal[i].bulgeb;
    bhm = (All.Time * AllGal[i].blackholeslope) + AllGal[i].blackholeb;
    dsr = (All.Time * AllGal[i].scaleslope) + AllGal[i].scaleb;
    J1 = (All.Time * AllGal[i].Jslope[0]) + AllGal[i].Jb[0];
    J2 = (All.Time * AllGal[i].Jslope[1]) + AllGal[i].Jb[1];
    J3 = (All.Time * AllGal[i].Jslope[2]) + AllGal[i].Jb[2];
    Vel1 = (All.Time * AllGal[i].Velslope[0]) + AllGal[i].Velb[0];
    Vel2 = (All.Time * AllGal[i].Velslope[1]) + AllGal[i].Velb[1];
    Vel3 = (All.Time * AllGal[i].Velslope[2]) + AllGal[i].Velb[2];
   
    /*if(ThisTask == 0)
    {
       printf("coldslope = %f, coldb = %f, scaleslope = %f, scaleb = %f, stellarslope = %f, stellarb = %f\n", AllGal[i].coldslope, AllGal[i].coldb, AllGal[i].scaleslope, AllGal[i].scaleb, AllGal[i].stellarslope, AllGal[i].stellarb);
       printf("J1 = %f, J2 = %f, J3 = %f\n", J1, J2, J3);
    }*/

    if(cg > 0.0)
    {
       AllGal[i].NewColdGas = cg;
    }

    if(sm > 0.0)
    {
       AllGal[i].NewStellarMass = sm;
    }

    if(bm > 0.0)
    {
      AllGal[i].NewBulgeMass = bm;
    }

    if(bhm > 0.0)
    {
      AllGal[i].NewBlackHoleMass = bhm;
    }

    if(dsr > 0.0)
    {
      AllGal[i].NewDiskScaleRadius = dsr;
    }

    
    AllGal[i].NewJ[0] = J1;
    AllGal[i].NewJ[1] = J2;  
    AllGal[i].NewJ[2] = J3;
    AllGal[i].NewGalVel[0] = Vel1;
    AllGal[i].NewGalVel[1] = Vel2;
    AllGal[i].NewGalVel[2] = Vel3; 

  }


  return;
}
