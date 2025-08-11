#include <mpi.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>   

#include "proto.h"
#include "allvars.h"

void damper()
{
   int i,j;
   //printf("Damp Particle\n");
   for(j=0;j<NumPart;j++)
   {
      if(P[j].Type == 2)
      {
        for(i=0;i<NumGalaxies;i++)
        {
           if(AllGal[i].GPID == P[j].ID)    // && AllGal[i].sub_len > 50000)
           {
              fict_part_damping(j, i);
              //printf("ThisTask fp: %d\n", ThisTask);
            //  printf("GPosx = %f, GPosy = %f, GPosz = %f, FPos[0] = %f, FPos[1] = %f, FPos[2] = %f\n", AllGal[i].Pos[0], AllGal[i].Pos[1], AllGal[i].Pos[2], P[j].Pos[0], P[j].Pos[1], P[j].Pos[2]);
           }
        }
      }
   }

return;
}

void fict_part_damping(int target, int gidx)
{
  double DFx, DFy, DFz;
  double kcon,Tscale;
  double dvx, dvy, dvz;
  double Accelx,Accely,Accelz;
  FILE *data;
  char filelen[1000];
 
  kcon = 0.0;
  Tscale = 0.0;
 
  DFx = 0.0;
  DFy = 0.0;
  DFz = 0.0;
  Accelx = P[target].g.GravAccel[0];
  Accely = P[target].g.GravAccel[1];
  Accelz = P[target].g.GravAccel[2];
  
//  if(snap_oldgal_count > 1)
 // {
   //  dvx = P[target].Vel[0] - AllGal[gidx].NewGalVel[0];
  //   dvy = P[target].Vel[1] - AllGal[gidx].NewGalVel[1];
  //   dvz = P[target].Vel[2] - AllGal[gidx].NewGalVel[2];
  //   printf("Extrap Vel:GPID = %u, PVelx = %f, PVely =%f, PVelz = %f, NVelx = %f, NVely = %f, NVelz = %f\n", AllGal[gidx].GPID, P[target].Vel[0], P[target].Vel[1], P[target].Vel[2], AllGal[gidx].NewGalVel[0], AllGal[gidx].NewGalVel[1], AllGal[gidx].NewGalVel[2]);
//  }
 // else
 // {
     dvx = P[target].Vel[0] - AllGal[gidx].Vel[0];
     dvy = P[target].Vel[1] - AllGal[gidx].Vel[1];
     dvz = P[target].Vel[2] - AllGal[gidx].Vel[2];
    // printf("NO Extrap Vel:GPID = %u, PVelx = %f, PVely =%f, PVelz = %f, NVelx = %f, NVely = %f, NVelz = %f\n", AllGal[gidx].GPID,P[target].Vel[0], P[target].Vel[1], P[target].Vel[2], AllGal[gidx].Vel[0], AllGal[gidx].Vel[1], AllGal[gidx].Vel[2]);

 // }


  Tscale = All.SofteningDisk/ (AllGal[gidx].veldisp * All.Time);
  kcon = 1.0 / Tscale;

  DFx = -1.0 * kcon * dvx;
  DFy = -1.0 * kcon * dvy;
  DFz = -1.0 * kcon * dvz;

  P[target].dampaccel[0] = DFx;
  P[target].dampaccel[1] = DFy;
  P[target].dampaccel[2] = DFz;

  //P[target].g.GravAccel[0] += DFx; 
  //P[target].g.GravAccel[1] += DFy;  
  //P[target].g.GravAccel[2] += DFz;  

  //printf(data, "%f %d %u %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", All.Time, All.SnapshotFileCount-1, AllGal[gidx].GPID, gidx, AllGal[gidx].Mvir, AllGal[gidx].Pos[0], AllGal[gidx].Pos[1], AllGal[gidx].Pos[2], AllGal[gidx].Vel[0], AllGal[gidx].Vel[1], AllGal[gidx].Vel[2], P[target].Pos[0], P[target].Pos[1], P[target].Pos[2], P[target].Vel[0], P[target].Vel[1], P[target].Vel[2], Accelx, Accely, Accelz, DFx, DFy, DFz, AllGal[gidx].veldisp);



 // printf("Galaxy = %d, kcon = %f, Tscale = %f, veldisp = %f, DFx = %f, DFy = %f, DFz = %f\n", gidx, kcon, Tscale, AllGal[gidx].veldisp, DFx, DFy, DFz); 
  
   if(All.SnapshotFileCount > 158 && All.SnapshotFileCount < 175)
   {
      sprintf(filelen, "%s/data/fp_data_%d.txt", All.OutputDir, ThisTask);
      data = fopen(filelen, "a+");
   //   fprintf(data, "%f\n", All.Time);
   fprintf(data, "%f %d %u %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", All.Time, All.SnapshotFileCount-1, AllGal[gidx].GPID, gidx, AllGal[gidx].Mvir, AllGal[gidx].Pos[0], AllGal[gidx].Pos[1], AllGal[gidx].Pos[2], AllGal[gidx].Vel[0], AllGal[gidx].Vel[1], AllGal[gidx].Vel[2], P[target].Pos[0], P[target].Pos[1], P[target].Pos[2], P[target].Vel[0], P[target].Vel[1], P[target].Vel[2], Accelx, Accely, Accelz, DFx, DFy, DFz, AllGal[gidx].veldisp, dvx,dvy,dvz);

    
   fclose(data);
   }
  
return;

}
