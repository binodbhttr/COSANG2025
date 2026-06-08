/*
** @author Krista McCord (2016), Revision Shahram Talei (2017)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>  
#include <mpi.h> 

#include "proto.h"
#include "allvars.h"
#include "hash.h"

static int *Gal_part_start = NULL;
static int *Gal_part_count = NULL;

void center_of_mass(const hash_t* Phash)
{
  int i, j, k, gidx;
  float pcmx, pcmy, pcmz, cmass, rcheck;
  int diff, itcount, pindex;
  int diff2;

  Gal_part_start = malloc(NumGalaxies * sizeof(int));
  Gal_part_count = malloc(NumGalaxies * sizeof(int));
  for(i = 0; i < NumGalaxies; i++)
  {
     Gal_part_start[i] = -1;
     Gal_part_count[i] = 0;
  }
  for(i = 0; i < counter3; i++)
  {
     int g = G_list[i];
     if(g >= 0 && g < NumGalaxies)
     {
        if(Gal_part_count[g] == 0)
        {
           Gal_part_start[g] = i;
        }
        Gal_part_count[g]++;
     }
  }

  malloc_CM();
  
  CM_Px = 0;
  CM_Py = 0;
  CM_Pz = 0;
  CM_M = 0;
  for(i=0; i<counter3; i++)
  {
     if((pindex = hash_lookup(Phash, P_list[i])) != HASH_INVALID)
     { 
        gidx = G_list[i];
        if(gidx >= 0 && gidx < NumGalaxies)
        {
           //Check if particle is within virial radius
           rcheck = pow((P[pindex].Pos[0] - AllGal[gidx].Pos[0]), 2.0) + pow((P[pindex].Pos[1] - AllGal[gidx].Pos[1]), 2.0) + pow((P[pindex].Pos[2] - AllGal[gidx].Pos[2]), 2.0);
        
           if(rcheck <= (AllGal[gidx].Rvir * AllGal[gidx].Rvir) && AllGal[gidx].sub_len >= 1000)
           {
              CM_Pxlist[gidx] += P[pindex].Mass * P[pindex].Pos[0];
              CM_Pylist[gidx] += P[pindex].Mass * P[pindex].Pos[1];
              CM_Pzlist[gidx] += P[pindex].Mass * P[pindex].Pos[2];
              CM_Mlist[gidx] += P[pindex].Mass;
           }
        }
     }
  }

  for(j=0; j<NumGalaxies; j++)
  {
     pcmx=0;
     pcmy=0;
     pcmz=0;
     cmass=0;
     MPI_Allreduce(&CM_Pxlist[j], &pcmx, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
     MPI_Allreduce(&CM_Pylist[j], &pcmy, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
     MPI_Allreduce(&CM_Pzlist[j], &pcmz, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
     MPI_Allreduce(&CM_Mlist[j], &cmass, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
     float cmx = pcmx / cmass;
     float cmy = pcmy / cmass;
     float cmz = pcmz / cmass;
     if(AllGal[j].sub_len >= 1000)
     {
        if (cmx != cmx)
           AllGal[j].CM_Pos[0] = AllGal[j].Pos[0];
        else
           AllGal[j].CM_Pos[0] = cmx;
        if (cmy != cmy)
           AllGal[j].CM_Pos[1] = AllGal[j].Pos[1];
        else
           AllGal[j].CM_Pos[1] = cmy;
        if (cmz != cmz)
           AllGal[j].CM_Pos[2] = AllGal[j].Pos[2];
        else
           AllGal[j].CM_Pos[2] = cmz;
     }
  }

  free(CM_Pxlist);
  free(CM_Pylist);
  free(CM_Pzlist);
  free(CM_Mlist);
  CM_Pxlist = NULL;
  CM_Pylist = NULL;
  CM_Pzlist = NULL;
  CM_Mlist = NULL;

  for(k=0; k<NumGalaxies; k++)
  {
     itcount = 1;
     diff = 1;
     if(AllGal[k].sub_len >= 1000)
     {
        do {
           diff = cm_iterate(Phash, k, itcount);
           itcount++;
           MPI_Allreduce(&diff, &diff2, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        } while(diff2 != 0 && itcount < 3);
     }
  }

  free(Gal_part_start);
  free(Gal_part_count);
  Gal_part_start = NULL;
  Gal_part_count = NULL;

  if(ThisTask == 0)
     printf("Finished CoM calculation\n");

  return;
}

void malloc_CM()
{
    CM_Pxlist = (float *)calloc(NumGalaxies, sizeof(float));
    CM_Pylist = (float *)calloc(NumGalaxies, sizeof(float));
    CM_Pzlist = (float *)calloc(NumGalaxies, sizeof(float));
    CM_Mlist = (float *)calloc(NumGalaxies, sizeof(float));
  return;
}

void malloc_CM_iterate()
{
    CM_Pxlist = (float *)calloc(1, sizeof(float));
    CM_Pylist = (float *)calloc(1, sizeof(float));
    CM_Pzlist = (float *)calloc(1, sizeof(float));
    CM_Mlist = (float *)calloc(1, sizeof(float));
  return;
}

int cm_iterate(const hash_t* Phash, int galid, int itcount)
{
  int j, partidx;
  float rsmall, rcheck, diffx, diffy, diffz;
  float cmx_orig, cmy_orig, cmz_orig;
  float pcmx, pcmy, pcmz, cmass;
  int diff;

  rsmall = AllGal[galid].Rvir;
  for(j=0; j<itcount; j++)
  {
      rsmall -= 0.2 * rsmall;
  }
  
  CM_Px = 0;
  CM_Py = 0;
  CM_Pz = 0;
  CM_M = 0;

  int start_idx = Gal_part_start[galid];
  int count = Gal_part_count[galid];
  if(start_idx >= 0 && count > 0)
  {
     for(j = 0; j < count; j++)
     {
        int list_idx = start_idx + j;
        if((partidx = hash_lookup(Phash, P_list[list_idx])) != HASH_INVALID)
        {
           rcheck = pow((P[partidx].Pos[0] - AllGal[galid].CM_Pos[0]), 2.0) + pow((P[partidx].Pos[1] - AllGal[galid].CM_Pos[1]), 2.0) + pow((P[partidx].Pos[2] - AllGal[galid].CM_Pos[2]), 2.0);
           if(rcheck <= (rsmall * rsmall))
           {
              CM_Px += P[partidx].Mass * P[partidx].Pos[0];
              CM_Py += P[partidx].Mass * P[partidx].Pos[1];
              CM_Pz += P[partidx].Mass * P[partidx].Pos[2];
              CM_M += P[partidx].Mass;
           }
        }
     }
  }

  MPI_Allreduce(&CM_Px, &pcmx, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&CM_Py, &pcmy, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&CM_Pz, &pcmz, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&CM_M, &cmass, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  
  if (AllGal[galid].CM_Pos[0] != AllGal[galid].CM_Pos[0])
     cmx_orig = AllGal[galid].Pos[0];
  else
     cmx_orig = AllGal[galid].CM_Pos[0];
  if (AllGal[galid].CM_Pos[1] != AllGal[galid].CM_Pos[1])
     cmy_orig = AllGal[galid].Pos[1];
  else
     cmy_orig = AllGal[galid].CM_Pos[1];
  if (AllGal[galid].CM_Pos[2] != AllGal[galid].CM_Pos[2])
     cmz_orig = AllGal[galid].Pos[2];
  else
     cmz_orig = AllGal[galid].CM_Pos[2];

  if(cmass != 0.000)
  {
     AllGal[galid].CM_Pos[0] = pcmx / cmass;
     AllGal[galid].CM_Pos[1] = pcmy / cmass;
     AllGal[galid].CM_Pos[2] = pcmz / cmass;
  
     diffx = fabs(AllGal[galid].CM_Pos[0] - cmx_orig);
     diffy = fabs(AllGal[galid].CM_Pos[1] - cmy_orig);
     diffz = fabs(AllGal[galid].CM_Pos[2] - cmz_orig);
     if(diffx <= 0.01 && diffy <= 0.01 && diffz <= 0.01)
     {
        diff = 0;
     }
     else
     {
        diff = 1;
     }
  }
  else
  {
     diff = 0;
  }
  return diff;
}

