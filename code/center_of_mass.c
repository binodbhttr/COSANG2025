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

void center_of_mass(const hash_t* Ghash,const hash_t* Phash)
{
  int i, j, k, gidx, processor;
  float pcmx, pcmy, pcmz, cmass, rcheck;
  int diff, itcount, pindex;
int diff2;
  malloc_CM();
  //if(ThisTask == 0)
    //printf("Find C0M\n");
  //find particles in each processor and calculate center of mass in every processor and send to main processor
	 //Point to point parallel 
  // search particles in processor then if its part of galaxy then calc local center of mass
	//int world_size;//Maybe NTask is defined to do this, I am not sure
	//MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	//printf("world size=%d\n",world_size); this works
	//printf("counter3=%d,NumPart=%d,NTask=%d\n",counter3,NumPart,ThisTask);
	CM_Px = 0;
	CM_Py = 0;
	CM_Pz = 0;
	CM_M = 0;
	for(i=0;i<counter3;i++)
	{
		if((pindex = hash_lookup(Phash,P_list[i])) != HASH_INVALID)
		{ 
			if((gidx = hash_lookup(Ghash, P[pindex].ID)) != HASH_INVALID)
			{
			//Check if particle is within virial radius
				rcheck = pow((P[pindex].Pos[0] - AllGal[gidx].Pos[0]),2.0) + pow((P[pindex].Pos[1] - AllGal[gidx].Pos[1]),2.0) + pow((P[pindex].Pos[2] - AllGal[gidx].Pos[2]),2.0);
        
			if(rcheck <= (AllGal[gidx].Rvir * AllGal[gidx].Rvir) && AllGal[gidx].sub_len >= 1000)
			{
				CM_Pxlist[gidx] += P[pindex].Mass * P[pindex].Pos[0];
				CM_Pylist[gidx] += P[pindex].Mass * P[pindex].Pos[1];
				CM_Pzlist[gidx] += P[pindex].Mass * P[pindex].Pos[2];
				CM_Mlist[gidx] += P[pindex].Mass;
			}
  //       printf("LISTS: rcheck = %f, Rvir = %f, CM_Pxlist= %f, CM_Mlist = %f\n", rcheck, AllGal[gidx].Rvir, CM_Pxlist[gidx], CM_Mlist[gidx]);
			}
		}
		///////////////
		
	}
	//collect codes. it is faster
	    //MPI_Allreduce(&CM_Px, &pcmx, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        //MPI_Allreduce(&CM_Py, &pcmy, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        //MPI_Allreduce(&CM_Pz, &pcmz, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        //MPI_Allreduce(&CM_M, &cmass, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	  ////
//	  if(ThisTask==0)
//		  printf("COM64\n");
 	 for(j=0;j<NumGalaxies;j++)
 	 {
		pcmx=0;
		pcmy=0;
		pcmz=0;
		cmass=0;
		MPI_Allreduce(&CM_Pxlist[j], &pcmx, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	        MPI_Allreduce(&CM_Pylist[j], &pcmy, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&CM_Pzlist[j], &pcmz, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&CM_Mlist[j], &cmass, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
		float cmx=pcmx/cmass;
		float cmy=pcmy/cmass;
		float cmz=pcmz/cmass;
	     if(AllGal[j].sub_len >= 1000)
	     {
		if (cmx !=cmx)
                	AllGal[j].CM_Pos[0] = AllGal[i].Pos[0];
        	else
                	AllGal[j].CM_Pos[0]=cmx;
        	if (cmy !=cmy)
                	AllGal[j].CM_Pos[1]=AllGal[i].Pos[1];
        	else
                	AllGal[j].CM_Pos[1]=cmy;
        	if (cmz !=cmz)
                	AllGal[j].CM_Pos[2]= AllGal[i].Pos[2];
        	else
			AllGal[j].CM_Pos[2]=cmz;

       // printf("CM 1st run:  CMPosx = %f, CMPosy = %f, CMPosz = %f, HaloxPos = %f, HaloPosy = %f, HaloPosz = %f\n", AllGal[j].CM_Pos[0], AllGal[j].CM_Pos[1], AllGal[j].CM_Pos[2], AllGal[j].Pos[0], AllGal[j].Pos[1], AllGal[j].Pos[2]); 
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
  //
  //Combine two if's after test
// if(ThisTask==0)
 //	printf("COM94\n");
 //if(ThisTask==0)
//{
  for(k=0; k<NumGalaxies; k++)
  {
    itcount = 1;
    diff = 1;
    if(AllGal[k].sub_len >= 1000)
    {
     do{
		diff = cm_iterate(Ghash, Phash, k, itcount);
		itcount++;
		 //printf("Gal = %d, diff in main = %d\n", k, diff);
		//this is to ignore iteration for this dynamics test
		MPI_Allreduce(&diff, &diff2, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		//printf("Gal = %d, diff_main = %d, diff2=\n", k, diff,diff2);
		}while(diff2 != 0 && itcount<3);
    }
  }
//}
  if(ThisTask == 0)
     printf("Finished CoM calculation\n");

  return;
}
/////////////////////////////////end main, start functions

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

int cm_iterate(const hash_t* hashtable,const hash_t* Phash, int galid, int itcount)
{
  int j, i, gidx,processor,partidx;
  float rsmall, rcheck, diffx, diffy, diffz;
  float cmx_orig, cmy_orig, cmz_orig;
  float pcmx, pcmy, pcmz, cmass;
  int diff;
  CM_Px = 0;
  CM_Py = 0;
  CM_Pz = 0;
  CM_M = 0;
  //malloc_CM_iterate();
  //if(ThisTask==0)
	 // printf("CoM Iterate started !\n");
  rsmall = AllGal[galid].Rvir;
  for(j=0; j<itcount;j++)
  {
      rsmall -= 0.2 * rsmall;
  }
  //S: put in cpu 0 to end-1
  // search particles in processor then if its part of galaxy then calc local center of mass
	//int world_size;//Maybe NTask is defined to to this, I am not sure
	//MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	CM_Px = 0;
	CM_Py = 0;
	CM_Pz = 0;
	CM_M = 0;
	for(j=0;j<counter3;j++)//counterk          
	{
		if((partidx = hash_lookup(Phash, P_list[j])) != HASH_INVALID)
		{
			for(i=0;i<NumPart;i++)
			{
				if((gidx = hash_lookup(hashtable, P[i].ID)) != HASH_INVALID)
				{
					if(gidx == galid)// galid is number if galaxy we are dealing with, it is input argument
					{
					//Check if particle is within shrinking sphere radius
						rcheck = pow((P[i].Pos[0] - AllGal[galid].CM_Pos[0]),2.0) + pow((P[i].Pos[1] - AllGal[galid].CM_Pos[1]),2.0) + pow((P[i].Pos[2] - AllGal[galid].CM_Pos[2]),2.0);
							//   printf("rcheck = %f, rsmall = %f, Rvir = %f\n", rcheck, rsmall, AllGal[galid].Rvir);
						if(rcheck <= (rsmall * rsmall) )
						{
							CM_Px += P[i].Mass * P[i].Pos[0];
							CM_Py += P[i].Mass * P[i].Pos[1];
							CM_Pz += P[i].Mass * P[i].Pos[2];
							CM_M += P[i].Mass;
						}
					}
				}
			}//
		}

	}//
  // First processor gathers and compare and return 
	  //collect codes. it is faster
	MPI_Allreduce(&CM_Px, &pcmx, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&CM_Py, &pcmy, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&CM_Pz, &pcmz, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&CM_M, &cmass, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	  ////	  
	  /////
	if (AllGal[galid].CM_Pos[0] != AllGal[galid].CM_Pos[0])
		cmx_orig = AllGal[i].Pos[0];
	else
		cmx_orig = AllGal[galid].CM_Pos[0];
	if (AllGal[galid].CM_Pos[1] !=AllGal[galid].CM_Pos[1])
                cmx_orig = AllGal[i].Pos[1];
        else
                cmx_orig = AllGal[galid].CM_Pos[1];
	if (AllGal[galid].CM_Pos[2] !=AllGal[galid].CM_Pos[2])
                cmx_orig = AllGal[i].Pos[2];
        else
                cmx_orig = AllGal[galid].CM_Pos[2];

//printf("cmx_orig=%f,cmy_orig=%f,cmz_orig=%f\n",cmx_orig,cmy_orig,cmz_orig); 
	if(cmass != 0.000)
	{
		AllGal[galid].CM_Pos[0] = pcmx / cmass;
		AllGal[galid].CM_Pos[1] = pcmy / cmass;
		AllGal[galid].CM_Pos[2] = pcmz / cmass;
  
		diffx = fabs(AllGal[galid].CM_Pos[0] - cmx_orig);
		diffy = fabs(AllGal[galid].CM_Pos[1] - cmy_orig);
		diffz = fabs(AllGal[galid].CM_Pos[2] - cmz_orig);
		//printf("diffx=%f, diffy=%f,diffz=%f\n",diffx,diffy,diffz);
			// this value is not generalized yet. So be careful about units.
			// for Mpc <=0.01 and for kpc <=10
		if(diffx <= 0.01 && diffy <= 0.01 && diffz <= 0.01)
		{
			diff = 0;
		}
		else
		{
		// should be one but for test
			diff = 1;
		}
	}
	else
	{
		diff = 0;
	}
	//printf("diff before return=%d\n",diff);
	return diff;
}

