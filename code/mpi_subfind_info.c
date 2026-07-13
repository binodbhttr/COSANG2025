/*
** @author Krista McCord, 2016
*/

#include <mpi.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>   

#include "proto.h"
#include "allvars.h"

void Get_subfind_halo_info()
{
    int i;
    int a;
    int b;
    int num;
    int shalo;
    int sublen, subgrnr, offset, grouplen;
    int tag1, tag2, tag3, tag4, tag5, tag6, tag7, tag8, tag9, tag10, tag11;
    float grouppos1, grouppos2, grouppos3, virrad;
    float jspin1, jspin2, jspin3, veldisp;
     
    tag1 = 1; 
    tag2 = 2;
    tag3 = 3;
    tag4 = 4;
    tag5 = 5;
    tag6 = 6;
    tag7 = 7;
    tag8 = 8;
    tag9 = 9;
    tag10 = 10;
    tag11 = 11;

    for(i=0; i<NumGalaxies; i++)
    {
      num = AllGal[i].subfnum;
      shalo = AllGal[i].subhaloidx;  //AllGal[i].MhaloIdx;   //AllGal[i].subhaloidx;
  

      if(ThisTask == 0)
      {
         if(num != 0)
         {
             

             MPI_Recv(&sublen, 1, MPI_INT, num, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
             MPI_Recv(&subgrnr, 1, MPI_INT, num, tag3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
             MPI_Recv(&jspin1, 1, MPI_FLOAT, num, tag4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
             MPI_Recv(&jspin2, 1, MPI_FLOAT, num, tag5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
             MPI_Recv(&jspin3, 1, MPI_FLOAT, num, tag6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
             MPI_Recv(&veldisp, 1, MPI_FLOAT, num, tag11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

             AllGal[i].sub_len = sublen;
             AllGal[i].sub_grnr = subgrnr;
             AllGal[i].J[0] = jspin1;
             AllGal[i].J[1] = jspin2;
             AllGal[i].J[2] = jspin3;
             AllGal[i].veldisp = veldisp;
           //  printf("MPI_TEST 1: galaxy = %d, ThisTask = %d, veldisp = %f\n", i, ThisTask, veldisp);     
            
         }
         else
         {
           
           AllGal[i].sub_len = subfind_tab[shalo].sub_len;
           AllGal[i].sub_grnr = subfind_tab[shalo].sub_grnr;
           AllGal[i].J[0] = subfind_tab[shalo].sub_spin[0];
           AllGal[i].J[1] = subfind_tab[shalo].sub_spin[1];
           AllGal[i].J[2] = subfind_tab[shalo].sub_spin[2];
           AllGal[i].veldisp = subfind_tab[shalo].sub_veldisp * All.Time;

          // printf("MPI_TEST 2: galaxy = %d, ThisTask = %d, sub_veldisp = %f, shalo = %d\n", i, ThisTask, subfind_tab[shalo].sub_veldisp, shalo);
           
         }    
    

      }

      if(ThisTask == num && num != 0)
      {
         

         sublen = subfind_tab[shalo].sub_len;
         subgrnr = subfind_tab[shalo].sub_grnr;
         jspin1 = subfind_tab[shalo].sub_spin[0];
         jspin2 = subfind_tab[shalo].sub_spin[1];
         jspin3 = subfind_tab[shalo].sub_spin[2];
         veldisp = subfind_tab[shalo].sub_veldisp * All.Time;
         
        // printf("MPI_TEST 3: galaxy = %d, ThisTask = %d, sub_veldisp = %f, shalo = %d\n", i, ThisTask, subfind_tab[shalo].sub_veldisp,shalo);

         MPI_Send(&sublen, 1, MPI_INT, 0, tag2, MPI_COMM_WORLD);
         MPI_Send(&subgrnr, 1, MPI_INT, 0, tag3, MPI_COMM_WORLD);
         MPI_Send(&jspin1, 1, MPI_FLOAT, 0, tag4, MPI_COMM_WORLD);
         MPI_Send(&jspin2, 1, MPI_FLOAT, 0, tag5, MPI_COMM_WORLD);
         MPI_Send(&jspin3, 1, MPI_FLOAT, 0, tag6, MPI_COMM_WORLD);
         MPI_Send(&veldisp, 1, MPI_FLOAT, 0, tag11, MPI_COMM_WORLD);
         
      }

     MPI_Barrier(MPI_COMM_WORLD);

     MPI_Bcast(&AllGal[i].sub_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
     MPI_Bcast(&AllGal[i].sub_grnr, 1, MPI_INT, 0, MPI_COMM_WORLD);
     MPI_Bcast(&AllGal[i].J[0], 1, MPI_FLOAT,0, MPI_COMM_WORLD);
     MPI_Bcast(&AllGal[i].J[1], 1, MPI_FLOAT,0, MPI_COMM_WORLD);
     MPI_Bcast(&AllGal[i].J[2], 1, MPI_FLOAT,0, MPI_COMM_WORLD);
     MPI_Bcast(&AllGal[i].halovir, 1, MPI_FLOAT,0,MPI_COMM_WORLD);
     MPI_Bcast(&AllGal[i].haloPos[0],1,MPI_FLOAT,0,MPI_COMM_WORLD);
     MPI_Bcast(&AllGal[i].haloPos[1],1,MPI_FLOAT,0,MPI_COMM_WORLD);
     MPI_Bcast(&AllGal[i].haloPos[2],1,MPI_FLOAT,0,MPI_COMM_WORLD);
     MPI_Bcast(&AllGal[i].veldisp, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
     //printf("MPI TEST: ThisTask = %d, galaxy = %d, Veldisp = %f\n", ThisTask, i, AllGal[i].veldisp);
    }
 


  if(ThisTask == 0)
  {
    int gcount;
 
    gcount = 0;

    for(a=0; a<NTask; a++)
    {
       gcount = gcount + NumGroups[a];
       SumNumgrps[a] = gcount; 
       
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Bcast(SumNumgrps, NTask, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(NumGroups, NTask, MPI_INT, 0, MPI_COMM_WORLD);
  
 

 int grpcounter;
 int grfile;
 grpcounter = 0;
 


 for(i=0; i<NumGalaxies; i++)
 {
    grfile = 0;
    if(ThisTask == 0)
    {
       if(AllGal[i].sub_grnr < SumNumgrps[0])
       {
          AllGal[i].poffset = subfind_tab[AllGal[i].sub_grnr].group_offset;
          AllGal[i].group_len = subfind_tab[AllGal[i].sub_grnr].group_len;
          AllGal[i].haloPos[0] = subfind_tab[AllGal[i].sub_grnr].group_pos[0];
          AllGal[i].haloPos[1] = subfind_tab[AllGal[i].sub_grnr].group_pos[1];
          AllGal[i].haloPos[2] = subfind_tab[AllGal[i].sub_grnr].group_pos[2];
          AllGal[i].halovir = subfind_tab[AllGal[i].sub_grnr].group_r_crit200;

      
       }
       else 
       {
          for(b=1; b<NTask; b++)
          {
            if(AllGal[i].sub_grnr < SumNumgrps[b] && AllGal[i].sub_grnr >= SumNumgrps[b-1])
            {
               grfile = b;
   
            }
          }
       }
     }

    MPI_Bcast(&grfile, 1, MPI_INT, 0, MPI_COMM_WORLD);
   

    if(ThisTask == grfile && ThisTask != 0)
    {
   
       grpcounter = AllGal[i].sub_grnr - SumNumgrps[grfile-1];         
       offset = subfind_tab[grpcounter].group_offset;
       grouplen = subfind_tab[grpcounter].group_len;
       grouppos1 = subfind_tab[grpcounter].group_pos[0];
       grouppos2 = subfind_tab[grpcounter].group_pos[1];
       grouppos3 = subfind_tab[grpcounter].group_pos[2];
       virrad = subfind_tab[grpcounter].group_r_crit200;
       MPI_Send(&offset, 1, MPI_INT, 0, tag4, MPI_COMM_WORLD);
       MPI_Send(&grouplen, 1, MPI_INT, 0, tag5, MPI_COMM_WORLD);
       MPI_Send(&grouppos1, 1, MPI_FLOAT,0,tag7, MPI_COMM_WORLD);
       MPI_Send(&grouppos2, 1, MPI_FLOAT,0,tag8, MPI_COMM_WORLD);
       MPI_Send(&grouppos3, 1,MPI_FLOAT,0,tag9,MPI_COMM_WORLD);
       MPI_Send(&virrad, 1, MPI_FLOAT,0,tag10,MPI_COMM_WORLD);
   

    }


   if(ThisTask == 0  && grfile != 0)
   {
     MPI_Recv(&offset, 1, MPI_INT, grfile, tag4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     MPI_Recv(&grouplen, 1, MPI_INT, grfile, tag5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     MPI_Recv(&grouppos1, 1, MPI_FLOAT, grfile, tag7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     MPI_Recv(&grouppos2, 1, MPI_FLOAT, grfile, tag8, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     MPI_Recv(&grouppos3, 1, MPI_FLOAT, grfile, tag9, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     MPI_Recv(&virrad, 1, MPI_FLOAT, grfile, tag10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     AllGal[i].poffset = offset;
     AllGal[i].group_len = grouplen;
     AllGal[i].haloPos[0] = grouppos1;
     AllGal[i].haloPos[1] = grouppos2;
     AllGal[i].haloPos[2] = grouppos3;
     AllGal[i].halovir = virrad;

   

   }


   MPI_Bcast(&AllGal[i].poffset, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&AllGal[i].group_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&AllGal[i].halovir, 1, MPI_FLOAT,0,MPI_COMM_WORLD);
   MPI_Bcast(&AllGal[i].haloPos[0],1,MPI_FLOAT,0,MPI_COMM_WORLD);
   MPI_Bcast(&AllGal[i].haloPos[1],1,MPI_FLOAT,0,MPI_COMM_WORLD);
   MPI_Bcast(&AllGal[i].haloPos[2],1,MPI_FLOAT,0,MPI_COMM_WORLD);
    
   

 }

 
  return;
}


