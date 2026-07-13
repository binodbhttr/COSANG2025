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
#include "hash.h"

void fict_particle()
{
  if(snap_fict_count == 0)
  {
     first_fict_creation();
     snap_fict_count = 1;
  }
  else
  {
     clear_GPID();
     second_fict_creation();
  }

return;
}

void first_fict_creation()
{
   int i;

   clear_GPID();
 
   for(i=0; i<NumGalaxies; i++)
   {
      if(AllGal[i].sub_len >50000)
      {
         create_fict_particle(i, AllGal[i].mostbound);
      }
   }

   New_part_count();

return;
}

void second_fict_creation()
{
   fict_part_check();
   
   check_creation();

return;
}

void clear_GPID()
{
   int i;

   for(i=0; i<NumGalaxies; i++)
   {
      AllGal[i].GPID = 0;
   } 
  

return;
}

void fict_part_check()
{
   /*******Checks to see if a ficticous particle has already been created for a galaxy
   Only needs to run after AllGal is created after a snapshot is written**************/

   int i, a, b, ii, j, aa, kk, hh;
   int *galidx;
   unsigned int *fictID;
   unsigned int fID;
   int totmpart[NTask];
   int sumidx[NTask];
   int totfids[NTask];
   int sumfids[NTask];
   unsigned int *totfidp;
   int totgidrecv;
   int *totgidx;
   int fid_count;
   //hash_t* hashtable;

   kk = 0;
   hh = 0;
   totgidrecv = 0;
   for(j=0; j<NumPart; j++)
   {
      if(P[j].Type == 2)
      {
           for(b=0; b<NumGalaxies; b++)
          {
             if(fabs(AllGal[b].Pos[0] - P[j].Pos[0]) <= 0.060 && fabs(AllGal[b].Pos[1] - P[j].Pos[1]) <= 0.060 && fabs(AllGal[b].Pos[2] - P[j].Pos[2]) <= 0.060 && AllGal[b].sub_len > 40000)
             {
                kk++;
             }
          }
      }

   }

   galidx = malloc(kk * sizeof(int));
   fictID = malloc(kk * sizeof(unsigned int));
   

   for(j=0; j<NumPart; j++)
   {
      if(P[j].Type == 2)
      {
          for(b=0; b<NumGalaxies; b++)
          {
            // printf("Galaxy: ThisTask = %d, GPosx = %f, GPosy = %f, GPosz = %f\n", ThisTask, AllGal[b].Pos[0], AllGal[b].Pos[1], AllGal[b].Pos[2]);
            // printf("ThisTask = %d, Diff x = %f, Diff y = %f, Diff z = %f\n", ThisTask, fabs(AllGal[b].Pos[0] - P[j].Pos[0]),fabs(AllGal[b].Pos[1] - P[j].Pos[1]),fabs(AllGal[b].Pos[2] - P[j].Pos[2]));
             if(fabs(AllGal[b].Pos[0] - P[j].Pos[0]) <= 0.060 && fabs(AllGal[b].Pos[1] - P[j].Pos[1]) <= 0.060 && fabs(AllGal[b].Pos[2] - P[j].Pos[2]) <= 0.060 && AllGal[b].sub_len > 40000)
             {
               
                  galidx[hh] = b;
                  fictID[hh] = P[j].ID;
                  hh++;
 
                  printf("GPosx = %f, GPosy = %f, GPosz = %f, Posx = %f, Posy = %f, Posz = %f\n", AllGal[b].Pos[0], AllGal[b].Pos[1], AllGal[b].Pos[2], P[j].Pos[0], P[j].Pos[1], P[j].Pos[2]);
               
             } 
                        
          }      
   
      }

   }


MPI_Allgather(&kk, 1, MPI_INT, totfids, 1, MPI_INT, MPI_COMM_WORLD);


for(j=0; j<NTask; j++)
{
    if(j == 0)
    {
       sumfids[j] = 0;
    }
    else
    {
       sumfids[j] = sumfids[j-1] + totfids[j-1];
    }
}

for(i=0; i<NTask; i++)
{
    totgidrecv = totgidrecv + totfids[i];
}

if((totfidp = malloc((totgidrecv) * sizeof(unsigned int))) == NULL)
{
      perror("Failed to allocate for totposx...");
      return;
}


if((totgidx = malloc((totgidrecv) * sizeof(int))) == NULL)
{
      perror("Failed to allocate for totgalidx...");
      return;
}


MPI_Allgatherv(fictID, kk, MPI_UNSIGNED, totfidp, totfids, sumfids, MPI_UNSIGNED, MPI_COMM_WORLD);
MPI_Allgatherv(galidx, kk, MPI_INT, totgidx, totfids, sumfids , MPI_INT, MPI_COMM_WORLD);

 MPI_Barrier(MPI_COMM_WORLD);
/***Fill AllGal.GPID at the same time on every processor for ficticous particle IDs that already exist***/
 
   for(i=0;i<NumGalaxies;i++)
   { 
      
       fid_count = 0;
       for(aa=0;aa<totgidrecv;aa++)
       {
          
          if(totgidx[aa] == i && fid_count==0)
          {
               AllGal[i].GPID = totfidp[aa];
               //printf("FID Count: GPID = %u\n", AllGal[i].GPID);
               fid_count = 1;
          }
       }      
    }

free(galidx);
free(fictID);
free(totgidx);
free(totfidp);
galidx = NULL;
fictID = NULL;
totgidx = NULL;
totfidp = NULL;


return;
}


void check_creation()
{
   int i;

   for(i=0; i<NumGalaxies; i++)
   {
       
       if(AllGal[i].GPID == 0 && AllGal[i].sub_len > 50000)
       {
          create_fict_particle(i, AllGal[i].mostbound);
       }
   }
   
   New_part_count();
}

void create_fict_particle(int gidx, unsigned int mbid)
{

   int a,i, j;
   int pindex;
   int procnum;
   int totmpart[NTask];
   int totalreceive;
   int countproc;
   int tag1;
   int galaxy_idx;
   size_t bytes;
   unsigned int GPID;
   tag1 = 0;
   countproc = 0;
   totalreceive = 0;
 
   procnum = -1;


   New_part_count();

   for(a=0; a<NumPart; a++)
   {

     if(P[a].ID == mbid)
     {
         procnum = ThisTask;
         pindex = a;
     

         if(ThisTask == procnum && procnum >= 0)
         {
           
             if(ThisTask == 0)
             {
               countproc = 1;
             }
             P[NumPart]=P[pindex];
             P[NumPart].ID = All.TotNumPart;  
             P[NumPart].Type = 2;
             P[NumPart].Mass = 0.0;
             P[NumPart].InitMass = 0.0;
             //P[NumPart].Pos[0] = AllGal[gidx].Pos[0];
            // P[NumPart].Pos[1] = AllGal[gidx].Pos[1];
             //P[NumPart].Pos[2] = AllGal[gidx].Pos[2];
             P[NumPart].Vel[0] = AllGal[gidx].Vel[0]; 
             P[NumPart].Vel[1] = AllGal[gidx].Vel[1];
             P[NumPart].Vel[2] = AllGal[gidx].Vel[2];        
             printf("Create Fict Particle: GPosx = %f, GPosy=%f, GPosz=%f, Posx = %f, Posy=%f, Posz =%f\n", AllGal[gidx].Pos[0], AllGal[gidx].Pos[1], AllGal[gidx].Pos[2], P[NumPart].Pos[0], P[NumPart].Pos[1], P[NumPart].Pos[2]);
             if(ThisTask == 0)
             {
              
                   AllGal[gidx].GPID = All.TotNumPart;//+ gidx;
               
             }

             NumPart = NumPart + 1;
            
             GPID = All.TotNumPart;
            
             
             if(ThisTask != 0)
             {
               
                MPI_Send(&GPID, 1, MPI_UNSIGNED, 0, tag1, MPI_COMM_WORLD);
             }

         }
     }

    
   }
   
  
   if(ThisTask == 0)
   {
      if(countproc != 1)
      {
         
         MPI_Recv(&GPID, 1, MPI_UNSIGNED, MPI_ANY_SOURCE, tag1, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);
         AllGal[gidx].GPID = GPID;
      }
   }

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Bcast(&AllGal[gidx].GPID, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

   //printf("GPID = %u \n", AllGal[gidx].GPID);
   
 return;
}

void New_part_count()
{

   int totmpart[NTask];
   long long totalreceive;
   int i;

   totalreceive = 0;

  
   MPI_Allgather(&NumPart, 1, MPI_INT, totmpart, 1, MPI_INT, MPI_COMM_WORLD);

   
   for(i=0; i<NTask; i++)
   {
     
      totalreceive = totalreceive + totmpart[i];
     
   } 

   All.TotNumPart = totalreceive;

   if(ThisTask == 0)
   {
      printf("TotNumPart = %lld\n", All.TotNumPart);
   }
return;
}
