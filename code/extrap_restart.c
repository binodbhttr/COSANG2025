/*
** @author Krista McCord, 2016
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>  
#include <mpi.h> 
#include <sys/stat.h>

#include "proto.h"
#include "allvars.h"


void OldGal_write(){

   FILE *data;
   char filelen[1000];
   char bufdir[1000];
   int i;

   for(i=0;i<OldNumGalaxies;i++)
   {
      OldGal[i].oldGalVel[0] = OldGal[i].oldGalVel[0] / OldGal[i].timestep;
      OldGal[i].oldGalVel[1] = OldGal[i].oldGalVel[1] / OldGal[i].timestep;
      OldGal[i].oldGalVel[2] = OldGal[i].oldGalVel[2] / OldGal[i].timestep;
   }


   sprintf(bufdir, "%s/extrap_restart", All.OutputDir);
   mkdir(bufdir, 02755);


   sprintf(filelen, "%s/extrap_restart/extrap_%d", All.OutputDir, All.SnapshotFileCount-2);
   data = fopen(filelen, "w");

   fwrite(&OldNumGalaxies, sizeof(int), 1, data);
   fwrite(OldGal, sizeof(struct Old_Gal_info), OldNumGalaxies, data);

   fclose(data);


return;
}

void read_extrap_restart(){
   int totgals;
   int i;

   totgals = 0;

   if(ThisTask == 0)
   {
     totgals = read_extrap_header();
   }

   MPI_Bcast(&totgals, 1, MPI_INT, 0, MPI_COMM_WORLD);

   malOldgal_extrap(totgals);

   if(ThisTask == 0)
   {
      read_extrap(); 
   }
   
   for(i=0; i<totgals; i++)
   {
      MPI_Bcast(&OldGal[i].oldColdGas, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&OldGal[i].oldStellarMass, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&OldGal[i].oldBulgeMass, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&OldGal[i].oldBlackHoleMass, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&OldGal[i].oldDiskScaleRadius, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&OldGal[i].timestep, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&OldGal[i].oldmostbound, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&OldGal[i].oldJ[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&OldGal[i].oldJ[1], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&OldGal[i].oldJ[2], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&OldGal[i].oldGalPos[0], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&OldGal[i].oldGalPos[1], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&OldGal[i].oldGalPos[2], 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&OldGal[i].oldsub_len, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

   }

   //printf("OldGalPos x = %f, y = %f, z = %f \n", OldGal[1].oldGalPos[0], OldGal[1].oldGalPos[1], OldGal[1].oldGalPos[2]);
   //printf("timestep = %f, oldmostbound = %u, oldcoldgas = %f\n", OldGal[0].timestep, OldGal[0].oldmostbound, OldGal[0].oldColdGas);

   return;
}


int read_extrap_header(){

    int Ntrees;
    int NtotGal;
    int totmal; 
    int i;
    FILE *fd;
    char file1[1000];
    totmal = 0;
     
    sprintf(file1, "%s/extrap_restart/extrap_%d", All.OutputDir, All.SnapshotFileCount-2);
            
    //Open extrap restart file
        
    fd = fopen(file1, "rb");
        
    if(NULL == fd)
    {
       printf("Cannot open Extrap Restart file");
       return(-1);
    }
      
    //Read in the header
         
    fread(&NtotGal, sizeof(int), 1, fd);     
         
    fclose(fd);
   
    printf("Number of Old Galaxies NtotGal = %d\n", NtotGal);
 
    return NtotGal;
}

void malOldgal_extrap(int totoldgals)
{   
    OldNumGalaxies = totoldgals;

    //Designate size of file for the array of structs
    if((OldGal = (struct Old_Gal_info*)malloc(OldNumGalaxies * sizeof(struct Old_Gal_info))) == NULL)
    {
        printf("Failed to allocate for OldGal...");
        return;
    }
  return;
}

void read_extrap(){

   int NtotGal,i;    
   char file1[1000];
   FILE *fd;
   
   sprintf(file1, "%s/extrap_restart/extrap_%d", All.OutputDir, All.SnapshotFileCount-2);
          
   printf("Reading extrap restart file = %s\n", file1);
     
   fd = fopen(file1, "rb");
      
   if(NULL == fd)
   {
      printf("Cannot open extrap restart file for reading");
      return;
   }
   
   //Read in the header
                                                                                                                                                          
   fread(&NtotGal, sizeof(int), 1, fd);
        
    
   //Read in info for each galaxy in file
   for(i=0;i<NtotGal;i++)
   {  
     fread(&OldGal[i], sizeof(struct Old_Gal_info), 1, fd);
   }
                                                                                                                                                      
   fclose(fd);
                                                                                                                                                      
      
   return;
}   



