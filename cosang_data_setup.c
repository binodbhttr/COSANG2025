/*
** @author Krista McCord, 2016
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include <assert.h>


#include "allvars.h"
#include "proto.h"
#include "core_proto.h"


void cosang_data_setup(){

   int b, c, d, j, q;
        
   if(All.SnapshotFileCount > 2 && TotHalos != 0)
   {
      MPI_Barrier(MPI_COMM_WORLD);

                   
#ifdef profiler
      struct ProfileInstance profile_instance;
      struct ProfileInstance subf_instance;
      struct ProfileInstance plist_instance;

      enum Bucket BucketNames = Galax;
      enum Bucket BucketSubread = Subread;
      enum Bucket Bucketplist = plist;

      profile_instance_start(&profile_instance, BucketNames, __FUNCTION__, __FILE__, __LINE__);
#endif

#ifdef Extrap
      galsearchcount = 0;
      snap_oldgal_count = snap_oldgal_count + 1;
          
      if(snap_oldgal_count > 1)
      {            
         free(OldGal);
         OldGal = NULL;
  
         OldGal_create();
      }

#endif
          
      if(Particle_List_tracker == 1)
      {
         free(P_list);
         free(G_list);
         free(O_list);
         P_list = NULL;
         G_list = NULL;
         O_list = NULL;
      } 

      sage_setup();

#ifdef profiler
      profile_instance_end(&profile_instance, BucketNames);
      profile_instance_start(&subf_instance, BucketSubread, __FUNCTION__, __FILE__, __LINE__);
#endif
          
      MPI_Barrier(MPI_COMM_WORLD);
      pidcounter = 0;
          
      if(NumGalaxies != 0)
      {       
        for(j=0; j<NTask; j++)
        {
          if(ThisTask == j)
          {      
             single_subf_read(j);  //Each processor reads in 1 subfind_tab file
          }
        }

        if(ThisTask == 0)
        {
          printf("Get halo info from subfind\n");
        }
         
        Get_subfind_halo_info(); 
          
        assert(SubFindHeader != NULL);
 
        free(SubFindHeader);

        assert(subfind_tab != NULL);

        free(subfind_tab);  
        SubFindHeader = NULL;
        subfind_tab = NULL;


        //Fill OldGal struct
#ifdef Extrap
        if(snap_oldgal_count == 1)
        {          
          malOldgal();

          First_fill();
        }
#endif

          
//Calculate rotation angles for the disk orientation from the halo angular momentum
             
        for(j=0; j<NumGalaxies; j++)
        {
           rotation_angles(j);
        }
         
  
        if(ThisTask == 0) 
        {
           printf("Read subfind ids\n");
        }
   
        read_subf_ids();  //Each processor reads in 1 subfind_id file
           
#ifdef profiler
       profile_instance_end(&subf_instance, BucketSubread);
       profile_instance_start(&plist_instance, Bucketplist, __FUNCTION__, __FILE__, __LINE__);
#endif
       MPI_Barrier(MPI_COMM_WORLD);
   
       if(ThisTask == 0)
       {
          printf("Start creating particle lists\n");
       }
        
       create_plist();
            
       MPI_Barrier(MPI_COMM_WORLD);

       assert(Subf_pids != NULL);

       free(Subf_pids);

       assert(Subf_pids_header != NULL);

       free(Subf_pids_header);
       Subf_pids = NULL;
       Subf_pids_header = NULL; 

       Particle_List_tracker = 1;
   
#ifdef profiler
       profile_instance_end(&plist_instance, Bucketplist);
#endif
       if(ThisTask == 0)
       {
          printf("Finished with particle lists\n");
       }
   
      }
         
   }


return;

}


void run_trees_sage(){

    int i;

    if(ThisTask == 0 && All.SnapshotFileCount > 2)
    {
       for(i=All.SnapshotFileCount-3; i<All.SnapshotFileCount-1; i++)
       {
           basetree(i);
       }

       halotree();
       if(TotHalos !=0)
       { 
          int a;
          sage(Sageparam);
         
          free(ZZ);
          ZZ = NULL;
         
          free(AA);      
          AA = NULL;
          
          free(Age);
          Age = NULL;
       }
         
    }
     
    MPI_Bcast(&TotHalos, 1, MPI_INT, 0, MPI_COMM_WORLD); 
    MPI_Barrier(MPI_COMM_WORLD);

    return;

}

void sage_setup(){

    int c;

    if(ThisTask == 0)
    {     
      Galaxies();
                
      free(sagepath);
      sagepath = NULL;
    }
             
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(&NumGalaxies, 1, MPI_INT, 0, MPI_COMM_WORLD);

    printf("NumGalaxies = %d\n", NumGalaxies);

    if(ThisTask != 0)
    {
       malallgal(NumGalaxies);
    }

    MPI_Barrier(MPI_COMM_WORLD);
     
    for(c = 0; c < NumGalaxies; c++)
    {
       MPI_Bcast(&AllGal[c].Mvir, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
       MPI_Bcast(&AllGal[c].Pos[0], 1, MPI_FLOAT,0, MPI_COMM_WORLD);
       MPI_Bcast(&AllGal[c].Pos[1], 1, MPI_FLOAT,0, MPI_COMM_WORLD);
       MPI_Bcast(&AllGal[c].Pos[2], 1, MPI_FLOAT,0, MPI_COMM_WORLD);
       MPI_Bcast(&AllGal[c].Rvir, 1, MPI_FLOAT,0,MPI_COMM_WORLD);
       MPI_Bcast(&AllGal[c].Vel[0],1,MPI_FLOAT,0,MPI_COMM_WORLD);
       MPI_Bcast(&AllGal[c].Vel[1],1,MPI_FLOAT,0,MPI_COMM_WORLD);
       MPI_Bcast(&AllGal[c].Vel[2],1,MPI_FLOAT,0,MPI_COMM_WORLD);
       MPI_Bcast(&AllGal[c].DiskScaleRadius, 1, MPI_FLOAT,0,MPI_COMM_WORLD);
       MPI_Bcast(&AllGal[c].BulgeMass, 1, MPI_FLOAT,0,MPI_COMM_WORLD);
       MPI_Bcast(&AllGal[c].BlackHoleMass, 1, MPI_FLOAT,0,MPI_COMM_WORLD);
       MPI_Bcast(&AllGal[c].StellarMass, 1, MPI_FLOAT,0,MPI_COMM_WORLD);
       MPI_Bcast(&AllGal[c].ColdGas, 1, MPI_FLOAT,0,MPI_COMM_WORLD);
       MPI_Bcast(&AllGal[c].mostbound, 1, MPI_LONG_LONG_INT,0,MPI_COMM_WORLD);
       MPI_Bcast(&AllGal[c].subfnum, 1, MPI_INT,0,MPI_COMM_WORLD);
       MPI_Bcast(&AllGal[c].subhaloidx, 1, MPI_INT,0,MPI_COMM_WORLD);
       MPI_Bcast(&AllGal[c].MhaloIdx, 1, MPI_INT,0,MPI_COMM_WORLD);

    }
       

   
return;

}
