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


void restart_cosang(){
    
     int j;
     
     sage_setup();

     
     MPI_Barrier(MPI_COMM_WORLD);
     pidcounter = 0;
       
     if(NumGalaxies != 0)
     {
       for(j=0; j<NTask; j++)
       {
          if(ThisTask == j)
          {
             single_subf_read(j);
          }
       }
          
       Get_subfind_halo_info();
          
       free(SubFindHeader);
       free(subfind_tab);
       SubFindHeader = NULL;
       subfind_tab = NULL;

       read_subf_ids();
          
       MPI_Barrier(MPI_COMM_WORLD);
         
       create_plist();
  
       MPI_Barrier(MPI_COMM_WORLD);

             
       for(j=0; j<NumGalaxies; j++)
       {
          rotation_angles(j);      
       }

       Particle_List_tracker = 1;

       free(Subf_pids);
       free(Subf_pids_header);
       Subf_pids = NULL;
       Subf_pids_header = NULL;
     }


  return;

}
