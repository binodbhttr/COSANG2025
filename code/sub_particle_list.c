/*
** @author Krista McCord, 2016
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>   

#include "proto.h"
#include "allvars.h"



/*void per_pid_list()
{
  
  unsigned n_local_idxs, *local_idxs;
  int ii, jj, kk;

  n_local_idxs =  totpartids(); 
   
  local_idxs = ALLOC(unsigned, n_local_idxs);

  for( ii = 0, kk = 0; ii < NumGalaxies; ++ii )
   {
      for( jj = AllGal[ii].poffset; jj < AllGal[ii].poffset + AllGal[ii].group_len; ++jj, ++kk )
         local_idxs[kk] = jj;
     // halos_displs[ii + 1] = halos_displs[ii] + halos_data[2*ii + 1] - halos_data[2*ii];
   }


  permute(subtotids, n_local_idxs, local_idxs, &Subf_pids.part_ids, MPI_UNSIGNED, MPI_COMM_WORLD ); 

  

}*/

void create_pid_list(int galindex)
{
  int a;
  int idfile;
  int startidx;
  int tag1, tag2, tag3, tag4, tag5, tag6, tag7;
  int numread;
  int numleft;
  int search;
 
  int idfile2;
  //unsigned int *Temp_pid_list1, *Temp_pid_list2;
  idfile2 = 0;
  numread = 0;
  numleft = 0;
  tag1 = 1;
  tag2 = 2;
  tag3 = 3;
  tag4 = 4;
  tag5 = 5;
  tag6 = 6;
  tag7 = 7;
 

  /*  if((Temp_pid_list = (unsigned int *)malloc((AllGal[galindex].group_len + 1) * sizeof(unsigned int))) == NULL)
    {
          perror("Failed to allocate for Temp particle list...");
          return;
    }*/

 // printf("Finished allocating memory for temp list\n");
 // printf("pid test = %d\n", Subf_pids[0].part_ids);


  if(ThisTask == 0)
  {  
     int i;
     int c;
     int readoffset;
     
     readoffset = 0;
     c = 0;
     search = 0;
     //int b;
     //int idfile;
     printf("poffset = %d, galindex = %d\n", AllGal[galindex].poffset, galindex);
     if(AllGal[galindex].poffset < sumNumids[0])
     {
       idfile = 0;

       for(a=0; a<AllGal[galindex].group_len; a++)
       {
          particle_list[a + pidcounter] = Subf_pids[AllGal[galindex].poffset + a].part_ids;
       }
       search = 1;

       printf("NumGalaxies =%d, galindex =%d\n", NumGalaxies, galindex);
      
       //pidcounter = pidcounter + AllGal[galindex].group_len;
     }
     else
     {
      for(i=1; i<NTask; i++)
      {

           if(AllGal[galindex].poffset >= sumNumids[i-1] && AllGal[galindex].poffset < sumNumids[i] && AllGal[galindex].poffset + AllGal[galindex].group_len < sumNumids[i])
        {
           idfile = i;

           startidx = abs(AllGal[galindex].poffset - sumNumids[i-1]);
           printf("Send startidx = %d, galindex = %d, file = %d, sumNumids = %d\n", startidx, galindex, idfile, sumNumids[i-1]);
           MPI_Send(&startidx, 1, MPI_INT, idfile, tag1, MPI_COMM_WORLD);

           search = 1;
        }

       /* if(AllGal[galindex].poffset >= sumNumids[i-1] && AllGal[galindex].poffset < sumNumids[i] && sumNumids[i] - AllGal[galindex].poffset + AllGal[galindex].group_len < Numids[i])
        {
           idfile = i;
           
           startidx = sumNumids[i] - AllGal[galindex].poffset;
           printf("Send startidx = %d, galindex = %d\n", startidx, galindex);
           MPI_Send(&startidx, 1, MPI_INT, idfile, tag1, MPI_COMM_WORLD);
       
           search = 1; 
        }*/
        
      }
     

     if(search == 0)
     {
     //  printf("Test c = %d\n", c);
       for(c=0; c<NTask-1; c++)
       {
         // printf("Test loop\n");
         // printf("poffset %d, sumNumids %d, Num to read %d, Numids %d, c %d \n", AllGal[galindex].poffset, sumNumids[c],  sumNumids[c] - AllGal[galindex].poffset + AllGal[galindex].group_len,  Numids[c], c);  
         // if(AllGal[galindex].poffset <= sumNumids[c] && sumNumids[c] - AllGal[galindex].poffset + AllGal[galindex].group_len >= Numids[c] && readoffset < AllGal[galindex].group_len)
         
          if(AllGal[galindex].poffset <= sumNumids[c] && AllGal[galindex].poffset + AllGal[galindex].group_len >= sumNumids[c] && readoffset < AllGal[galindex].group_len)

          {
               idfile = c;
               idfile2 = c+1;

               if(c!=0)
               {

                   startidx = abs(AllGal[galindex].poffset - sumNumids[c-1]);
          //         printf("Send 1st startidx = %d, galindex = %d\n", startidx, galindex);
               }
               else
               {
                   startidx = AllGal[galindex].poffset;
               }

               numread = abs(Numids[c] - startidx);
               MPI_Send(&startidx, 1, MPI_INT, idfile, tag3, MPI_COMM_WORLD);
               MPI_Send(&numread, 1, MPI_INT, idfile, tag4, MPI_COMM_WORLD);
           
               numleft = abs(AllGal[galindex].group_len - numread);
               readoffset = numread + numleft;
               MPI_Send(&numleft, 1, MPI_INT, idfile2, tag5, MPI_COMM_WORLD);    
          }
       }
     }
   }
   /* if(search = 1)
    {

       MPI_Recv(&Temp_pid_list, AllGal[galindex].group_len, MPI_UNSIGNED, idfile, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     
       for(b=0; b<AllGal[galindex].group_len; b++)
       {
          printf("Temp_pid_list = %d\n", Temp_pid_list[b]);
          particle_list[b + pidcounter] = Temp_pid_list[b];
       }
    }*/
   }

  //MPI_Barrier(MPI_COMM_WORLD);
  //printf("Bcast idfile\n");

   MPI_Bcast(&idfile, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&search, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&idfile2, 1, MPI_INT, 0, MPI_COMM_WORLD);

   printf("search %d, idfile2 %d, idfile %d galindex %d\n", search, idfile2, idfile, galindex);
  
   if(search == 1 && idfile != 0 && ThisTask == idfile)
   {
     printf("Before Receive: proc %d, galindex = %d\n", idfile, galindex);
     MPI_Recv(&startidx, 1, MPI_INT, 0, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     printf("Receive startidx = %d, proc %d\n", startidx, idfile);

     malloc_temp_pid(AllGal[galindex].group_len);

     for(a=0; a<AllGal[galindex].group_len; a++)
     {
        Temp_pid_list[a] = Subf_pids[startidx + a].part_ids;
      //  printf("Subf_pids = %d, group_len = %d\n", Subf_pids[startidx +a].part_ids, AllGal[galindex].group_len);
     }

     MPI_Send(Temp_pid_list, AllGal[galindex].group_len, MPI_UNSIGNED, 0, tag2, MPI_COMM_WORLD);

    // free(Temp_pid_list);
   }
  else if(search == 0 && ThisTask == idfile)
  {
     //int testarr; //* testarr[5] = {0,1,2,3,4};

     MPI_Recv(&startidx, 1, MPI_INT, 0, tag3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     MPI_Recv(&numread, 1, MPI_INT, 0, tag4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

     malloc_temp_list1(numread);

     for(a=0; a<numread; a++)
     {
       Temp_pid_list1[a] = Subf_pids[startidx + a].part_ids;
       //printf("Subf_pids = %d, numread = %d\n", Temp_pid_list1[a], numread);
     }

     //testarr = 5;
     //MPI_Send(&testarr, 1, MPI_INT, 0, tag2, MPI_COMM_WORLD);

     MPI_Send(Temp_pid_list1, numread, MPI_UNSIGNED, 0, tag6, MPI_COMM_WORLD);

     //free(Temp_pid_list1);
  }
  
  if(search == 0 && ThisTask == idfile2 && idfile2 != 0)
  {
    MPI_Recv(&numleft, 1, MPI_INT, 0, tag5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   
    malloc_temp_list2(numleft);

    for(a=0; a<numleft; a++)
    {
      Temp_pid_list2[a] = Subf_pids[a].part_ids;
      //printf("Subf_pids = %d, numleft = %d\n", Temp_pid_list2[a], numleft);
    }

    MPI_Send(Temp_pid_list2, numleft, MPI_UNSIGNED, 0, tag7, MPI_COMM_WORLD);
    //printf("1. Sending temp pid list\n");

   // free(Temp_pid_list2);
  }  

   MPI_Barrier(MPI_COMM_WORLD);

 
   if(ThisTask == 0 && search == 1 && idfile != 0)
   {
     int b;
    
     
     malloc_temp_pid(AllGal[galindex].group_len);

     printf("Receive Temp pid list\n");
     MPI_Recv(Temp_pid_list, AllGal[galindex].group_len, MPI_UNSIGNED, idfile, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

     for(b=0; b<AllGal[galindex].group_len; b++)
     {
      // printf("Temp_pid_list = %d\n", Temp_pid_list[b]);
       particle_list[b + pidcounter] = Temp_pid_list[b];
     }

   // free(Temp_pid_list);
   }
  

   else if(search == 0 && ThisTask == 0 && idfile2 != 0)
   {
     int b,i;
     int d;
     int testarr;     
     malloc_temp_list1(numread);
     malloc_temp_list2(numleft);

    // printf("1. Receiving pid list idfile %d\n", idfile);
   
    // MPI_Recv(&testarr, 1, MPI_INT, idfile, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     
    // printf("testarr = %d\n", testarr);
     /*for(i=0; i<5; i++)
     {
       printf( "testarr = %d\n", testarr[i]);
     }*/
 
     MPI_Recv(Temp_pid_list1, numread, MPI_UNSIGNED, idfile, tag6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     printf("Receiving list numread = %d, Temp_pid_list %d, numleft = %d, galindex = %d\n", numread, Temp_pid_list1[0], numleft, galindex);
     for(b=0; b<numread; b++)
     {
      // printf(" b %d\n", b);
       particle_list[b + pidcounter] = Temp_pid_list1[b];
      // printf("part list 1st = %d\n", Temp_pid_list1[b]);
     }

     MPI_Recv(Temp_pid_list2, numleft, MPI_UNSIGNED, idfile2, tag7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

     for(d=0; d<numleft; d++)
     {
       //printf("d %d\n", d);
       particle_list[d+pidcounter+numread] = Temp_pid_list2[d];
       //printf("part list 2nd = %d\n", Temp_pid_list2[d]);
     }

    //free(Temp_pid_list1);
   // free(Temp_pid_list2);
   }
  
  

  pidcounter = pidcounter + AllGal[galindex].group_len;
  
 if(ThisTask == 0)
 {
   printf("pidcounter %d, galindex = %d\n", pidcounter, galindex);  
 }

  Num_Part_list[galindex] = AllGal[galindex].group_len;
 
/* if(ThisTask == 0)
  {
    printf("pidcounter = %d, Part_list = %d, galindex = %d\n", pidcounter, Num_Part_list[galindex], galindex);
  }*/

 MPI_Barrier(MPI_COMM_WORLD);

  if(search == 0 )
  {
    free(Temp_pid_list1);
    free(Temp_pid_list2);

    Temp_pid_list1 = NULL;
    Temp_pid_list2 = NULL;
  }

  if(search == 1 && idfile != 0)
  {
    free(Temp_pid_list);  
  
    Temp_pid_list = NULL;
  }

  //free(Temp_pid_list);
  return;
}


void create_pid_arrays()
{
  int totalids;

  totalids = 0;

  totalids = tot_partids();
  //printf("totalids = %d\n", totalids);
  malloc_pidlist(totalids);

  return;
}


void malloc_temp_pid(int len)
{
  if((Temp_pid_list = (unsigned int *)malloc(len * sizeof(unsigned int))) == NULL)
  {
    //perror("Failed to allocate for temp pid list");
    printf("Failed to allocate for temp pid list");
    return;
  }

  return;
}

void malloc_temp_list1(int len)
{
  if((Temp_pid_list1 = (unsigned int *)malloc(len * sizeof(unsigned int))) == NULL)
  {
    //perror("Failed to allocate for temp pid list 1");
    printf("Failed to allocate for temp pid list 1");
    return;
  }

  return;
}

void malloc_temp_list2(int len)
{
  if((Temp_pid_list2 = (unsigned int *)malloc(len * sizeof(unsigned int))) == NULL)
  {
    // perror("Failed to allocate for temp pid list 2");
    printf("Failed to allocate for temp pid list2");
    return;
  }

  return;
}

void malloc_pidlist(int totids)
{
    printf("totids = %d, NumGalaxies = %d\n", totids, NumGalaxies);

    if((particle_list = (unsigned int *)malloc((totids) * sizeof(unsigned int))) == NULL)
    {
          perror("Failed to allocate for particle list...");
          return;
    }

    if((Num_Part_list = (int *)malloc(NumGalaxies * sizeof(int))) == NULL)
    {
          perror("Failed to allocate Num_Part_list...");
          return;
    }

  return;
}


int tot_partids()
{
   int totids;
   int i;

   totids = 0;

   //printf("NumGalaxies = %d\n", NumGalaxies);

   for(i=0; i<NumGalaxies; i++)
    {
       totids = totids + AllGal[i].group_len;
       //printf("totids = %d, grouplen = %d, i = %d \n", totids, AllGal[i].group_len, i);
    }

  return totids;

}
