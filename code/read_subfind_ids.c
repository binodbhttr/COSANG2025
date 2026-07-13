/*
** @author Krista McCord, 2016
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>   

#include "proto.h"
#include "allvars.h"

void read_subf_ids()
{
   unsigned int tots;
   int totids;
   int tag1;
   int idfilenum;

   idfilenum = ThisTask;
   tag1 =1;
 
 
   malloc_subfid_header();
   
   tots = read_subfid_header(idfilenum);


   malloc_subf_ids(tots);
  

   load_numids(idfilenum);

     if((Numids = (int *)malloc(NTask * sizeof(int))) == NULL)
     {
          perror("Failed to allocate for Array for number of ids in each file...");
          return;
     }

     if((sumNumids = (int *)malloc(NTask * sizeof(int))) == NULL)
     {
          perror("Failed to allocate for Array for number of ids in each file...");
          return;
     }



   if(ThisTask == 0 )
   {
     int i;
     int sumcount;
     int a;
     sumcount = 0;
     


     Numids[0] = Subf_pids_header[0].idnids;


     for(i=1; i<NTask; i++)
     {
    

        MPI_Recv(&totids, 1, MPI_INT, i, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
 
            
 
        Numids[i] = totids;
     }

     for(a=0; a<NTask; a++)
     {
       sumNumids[a] = Numids[a] + sumcount;
       sumcount = Numids[a] + sumcount; 

      
     }

   }


   if(ThisTask != 0)
   { 
      totids = Subf_pids_header[0].idnids;
      
      MPI_Send(&totids, 1, MPI_INT, 0, tag1, MPI_COMM_WORLD);
   }      


  MPI_Barrier(MPI_COMM_WORLD);
  
  MPI_Bcast(sumNumids, NTask, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(Numids, NTask, MPI_INT, 0, MPI_COMM_WORLD);
   
  return;
}


void malloc_subfid_header()
{
    if((Subf_pids_header = (struct SubFind_Ids_Header*)malloc(sizeof(struct SubFind_Ids_Header))) == NULL)
    {
          perror("Failed to allocate for ID header...");
          return;
    }

}


void malloc_subf_ids(unsigned int totmal)
{


   
    if((Subf_pids = (struct Subfind_ids*)malloc(totmal * sizeof(struct Subfind_ids))) == NULL)
    {
          perror("Failed to allocate for Subfind IDs...");
          return;
    }

    return;

}


unsigned int read_subfid_header(int idfilenum)
{
    int idngroups;
    int idtotngroups;
    int idnids;
    int idtotnids;
    int idntask;
    int idoffset;
    int tots;
    char fname[10000];
    int num;

    num = All.SnapshotFileCount-1; 

    tots = 0;
    FILE *sf;
    
    sprintf(fname, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "subhalo_ids", num, idfilenum);
    sf = fopen(fname, "rb");
    

    //Read in the header
    fread(&Subf_pids_header[0].idngroups, sizeof(unsigned int), 1, sf);
    fread(&Subf_pids_header[0].idtotngroups, sizeof(unsigned int), 1, sf);
    fread(&Subf_pids_header[0].idnids, sizeof(unsigned int), 1, sf);
    fread(&Subf_pids_header[0].idtotnids, sizeof(long long), 1, sf);
    fread(&Subf_pids_header[0].idntask, sizeof(unsigned int), 1, sf);
    fread(&Subf_pids_header[0].idoffset, sizeof(unsigned int), 1, sf);
   
    fclose(sf);

    tots = Subf_pids_header[0].idnids;

    return tots;

}


void load_numids(int idfilenum)
{
    char fname[10000];

    int num;
    int offset;
    int i;
    int a;
    num = All.SnapshotFileCount-1; 
    offset = 0;

    FILE *sf;
    
    
    sprintf(fname, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "subhalo_ids", num, idfilenum);
    sf = fopen(fname, "rb");
 

    //Read in the header
    fseek(sf, sizeof(unsigned int), SEEK_CUR);
    fseek(sf, sizeof(unsigned int), SEEK_CUR);
    fseek(sf, sizeof(unsigned int), SEEK_CUR);
    fseek(sf, sizeof(long long), SEEK_CUR);
    fseek(sf, sizeof(unsigned int), SEEK_CUR);
    fseek(sf, sizeof(unsigned int), SEEK_CUR);
            
    for(a=0; a<Subf_pids_header[0].idnids; a++)
    {
        fread(&Subf_pids[a].part_ids,sizeof(unsigned int), 1, sf);
           
             
    }

   
        fclose(sf);
     

    return;
}

void determine_ids(int galindex)
{
    char fname[5120];
    char fname2[10000];
    int num;
    int offset;
    int i, a, b, c, d, g;
    num = All.SnapshotFileCount-1;
    offset = 0;
    int idcount;
    int endid;
     int goffset;
    idcount = 0;
    endid = 0;
    goffset = 0;
    FILE *sf;
    FILE *hf;
    
    for(i=0; i<NTask; i++)
    {
              
 //Is the particle offset and number of particles within the range of the number of ids in this specific file
        if(offset < AllGal[galindex].poffset + AllGal[galindex].group_len &&  AllGal[galindex].poffset + AllGal[galindex].group_len <= Subf_pids_header[i].idnids + offset)
        {
            sprintf(fname, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "subhalo_ids", num, i);
      
            sf = fopen(fname, "rb");
        
           
            //Skip the header
            fseek(sf, sizeof(unsigned int), SEEK_CUR);
            fseek(sf, sizeof(unsigned int), SEEK_CUR);
            fseek(sf, sizeof(unsigned int), SEEK_CUR);
            fseek(sf, sizeof(long long), SEEK_CUR);
            fseek(sf, sizeof(unsigned int), SEEK_CUR);
            fseek(sf, sizeof(unsigned int), SEEK_CUR);
            
            for(a=0; a<AllGal[galindex].poffset-offset; a++)
            {
              
               fseek(sf, sizeof(unsigned int), SEEK_CUR);
            }
          
            //for loop is the range of lines to read in the file
            for(b=AllGal[galindex].poffset-offset; b<AllGal[galindex].poffset - offset + AllGal[galindex].group_len; b++)
            {
               fread(&Subf_pids[goffset].part_ids,sizeof(unsigned int), 1, sf);
             
               goffset++;
             }

          
	    fclose(sf);
         }
        offset = offset + Subf_pids_header[i].idnids;
        
             
 }
    
    return;

}
