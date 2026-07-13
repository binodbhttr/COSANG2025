/*
** @author Krista McCord, 2016
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"


void Galaxies()
{
   int totmal;
   int tots;

   printf("reading sage file paths\n");
   sagefiles();
   
   
   totmal = read_sageheader();
   NumGalaxies = totmal;
   mal_gal(totmal);
    
   printf("Reading all Sage files\n");
   read_sage();
   

   malallgal(totmal);
  
   read_trees(totmal);
  

   fill_gal(totmal);
  
   

   free(treepath);
   
               
}

int read_sageheader()
{

    int Ntrees;
    int NtotGal;
    int totmal; 
    int i;
    FILE *fd;
    char file1[1000];
    totmal = 0;
  
    for(i=0; i<NumSfiles; i++)
    {
         sprintf(file1, "%s", sagepath[i].paths);
            

         //Open Sage file
        
         fd = fopen(file1, "rb");
        
         if(NULL == fd)
         {
            printf("Cannot open sage file\n");
            return(-1);
         }
      
         //Read in the header
         
         fread(&Ntrees, sizeof(int), 1, fd);
         fread(&NtotGal, sizeof(int), 1, fd);
         
         totmal = totmal + NtotGal;
         fclose(fd);
    }
  return totmal;
}

void mal_gal(int totmal)
{   
    //Designate size of file for the array of structs
    if((SageOutput = (struct sage_galaxies*)malloc(totmal * sizeof(struct sage_galaxies))) == NULL)
    {
        printf("Failed to allocate for sage...");
        return;
    }
  return;
}

void read_sage()
{  
   int Ntrees;
   int NtotGal;    
   int *galpertree;
   char file1[1000];
   int offset; 
   int i;
   FILE *fd;
   offset = 0;
   
   for(i=0; i<NumSfiles; i++)
   {
      sprintf(file1, "%s", sagepath[i].paths);
          
      printf("file path in read_sage = %s\n", file1);
      //Open Sage file
      fd = fopen(file1, "rb");
      
      if(NULL == fd)
      {
         printf("Cannot open sage for reading");
         return;
      }
   
      //Read in the header
                                                                                                                                                          
      fread(&Ntrees, sizeof(int), 1, fd);
      fread(&NtotGal, sizeof(int), 1, fd);
        
      //Part of the header allocate the needed array size
      galpertree = (int*)malloc(Ntrees*sizeof(int));
    
      //Read in the rest of the header into an array
      fread(&galpertree[0], sizeof(int), Ntrees, fd);
    
      //Read in info for each galaxy in file
      fread(&SageOutput[offset], sizeof(struct sage_galaxies), NtotGal, fd);
                                                                                                                                                      
      fclose(fd);
                                                                                                                                                      
      offset = offset + NtotGal;
                                                                                                                                                    }
      free(galpertree);  
      return;
}   


void malallgal(int totmal)
{   
    //Designate size of file for the array of structs
    if((AllGal = (struct All_Gal_info*)malloc(totmal * sizeof(struct All_Gal_info))) == NULL)
    {
        printf("Failed to allocate for AllGal...");
        return;
    }
  return;
}

void fill_gal(int totmal)
{
    int i;
    int a;

    for(i=0; i<totmal; i++)
    {
       AllGal[i].Mvir = SageOutput[i].Mvir;

       for(a=0; a<3; a++)
       {
            AllGal[i].Pos[a] = SageOutput[i].Pos[a];
            AllGal[i].Vel[a] = SageOutput[i].Vel[a] * All.Time;
       }
   
       AllGal[i].Rvir = SageOutput[i].Rvir;
       
    }
  free(SageOutput);
  return;
}
