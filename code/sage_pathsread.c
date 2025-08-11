/*
** @author Krista McCord, 2016
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"



void sagefiles()
{
    int numfiles;
    
    numfiles = malloc_sagenames();
    NumSfiles = numfiles;
    printf("NumSfiles in function = %d\n", NumSfiles);
    mal_sagestrings(numfiles);
    read_fnames();
   
    return;
}
    


void read_fnames()
{
    int i=0;
    char line[10000];
    char file1[10000];
    FILE *fd; 

    sprintf(file1, "%s/sage_out/sagefile_%03d.txt", All.OutputDir, All.SnapshotFileCount-1);     
         
    fd = fopen(file1, "r");
	if(NULL == fd)
         {
            printf("Cannot open sage file in read\n");
            return(-1);
         }

    while(EOF != fscanf(fd, "%s%*[^\n]", line))
    {
       strcpy(sagepath[i].paths, line);
      
       i++;
    }

       
   fclose(fd);
   return; 

}

int malloc_sagenames()
{
    int i=0;
    int n=0;
    char line[10000];
    char file1[10000];
    FILE *fd; 

    sprintf(file1, "%s/sage_out/sagefile_%03d.txt", All.OutputDir, All.SnapshotFileCount-1);     
         
    fd = fopen(file1, "r");
	if(NULL == fd)
         {
            printf("Cannot open sage file in malloc\n");
            return(-1);
         }

   
    while(fgets(line, sizeof(line), fd)!=NULL) 
    {
       n += 1;
    }

    fclose(fd);
    return n;
}


void mal_sagestrings(int numfiles)
{

    sagepath = (struct Sage_Path_Names*)malloc(numfiles * sizeof(struct Sage_Path_Names));

    return;
}



