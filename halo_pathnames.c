/*
** @author Krista McCord, 2016
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

void Treefiles()
{
    int numfiles;
   
    numfiles = malloc_treenames();
    NumTreefiles = numfiles;
   
    mal_treestrings(numfiles);
    read_ftreenames();
 
    return;
}

void read_ftreenames()
{
    int i=0;
    char line[1000];
    char file1[1000];
    FILE *fd; 
    
    sprintf(file1, "%s/treedata/Treefile_%03d.txt", All.OutputDir, All.SnapshotFileCount-1);     
         
    fd = fopen(file1, "r");

    while(EOF != fscanf(fd, "%s%*[^\n]", line))
    {
       strcpy(treepath[i].tpaths, line);
       i++;
    }

       
   fclose(fd);
   return; 

}

int malloc_treenames()
{
    int i=0;
    int n=0;
    char line[1000];
    char file1[1000];
    FILE *fd; 
    
    sprintf(file1, "%s/treedata/Treefile_%03d.txt", All.OutputDir, All.SnapshotFileCount-1);      
    fd = fopen(file1, "r");

    while(fgets(line, sizeof(line), fd)!=NULL) 
    {
       n += 1;
    }

    fclose(fd);
    return n;
}


void mal_treestrings(int numfiles)
{

    treepath = (struct Tree_Path_Names*)malloc(numfiles * sizeof(struct Tree_Path_Names));

    return;
}



