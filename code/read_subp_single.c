/*
** @author Krista McCord, 2016
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>   

#include "proto.h"
#include "allvars.h"



void single_subf_read(int subfile) 
{
    int tots;
    int shalo;
    
    tots = sub_readh(subfile);

    malsub(tots);
   
    read_subfind(subfile);
    
    subtotids = SubFindHeader[0].totnids;

    return;
}

int sub_readh(int subfilenum)
{
    int ngroups;
    int totngroups;
    int nids;
    int totnids;
    int ntask;
    int nsubs;
    int totnsubs;
    int tots;
    char fname[10000];
    int num;

    num = All.SnapshotFileCount-1; 
   
    FILE *sf;

    tots = 0;

    sprintf(fname, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "subhalo_tab", num, subfilenum);
   
    sf = fopen(fname, "rb");

    if(NULL == sf)
    {
       perror("Cannot open file to read header\n");
       return;
    }

    //Read in the header
    fread(&ngroups, sizeof(unsigned int), 1, sf);
    fread(&totngroups, sizeof(unsigned int), 1, sf);
    fread(&nids, sizeof(unsigned int), 1, sf);
    fread(&totnids, sizeof(long long), 1, sf);
    fread(&ntask, sizeof(unsigned int), 1, sf);
    fread(&nsubs, sizeof(unsigned int), 1, sf);
    fread(&totnsubs, sizeof(unsigned int), 1, sf);

    fclose(sf);

    tots = ngroups + nsubs;
   
    return tots;

}

void malsub(int tots)
{

    if((SubFindHeader = (struct SubFind_Header_data*)malloc(sizeof(struct SubFind_Header_data))) == NULL)
    {
          perror("Failed to allocate for subfind header...");
          return;
    }

    if((subfind_tab = (struct Subfind_data_tab*)malloc((tots+1) * sizeof(struct Subfind_data_tab))) == NULL)
    {
          perror("Failed to allocate for subfind...");
          return;
    }

    return;
}





void read_subfind(int subfnr)
{
    int i=0;
    int gcount;
    int ngroups;
    int ngroups2;
    int nsubs;
    int N;
    int a;
    FILE *sf;
    gcount = 0;
    N = 0;
    char fname[10000];
    int num;
    int tag1;
    ngroups = 0;
    ngroups2 = 0;
    gcount = 0;
    nsubs = 0;
    tag1 = 1;
   
    //Open Subfind file

    num = All.SnapshotFileCount-1; 


    if((NumGroups = (int *)malloc(NTask * sizeof(int))) == NULL)
    {
       
       printf("Failed to allocate for NumGroups");
       return;
    }
 

    if((SumNumgrps = (int *)malloc(NTask * sizeof(int))) == NULL)
    {
       
       printf("Failed to allocate for SumNumGrps");
       return;
    }


    sprintf(fname, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "subhalo_tab", num, subfnr);
   

    sf = fopen(fname, "rb");


    //Read in the header
    fread(&SubFindHeader[0].ngroups, sizeof(unsigned int), 1, sf);
    fread(&SubFindHeader[0].totngroups, sizeof(unsigned int), 1, sf);
    fread(&SubFindHeader[0].nids, sizeof(unsigned int), 1, sf);
    fread(&SubFindHeader[0].totnids, sizeof(long long), 1, sf);
    fread(&SubFindHeader[0].ntask, sizeof(unsigned int), 1, sf);
    fread(&SubFindHeader[0].nsubs, sizeof(unsigned int), 1, sf);
    fread(&SubFindHeader[0].totnsubs, sizeof(unsigned int), 1, sf);

            
    N = SubFindHeader[0].ngroups + SubFindHeader[0].nsubs;
    ngroups = SubFindHeader[0].ngroups;
   

    if(ThisTask == 0) 
    {
       NumGroups[0] = ngroups;
      
    }

    for(a=1; a<NTask; a++)
    {
       if(ThisTask == a && ThisTask !=0)
       {
          MPI_Send(&ngroups, 1, MPI_INT, 0, tag1, MPI_COMM_WORLD);
       }

       if(ThisTask == 0)
       {
 
           MPI_Recv(&ngroups2,1, MPI_INT, a, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
           NumGroups[a] = ngroups2;
        
       }
     }
    
   
     nsubs = SubFindHeader[0].nsubs;
    
    //Designate size of file for the array of structs
    
    for(i=0; i<ngroups; i++)
    {
       fread(&subfind_tab[i].group_len, sizeof(unsigned int), 1, sf);
   
    }

    for(i=0; i<ngroups; i++)
    {
       fread(&subfind_tab[i].group_offset, sizeof(unsigned int), 1,sf);
    }

    for(i=0; i<ngroups; i++)
    {
       fread(&subfind_tab[i].group_mass, sizeof(float), 1,sf);
    }

    
    for(i=0; i<ngroups; i++)
    {
        for(a=0; a<3; a++)
        {
       
              fread(&subfind_tab[i].group_pos[a], sizeof(float), 1, sf);
        }
    }

    for(i=0; i<ngroups; i++)
    {
       fread(&subfind_tab[i].group_m_mean200, sizeof(float), 1, sf);
    }

    for(i=0; i<ngroups; i++)
    {
       fread(&subfind_tab[i].group_r_mean200, sizeof(float), 1, sf);
    }

    for(i=0; i<ngroups; i++)
    {
       fread(&subfind_tab[i].group_m_crit200, sizeof(float), 1, sf);
    }

    for(i=0; i<ngroups; i++)
    {
       fread(&subfind_tab[i].group_r_crit200, sizeof(float), 1, sf);
    }

    for(i=0; i<ngroups; i++)
    {
       fread(&subfind_tab[i].group_m_tophat200, sizeof(float), 1, sf);
    }

    for(i=0; i<ngroups; i++)
    {
       fread(&subfind_tab[i].group_r_tophat200, sizeof(float), 1, sf);
    }

    for(i=0; i<ngroups; i++)
    {
       fread(&subfind_tab[i].group_contamination_count, sizeof(unsigned int), 1, sf);
    }

    for(i=0; i<ngroups; i++)
    {
       fread(&subfind_tab[i].group_contamination_mass, sizeof(float), 1, sf);
    }

    for(i=0; i<ngroups; i++)
    {
       fread(&subfind_tab[i].group_nsubs, sizeof(unsigned int), 1, sf);
    }

    for(i=0; i<ngroups; i++)
    {
       fread(&subfind_tab[i].group_firstsub, sizeof(unsigned int), 1, sf);
    }

    for(i=0; i<nsubs; i++)
    {
       fread(&subfind_tab[i].sub_len, sizeof(unsigned int), 1, sf);
    }
   
    for(i=0; i<nsubs; i++)
    {
       fread(&subfind_tab[i].sub_offset, sizeof(unsigned int), 1, sf);
    }

    for(i=0; i<nsubs; i++)
    {
       fread(&subfind_tab[i].sub_parent, sizeof(unsigned int), 1, sf);
    }

    for(i=0; i<nsubs; i++)
    {
       fread(&subfind_tab[i].sub_mass, sizeof(float), 1, sf);
    }

    for(i=0; i<nsubs; i++)
    {
       for(a=0; a<3; a++)
       {
           fread(&subfind_tab[i].sub_pos[a], sizeof(float), 1, sf);
       }
    }

    for(i=0; i<nsubs; i++)
    {
       for(a=0; a<3; a++)
       {
           fread(&subfind_tab[i].sub_vel[a], sizeof(float), 1, sf);
       }
    }

    for(i=0; i<nsubs; i++)
    {
       for(a=0; a<3; a++)
       {
           fread(&subfind_tab[i].sub_cm[a], sizeof(float), 1, sf);
       }
    }

    for(i=0; i<nsubs; i++)
    {
       for(a=0; a<3; a++)
       {
           fread(&subfind_tab[i].sub_spin[a], sizeof(float), 1, sf);
       }
    }

    for(i=0; i<nsubs; i++)
    {
       fread(&subfind_tab[i].sub_veldisp, sizeof(float), 1, sf);
    }

    for(i=0; i<nsubs; i++)
    {
       fread(&subfind_tab[i].sub_vmax, sizeof(float), 1, sf);
    }

    for(i=0; i<nsubs; i++)
    {
       fread(&subfind_tab[i].sub_vmaxrad, sizeof(float), 1, sf);
    }

    for(i=0; i<nsubs; i++)
    {
       fread(&subfind_tab[i].sub_halfmassrad, sizeof(float), 1, sf);
    }

    for(i=0; i<nsubs; i++)
    {
       fread(&subfind_tab[i].sub_id_mostbound, sizeof(unsigned int), 1, sf);
    }

    for(i=0; i<nsubs; i++)
    {
       fread(&subfind_tab[i].sub_grnr, sizeof(unsigned int), 1, sf);
    }
 

    fclose(sf);
    
    
    return;
}


