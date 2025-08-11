/*
** @author Krista McCord, 2016
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>   

#include "proto.h"
#include "allvars.h"



void read_subf_multiple(int galindex)
{
   int tots;
   int i;

   tots = sub_readh_multiple();
   malsub_multiple(tots);

   printf("Reading all Subfind files\n");
   read_subfind_multiple_files();
  
   AllGal[galindex].group_firstsub = subf_multi[AllGal[galindex].subhaloidx].group_firstsub;  //gets the first subhalo in the group
   AllGal[galindex].group_len = subf_multi[AllGal[galindex].subhaloidx].group_len;  //gets the number of particles in the group
   AllGal[galindex].poffset = subf_multi[AllGal[galindex].subhaloidx].group_offset;
   AllGal[galindex].firstlen = subf_multi[AllGal[galindex].group_firstsub].sub_len; //number of particle in first sub halo
   //First halo should be the parent halo containing all subhalos.
   //


   
  
   //Find the offset for the halo in the id file.  This makes sure that the lengths match to get the correct offset.


   free(SubHead_all);
   free(subf_multi);
   SubHead_all = NULL;
   subf_multi = NULL; 
}

int sub_readh_multiple()
{
    int ngroups;
    int totngroups;
    int nids;
    int totnids;
    int ntask;
    int nsubs;
    int totnsubs;
    int tots;
    char fname[100];
    int num;

    num = All.SnapshotFileCount-1;


    FILE *sf;
    
    sprintf(fname, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "subhalo_tab", num, 0);
    sf = fopen(fname, "rb");

    //Read in the header
    fread(&ngroups, sizeof(unsigned int), 1, sf);
    fread(&totngroups, sizeof(unsigned int), 1, sf);
    fread(&nids, sizeof(unsigned int), 1, sf);
    fread(&totnids, sizeof(long long), 1, sf);
    fread(&ntask, sizeof(unsigned int), 1, sf);
    fread(&nsubs, sizeof(unsigned int), 1, sf);
    fread(&totnsubs, sizeof(unsigned int), 1, sf);

    fclose(sf);

    tots = totngroups + totnsubs;

    return tots;
    


}

void malsub_multiple(int tots)
{

    if((SubHead_all = (struct SubFind_Header_multiple*)malloc(NTask * sizeof(struct SubFind_Header_multiple))) == NULL)
    {
          perror("Failed to allocate for subfind header...");
          return;
    }

    if((subf_multi = (struct Subfind_multi_tab*)malloc(tots * sizeof(struct Subfind_multi_tab))) == NULL)
    {
          perror("Failed to allocate for subfind...");
          return;
    }

    return;
}





void read_subfind_multiple_files()
{
    int i=0;
    int gcount;
    int ngroups;
    int nsubs;
    int N;
    int a;
    FILE *sf;
    gcount = 0;
    N = 0;
    char fname[100];
    int num;
    int offset, offset2, offset3;
    
    num = All.SnapshotFileCount-1;

    

    offset = 0;
    offset2 = 0;
    offset3 = 0;
    
    for(a=0; a<NTask; a++)
    {
       //Open Subfind file
       sprintf(fname, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "subhalo_tab", num, a);
       sf = fopen(fname, "rb");


       //Read in the header
   
       fread(&SubHead_all[offset], sizeof(struct SubFind_Header_multiple), 1, sf);
            
       N = SubHead_all[offset].ngroups + SubHead_all[offset].nsubs;
       ngroups = SubHead_all[offset].ngroups;
       nsubs = SubHead_all[offset].nsubs;
    
       offset = offset + 1;

       //Designate size of file for the array of structs
       //struct subfind subf[nsubs + ngroups];
    
       for(i=0; i<ngroups; i++)
       {
          fread(&subf_multi[i+offset2].group_len, sizeof(unsigned int), 1, sf);
       }

       for(i=0; i<ngroups; i++)
       {
          fread(&subf_multi[i+offset2].group_offset, sizeof(unsigned int), 1,sf);
       }

       for(i=0; i<ngroups; i++)
       {
          fread(&subf_multi[i+offset2].group_mass, sizeof(float), 1,sf);
       }


    
       for(i=0; i<ngroups; i++)
       {
          for(a=0; a<3; a++)
          {
       
              fread(&subf_multi[i+offset2].group_pos[a], sizeof(float), 1, sf);
          }
       }

       for(i=0; i<ngroups; i++)
       {
          fread(&subf_multi[i+offset2].group_m_mean200, sizeof(float), 1, sf);
       }

       for(i=0; i<ngroups; i++)
       {
          fread(&subf_multi[i+offset2].group_r_mean200, sizeof(float), 1, sf);
       }

       for(i=0; i<ngroups; i++)
       {
          fread(&subf_multi[i+offset2].group_m_crit200, sizeof(float), 1, sf);
       }

       for(i=0; i<ngroups; i++)
       {
          fread(&subf_multi[i+offset2].group_r_crit200, sizeof(float), 1, sf);
       }

       for(i=0; i<ngroups; i++)
       {
          fread(&subf_multi[i+offset2].group_m_tophat200, sizeof(float), 1, sf);
       }

       for(i=0; i<ngroups; i++)
       {
          fread(&subf_multi[i+offset2].group_r_tophat200, sizeof(float), 1, sf);
       }

       for(i=0; i<ngroups; i++)
       {
          fread(&subf_multi[i+offset2].group_contamination_count, sizeof(unsigned int), 1, sf);
       }

       for(i=0; i<ngroups; i++)
       {
          fread(&subf_multi[i+offset2].group_contamination_mass, sizeof(float), 1, sf);
       }

       for(i=0; i<ngroups; i++)
       {
          fread(&subf_multi[i+offset2].group_nsubs, sizeof(unsigned int), 1, sf);
       }

       for(i=0; i<ngroups; i++)
       {
          fread(&subf_multi[i+offset2].group_firstsub, sizeof(unsigned int), 1, sf);
       }

       for(i=0; i<nsubs; i++)
       {
          fread(&subf_multi[i+offset3].sub_len, sizeof(unsigned int), 1, sf);
       }
   
       for(i=0; i<nsubs; i++)
       {
          fread(&subf_multi[i+offset3].sub_offset, sizeof(unsigned int), 1, sf);
       }

       for(i=0; i<nsubs; i++)
       {
          fread(&subf_multi[i+offset3].sub_parent, sizeof(unsigned int), 1, sf);
       }

       for(i=0; i<nsubs; i++)
       {
          fread(&subf_multi[i+offset3].sub_mass, sizeof(float), 1, sf);
       }

       for(i=0; i<nsubs; i++)
       {
          for(a=0; a<3; a++)
          {
             fread(&subf_multi[i+offset3].sub_pos[a], sizeof(float), 1, sf);
          }
       }

       for(i=0; i<nsubs; i++)
       {
          for(a=0; a<3; a++)
          {
             fread(&subf_multi[i+offset3].sub_vel[a], sizeof(float), 1, sf);
          }
       }

       for(i=0; i<nsubs; i++)
       {
          for(a=0; a<3; a++)
          {
             fread(&subf_multi[i+offset3].sub_cm[a], sizeof(float), 1, sf);
          }
       }

       for(i=0; i<nsubs; i++)
       {
          for(a=0; a<3; a++)
          {
             fread(&subf_multi[i+offset3].sub_spin[a], sizeof(float), 1, sf);
          }
       }

       for(i=0; i<nsubs; i++)
       {
          fread(&subf_multi[i+offset3].sub_veldisp, sizeof(float), 1, sf);
       }

       for(i=0; i<nsubs; i++)
       {
          fread(&subf_multi[i+offset3].sub_vmax, sizeof(float), 1, sf);
       }

       for(i=0; i<nsubs; i++)
       {
          fread(&subf_multi[i+offset3].sub_vmaxrad, sizeof(float), 1, sf);
       }

       for(i=0; i<nsubs; i++)
       {
          fread(&subf_multi[i+offset3].sub_halfmassrad, sizeof(float), 1, sf);
       }

       for(i=0; i<nsubs; i++)
       {
          fread(&subf_multi[i+offset3].sub_id_mostbound, sizeof(unsigned int), 1, sf);
       }

       for(i=0; i<nsubs; i++)
       {
          fread(&subf_multi[i+offset3].sub_grnr, sizeof(unsigned int), 1, sf);
       }
 
       offset2 = offset2 + ngroups;
       offset3 = offset3 + nsubs;
       fclose(sf);
    
    }    
    return;
}


