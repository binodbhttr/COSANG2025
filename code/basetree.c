#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>

#include "allvars.h"
#include "proto.h"
//#include "allvars_bt.h"
#include "proto_bt.h"


#define ALPHA (2.0/3)


#define WEIGHT_FAK (3.0)

#define int4bytes int
int4bytes blksize,swap=0;
#define SKIP  {my_fread(&blksize,sizeof(int),1,fd); swap_Nbyte_bt((char*)&blksize,1,4);}

void basetree(int snapnumber)
{

  
  WallTime = second_bt();
 
  FilesPerSnapshot = NTask;
  LastSnapShotNr = All.SnapshotFileCount-1;
  
  SnapshotNum = snapnumber;

  swap = 0;

#ifdef zoomsim

  totalpart = Ntype[1] + Ntype[5];
#else
  totalpart = Ntype[1];
#endif zoomsim
  printf("loading catalogues...\n"); fflush(stdout);

  load_subhalo_catalogue(SnapshotNum, &CatA);
  printf("CatA loaded\n"); fflush(stdout);
  time_readA=measure_time_bt();

  load_subhalo_catalogue(SnapshotNum + 1, &CatB);
  printf("CatB loaded\n"); fflush(stdout);
  time_readB=measure_time_bt();

  if(SnapshotNum+2 <= LastSnapShotNr)
    {
      load_subhalo_catalogue(SnapshotNum + 2, &CatC);
      printf("CatC loaded\n"); fflush(stdout);
      time_readC=measure_time_bt();
    }

  printf("done.\n"); fflush(stdout);
  time_read = time_readA + time_readB + time_readC;

  printf("preparing ID-to-halo tables...\n"); fflush(stdout);

  prepare_index_list(&CatB); 
  printf("index B done.\n"); fflush(stdout);
  time_prepB=measure_time_bt();

  if(SnapshotNum + 2 <= LastSnapShotNr)
    {
      prepare_index_list(&CatC);
      printf("index C done.\n"); fflush(stdout);
      time_prepC=measure_time_bt();
    }
  time_prep = time_prepA + time_prepB + time_prepC;
 

  printf("determine_descendants...\n");  fflush(stdout);
  determine_descendants(&CatA, &CatB, 0, SnapshotNum + 1);
  printf("desc AB done.\n"); fflush(stdout);
  time_findAB=measure_time_bt();


  if(SnapshotNum + 2 <= LastSnapShotNr)
    {
      determine_descendants(&CatB, &CatC, 0, SnapshotNum + 2);
      printf("desc BC done.\n"); fflush(stdout); 
      time_findBC=measure_time_bt();

      determine_descendants(&CatA, &CatC, 1, SnapshotNum + 2);  /* secondary descendant */
      printf("desc AC done.\n"); fflush(stdout);
      time_findAC=measure_time_bt();
    }

  printf("descendants done.\n"); fflush(stdout);
  time_find = time_findAB + time_findBA + time_findBC + time_findAC;


  if(SnapshotNum + 2 <= LastSnapShotNr)
    {
      printf("decide whether we should take secondary descendant...\n"); fflush(stdout);

      count_progenitors(&CatA, &CatB);
      printf("progcount AB done\n");  fflush(stdout);
      time_countAB=measure_time_bt();

      count_progenitors(&CatB, &CatC);
      printf("progcount BC done\n");  fflush(stdout);
      time_countBC=measure_time_bt();

      decide_upon_descendant();
      printf("decision made\n");  fflush(stdout);
      time_decide_desc=measure_time_bt();
    }

   time_count = time_countAB + time_countBC + time_countAB_back;
   time_decide = time_decide_desc + time_decide_back;

  printf("saving descendants...\n"); fflush(stdout);

  save_decendant_list();

  printf("saving done.\n"); fflush(stdout);
  time_write=measure_time_bt();


  time_total = time_read + time_prep + time_find + time_count + time_decide + time_write;

  printf("Timings:\n");
  printf("reading:             %9.2f sec  (%5.2f%%)\n",time_read,time_read/time_total*100);
  printf("   CatA:             %9.2f sec  (%5.2f%%)\n",time_readA,time_readA/time_total*100);
  printf("   CatB:             %9.2f sec  (%5.2f%%)\n",time_readB,time_readB/time_total*100);
  printf("   CatC:             %9.2f sec  (%5.2f%%)\n",time_readC,time_readC/time_total*100);
  printf("prepare index:       %9.2f sec  (%5.2f%%)\n",time_prep,time_prep/time_total*100);
  printf("   CatA:             %9.2f sec  (%5.2f%%)\n",time_prepA,time_prepA/time_total*100);
  printf("   CatB:             %9.2f sec  (%5.2f%%)\n",time_prepB,time_prepB/time_total*100);
  printf("   CatC:             %9.2f sec  (%5.2f%%)\n",time_prepC,time_prepC/time_total*100);
  printf("find decendants:     %9.2f sec  (%5.2f%%)\n",time_find,time_find/time_total*100);
  printf("   AB:               %9.2f sec  (%5.2f%%)\n",time_findAB,time_findAB/time_total*100);
  printf("   BA:               %9.2f sec  (%5.2f%%)\n",time_findBA,time_findBA/time_total*100);
  printf("   BC:               %9.2f sec  (%5.2f%%)\n",time_findBC,time_findBC/time_total*100);
  printf("   AC:               %9.2f sec  (%5.2f%%)\n",time_findAC,time_findAC/time_total*100);
  printf("count progenitors:   %9.2f sec  (%5.2f%%)\n",time_count,time_count/time_total*100);
  printf("decide:              %9.2f sec  (%5.2f%%)\n",time_decide,time_decide/time_total*100);
  printf("writing:             %9.2f sec  (%5.2f%%)\n",time_write,time_write/time_total*100);
  printf("---------------------------------------------\n");
  printf("total:               %9.2f sec \n",time_total);

  //free(&CatA);
  //*CatA = NULL;
  //free(&CatB);
  //*CatB = NULL;
  //free(&CatC);
 // *CatC = NULL;
 // free(&Cats);
 // Cats = NULL;
  return;
}

#define HALO_SIZE_INCREASE_FOR_SWITCHING 1.5


void decide_upon_descendant(void)
{
  int i, index_b, index_c;
  int count_b, count_c, count_w, count_n;
  double sumpart;

  count_b = count_c = count_w = count_n = 0;
  sumpart = 0.0;

  printf("Deciding decendants ...");
  fflush(stdout);

  for(i = 0; i < CatA.TotNsubhalos; i++)
    {
      index_b = CatA.Descendant[i].HaloIndex[0];
      index_c = CatA.Descendant[i].HaloIndex[1];

      if(index_b >= 0)
	count_b++;

      if(index_b >= 0 && index_c >= 0)
	{
	  if(CatB.CountProgenitors[index_b] > 1 && CatC.CountProgenitors[index_c] == 0)
	    {
	      CatB.CountProgenitors[index_b]--;
	      CatC.CountProgenitors[index_c]++;
	      CatA.Descendant[i].HaloIndex[0] =  CatA.Descendant[i].HaloIndex[1];
	      CatA.Descendant[i].SnapNum[0] =  CatA.Descendant[i].SnapNum[1];
	      count_c++;
	    }

	}

      if(index_b < 0 && index_c >= 0)
	{
	  CatA.Descendant[i].HaloIndex[0] =  CatA.Descendant[i].HaloIndex[1];
	  CatA.Descendant[i].SnapNum[0] =  CatA.Descendant[i].SnapNum[1];
	  CatC.CountProgenitors[index_c]++;
	  count_w++;
	}

      if(index_b < 0 && index_c < 0)
	{
	  /*
	    printf("len=%d\n", CatA.SubLen[i]);
	  */
	  sumpart += CatA.SubLen[i];
	  count_n++;
	}
    }
  printf("Out of %d primary descendants, %d have been rerouted to the secondary descendant.\n",
	 count_b, count_c);

  printf("Additionally, %d have been pointed to the secondary because they had no primary.\n",
	 count_w);

  printf("This leaves %d without descendant, of average size = %g particles.\n", count_n, sumpart/count_n);
  fflush(stdout);
}



void count_progenitors(struct halo_catalogue *catA, struct halo_catalogue *catB)
{
  int i;

  printf("Setting count to zero ...\n");
  fflush(stdout);


  for(i = 0; i < catB->TotNsubhalos; i++)
    catB->CountProgenitors[i] = 0;

  printf("Counting ...");
  fflush(stdout);

  for(i = 0; i < catA->TotNsubhalos; i++)
    {
      if(catA->Descendant[i].HaloIndex[0] >= 0)
	{
	  catB->CountProgenitors[catA->Descendant[i].HaloIndex[0]]++;
	}
    }
}


struct cand_data
{
    int haloindex;
    float weight;
};

struct index_table
{
  unsigned long long index;
  unsigned long long ID;
};


void determine_descendants(struct halo_catalogue *catA, struct halo_catalogue *catB, int entry, int snapnum)
{
  int i, j, ndiff, ncand, haloB, prev, maxlen;
  unsigned long long id;
  float weightmax;
  int halomax;
  struct cand_data *candlist, *difflist;

  maxlen = 0;

  printf("Finding maximum length ...\n");
  fflush(stdout);

  for(i = 0; i < catA->TotNsubhalos; i++)
    if(catA->SubLen[i] > maxlen)
      maxlen = catA->SubLen[i];

  printf("Determin decendants ...");
  fflush(stdout);

  {
    candlist = mymalloc_bt(maxlen * sizeof(struct cand_data));
    difflist = mymalloc_bt(maxlen * sizeof(struct cand_data));

    for(i = 0; i < catA->TotNsubhalos; i++)
      {
	ncand = 0;
	for(j = 0; j < catA->SubLen[i]; j++)
	  {
	    id = catA->IdList[catA->SubOffset[i] + j] - 1;

	    if(id>=0 && id<totalpart) 
	      {
		haloB = catB->IdToHalo[id];
	    
		if(haloB >= 0)
		  {
		    candlist[ncand].haloindex = haloB;
		    candlist[ncand].weight = 1.0/pow(j+1, ALPHA);
		    ncand++;
		  }
	      }
      	    else
	      {
		printf("bummer! i = %d, j = %d, id = %d, totalpart = %d\n", i, j, id, totalpart); //exit(4);
	      }
	  }
	
	qsort(candlist, ncand, sizeof(struct cand_data), sort_candlist);
	
	for(j = 0, ndiff = 0, prev = -1; j < ncand; j++)
	  {
	    if(candlist[j].haloindex != prev)
	      {
		ndiff++;
		difflist[ndiff-1].haloindex = candlist[j].haloindex;
		difflist[ndiff-1].weight = 0;
	      }
	    difflist[ndiff-1].weight += candlist[j].weight;
	    prev = candlist[j].haloindex;
	  }
	
	weightmax = 0;
	halomax = -1;
	
	for(j = 0; j < ndiff; j++)
	  {
	    if(difflist[j].weight > weightmax)
	      {
		weightmax = difflist[j].weight;
		halomax = difflist[j].haloindex;
	      }
	  }
	

	if(ndiff > 0 && halomax >=0)
	  {
	    catA->Descendant[i].HaloIndex[entry] = halomax;
	    catA->Descendant[i].SnapNum[entry] = snapnum;

	  }
	else
	  {
	    catA->Descendant[i].HaloIndex[entry] = -1;
	    catA->Descendant[i].SnapNum[entry] = -1;

	  }
      }
  }

 free(candlist);
 candlist = NULL;
 free(difflist);
 difflist = NULL;
}



int sort_candlist(const void *a, const void *b)
{
  if(((struct cand_data *)a)->haloindex < ((struct cand_data *)b)->haloindex)
    return -1;

  if(((struct cand_data *)a)->haloindex > ((struct cand_data *)b)->haloindex)
    return +1;

  return 0;
}

int sort_indextable(const void *a, const void *b)
{
  if(((struct index_table *)a)->ID < ((struct index_table *)b)->ID)
    return -1;

  if(((struct index_table *)a)->ID > ((struct index_table *)b)->ID)
    return +1;

  return 0;
}

void prepare_index_list(struct halo_catalogue *cat)
{
  unsigned long long id;
  int i, j;

  cat->IdToHalo = mymalloc_bt(sizeof(int) * totalpart);

  printf("Reset index ...\n");
  fflush(stdout);

  for(id = 0; id < totalpart; id++)
    cat->IdToHalo[id] = -1;

  printf("Fill index list ...\n");
  fflush(stdout);

  for(i = 0; i < cat->TotNsubhalos; i++)
    {
      for(j = 0; j < cat->SubLen[i]; j++)
	{
	  
        //  if (cat->SubOffset[i]+j > cat->TotNids) 
	//s{
	      // printf("Error Here %d %d %d %d\n",cat->SubOffset[i],j, cat->TotNids,  cat->TotNsubhalos);
		//fflush(stdout);
         // }
	  id = cat->IdList[cat->SubOffset[i]+j] - 1;
	  
          
	  if(id>=0 && id<totalpart)
	    cat->IdToHalo[id] = i;
	  else
	    {
	      printf("bummer! i=%d j=%d, (%lld %lld)\n", i, j,id,totalpart);// exit(1);
	    }
	}
    }

}


/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to swap ENDIAN -------------------------------*/
/*-------- char *data:    Pointer to the data ---------------------------------*/
/*-------- int n:         Number of elements to swap --------------------------*/
/*-------- int m:         Size of single element to swap ----------------------*/
/*--------                int,float = 4 ---------------------------------------*/
/*--------                double    = 8 ---------------------------------------*/
/*-----------------------------------------------------------------------------*/
void swap_Nbyte_bt(char *data,int n,int m)
{
  int i,j;
  char old_data[16];

  if(swap>0)
    {
      for(j=0;j<n;j++)
	{
          memcpy(&old_data[0],&data[j*m],m);
          for(i=0;i<m;i++)
            {
              data[j*m+i]=old_data[m-i-1];
	    }
	}
    }
}



void load_subhalo_catalogue(int num, struct halo_catalogue *cat)
{
  int i, ngroups, nids, nFiles, nsubhalos, subcount, filenr;
  long long idcount;
  char buf[1000];
  FILE *fd;
  int *suboffset_temp;


  unsigned int *IdListtmp;

  cat->TotNgroups = 0;
  cat->TotNids = 0;
  cat->TotNsubhalos = 0;



  for(i = 0; i < FilesPerSnapshot; i++)
    {

      sprintf(buf, "%sgroups_%03d/subhalo_tab_%03d.%d", All.OutputDir, num, num, i);

      if(!(fd = fopen(buf, "r")))
	{
	  printf("Tried to sum up total halos, can't open file `%s'\n", buf);
	  exit(1);
	}

      fread(&ngroups, sizeof(int), 1, fd);

      fread(&cat->TotNgroups, sizeof(int), 1, fd);

      fread(&nids, sizeof(int), 1, fd);

      fseek(fd, sizeof(long long), SEEK_CUR);      /*- Skip TotNids -*/

      fread(&nFiles, sizeof(int), 1, fd);
      fread(&nsubhalos, sizeof(int), 1, fd);

      fseek(fd, sizeof(int), SEEK_CUR);      /*- Skip TotNsubhalos -*/

      swap_Nbyte_bt((char*)&ngroups,1,4);
      swap_Nbyte_bt((char*)&nids,1,4);
      swap_Nbyte_bt((char*)&cat->TotNgroups,1,4);
      swap_Nbyte_bt((char*)&nFiles,1,4);
      swap_Nbyte_bt((char*)&nsubhalos,1,4);

      cat->TotNids += nids;
      cat->TotNsubhalos += nsubhalos;
      fclose(fd);

      fflush(stdout);
    }

  cat->IdList = mymalloc_bt(sizeof(unsigned long long) * cat->TotNids);
  cat->SubLen = mymalloc_bt(sizeof(int) * cat->TotNsubhalos);
  cat->SubOffset = mymalloc_bt(sizeof(long long) * cat->TotNsubhalos);
  cat->SubParentHalo = mymalloc_bt(sizeof(int) * cat->TotNsubhalos);
  cat->Descendant = mymalloc_bt(sizeof(struct descendant_data) * cat->TotNsubhalos);
  cat->CountProgenitors = mymalloc_bt(sizeof(int) * cat->TotNsubhalos);

  subcount = 0;
  idcount = 0;


  printf("Loading  %sgroups_%03d/subhalo_tab_XXX\n",  All.OutputDir, num);

  printf("cat->TotNsubhalos = %d\n", cat->TotNsubhalos);
  printf("cat->TotNids = %d%09d\n\n", 
	 (int) (cat->TotNids / 1000000000),
	 (int) (cat->TotNids % 1000000000));
  fflush(stdout);

  for(filenr = 0; filenr < FilesPerSnapshot; filenr++)
    {


      sprintf(buf, "%sgroups_%03d/subhalo_tab_%03d.%d", All.OutputDir, num, num, filenr);

      if(!(fd = fopen(buf, "r")))
	{
	  printf("Tried to read subfind output, can't open file `%s'\n", buf);
	  exit(1);
	}

      fread(&ngroups, sizeof(int), 1, fd);


      fread(&cat->TotNgroups, sizeof(int), 1, fd);

      fread(&nids, sizeof(int), 1, fd);

      fseek(fd, sizeof(long long), SEEK_CUR);      /*- Skip TotNids -*/


      fread(&nFiles, sizeof(int), 1, fd);
      fread(&nsubhalos, sizeof(int), 1, fd);

      fseek(fd, sizeof(int), SEEK_CUR);      /*- Skip TotNsubhalos -*/

      swap_Nbyte_bt((char*)&ngroups,1,4);
      swap_Nbyte_bt((char*)&nids,1,4);
      swap_Nbyte_bt((char*)&cat->TotNgroups,1,4);
      swap_Nbyte_bt((char*)&nFiles,1,4);
      swap_Nbyte_bt((char*)&nsubhalos,1,4);



      fseek(fd, sizeof(int) * ngroups, SEEK_CUR);       /*- Skip HaloLen -*/

      fseek(fd, sizeof(int) * ngroups, SEEK_CUR);       /*- Skip HaloMemberID -*/
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);     /*- Skip HaloMass -*/
      fseek(fd, sizeof(float) * 3 * ngroups, SEEK_CUR); /*- Skip HaloPos -*/
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);     /*- Skip Halo_M_Mean200 -*/
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);     /*- Skip Halo_R_Mean200 -*/
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);     /*- Skip Halo_M_Crit200 -*/
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);     /*- Skip Halo_R_Crit200 -*/
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);     /*- Skip Halo_M_TopHat200 -*/
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);     /*- Skip Halo_R_TopHat200 -*/


      fseek(fd, sizeof(int) * ngroups, SEEK_CUR);       /*- Skip HaloContCount -*/
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);     /*- Skip HaloContermination -*/


      /* skip NsubPerHalo and FirstSubOfHalo */
      fseek(fd, sizeof(int) * ngroups, SEEK_CUR);
      fseek(fd, sizeof(int) * ngroups, SEEK_CUR);
      /*
         fread(NsubPerHalo, sizeof(int), Ngroups, fd);
         fread(FirstSubOfHalo, sizeof(int), Ngroups, fd);
       */

      fread(&cat->SubLen[subcount], sizeof(int), nsubhalos, fd);

      suboffset_temp = mymalloc_bt(sizeof(int) * nsubhalos);
      fread(suboffset_temp, sizeof(int), nsubhalos, fd);
      
     /* for (i = 0; i < nsubhalos; i++) {
	printf("%d %d\n",suboffset_temp[i],nsubhalos);
        fflush(stdout);
	}*/  //DEBUG loop

      fread(&cat->SubParentHalo[subcount], sizeof(int), nsubhalos, fd);

      fclose(fd);

      if(swap>0)
	{
	  printf("Swapping Len/Offset/ParHal .,,\n");
	  fflush(stdout);

	  swap_Nbyte_bt((char*)&cat->SubLen[subcount],nsubhalos,4);
	  printf(" .,,\n");
	  fflush(stdout);

	  swap_Nbyte_bt((char*)suboffset_temp,nsubhalos,4);
	  printf(" .,,\n");
	  fflush(stdout);


	  swap_Nbyte_bt((char*)&cat->SubParentHalo[subcount],nsubhalos,4);


	  printf("Finished swapping .,,\n");
	  fflush(stdout);
	}

      //sprintf(buf, "%s/groups_%03d/group_ids_%03d.%d", All.OutputDir, num, num, filenr);
      sprintf(buf, "%s/groups_%03d/subhalo_ids_%03d.%d", All.OutputDir, num, num, filenr);

      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s'\n", buf);
	  exit(1);
	}
      else

      fread(&ngroups, sizeof(int), 1, fd);


      fread(&cat->TotNgroups, sizeof(int), 1, fd);

      fread(&nids, sizeof(int), 1, fd);

      fseek(fd, sizeof(long long), SEEK_CUR);      /*- Skip TotNids -*/

      fread(&nFiles, sizeof(int), 1, fd);

      fseek(fd, sizeof(int), SEEK_CUR);      /*- Skip TotNsubhalos -*/


      if(swap>0)
	printf("Swapping header ,,, %d %d %d %d\n",ngroups,nids,cat->TotNgroups,nFiles);
      fflush(stdout);
      swap_Nbyte_bt((char*)&ngroups,1,4);
      swap_Nbyte_bt((char*)&nids,1,4);
      swap_Nbyte_bt((char*)&cat->TotNgroups,1,4);
      swap_Nbyte_bt((char*)&nFiles,1,4);
      if(swap>0)
	printf("Finished swapping header ,,, %d %d %d %d\n",ngroups,nids,cat->TotNgroups,nFiles);
      fflush(stdout);


      IdListtmp = mymalloc_bt(sizeof(int) * nids);
      fread(IdListtmp, sizeof(int), nids, fd);
     // printf("id_0=%d id_n-1=%d\n",IdListtmp[0],IdListtmp[nids-1]);
      for (i=0;i<nids;i++)
	cat->IdList[idcount+i]=IdListtmp[i];
      free(IdListtmp);
      IdListtmp = NULL;
      if(swap>0)
	printf("Swapping IdList .,,\n");
      fflush(stdout);
      swap_Nbyte_bt((char*)&cat->IdList[idcount],nids,8);
      if(swap>0)
	printf("finished swapping .,,\n");
      fflush(stdout);



      fclose(fd);

      for(i = 0; i < nsubhalos; i++)
	{
	   cat->SubOffset[subcount + i] = suboffset_temp[i];
	}

      free(suboffset_temp);
      suboffset_temp = NULL;
      for(i = 0; i < nids; i++)
	cat->IdList[idcount + i] &= ((((long long) 1) << 34) - 1);


      subcount += nsubhalos;
      idcount += nids;

      fflush(stdout);
    }

    fflush(stdout);

}





void save_decendant_list(void)
{
  int i, *data;
  char buf[1000];
  FILE *fd;


  sprintf(buf, "%s/groups_%03d/sub_desc_%03d", All.OutputDir, SnapshotNum, SnapshotNum);

  if(!(fd = fopen(buf, "w")))
    {
      printf("can't open file `%s'\n", buf);
      exit(1);
    }
  
  fwrite(&CatA.TotNsubhalos, sizeof(int), 1, fd);

  data= mymalloc_bt(sizeof(int) * CatA.TotNsubhalos);

  for(i=0; i < CatA.TotNsubhalos; i++)
    data[i] = CatA.Descendant[i].HaloIndex[0];

  fwrite(data, sizeof(int), CatA.TotNsubhalos, fd);

  for(i=0; i < CatA.TotNsubhalos; i++)
    data[i] = CatA.Descendant[i].SnapNum[0];

  fwrite(data, sizeof(int), CatA.TotNsubhalos, fd);

  free(data);
  data = NULL;
  fclose(fd);
}


