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

//#define MY_DEBUG
#define int4bytes int
//On Anvil got error with blksize (repeated definition, first def in basetree,c. So let's remove it/
int4bytes swap_ht=0;
//int4bytes blksize,swap_ht=0;
#define SKIP  {my_fread(&blksize,sizeof(int),1,fd); swap_Nbyte_ht((char*)&blksize,1,4);}

struct halo_data
{
  int Descendant;
  int FirstProgenitor;
  int NextProgenitor;
  int FirstHaloInFOFgroup;
  int NextHaloInFOFgroup;

  /* properties of halo */

  int Len;
  float M_Mean200, M_Crit200, M_TopHat;
  float Pos[3];
  float Vel[3];
  float VelDisp;
  float Vmax;
  float Spin[3];
  unsigned long long MostBoundID;


  /* original position in subfind output */

  int SnapNum, FileNr, SubhaloIndex;
  float SubhalfMass;
#ifdef SAVE_MASS_TAB
  float SubMassTab[6];
#endif
}
 *Halo, *HaloList;

struct halo_aux_data
{
  int UsedFlag;
  int HaloFlag;
  int TargetIndex;
  int Origin;
}
 *HaloAux;


int CountUsed, CountSumUsed;

void walk_it(int i, int flag)
{
  int p;

  HaloAux[i].UsedFlag = 1;
  HaloAux[i].TargetIndex = CountUsed;
  HaloAux[CountUsed].Origin = i;

  if(flag == 1)
    HaloList[CountUsed] = Halo[i];

  CountUsed++;

  if(Halo[i].Descendant >= 0)
    {
      if(HaloAux[Halo[i].Descendant].UsedFlag == 0)
	walk_it(Halo[i].Descendant, flag);
    }

  p = Halo[i].FirstProgenitor;
  while(p >= 0)
    {
      if(HaloAux[p].UsedFlag == 0)
	walk_it(p, flag);

      p = Halo[p].NextProgenitor;
    }

  p = Halo[i].FirstHaloInFOFgroup;
  if(HaloAux[p].HaloFlag == 0)
    {
      HaloAux[p].HaloFlag = 1;
      while(p >= 0)
	{
	  if(HaloAux[p].UsedFlag == 0)
	    walk_it(p, flag);
	  p = Halo[p].NextHaloInFOFgroup;
	}
    }
}





void halotree()
{
  int num, count;
  count = 0;
  FirstSnapShotNr = 0;
  
  fileset = LastSnapShotNr;

  Cats = mymalloc_ht(sizeof(struct halo_catalogue) * (LastSnapShotNr + 1));
  FirstHaloInSnap = mymalloc_ht(sizeof(int) * (LastSnapShotNr + 1));

  count_halos();
  if(TotHalos != 0)
  {
  printf("TotHalos = %d\n", TotHalos);
  printf("Want to allocate %g GB\n",
	 ((double) TotHalos) * (sizeof(struct halo_data) +
				sizeof(struct halo_aux_data)) / (1024.0 * 1024.0 * 1024.0));

  Halo = malloc(TotHalos * sizeof(struct halo_data));
  HaloAux = malloc(TotHalos * sizeof(struct halo_aux_data));


  printf("loading halo catalogues and descendant tree files...\n");
  fflush(stdout);

  for(num = LastSnapShotNr, count = 0; num >= FirstSnapShotNr; num--)
    {
      FirstHaloInSnap[num] = count;

      load_subhalo_catalogue_ht(num);

      count += Cats[num].TotNsubhalos;
    }

  printf("done.\n");


  set_progenitor_pointers();

  printf("progenitor pointers done.\n");
  fflush(stdout);

  generate_trees();
  }

  reset_vars_ht();
  free(Cats);
  Cats = NULL;
  free(Halo);
  Halo = NULL;
  free(HaloAux);
  HaloAux = NULL;
  free(FirstHaloInSnap);
  FirstHaloInSnap = NULL;
  return;
}



void generate_trees(void)
{
  int i, k, maxhalos, filenr, treenr;
  int *NtreesPerFile, *NhalosPerFile, *npertree;
  char buf[500];
  FILE *fd;

  NtreesPerFile = mymalloc_ht(sizeof(int) * FilesPerSnapshot);
  NhalosPerFile = mymalloc_ht(sizeof(int) * FilesPerSnapshot);

  for(i = 0; i < FilesPerSnapshot; i++)
    {
      NtreesPerFile[i] = 0;
      NhalosPerFile[i] = 0;
    }

  CountSumUsed = 0;

  for(i = 0; i < TotHalos; i++)
    HaloAux[i].UsedFlag = HaloAux[i].HaloFlag = 0;

  maxhalos = 0;

  for(filenr = 0; filenr < FilesPerSnapshot; filenr++)
    for(i = 0; i < Cats[LastSnapShotNr].TotNsubhalos; i++)
      {
	if(HaloAux[i].UsedFlag == 0)
	  {
	    if(filenr == whichfile(Halo[i].Pos))
	      {
		CountUsed = 0;

		walk_it(i, 0);

		NtreesPerFile[filenr] += 1;
		NhalosPerFile[filenr] += CountUsed;

		if(CountUsed > maxhalos)
		  maxhalos = CountUsed;

		CountSumUsed += CountUsed;
	      }
	  }
      }

  printf("TotHalos=%d   Used=%d   maxhalos=%d\n", TotHalos, CountSumUsed, maxhalos);
  fflush(stdout);

  for(i = 0; i < TotHalos; i++)
    HaloAux[i].UsedFlag = HaloAux[i].HaloFlag = 0;


  HaloList = mymalloc_ht(maxhalos * sizeof(struct halo_data));

  sprintf(buf, "%s/treedata", All.OutputDir);
  mkdir(buf, 02755);


  for(filenr = 0; filenr < FilesPerSnapshot; filenr++)
    {
      sprintf(buf, "%s/treedata/trees_%03d.%d", All.OutputDir, LastSnapShotNr, filenr);

      printf("starting: %s\n", buf);
      fflush(stdout);

      if(!(fd = fopen(buf, "w")))
	{
	  printf("can't open file `%s'\n", buf);
	  exit(1);
	}

      fwrite(&NtreesPerFile[filenr], 1, sizeof(int), fd);
      fwrite(&NhalosPerFile[filenr], 1, sizeof(int), fd);

      fseek(fd, NtreesPerFile[filenr] * sizeof(int), SEEK_CUR);

      npertree = mymalloc_ht(NtreesPerFile[filenr] * sizeof(int));
      for(i = 0; i < NtreesPerFile[filenr]; i++)
	npertree[i] = 0;

      treenr = 0;

      for(i = 0; i < Cats[LastSnapShotNr].TotNsubhalos; i++)
	{
	  if(HaloAux[i].UsedFlag == 0)
	    {
	      if(filenr == whichfile(Halo[i].Pos))
		{
		  CountUsed = 0;

		  walk_it(i, 1);

		  for(k = 0; k < CountUsed; k++)
		    {
		      if(HaloList[k].Descendant >= 0)
			HaloList[k].Descendant = HaloAux[HaloList[k].Descendant].TargetIndex;

		      if(HaloList[k].FirstProgenitor >= 0)
			HaloList[k].FirstProgenitor = HaloAux[HaloList[k].FirstProgenitor].TargetIndex;

		      if(HaloList[k].NextProgenitor >= 0)
			HaloList[k].NextProgenitor = HaloAux[HaloList[k].NextProgenitor].TargetIndex;

		      if(HaloList[k].FirstHaloInFOFgroup >= 0)
			HaloList[k].FirstHaloInFOFgroup =
			  HaloAux[HaloList[k].FirstHaloInFOFgroup].TargetIndex;

		      if(HaloList[k].NextHaloInFOFgroup >= 0)
			HaloList[k].NextHaloInFOFgroup = HaloAux[HaloList[k].NextHaloInFOFgroup].TargetIndex;
		    }

		  fwrite(HaloList, CountUsed, sizeof(struct halo_data), fd);

		  npertree[treenr] = CountUsed;
		  treenr++;
		}
	    }
	}

      fclose(fd);

      if(!(fd = fopen(buf, "r+")))
	{
	  printf("can't open file `%s'\n", buf);
	  exit(1);
	}
      fseek(fd, 2 * sizeof(int), SEEK_SET);
      fwrite(npertree, NtreesPerFile[filenr], sizeof(int), fd);
      fclose(fd);

      myfree_ht(npertree);
    }

  myfree_ht(HaloList);
  HaloList = NULL;
  printf("Saved=%d\n", CountSumUsed);
}



/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to swap ENDIAN -------------------------------*/
/*-------- char *data:    Pointer to the data ---------------------------------*/
/*-------- int n:         Number of elements to swap --------------------------*/
/*-------- int m:         Size of single element to swap ----------------------*/
/*--------                int,float = 4 ---------------------------------------*/
/*--------                double    = 8 ---------------------------------------*/
/*-----------------------------------------------------------------------------*/
void swap_Nbyte_ht(char *data,int n,int m)
{
  int i,j;
  char old_data[16];

  if(swap_ht>0)
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


void count_halos(void)
{
  int num, nFiles, nsubhalos, filenr;
  char buf[1000];
  FILE *fd;

  int ngroups, nids;

  TotHalos = 0;

  printf("Counting Halos: ");
  fflush(stdout);
  for(num = LastSnapShotNr; num >= FirstSnapShotNr; num--)
    {
      Cats[num].TotNsubhalos = 0;
      Cats[num].TotNgroups = 0;
      printf("%d ",num);
      fflush(stdout);
      for(filenr = 0; filenr < FilesPerSnapshot; filenr++)
	{
          sprintf(buf, "%s/groups_%03d/subhalo_tab_%03d.%d", All.OutputDir, num, num, filenr);
          if(!(fd = fopen(buf, "r")))
	    {
	      printf("can't open file `%s'\n", buf);
	      exit(1);
	    }

	  fread(&ngroups, sizeof(int), 1, fd);

	  fread(&Cats[num].TotNgroups, sizeof(int), 1, fd);

	  fread(&nids, sizeof(int), 1, fd);

	  fseek(fd, sizeof(long long), SEEK_CUR);      /*- Skip TotNids -*/

	  fread(&nFiles, sizeof(int), 1, fd);
	  fread(&nsubhalos, sizeof(int), 1, fd);

	  fseek(fd, sizeof(int), SEEK_CUR);      /*- Skip TotNsubhalos -*/
	  fclose(fd);

	  Cats[num].TotNsubhalos += nsubhalos;
	}

      TotHalos += Cats[num].TotNsubhalos;
    }

  printf("\n total number of halos=%d\n", TotHalos);

}





void load_subhalo_catalogue_ht(int num)
{
  int i, ngroups, nids, nFiles, nsubhalos, subcount;
  int groupcount, filenr, ncount;
  int subgr, gr, nh, sc, gr_nh;
  char buf[1000];
  FILE *fd;
  int *nsubPerHalo, *subLen, *descendant_haloindex, *descendant_snapnum, *filenrOfHalo, *subhaloindex;
  float *halo_M_Mean200, *halo_M_Crit200, *halo_M_TopHat;
  float *subpos, *subvel, *subveldisp, *subvmax, *subspin, *subhalfmass;
  unsigned long long *subMostBoundID;
  unsigned int *IdListtmp;

  printf("Catalogue num=%d (%d,%d)\n", num, Cats[num].TotNgroups, Cats[num].TotNsubhalos);

  nsubPerHalo = mymalloc_ht(sizeof(int) * Cats[num].TotNgroups);
  subLen = mymalloc_ht(sizeof(int) * Cats[num].TotNsubhalos);
  descendant_haloindex = mymalloc_ht(sizeof(int) * Cats[num].TotNsubhalos);
  descendant_snapnum = mymalloc_ht(sizeof(int) * Cats[num].TotNsubhalos);
  filenrOfHalo = mymalloc_ht(sizeof(int) * Cats[num].TotNsubhalos);
  subhaloindex = mymalloc_ht(sizeof(int) * Cats[num].TotNsubhalos);

  halo_M_Mean200 = mymalloc_ht(sizeof(float) * Cats[num].TotNgroups);
  halo_M_Crit200 = mymalloc_ht(sizeof(float) * Cats[num].TotNgroups);
  halo_M_TopHat = mymalloc_ht(sizeof(float) * Cats[num].TotNgroups);

  subpos = mymalloc_ht(3 * sizeof(float) * Cats[num].TotNsubhalos);
  subvel = mymalloc_ht(3 * sizeof(float) * Cats[num].TotNsubhalos);
  subveldisp = mymalloc_ht(sizeof(float) * Cats[num].TotNsubhalos);
  subvmax = mymalloc_ht(sizeof(float) * Cats[num].TotNsubhalos);
  subspin = mymalloc_ht(3 * sizeof(float) * Cats[num].TotNsubhalos);
  subMostBoundID = mymalloc_ht(sizeof(unsigned long long) * Cats[num].TotNsubhalos);
  subhalfmass = mymalloc_ht(sizeof(float) * Cats[num].TotNsubhalos);

  subcount = 0;
  groupcount = 0;

  for(filenr = 0; filenr < FilesPerSnapshot; filenr++)
    {

      sprintf(buf, "%s/groups_%03d/subhalo_tab_%03d.%d", All.OutputDir, num, num, filenr);

      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s'\n", buf);
	  exit(1);
	}


      /*----------------- Header: --------------------------------------------*/
      fread(&ngroups, sizeof(int), 1, fd);
      fread(&Cats[num].TotNgroups, sizeof(int), 1, fd);

      fread(&nids, sizeof(int), 1, fd);

      fseek(fd, sizeof(long long), SEEK_CUR);      /*- Skip TotNids -*/

      fread(&nFiles, sizeof(int), 1, fd);
      fread(&nsubhalos, sizeof(int), 1, fd);

      fseek(fd, sizeof(int), SEEK_CUR);      /*- Skip TotNsubhalos -*/


      /*----------------- Halos: --------------------------------------------*/


      fseek(fd, sizeof(int) * ngroups, SEEK_CUR);       /*- Skip HaloLen -*/
      fseek(fd, sizeof(int) * ngroups, SEEK_CUR);       /*- Skip HaloMemberID -*/
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);     /*- Skip HaloMass -*/
      fseek(fd, sizeof(float) * 3 * ngroups, SEEK_CUR); /*- Skip HaloPos -*/
      fread(&halo_M_Mean200[groupcount], sizeof(float), ngroups, fd);
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skipp R200 */
      fread(&halo_M_Crit200[groupcount], sizeof(float), ngroups, fd);
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skipp R200 */
      fread(&halo_M_TopHat[groupcount], sizeof(float), ngroups, fd);
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skipp R200 */
      fseek(fd, sizeof(int) * ngroups, SEEK_CUR);       /*- Skip HaloContCount -*/
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);     /*- Skip HaloContermination -*/

      fread(&nsubPerHalo[groupcount], sizeof(int), ngroups, fd);

      /*  fread(firstSubOfHalo[groupcount], sizeof(int), ngroups, fd); */
      fseek(fd, sizeof(int) * ngroups, SEEK_CUR);


      /*----------------- SubHalos: --------------------------------------------*/

      fread(&subLen[subcount], sizeof(int), nsubhalos, fd);  /*- read SubLen -*/

      /*
         fread(&subOffset[subcount], sizeof(int), nsubhalos, fd);
         fread(&subParentHalo[subcount], sizeof(int), nsubhalos, fd);
       */
      fseek(fd, sizeof(int) * nsubhalos, SEEK_CUR);          /*- Skip Offset -*/

      fseek(fd, sizeof(int) * nsubhalos, SEEK_CUR);          /*- Skip GrNr -*/

      fseek(fd, sizeof(float) * nsubhalos, SEEK_CUR);	/* skipp TMass */

      fread(&subpos[3 * subcount], 3 * sizeof(float), nsubhalos, fd);
      fread(&subvel[3 * subcount], 3 * sizeof(float), nsubhalos, fd);

      fseek(fd, sizeof(float) * 3 * nsubhalos , SEEK_CUR);	/* skipp SubCM */
      fread(&subspin[3 * subcount], 3 * sizeof(float), nsubhalos, fd);

      fread(&subveldisp[subcount], sizeof(float), nsubhalos, fd);
      fread(&subvmax[subcount], sizeof(float), nsubhalos, fd);

      fseek(fd, sizeof(float) * nsubhalos , SEEK_CUR);	/* skipp SubVmaxRad */
      fread(&subhalfmass[subcount], sizeof(float), nsubhalos, fd);




// ORIGINAL
// fread(&subMostBoundID[subcount], sizeof(long long), nsubhalos, fd);
// CHANGED TO

      //fread(&subMostBoundID[subcount], sizeof(long long), nsubhalos, fd);
	  IdListtmp = mymalloc_ht(sizeof(unsigned int) * nsubhalos);
	 // n = read_gadget_uint(IdListtmp,"MBID", fd);
      	  fread(IdListtmp, sizeof(unsigned int), nsubhalos, fd);
	  for (i=0;i<nsubhalos;i++) 
	    subMostBoundID[subcount+i]=IdListtmp[i];
	  myfree_ht(IdListtmp);
          IdListtmp = NULL;
// END CHANGE

      fseek(fd, sizeof(int) * nsubhalos , SEEK_CUR);	/* skipp SubGroupNumber */

      fclose(fd);


      for(subgr = 0; subgr < nsubhalos; subgr++)
	filenrOfHalo[subcount + subgr] = filenr;

      for(subgr = 0; subgr < nsubhalos; subgr++)
	subhaloindex[subcount + subgr] = subgr;

      subcount += nsubhalos;
      groupcount += ngroups;

    }


  if(num < LastSnapShotNr)
    {
      sprintf(buf, "%s/groups_%03d/sub_desc_%03d", All.OutputDir, num, num);

      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s'\n", buf);
	  exit(1);
	}

      fread(&ncount, sizeof(int), 1, fd);
      fread(descendant_haloindex, sizeof(int), Cats[num].TotNsubhalos, fd);
      fread(descendant_snapnum, sizeof(int), Cats[num].TotNsubhalos, fd);

      fclose(fd);
    }

  nh = FirstHaloInSnap[num];
  sc = 0;

  for(gr = 0; gr < Cats[num].TotNgroups; gr++)
    {
      for(subgr = 0, gr_nh = nh; subgr < nsubPerHalo[gr]; subgr++, sc++, nh++)
	{
	  /*
          if(subgr>=Cats[num].TotNsubhalos || gr>=Cats[num].TotNgroups || nh >=TotHalos)
             printf("Local sub: %d/%d halo: %d/%d, Global: %d/%d\n",subgr,Cats[num].TotNsubhalos,gr,Cats[num].TotNgroups,nh,TotHalos);
	  */
	  Halo[nh].FirstHaloInFOFgroup = gr_nh;
	  if(subgr == nsubPerHalo[gr] - 1)
	    Halo[nh].NextHaloInFOFgroup = -1;
	  else
	    Halo[nh].NextHaloInFOFgroup = nh + 1;

	  if(num < LastSnapShotNr)
	    {
	      if(descendant_haloindex[sc] >= 0)
		Halo[nh].Descendant = FirstHaloInSnap[descendant_snapnum[sc]] + descendant_haloindex[sc];
	      else
		Halo[nh].Descendant = -1;
	    }
	  else
	    Halo[nh].Descendant = -1;

	  Halo[nh].FirstProgenitor = -1;
	  Halo[nh].NextProgenitor = -1;

	  /* assign properties */

	  Halo[nh].Len = subLen[sc];

	  if(subgr == 0)
	    {
	      Halo[nh].M_Mean200 = halo_M_Mean200[gr];
	      Halo[nh].M_Crit200 = halo_M_Crit200[gr];
	      Halo[nh].M_TopHat = halo_M_TopHat[gr];
	    }
	  else
	    {
	      Halo[nh].M_Mean200 = 0;
	      Halo[nh].M_Crit200 = 0;
	      Halo[nh].M_TopHat = 0;
	    }


	  for(i = 0; i < 3; i++)
	    {
	      Halo[nh].Pos[i] = subpos[3 * sc + i];
	      Halo[nh].Vel[i] = subvel[3 * sc + i];
	      Halo[nh].Spin[i] = subspin[3 * sc + i];
	    }
	  Halo[nh].VelDisp = subveldisp[sc];
	  Halo[nh].Vmax = subvmax[sc];
	  Halo[nh].MostBoundID = subMostBoundID[sc];


	  /* store position of halo in subfind output */

	  Halo[nh].SnapNum = num;
	  Halo[nh].FileNr = filenrOfHalo[sc];
	  Halo[nh].SubhaloIndex = subhaloindex[sc];
	  Halo[nh].SubhalfMass = subhalfmass[sc];

	  /* auxiliary stuff */

	  HaloAux[nh].UsedFlag = 0;
	}
    }


  for(gr = 0; gr < nh; gr++)
    {
      if(Halo[gr].NextHaloInFOFgroup == gr)
	{
	  printf("bummer! %d\n", gr);
	}
    }


  myfree_ht(subhalfmass);
  subhalfmass = NULL;
  myfree_ht(subMostBoundID);
  subMostBoundID = NULL;
  myfree_ht(subspin);
  subspin = NULL;
  myfree_ht(subvmax);
  subvmax = NULL;
  myfree_ht(subveldisp);
  subveldisp = NULL;
  myfree_ht(subvel);
  subvel = NULL;
  myfree_ht(subpos);
  subpos = NULL;

  myfree_ht(halo_M_TopHat);
  halo_M_TopHat = NULL;
  myfree_ht(halo_M_Crit200);
  halo_M_Crit200 = NULL;
  myfree_ht(halo_M_Mean200);
  halo_M_Mean200 = NULL;
  myfree_ht(subhaloindex);
  subhaloindex = NULL;
  myfree_ht(filenrOfHalo);
  filenrOfHalo = NULL;
  myfree_ht(descendant_snapnum);
  descendant_snapnum = NULL;
  myfree_ht(descendant_haloindex);
  descendant_haloindex = NULL;
  myfree_ht(subLen);
  subLen = NULL;
  myfree_ht(nsubPerHalo);
  nsubPerHalo = NULL;

}


void set_progenitor_pointers(void)
{
  int i, first, desc, k = 0;

  printf("Setting Progenitors, %d steps:",TotHalos/10000);

  for(i = 0; i < TotHalos; i++)
    {
      k++;
      if(k>10000)
	{
	  printf(".");
          k=0;
	}
      if((desc = Halo[i].Descendant) >= 0)
	{
	  if((first = Halo[desc].FirstProgenitor) >= 0)
	    {
	      if(Halo[i].Len >= Halo[first].Len)
		{
		  Halo[i].NextProgenitor = first;
		  Halo[desc].FirstProgenitor = i;
		}
	      else
		{
		  Halo[i].NextProgenitor = Halo[first].NextProgenitor;
		  Halo[first].NextProgenitor = i;
		}
	    }
	  else
	    {
	      Halo[desc].FirstProgenitor = i;
	    }
	}
    }
  printf("\n");
}
