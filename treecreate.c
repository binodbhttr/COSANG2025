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


/* ============================================================================
   [STACK SIZE OVERRIDE FIX] - Added by Binod on July 11, 2026.
   To bypass queue daemon pbs_mom's hard stack limits (which block setrlimit),
   we implement an iterative DFS traversal using an explicit heap-allocated 
   stack array (walk_stack).
   
   To switch back to the original recursive implementation, comment out the
   definition of USE_ITERATIVE_DFS below.
   ============================================================================ */
#define USE_ITERATIVE_DFS 1

int CountUsed, CountSumUsed;

#ifdef USE_ITERATIVE_DFS
struct stack_frame {
  int i;
  int state;
  int p;
};

struct stack_frame *walk_stack = NULL;

void print_walk_stack_trace(int top, int TotHalos)
{
  int k;
  printf("  Stack Traceback (depth=%d):\n", top);
  for(k = top; k >= 0 && k > top - 20; k--)
    {
      printf("    Frame[%d]: i=%d state=%d p=%d (FirstProg=%d NextProg=%d)\n", 
	     k, walk_stack[k].i, walk_stack[k].state, walk_stack[k].p,
	     (walk_stack[k].i >= 0 && walk_stack[k].i < TotHalos) ? Halo[walk_stack[k].i].FirstProgenitor : -2,
	     (walk_stack[k].i >= 0 && walk_stack[k].i < TotHalos) ? Halo[walk_stack[k].i].NextProgenitor : -2);
    }
  fflush(stdout);
}

void walk_it(int start_i, int flag)
{
  int top = 0;

  walk_stack[top].i = start_i;
  walk_stack[top].state = 0;
  walk_stack[top].p = -1;

  while(top >= 0)
    {
      if(top >= TotHalos)
	{
	  printf("Task %d: ERROR walk_it: Stack overflow top=%d >= TotHalos=%d\n", ThisTask, top, TotHalos);
	  fflush(stdout);
	  endrun(1115);
	}

      int i = walk_stack[top].i;
      int state = walk_stack[top].state;

      if(state == 0)
	{
	  if(CountUsed >= TotHalos)
	    {
	      printf("Task %d: ERROR walk_it: CountUsed=%d exceeds TotHalos=%d\n", ThisTask, CountUsed, TotHalos);
	      fflush(stdout);
	      endrun(1116);
	    }
	  HaloAux[i].UsedFlag = 1;
	  HaloAux[i].TargetIndex = CountUsed;
	  HaloAux[CountUsed].Origin = i;

	  if(flag == 1)
	    HaloList[CountUsed] = Halo[i];

	  CountUsed++;

	  walk_stack[top].state = 1;
	  if(Halo[i].Descendant >= 0)
	    {
	      int desc = Halo[i].Descendant;
	      if(desc >= TotHalos)
		{
		  printf("Task %d: ERROR walk_it: Descendant index %d is out of bounds [0, %d). i=%d, top=%d\n", ThisTask, desc, TotHalos, i, top);
		  print_walk_stack_trace(top, TotHalos);
		  fflush(stdout);
		  endrun(1111);
		}
	      if(HaloAux[desc].UsedFlag == 0)
		{
		  top++;
		  walk_stack[top].i = desc;
		  walk_stack[top].state = 0;
		  walk_stack[top].p = -1;
		  continue;
		}
	    }
	}
      else if(state == 1)
	{
	  walk_stack[top].state = 2;
	  walk_stack[top].p = Halo[i].FirstProgenitor;
	}

      if(walk_stack[top].state == 2)
	{
	  int p = walk_stack[top].p;
	  while(p >= 0)
	    {
	      if(p >= TotHalos)
		{
		  printf("Task %d: ERROR walk_it: Progenitor index %d is out of bounds [0, %d). i=%d, top=%d\n", ThisTask, p, TotHalos, i, top);
		  print_walk_stack_trace(top, TotHalos);
		  fflush(stdout);
		  endrun(1112);
		}
	      if(HaloAux[p].UsedFlag == 0)
		{
		  walk_stack[top].p = Halo[p].NextProgenitor;
		  top++;
		  walk_stack[top].i = p;
		  walk_stack[top].state = 0;
		  walk_stack[top].p = -1;
		  break;
		}
	      p = Halo[p].NextProgenitor;
	    }
	  if(walk_stack[top].state == 0)
	    continue;

	  walk_stack[top].state = 3;
	  int first_fof = Halo[i].FirstHaloInFOFgroup;
	  if(first_fof >= 0)
	    {
	      if(first_fof >= TotHalos)
		{
		  printf("Task %d: ERROR walk_it: first_fof index %d is out of bounds [0, %d). i=%d, top=%d\n", ThisTask, first_fof, TotHalos, i, top);
		  print_walk_stack_trace(top, TotHalos);
		  fflush(stdout);
		  endrun(1113);
		}
	      if(HaloAux[first_fof].HaloFlag == 0)
		{
		  HaloAux[first_fof].HaloFlag = 1;
		  walk_stack[top].p = first_fof;
		}
	      else
		{
		  walk_stack[top].p = -1;
		}
	    }
	  else
	    {
	      walk_stack[top].p = -1;
	    }
	}

      if(walk_stack[top].state == 3)
	{
	  int p = walk_stack[top].p;
	  while(p >= 0)
	    {
	      if(p >= TotHalos)
		{
		  printf("Task %d: ERROR walk_it: FOF index %d is out of bounds [0, %d). i=%d, top=%d\n", ThisTask, p, TotHalos, i, top);
		  print_walk_stack_trace(top, TotHalos);
		  fflush(stdout);
		  endrun(1114);
		}
	      if(HaloAux[p].UsedFlag == 0)
		{
		  walk_stack[top].p = Halo[p].NextHaloInFOFgroup;
		  top++;
		  walk_stack[top].i = p;
		  walk_stack[top].state = 0;
		  walk_stack[top].p = -1;
		  break;
		}
	      p = Halo[p].NextHaloInFOFgroup;
	    }
	  if(walk_stack[top].state == 0)
	    continue;

	  top--;
	}
    }

}
#else
/* Original recursive implementation of walk_it */
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
#endif






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

   if (Halo == NULL || HaloAux == NULL)
     {
       printf("Failed to allocate Halo or HaloAux in halotree()\n");
       exit(1);
     }

   /* ============================================================================
      [HALO MEMORY INITIALIZATION FIX] - Added by Antigravity on July 11, 2026.
      Because some snapshots contain subhalos that do not belong to any FOF groups
      (e.g., filtered out small groups), the catalog loading loop over FOF groups
      can finish early, leaving "gap" halos at the end of snapshot ranges
      unpopulated. Since malloc doesn't zero memory, these gap halos contain
      random heap garbage, which corrupts progenitor linking and causes crashes.
      
      We explicitly initialize the entire Halo array to safe default values here
      so that any unpopulated gap halos are cleanly disconnected (pointers = -1).
      ============================================================================ */
   int init_idx;
   for(init_idx = 0; init_idx < TotHalos; init_idx++)
     {
       Halo[init_idx].Descendant = -1;
       Halo[init_idx].FirstProgenitor = -1;
       Halo[init_idx].NextProgenitor = -1;
       Halo[init_idx].FirstHaloInFOFgroup = -1;
       Halo[init_idx].NextHaloInFOFgroup = -1;
       Halo[init_idx].Len = 0;
       
       Halo[init_idx].M_Mean200 = 0.0;
       Halo[init_idx].M_Crit200 = 0.0;
       Halo[init_idx].M_TopHat = 0.0;
       Halo[init_idx].VelDisp = 0.0;
       Halo[init_idx].Vmax = 0.0;
       Halo[init_idx].MostBoundID = 0;
       
       Halo[init_idx].SnapNum = -1;
       Halo[init_idx].FileNr = -1;
       Halo[init_idx].SubhaloIndex = -1;
       Halo[init_idx].SubhalfMass = 0.0;
     }

   memset(HaloAux, 0, TotHalos * sizeof(struct halo_aux_data));
#ifdef USE_ITERATIVE_DFS
  walk_stack = malloc(TotHalos * sizeof(struct stack_frame));
  if(walk_stack == NULL)
    {
      printf("Failed to allocate walk_stack in halotree()\n");
      exit(1);
    }
#endif


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

#ifdef USE_ITERATIVE_DFS
  if(walk_stack)
    {
      free(walk_stack);
      walk_stack = NULL;
    }
#endif

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
	/* ============================================================================
	   [GAP ROOT SAFETY CHECK]
	   Added on July 12, 2026, 21:47:12-07:00.
	   
	   Do not select unpopulated gap subhalos (FirstHaloInFOFgroup < 0) 
	   as tree roots, to prevent SAGE from constructing empty trees.
     // if(HaloAux[i].UsedFlag == 0)
	   ============================================================================ */
	if(HaloAux[i].UsedFlag == 0 && Halo[i].FirstHaloInFOFgroup >= 0)
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
	  /* ============================================================================
	     [GAP ROOT SAFETY CHECK]
	     Added on July 12, 2026, 21:47:12-07:00.
	     
	     Do not select unpopulated gap subhalos (FirstHaloInFOFgroup < 0) 
	     as tree roots, to prevent SAGE from constructing empty trees.
	     ============================================================================ */
	  if(HaloAux[i].UsedFlag == 0 && Halo[i].FirstHaloInFOFgroup >= 0)
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
	      if(descendant_haloindex[sc] >= 0 && descendant_snapnum[sc] >= 0 && descendant_snapnum[sc] <= LastSnapShotNr)
		{
		  int desc_idx = FirstHaloInSnap[descendant_snapnum[sc]] + descendant_haloindex[sc];
		  if(desc_idx >= 0 && desc_idx < TotHalos)
		    {
		      /* ============================================================================
			 [GAP DESCENDANT SAFETY CHECK]
			 Added on July 12, 2026, 21:30:29-07:00.
			 
			 SAGE pre-initializes unpopulated gap halos (subhalos listed in catalogs 
			 but filtered out of FOF groups) with FirstHaloInFOFgroup = -1. 
			 If an active halo lists one of these unpopulated gap halos as its descendant, 
			 we must NOT link them. Doing so would cause SAGE to walk into the gap halo, 
			 load FirstHaloInFOFgroup = -1, and crash with ABORT(54) in evolve_galaxies().
			 
			 We verify that the descendant halo has a valid FirstHaloInFOFgroup index (>= 0). 
			 If it does not, we set Descendant = -1, treating the progenitor as disrupted/merged, 
			 which is clean and safe.
			 ============================================================================ */
		      if(Halo[desc_idx].FirstHaloInFOFgroup >= 0)
			{
			  Halo[nh].Descendant = desc_idx;
			}
		      else
			{
			  Halo[nh].Descendant = -1;
			}
		    }
		  else
		    {
		      printf("Task %d: load_subhalo_catalogue_ht ERROR: nh=%d desc=%d out of bounds [0, %d). num=%d descendant_snapnum=%d descendant_haloindex=%d\n", 
			     ThisTask, nh, desc_idx, TotHalos, num, descendant_snapnum[sc], descendant_haloindex[sc]);
		      fflush(stdout);
		      Halo[nh].Descendant = -1;
		    }
		}
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
	  if(desc >= TotHalos)
	    {
	      printf("\nTask %d: WARNING descendant index %d is out of bounds [0, %d). i=%d. Setting Descendant to -1\n", ThisTask, desc, TotHalos, i);
	      fflush(stdout);
	      Halo[i].Descendant = -1;
	      continue;
	    }
	  if((first = Halo[desc].FirstProgenitor) >= 0)
	    {
	      if(first >= TotHalos)
		{
		  printf("\nTask %d: WARNING FirstProgenitor index %d is out of bounds [0, %d). desc=%d. Skipping descendant update\n", ThisTask, first, TotHalos, desc);
		  fflush(stdout);
		  continue;
		}
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
