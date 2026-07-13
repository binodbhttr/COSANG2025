#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>
#include <mpi.h>

#include "allvars.h"
#include "core_proto.h"
#include "proto.h"

char bufz0[1000];
int exitfail = 1;

struct sigaction saveaction_XCPU;
volatile sig_atomic_t gotXCPU = 0;



void termination_handler(int signum)
{
  gotXCPU = 1;
  sigaction(SIGXCPU, &saveaction_XCPU, NULL);
  if(saveaction_XCPU.sa_handler != NULL)
    (*saveaction_XCPU.sa_handler) (signum);
}



void myexit(int signum)
{
  //printf("Task: %d\tnode: %s\tis exiting.\n\n\n", ThisTask, ThisNode);
  exit(signum);
}



void bye()
{
  //MPI_Finalize();
  //free(ThisNode);

  if(exitfail)
  {
    unlink(bufz0);

    if(ThisTask == 0 && gotXCPU == 1)
      printf("Received XCPU, exiting. But we'll be back.\n");
  }
}



void sage(char *fname)
{
  int filenr, tree, halonr;
 // struct sigaction current_XCPU;

  struct stat filestatus;
  FILE *fd;
  time_t start, current;
  FILE *treep;
  char filelen[1000];
  char bufdir[1000];

  sprintf(bufdir, "%s/sage_out", OutputDir_s);
  mkdir(bufdir, 02755);

  MAXSNAPS = fileset + 1;
  //MPI_Init(&argc, &argv);
  //MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  //MPI_Comm_size(MPI_COMM_WORLD, &NTask);

/*  ThisNode = malloc(MPI_MAX_PROCESSOR_NAME * sizeof(char));

  MPI_Get_processor_name(ThisNode, &nodeNameLen);
  if (nodeNameLen >= MPI_MAX_PROCESSOR_NAME) 
  {
    printf("Node name string not long enough!...\n");
    ABORT(701);
  }*/

  /*if(argc != 2)
  {
    printf("\n  usage: L-Galaxies <parameterfile>\n\n");
    ABORT(1);
  }*/

  //atexit(bye);

 // sigaction(SIGXCPU, NULL, &saveaction_XCPU);
 // current_XCPU = saveaction_XCPU;
 // current_XCPU.sa_handler = termination_handler;
 // sigaction(SIGXCPU, &current_XCPU, NULL);

//  printf("MAXSNAPS = %d\n", MAXSNAPS);
 
//  memset(metallicities_sage, 0.0, 8*sizeof(double));
  ZZ = malloc(MAXSNAPS * sizeof(double));//calloc(MAXSNAPS, sizeof(double));  //malloc(sizeof(double) * MAXSNAPS);
  AA = malloc(MAXSNAPS * sizeof(double));//calloc(MAXSNAPS, sizeof(double));  //malloc(sizeof(double) * MAXSNAPS);
  Age = malloc(MAXSNAPS * sizeof(double)); //calloc(MAXSNAPS, sizeof(double)); //malloc(sizeof(double) * MAXSNAPS);
 // printf("Malloc size = %d\n", MAXSNAPS * sizeof(double));
  //read_parameter_file_s(fname);
  init_s();
  printf("After init_s()\n");
  /* a small delay so that processors dont use the same file */
  time(&start);
  do
    time(&current);
  while(difftime(current, start) < 5.0 * ThisTask);
  printf("After times\n");
  for(filenr = FirstFile; filenr <= LastFile; filenr++)
  {
    // if(filenr = FirstFile)
    // {
         printf("Main: filenr = %d\n", filenr);
         printf("LastSnapShotNr = %d\n", LastSnapShotNr);
         sprintf(filelen, "%s/treedata/Treefile_%03d.txt", SimulationDir_s, LastSnapShotNr);
         treep = fopen(filelen, "w+");
         fprintf(treep, "%s/treedata/trees_%03d.%d\n", SimulationDir_s, LastSnapShotNr, filenr);
         fclose(treep);
    // }

       sprintf(bufz0, "%s/treedata/trees_%03d.%d", SimulationDir_s, LastSnapShotNr, filenr);
    if(!(fd = fopen(bufz0, "r")))
      continue;  // tree file does not exist, move along
    else
      fclose(fd);

    sprintf(bufz0, "%s/sage_out/%s_z%1.3f_%d", OutputDir_s, FileNameGalaxies, ZZ[ListOutputSnaps[0]], filenr);
    if(stat(bufz0, &filestatus) == 0)	 // seems to exist, move along
      continue;

    if((fd = fopen(bufz0, "w")))
      fclose(fd);

    load_tree_table(filenr);

    for(tree = 0; tree < Ntrees; tree++)
    {
      
    //  if(gotXCPU)
  //      ABORT(5);
     // printf("tree = %d\n", tree);
      if(tree % 10000 == 0)
      {
        // printf("\ttask: %d\tnode: %s\tfile: %i\ttree: %i of %i\n", ThisTask, ThisNode, filenr, tree, Ntrees);
       
         printf("\ttask: %d\tfile: %i\ttree: %i of %i\n", ThisTask, filenr, tree, Ntrees);
         fflush(stdout);
      }
   //   printf("load tree\n");
      load_tree(filenr, tree);
     // printf("gsl_rng_set\n");
      gsl_rng_set(random_generator, filenr * 100000 + tree);
      NumGals = 0;
      GalaxyCounter = 0;
      //printf("construct galaxsies\n");
      for(halonr = 0; halonr < TreeNHalos[tree]; halonr++)
        if(HaloAux_s[halonr].DoneFlag == 0)
        construct_galaxies(halonr);

      //printf("save galaxies\n");
      save_galaxies(filenr, tree);
     // printf("free_galaxies\n");
      free_galaxies_and_tree();
    }
   // printf("finalize galaxies\n");
    finalize_galaxy_file(filenr);
   // printf("free table\n");
    free_tree_table();

    printf("\ndone file %d\n\n", filenr);
  }

  reset_vars();
  memset(bufz0, 0, sizeof(bufz0));
  //memset(metallicities_sage, 0.0, 8 * sizeof(double));
//  free(ZZ);
 // ZZ = NULL;
 // free(AA);
 // AA=NULL;
 // free(Age);
 // Age = NULL;
  exitfail = 0;

  printf("sage has finished\n");
  return;
}

