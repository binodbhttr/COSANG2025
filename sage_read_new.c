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

   printf("reading sage file paths\n");
   sagefiles();

   /* Basic sanity check: must have at least one SAGE file */
   if (NumSfiles <= 0)
   {
      printf("ERROR: No SAGE files found (NumSfiles = %d). Aborting Galaxies().\n", NumSfiles);
      return;
   }

   totmal = read_sageheader();

   if (totmal <= 0)
   {
      printf("ERROR: Failed to read SAGE header or no galaxies found (totmal = %d). Aborting Galaxies().\n", totmal);
      return;
   }

   NumGalaxies = totmal;

   mal_gal(totmal);
   if (SageOutput == NULL)
   {
      printf("ERROR: SageOutput allocation failed in mal_gal(). Aborting Galaxies().\n");
      return;
   }

   printf("Reading all Sage files\n");
   read_sage();

   malallgal(totmal);
   if (AllGal == NULL)
   {
      printf("ERROR: AllGal allocation failed in malallgal(). Aborting Galaxies().\n");
      return;
   }

   read_trees(totmal);
   fill_gal(totmal);

   /* Only free treepath if it is non-NULL, and then NULL it out */
   if (treepath != NULL)
   {
      free(treepath);
      treepath = NULL;
   }
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

    /* Guard against invalid NumSfiles */
    if (NumSfiles <= 0)
    {
        printf("ERROR: read_sageheader called with NumSfiles = %d.\n", NumSfiles);
        return -1;
    }

    for(i = 0; i < NumSfiles; i++)
    {
        sprintf(file1, "%s", sagepath[i].paths);

        /* Open Sage file */
        fd = fopen(file1, "rb");
        if (fd == NULL)
        {
            printf("Cannot open sage file: %s\n", file1);
            return -1;
        }

        /* Read in the header */
        if (fread(&Ntrees, sizeof(int), 1, fd) != 1)
        {
            printf("ERROR: Failed to read Ntrees from %s\n", file1);
            fclose(fd);
            return -1;
        }

        if (fread(&NtotGal, sizeof(int), 1, fd) != 1)
        {
            printf("ERROR: Failed to read NtotGal from %s\n", file1);
            fclose(fd);
            return -1;
        }

        totmal += NtotGal;

        fclose(fd);
    }

    return totmal;
}


void mal_gal(int totmal)
{
    /* Designate size of file for the array of structs */
    if (totmal <= 0)
    {
        printf("ERROR: mal_gal called with non-positive totmal = %d.\n", totmal);
        SageOutput = NULL;
        return;
    }

    SageOutput = (struct sage_galaxies *) malloc(totmal * sizeof(struct sage_galaxies));
    if (SageOutput == NULL)
    {
        printf("Failed to allocate for SageOutput (sage galaxies)...\n");
        return;
    }

    return;
}


void read_sage()
{
   int Ntrees;
   int NtotGal;
   int *galpertree = NULL;
   char file1[1000];
   int offset;
   int i;
   FILE *fd;

   offset = 0;

   if (NumSfiles <= 0)
   {
      printf("ERROR: read_sage called with NumSfiles = %d\n", NumSfiles);
      return;
   }

   if (SageOutput == NULL)
   {
      printf("ERROR: read_sage called but SageOutput is NULL.\n");
      return;
   }

   for(i = 0; i < NumSfiles; i++)
   {
      sprintf(file1, "%s", sagepath[i].paths);
      printf("file path in read_sage = %s\n", file1);

      /* Open Sage file */
      fd = fopen(file1, "rb");
      if (fd == NULL)
      {
         printf("Cannot open sage for reading: %s\n", file1);
         /* Clean up any allocated buffer before returning */
         if (galpertree != NULL)
         {
            free(galpertree);
            galpertree = NULL;
         }
         return;
      }

      /* Read in the header */
      if (fread(&Ntrees, sizeof(int), 1, fd) != 1 ||
          fread(&NtotGal, sizeof(int), 1, fd) != 1)
      {
         printf("ERROR: Failed to read SAGE header from %s\n", file1);
         fclose(fd);
         if (galpertree != NULL)
         {
            free(galpertree);
            galpertree = NULL;
         }
         return;
      }

      if (Ntrees <= 0 || NtotGal < 0)
      {
         printf("ERROR: Invalid header values in %s (Ntrees=%d, NtotGal=%d)\n", file1, Ntrees, NtotGal);
         fclose(fd);
         if (galpertree != NULL)
         {
            free(galpertree);
            galpertree = NULL;
         }
         return;
      }

      /* Allocate the needed array size for galpertree */
      galpertree = (int *) malloc(Ntrees * sizeof(int));
      if (galpertree == NULL)
      {
         printf("ERROR: malloc failed for galpertree in read_sage\n");
         fclose(fd);
         return;
      }

      /* Read in the rest of the header into an array */
      if (fread(galpertree, sizeof(int), Ntrees, fd) != (size_t) Ntrees)
      {
         printf("ERROR: Failed to read galpertree from %s\n", file1);
         fclose(fd);
         free(galpertree);
         galpertree = NULL;
         return;
      }

      /* Read in info for each galaxy in file */
      if (fread(&SageOutput[offset], sizeof(struct sage_galaxies), NtotGal, fd) != (size_t) NtotGal)
      {
         printf("ERROR: Failed to read SageOutput galaxy data from %s\n", file1);
         fclose(fd);
         free(galpertree);
         galpertree = NULL;
         return;
      }

      fclose(fd);

      offset += NtotGal;

      /* We only need galpertree within this iteration, so free it here */
      free(galpertree);
      galpertree = NULL;
   }

   return;
}


void malallgal(int totmal)
{
    /* Designate size of file for the array of structs */
    if (totmal <= 0)
    {
        printf("ERROR: malallgal called with non-positive totmal = %d.\n", totmal);
        AllGal = NULL;
        return;
    }

    AllGal = (struct All_Gal_info *) malloc(totmal * sizeof(struct All_Gal_info));
    if (AllGal == NULL)
    {
        printf("Failed to allocate for AllGal...\n");
        return;
    }

    return;
}


void fill_gal(int totmal)
{
    int i;
    int a;

    if (AllGal == NULL || SageOutput == NULL)
    {
        printf("ERROR: fill_gal called with NULL AllGal or SageOutput.\n");
        return;
    }

    for(i = 0; i < totmal; i++)
    {
       AllGal[i].Mvir = SageOutput[i].Mvir;

       for(a = 0; a < 3; a++)
       {
            AllGal[i].Pos[a] = SageOutput[i].Pos[a];
            AllGal[i].Vel[a] = SageOutput[i].Vel[a] * All.Time;
       }

       AllGal[i].Rvir = SageOutput[i].Rvir;
    }

    /* Free SageOutput now that we've copied everything we need */
    free(SageOutput);
    SageOutput = NULL;

    return;
}
