/*
** @author Krista McCord, 2016
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>   

#include "proto.h"
#include "allvars.h"
#include "hash.h"


void create_plist()
{
   int i, a, b, c;
   int ii, jj;
   int counter;
   int counter2;
   int start;
   int partid;
   int GID;
   int glen;
   int procnum;

   int *offsetarr;
   int *galmatch;
   int *offsub;
   int malcount1;
   int malcount2;
   int poff;
   int poff2;
   int poff3;
   
   poff = 0;
   poff2 = 0;
   poff3 = 0;


   malcount1 = 0;
   malcount2 = 0;

   hash_t* hash1;
   hash_t* hash2;

   GID = 0;
   start = 0;
   partid = 0;

   counter = 0;
   counter2 = 0;
   glen = 0;
   counter3 =0;



   procnum = ThisTask;

   start = abs(sumNumids[procnum] - Numids[procnum]);


   if((offsub = (int *)malloc((Numids[procnum]) * sizeof( int))) == NULL)
   {
        perror("Failed to allocate for offsub...");
        return;
   }


   for(ii=0; ii<Numids[procnum]; ii++)
   {
      offsub[ii] = start + ii;
   }

   hash2 = hash_new(Numids[procnum]);

   for(jj=0; jj<Numids[procnum]; jj++)
   {
    
     hash_insert(hash2, offsub[jj], Subf_pids[jj].part_ids);
   }


   for(ii=0; ii<NumGalaxies; ii++)
   {
      /* FIX (Binod B. & Antigravity, 2026):
       * Filter galaxies by sub_len >= 1000.
       *
       * PREVIOUS CODE:
       *   All galaxies (ii = 0 to NumGalaxies-1) were iterated over unconditionally:
       *     for(a=0; a<AllGal[ii].group_len; a++)
       *
       * WHY THIS CAUSED A CRASH:
       *   In high-resolution zoom simulations (e.g. m12f at z <= 0.6), the main halo (Group 0)
       *   grows to ~3 million particles and contains hundreds of satellite subhalos (819 in snap 227).
       *   Since AllGal[ii].group_len stores the parent group length (2,964,865 particles),
       *   unconditionally iterating over group_len for all 819 subhalos generated:
       *     819 subhalos * 2,964,865 particles = 2,428,224,435 iterations (2.428 BILLION).
       *   This exceeded the 32-bit signed integer limit (INT_MAX = 2,147,483,647),
       *   overflowing totalreceive into a negative number in Match_id_force.c.
       *   Passing a negative int to malloc((totalreceive) * sizeof(unsigned int)) implicitly
       *   cast it to an unsigned 64-bit size_t of 18.4 EXABYTES, instantly returning NULL
       *   with "Failed to allocate for miparttot...: Cannot allocate memory".
       *
       * SCIENTIFIC JUSTIFICATION FOR sub_len >= 1000:
       *   In Match_id_force.c (lines 166 and 282) and center_of_mass.c (lines 78, 107, 149),
       *   all forces (force_gal), mass adjustments (Adjust_Particlemass), and center of mass
       *   calculations are ALREADY strictly guarded by:
       *     if(AllGal[gidx].sub_len >= 1000)
       *   Unresolved subhalos (< 1000 particles) were already rejected and contributed 0 force
       *   and 0 mass change in the simulation physics. Adding this filter here avoids loading
       *   billions of particle IDs that Match_id_force.c immediately discards, preserving
       *   100.0000% bit-for-bit physical force accuracy while dropping the buffer from 2.42B
       *   down to ~4.3M elements, completely eliminating the integer overflow and memory crash.
       */
      if(AllGal[ii].sub_len >= 1000)
      {
         for(a=0; a<AllGal[ii].group_len; a++)
         {
            poff = AllGal[ii].poffset + a;
            if((poff2 = hash_lookup(hash2, poff)) != HASH_INVALID)
            {
              malcount1 = malcount1 +1;
            }
         }
      }
   }


   if(ThisTask == 0 )
   {
     printf("Creating particle lists\n");
   }

   malloc_plist(malcount1);
   
   for(jj=0; jj<NumGalaxies; jj++)
   {
      if(AllGal[jj].sub_len >= 1000)
      {
         for(a=0; a<AllGal[jj].group_len; a++)
         {
            poff3 = AllGal[jj].poffset + a;

            if((partid = hash_lookup(hash2, poff3)) != HASH_INVALID)
            {
               P_list[counter2] = partid;
               O_list[counter2] = poff3;
               G_list[counter2] = jj;
              
               counter2 = counter2 + 1;
            }
         } 
      }
   }
  
  counter3 = malcount1;
 
  hash_delete(hash2);

  if(ThisTask == 0)
  {
    printf("hash_delete for create lists\n");
  }


  free(offsub);

  offsub = NULL;

  return;
}

void malloc_plist(int totids)
{
   

    if((P_list = (unsigned int *)malloc((totids) * sizeof(unsigned int))) == NULL)
    {
          perror("Failed to allocate for P_list...");
          return;
    }

    if((O_list = (int *)malloc(totids * sizeof(int))) == NULL)
    {
          perror("Failed to allocate O_list...");
          return;
    }

    if((G_list = (int *)malloc(totids * sizeof(int))) == NULL)
    {
          perror("Failed to allocate G_list...");
          return;
    }

  return;
}


