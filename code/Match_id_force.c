/*
** @author Krista McCord, 2016
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>  
#include <mpi.h> 

#include "proto.h"
#include "allvars.h"
#include "hash.h"

void Match_ids()
{ 

#ifdef profiler
    struct ProfileInstance profile_instance;
    struct ProfileInstance hash_instance;
    struct ProfileInstance ppert_instance;
    struct ProfileInstance missing_instance;
    struct ProfileInstance mipart1_instance;
    struct ProfileInstance allgath_instance;
    struct ProfileInstance forceadd_instance;

    enum Bucket BucketNames = Matchids;
    enum Bucket Buckethash = createhash;
    enum Bucket Bucketperturb = perturb;
    enum Bucket Bucketmissing = missingpart;
    enum Bucket Bucketmiss1 = misspart1;
    enum Bucket Bucketallgath = allgath;
    enum Bucket Bucketforce = forceadd;

    profile_instance_start(&profile_instance, BucketNames, __FUNCTION__, __FILE__, __LINE__);
    
#endif

    int i, b, ii, j, kk;
    int partidx;
    int gidx;
    int pcounter;
    int misscount;
    int countm;
    int totmpart[NTask];
    int sumidx[NTask];  
    unsigned int *mipart;
    int *migidx;
    int totalreceive, totgalrecv;
    unsigned int *temppart;
    int *tempgidx;
    int tag1, tag2, tag3, tag4, tag5, tag6;
    
    int *misgidtot;
    unsigned int *mispidtot;
    tag1 = 1; 
    tag2 = 2;
    tag3 = 3;
    tag4 = 4;
    tag5 = 5;
    tag6 = 6;
    totalreceive = 0;
    totgalrecv = 0;
    countm = 0;
    misscount = 0;
    partidx = 0;
    gidx = 0;
    pcounter = 0;
    kk = 0;
    hash_t* hashtable;
    hash_t* hash3;
    

#ifdef profiler
    profile_instance_start(&hash_instance, Buckethash, __FUNCTION__, __FILE__, __LINE__);
#endif

    hashtable = hash_new(NumPart);
    //printf("ThisTask = %d, NumPart = %d\n", ThisTask, NumPart);
    for(b=0; b<NumPart; b++)
    {
       hash_insert(hashtable, P[b].ID, b);
    }


   hash3 = hash_new(counter3);

   for(ii=0; ii<counter3; ii++)
   {
     hash_insert(hash3, P_list[ii], G_list[ii]);
   }


#ifdef profiler
   profile_instance_end(&hash_instance, Buckethash);
#endif

   for(i=0;i<NumGalaxies; i++)
   {
      AllGal[i].CM_Pos[0] = AllGal[i].Pos[0];
      AllGal[i].CM_Pos[1] = AllGal[i].Pos[1];
      AllGal[i].CM_Pos[2] = AllGal[i].Pos[2];
   }
/////////////////////////////////////////////////////
   center_of_mass(hash3, hashtable);
   MPI_Barrier(MPI_COMM_WORLD);


/******************************************************************************************************/
#ifdef Extrap
   if(snap_oldgal_count > 1 && galsearchcount == 0)
   {
      
      Gal_search();

      galsearchcount = 1;
 
     
   }
#endif
  
   
#ifdef profiler
   profile_instance_start(&missing_instance, Bucketmissing, __FUNCTION__, __FILE__, __LINE__); 
#endif


//counter3 is the size of P_list, G_list and O_list

   for(i=0; i<counter3; i++)
   {
      if((partidx = hash_lookup(hashtable, P_list[i])) != HASH_INVALID)
        {
           //do nothing
        }
       else
       {

          misscount = misscount + 1;
       }
   }
 
    if((mipart = (unsigned int *)malloc((misscount) * sizeof(unsigned int))) == NULL)
    {
          perror("Failed to allocate for mipart...");
          return;
    }

    if((migidx = (int *)malloc((misscount) * sizeof(int))) == NULL)
    {
          perror("Failed to allocate for migidx...");
          return;
    }
	 //test printf
    //printf("MID1\n");
   for(i=0; i<counter3; i++)
   {
      if((partidx = hash_lookup(hashtable, P_list[i])) != HASH_INVALID)
        {
           if((gidx = hash_lookup(hash3, P[partidx].ID)) != HASH_INVALID)
           {
             if(AllGal[gidx].sub_len >= 1000)
             {
#ifdef Mass_Adjust
                Adjust_Particlemass(gidx, partidx);
#endif
            
		//test printf
		//printf("MID1\n");
                force_gal(gidx, partidx);
             }      
                 
           }
        }
        else
        {
          mipart[countm] = P_list[i];
          migidx[countm] = G_list[i];
          countm = countm + 1;
        }
   }

   hash_delete(hash3);


#ifdef profiler
   profile_instance_end(&missing_instance, Bucketmissing);
   profile_instance_start(&ppert_instance, Bucketperturb, __FUNCTION__, __FILE__, __LINE__);
#endif


   MPI_Allgather(&countm, 1, MPI_INT, totmpart, 1, MPI_INT, MPI_COMM_WORLD);


   for(i=0; i<NTask; i++)
   {
     
      totalreceive = totalreceive + totmpart[i];
   } 


   
#ifdef profiler
   profile_instance_start(&mipart1_instance, Bucketmiss1, __FUNCTION__, __FILE__, __LINE__);

#endif
    if((temppart = (unsigned int *)malloc((countm) * sizeof(unsigned int))) == NULL)
    {
          perror("Failed to allocate for mipart...");
          return;
    }

    if((tempgidx = (int *)malloc((countm) * sizeof(int))) == NULL)
    {
          perror("Failed to allocate for migidx...");
          return;
    }
    
    if((mispidtot = (unsigned int *)malloc((totalreceive) * sizeof(unsigned int))) == NULL)
    {
          perror("Failed to allocate for miparttot...");
          return;
    }

    if((misgidtot = (int *)malloc((totalreceive) * sizeof(int))) == NULL)
    {
          perror("Failed to allocate for migidxtot...");
          return;
    }

        for(j=0; j<countm; j++)
        {
           temppart[j] = mipart[j];
           tempgidx[j] = migidx[j];
        }

  

     for(j=0; j<NTask; j++)
     {
        if(j == 0)
        {
           sumidx[j] = 0;         
        }
        else
        {
           sumidx[j] = sumidx[j-1] + totmpart[j-1];         
        }
     }
#ifdef profiler
    profile_instance_end(&mipart1_instance, Bucketmiss1);
    profile_instance_start(&allgath_instance, Bucketallgath, __FUNCTION__, __FILE__, __LINE__);
    if(ThisTask <= 63)
    {
       printf(" counts sent = %d, totalreceive = %d, counter3 = %d\n", countm, totalreceive, counter3);
    }
#endif    
     MPI_Allgatherv(temppart, countm, MPI_UNSIGNED, mispidtot, totmpart, sumidx , MPI_UNSIGNED, MPI_COMM_WORLD);
     MPI_Allgatherv(tempgidx, countm, MPI_INT, misgidtot, totmpart, sumidx , MPI_INT, MPI_COMM_WORLD);

#ifdef profiler
     profile_instance_end(&allgath_instance, Bucketallgath);
     profile_instance_start(&forceadd_instance, Bucketforce, __FUNCTION__, __FILE__, __LINE__);
#endif

	 //test printf
	//printf("MID2\n");
     for(i=0; i<totalreceive; i++)
     {
      if((partidx = hash_lookup(hashtable, mispidtot[i])) != HASH_INVALID)
        {
           gidx = misgidtot[i];
           if(AllGal[gidx].sub_len >= 1000)
           {
#ifdef Mass_Adjust
             Adjust_Particlemass(gidx, partidx);
#endif
		//test printf
		//printf("MID2\n");           
             force_gal(gidx, partidx); 
           }
        }
     }
#ifdef profiler
     profile_instance_end(&forceadd_instance, Bucketforce);
#endif
     

     MPI_Barrier(MPI_COMM_WORLD);

     free(temppart);
     free(tempgidx);
     free(mispidtot);
     free(misgidtot);
     temppart = NULL;
     tempgidx = NULL;
     mispidtot = NULL;
     misgidtot = NULL;


#ifdef profiler
     profile_instance_end(&ppert_instance, Bucketperturb);
#endif
     if(ThisTask == 0)
     {
        printf("DelHashTable\n");
     }
     hash_delete(hashtable);
    // printf("After deleting hash table\n");
     free(mipart);
     free(migidx);
     mipart = NULL;
     migidx = NULL;     


/*#ifdef zoomsim
    for(i = 0; i<NumGalaxies; i++)
    {
       for(b=0; b<NumPart; b++)
       {
#ifdef Mass_Adjusti
      
           Adjust_Particlemass(i, b);
#endif
           force_gal(i, b);
       }

    }


#endif*/


#ifdef profiler
      profile_instance_end(&profile_instance, BucketNames);
#endif
     return;
}


void Adjust_Particlemass(int gal, int target)
{
   double TotGalMass, TotMassHalo, Tpartmass;
   double px, py, pz, r, dm, fracmass;
   int Parthalo;
   double partmass, masspart2;

   px = P[target].Pos[0] - AllGal[gal].CM_Pos[0];
   py = P[target].Pos[1] - AllGal[gal].CM_Pos[1];
   pz = P[target].Pos[2] - AllGal[gal].CM_Pos[2];

   Parthalo = AllGal[gal].group_len;
   partmass = P[target].InitMass;
 
 

   r = px*px + py*py + pz*pz;

#ifdef Extrap
   if(snap_oldgal_count > 1)
   {
      TotGalMass = AllGal[gal].NewBlackHoleMass + AllGal[gal].NewBulgeMass + (AllGal[gal].NewStellarMass - AllGal[gal].NewBulgeMass) + AllGal[gal].NewColdGas; 
   }
   else
   {
      TotGalMass = AllGal[gal].BlackHoleMass + AllGal[gal].BulgeMass + (AllGal[gal].StellarMass - AllGal[gal].BulgeMass) + AllGal[gal].ColdGas; 
   }
#else

     TotGalMass = AllGal[gal].BlackHoleMass + AllGal[gal].BulgeMass + (AllGal[gal].StellarMass - AllGal[gal].BulgeMass) + AllGal[gal].ColdGas;
  
#endif

   fracmass = partmass / AllGal[gal].Mvir;

   dm = fracmass * TotGalMass;


   if(r <= ((AllGal[gal].Rvir / All.Time) * ( AllGal[gal].Rvir / All.Time)))
   {

      P[target].Mass = P[target].Mass - dm;
   }
 
return;

}  
