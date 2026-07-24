/*
** @author Krista McCord, 2016
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

void profile_instance_start(struct ProfileInstance *profile_instance, int bucketname, char *func_name, char *file_name, int line_number)
{
	profile_instance->time = clock();
	profile_instance->bucketname = bucketname;
	if(!gProfileData[bucketname].num_instances)
	{
		//if this is the first time running, we need to fill some standard stuff in.
		//fill in the func name, file, and line into the global array.
	        struct ProfileData* data = &gProfileData[bucketname];
                data->procRank = ThisTask;
                data->funcname = func_name;
                data->filename = file_name;
                data->line = line_number;
	}
}

void profile_instance_end(struct ProfileInstance *profile_instance, int bucketname)
{
	profile_instance->timediff = ((double) (clock() - profile_instance->time)) / CLOCKS_PER_SEC;
	struct ProfileData* data = &gProfileData[bucketname];
        if(!gProfileData[bucketname].num_instances)
	{
            data->min_time=profile_instance->timediff;
            data->max_time=profile_instance->timediff;
        }
	
        if(gProfileData[bucketname].min_time > profile_instance->timediff)
        {
            data->min_time = profile_instance->timediff;
        }

        if(gProfileData[bucketname].max_time < profile_instance->timediff)
        {
            data->max_time = profile_instance->timediff;
        }

	data->num_instances++;
	data->total_time += profile_instance->timediff;

        if(ThisTask == 0)
        {
          printf("total_time = %f, bucketname = %d, max_time = %f\n", gProfileData[bucketname].total_time, bucketname, gProfileData[bucketname].max_time);
        }
	
}

void write_profile_data()
{
        char timefile[1000];
        FILE *timef; 
        enum Bucket BucketNames = NUM_BUCKETS;
        int i;
        //printf("%d %f, %f, %f, %d, %s\n", gProfileData[0].procRank, gProfileData[0].min_time, gProfileData[0].max_time, gProfileData[0].total_time, gProfileData[0].num_instances, gProfileData[0].funcname);
 
        sprintf(timefile, "%s/profile_data/PROC_%d.txt", All.OutputDir, ThisTask);
        timef = fopen(timefile, "a");
        for(i=0; i<BucketNames; i++)
        {
            fprintf(timef, "%d %f %f %f %d %s %s %d \n",  gProfileData[i].procRank, gProfileData[i].min_time, gProfileData[i].max_time, gProfileData[i].total_time, gProfileData[i].num_instances, gProfileData[i].funcname, gProfileData[i].filename, gProfileData[i].line);
        } 
	//fwrite(&gProfileData[0], sizeof(gProfileData),  BucketNames, timef);
  
        fclose(timef);
}

