/*
** @author Krista McCord, 2016
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>

#include "allvars.h"
#include "proto.h"


int create_lock_file()
{
  //struct flock lock, savelock;
  //int fd;
  FILE *fd;
  fd = fopen("/lustre/projects/p004_swin/kmccord/cosang_run1/lock_file.txt", "w"); 
  /*fd = open("lock_file.txt", O_RDONLY | O_CREAT);
  lock.l_type = F_RDLCK;
  lock.l_start = 0;
  lock.l_whence = SEEK_SET;
  lock.l_len = 0;
  savelock = lock;
  fcntl(fd, F_GETLK, &lock);
  if (lock.l_type == F_WRLCK)
  {
      printf("File is write-locked by process %ld.\n", lock.l_pid);
      exit(1);
  }
  fcntl(fd, F_SETLK, &savelock);*/
  fclose(fd);
}

int lock_file_exists()
{
   if( (access( "/lustre/projects/p004_swin/kmccord/cosang_run1/lock_file.txt", F_OK )) != -1 )
   {
      printf( "File Lock file exists\n" );
      return 1;
   }

return 0;
}


