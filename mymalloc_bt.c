#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"
//#include "allvars_bt.h"
#include "proto_bt.h"



void *mymalloc_bt(size_t n)
{
  void *p;

  if(!(p = malloc(n)))
    {
      if(n)
	{
	  printf("Failed to allocate memory for %u bytes.\n", (int) n);
	  exit(2);
	}
    }

  return p;
}

