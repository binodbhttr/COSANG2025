#include<time.h>
#include<math.h>

#include "allvars.h"
#include "proto.h"
//#include "allvars_bt.h"
#include "proto_bt.h"

/* returns the number of cpu-ticks in seconds that
 * have elapsed. (or the wall-clock time)
 */

double second_bt(void)
{

  return ((double) clock()) / CLOCKS_PER_SEC;


  /* note: on AIX and presumably many other 32bit systems,
   * clock() has only a resolution of 10ms=0.01sec
   */
}

double measure_time_bt(void)	/* strategy: call this at end of functions to account for time in this function, and before another (nontrivial) function is called */
{
  double t, dt;

  t = second_bt();
  dt = t - WallTime;
  WallTime = t;

  return dt;
}

/* returns the time difference between two measurements
 * obtained with second(). The routine takes care of the
 * possible overflow of the tick counter on 32bit systems.
 */
double timediff_bt(double t0, double t1)
{
  double dt;

  dt = t1 - t0;

  if(dt < 0)			/* overflow has occured (for systems with 32bit tick counter) */
    {

      dt = t1 + pow(2, 32) / CLOCKS_PER_SEC - t0;

    }

  return dt;
}

