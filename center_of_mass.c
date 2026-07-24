/*
** @author Krista McCord (2016), Revision Shahram Talei (2017)
** Corrected & cleaned Binod Bhattarai (2025)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "proto.h"
#include "allvars.h"
#include "hash.h"

/* Assuming the following globals exist:
   - int ThisTask, NumGalaxies, NumPart, counter3;
   - float CM_Px, CM_Py, CM_Pz, CM_M;  // iteration accumulators
   - float *CM_Pxlist, *CM_Pylist, *CM_Pzlist, *CM_Mlist; // per-galaxy accumulators
   - struct { float Pos[3], CM_Pos[3], Rvir; int sub_len; } *AllGal;
   - struct { long long ID; float Pos[3]; float Mass; } *P;
   - long long *P_list;
*/

void center_of_mass(const hash_t* Ghash, const hash_t* Phash)
{
    /*
	Function to compute the center of mass (CoM) positions for all galaxies
	using mass-weighted averages of particle positions within the virial radius (Rvir).
	The calculation is performed in parallel across MPI tasks, with an initial
	accumulation followed by a shrinking-sphere refinement for better accuracy.

	Arguments:
	- const hash_t* Ghash: Hash table mapping galaxy IDs to their indices in AllGal.
	- const hash_t* Phash: Hash table mapping particle IDs to their indices in P.

	Returns:
	- void: The function updates the CM_Pos field of each galaxy in the AllGal array.
	
	*/

	int i, j, k;
    float rcheck;

    malloc_CM(); /* Allocate per-galaxy accumulators */

    /* Accumulate local (per-task) mass-weighted sums inside Rvir for sufficiently
       resolved subhalos */
    for (i = 0; i < counter3; i++) {
        int pindex = hash_lookup(Phash, P_list[i]);
        if (pindex == HASH_INVALID) continue;

        int gidx = hash_lookup(Ghash, P[pindex].ID);
        if (gidx == HASH_INVALID) continue;
        if (gidx < 0 || gidx >= NumGalaxies) continue; /* safety check */ 
        /*  In the old code this was not checked, which could lead to out-of-bounds access.  
        gidx is whatever value we pulled out of Ghash for a particle’s ID, so it’s only meaningful if 
        it coincides with a valid slot in AllGal[]. In the normal run every call to hash_insert() for 
        Ghash uses the galaxy’s index 
        (0 … NumGalaxies-1) as the stored value, but nothing in the hash layer enforces that. 
        If the bookkeeping gets out of sync—e.g. we rebuild AllGal with a smaller NumGalaxies after filtering galaxies 
        but continue to reuse a hash table created before the trim (the code in cosang_data_setup.c 
        (lines 214-236) only reports/broadcasts the new count; existing hash entries still carry the old indices),
        Or if we restart from a snapshot where the galaxy list shrank yet accidentally reuse a persistent Ghash,
        or the hash gets corrupted because P_list/G_list were filled from bad Subfind offsets 
        (create_plists.c (lines 85-128)) and we inserted garbage indices,
        then hash_lookup() can hand back any unsigned int it happens to have stored—including 
        values >= NumGalaxies. The guard in center_of_mass.c line 42 bails out of that particle 
        when it detects such a stray index, preventing an out-of-bounds access to AllGal[gidx] or the CM 
        accumulator arrays.
        */
        float dx = P[pindex].Pos[0] - AllGal[gidx].Pos[0]; /* vectorized diff */
        float dy = P[pindex].Pos[1] - AllGal[gidx].Pos[1];
        float dz = P[pindex].Pos[2] - AllGal[gidx].Pos[2];
        rcheck   = dx*dx + dy*dy + dz*dz; /* Check if the particle is within the virial radius */

        if (rcheck <= AllGal[gidx].Rvir * AllGal[gidx].Rvir &&
            AllGal[gidx].sub_len >= 1000) /* Only consider well-resolved subhalos */
        {
            CM_Pxlist[gidx] += P[pindex].Mass * P[pindex].Pos[0];
            CM_Pylist[gidx] += P[pindex].Mass * P[pindex].Pos[1];
            CM_Pzlist[gidx] += P[pindex].Mass * P[pindex].Pos[2];
            CM_Mlist[gidx]  += P[pindex].Mass;
        }
        //printf("LISTS: rcheck = %f, Rvir = %f, CM_Pxlist= %f, CM_Mlist = %f\n", rcheck, AllGal[gidx].Rvir, CM_Pxlist[gidx], CM_Mlist[gidx]);

    }

    /* Vectorized reduction of per-galaxy accumulators across tasks */
    double *send_buf = (double *)malloc(NumGalaxies * 4 * sizeof(double));
    double *recv_buf = (double *)malloc(NumGalaxies * 4 * sizeof(double));
    if (!send_buf || !recv_buf) {
        printf("Task %d: Failed to allocate MPI reduction buffers in center_of_mass\n", ThisTask);
        ABORT(113);
    }

    for (j = 0; j < NumGalaxies; j++) {
        send_buf[4 * j + 0] = (double)CM_Pxlist[j];
        send_buf[4 * j + 1] = (double)CM_Pylist[j];
        send_buf[4 * j + 2] = (double)CM_Pzlist[j];
        send_buf[4 * j + 3] = (double)CM_Mlist[j];
    }

    MPI_Allreduce(send_buf, recv_buf, NumGalaxies * 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (j = 0; j < NumGalaxies; j++) {
        if (AllGal[j].sub_len >= 1000) {
            double cmass = recv_buf[4 * j + 3];
            if (cmass > 0.0) {
                double cmx = recv_buf[4 * j + 0] / cmass;
                double cmy = recv_buf[4 * j + 1] / cmass;
                double cmz = recv_buf[4 * j + 2] / cmass;

                AllGal[j].CM_Pos[0] = isnan(cmx) ? AllGal[j].Pos[0] : (float)cmx;
                AllGal[j].CM_Pos[1] = isnan(cmy) ? AllGal[j].Pos[1] : (float)cmy;
                AllGal[j].CM_Pos[2] = isnan(cmz) ? AllGal[j].Pos[2] : (float)cmz;
            } else {
                /* No mass contributed on any rank: fall back to halo position */
                AllGal[j].CM_Pos[0] = AllGal[j].Pos[0];
                AllGal[j].CM_Pos[1] = AllGal[j].Pos[1];
                AllGal[j].CM_Pos[2] = AllGal[j].Pos[2];
            }
        } else {
            /* For unresolved halos, fall back to default position */
            AllGal[j].CM_Pos[0] = AllGal[j].Pos[0];
            AllGal[j].CM_Pos[1] = AllGal[j].Pos[1];
            AllGal[j].CM_Pos[2] = AllGal[j].Pos[2];
        }
    }

    free(send_buf);
    free(recv_buf);

    /* Free the lists (not needed beyond initial CM) */
    free(CM_Pxlist); CM_Pxlist = NULL;
    free(CM_Pylist); CM_Pylist = NULL;
    free(CM_Pzlist); CM_Pzlist = NULL;
    free(CM_Mlist);  CM_Mlist  = NULL;

    /* Shrinking-sphere refinement (fully vectorized) */
    int *active_refine = (int *)malloc(NumGalaxies * sizeof(int));
    if (!active_refine) {
        printf("Task %d: Failed to allocate active_refine buffer\n", ThisTask);
        ABORT(114);
    }

    int any_active = 0;
    for (k = 0; k < NumGalaxies; k++) {
        if (AllGal[k].sub_len >= 1000) {
            active_refine[k] = 1;
            any_active = 1;
        } else {
            active_refine[k] = 0;
        }
    }

    if (any_active) {
        double *iter_px = (double *)malloc(NumGalaxies * sizeof(double));
        double *iter_py = (double *)malloc(NumGalaxies * sizeof(double));
        double *iter_pz = (double *)malloc(NumGalaxies * sizeof(double));
        double *iter_m  = (double *)malloc(NumGalaxies * sizeof(double));
        double *iter_send = (double *)malloc(NumGalaxies * 4 * sizeof(double));
        double *iter_recv = (double *)malloc(NumGalaxies * 4 * sizeof(double));
        int *local_diff = (int *)malloc(NumGalaxies * sizeof(int));
        int *global_diff = (int *)malloc(NumGalaxies * sizeof(int));

        if (!iter_px || !iter_py || !iter_pz || !iter_m || !iter_send || !iter_recv || !local_diff || !global_diff) {
            printf("Task %d: Failed to allocate active refinement iteration buffers\n", ThisTask);
            ABORT(115);
        }

        int itcount;
        for (itcount = 1; itcount <= 2; itcount++) {
            memset(iter_px, 0, NumGalaxies * sizeof(double));
            memset(iter_py, 0, NumGalaxies * sizeof(double));
            memset(iter_pz, 0, NumGalaxies * sizeof(double));
            memset(iter_m,  0, NumGalaxies * sizeof(double));

            for (j = 0; j < counter3; j++) {
                int partidx = hash_lookup(Phash, P_list[j]);
                if (partidx == HASH_INVALID) continue;

                int gidx = hash_lookup(Ghash, P[partidx].ID);
                if (gidx == HASH_INVALID || gidx < 0 || gidx >= NumGalaxies) continue;
                if (!active_refine[gidx]) continue;

                float rsmall = AllGal[gidx].Rvir;
                int s;
                for (s = 0; s < itcount; s++) rsmall *= 0.8f;

                float dx = P[partidx].Pos[0] - AllGal[gidx].CM_Pos[0];
                float dy = P[partidx].Pos[1] - AllGal[gidx].CM_Pos[1];
                float dz = P[partidx].Pos[2] - AllGal[gidx].CM_Pos[2];
                float r2 = dx*dx + dy*dy + dz*dz;

                if (r2 <= rsmall * rsmall) {
                    iter_px[gidx] += P[partidx].Mass * P[partidx].Pos[0];
                    iter_py[gidx] += P[partidx].Mass * P[partidx].Pos[1];
                    iter_pz[gidx] += P[partidx].Mass * P[partidx].Pos[2];
                    iter_m[gidx]  += P[partidx].Mass;
                }
            }

            for (k = 0; k < NumGalaxies; k++) {
                iter_send[4 * k + 0] = iter_px[k];
                iter_send[4 * k + 1] = iter_py[k];
                iter_send[4 * k + 2] = iter_pz[k];
                iter_send[4 * k + 3] = iter_m[k];
            }

            MPI_Allreduce(iter_send, iter_recv, NumGalaxies * 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            memset(local_diff, 0, NumGalaxies * sizeof(int));
            for (k = 0; k < NumGalaxies; k++) {
                if (active_refine[k]) {
                    double cmass = iter_recv[4 * k + 3];
                    if (cmass > 0.0) {
                        double cmx = iter_recv[4 * k + 0] / cmass;
                        double cmy = iter_recv[4 * k + 1] / cmass;
                        double cmz = iter_recv[4 * k + 2] / cmass;

                        double dx = fabs(cmx - AllGal[k].CM_Pos[0]);
                        double dy = fabs(cmy - AllGal[k].CM_Pos[1]);
                        double dz = fabs(cmz - AllGal[k].CM_Pos[2]);

                        AllGal[k].CM_Pos[0] = (float)cmx;
                        AllGal[k].CM_Pos[1] = (float)cmy;
                        AllGal[k].CM_Pos[2] = (float)cmz;

                        if (dx > 0.01 || dy > 0.01 || dz > 0.01) {
                            local_diff[k] = 1;
                        }
                    }
                }
            }

            MPI_Allreduce(local_diff, global_diff, NumGalaxies, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

            int any_still_active = 0;
            for (k = 0; k < NumGalaxies; k++) {
                if (active_refine[k]) {
                    if (global_diff[k] == 0) {
                        active_refine[k] = 0;
                    } else {
                        any_still_active = 1;
                    }
                }
            }

            if (!any_still_active) break;
        }

        free(iter_px); free(iter_py); free(iter_pz); free(iter_m);
        free(iter_send); free(iter_recv);
        free(local_diff); free(global_diff);
    }

    free(active_refine);

    if (ThisTask == 0) {
        printf("Finished CoM calculation\n");
        fflush(stdout);
    }
}

/* ========================= Helpers ========================= */

void malloc_CM(void)

{
    CM_Pxlist = (float *)calloc(NumGalaxies, sizeof(float)); /*Allocate memory for per-galaxy accumulators */
    CM_Pylist = (float *)calloc(NumGalaxies, sizeof(float)); 
    CM_Pzlist = (float *)calloc(NumGalaxies, sizeof(float));
    CM_Mlist  = (float *)calloc(NumGalaxies, sizeof(float));
}

void malloc_CM_iterate(void)
{
    /* Kept for compatibility with proto.h; not used in this file’s flow */
    CM_Pxlist = (float *)calloc(1, sizeof(float));
    CM_Pylist = (float *)calloc(1, sizeof(float));
    CM_Pzlist = (float *)calloc(1, sizeof(float));
    CM_Mlist  = (float *)calloc(1, sizeof(float));
}

/* One shrinking-sphere iteration for a single galaxy */
int cm_iterate(const hash_t* hashtable, const hash_t* Phash, int galid, int itcount)
{
    int j;
    float rsmall = AllGal[galid].Rvir;
    for (j = 0; j < itcount; j++) rsmall *= 0.8f;  /* shrink by 20% per iter */

    /* local accumulators */
    CM_Px = 0.0f; CM_Py = 0.0f; CM_Pz = 0.0f; CM_M = 0.0f;

    /* Loop only over local particle list via Phash -> partidx (no O(NumPart) scan) */
    for (j = 0; j < counter3; j++) {
        int partidx = hash_lookup(Phash, P_list[j]);
        if (partidx == HASH_INVALID) continue;

        /*
        /* Compared to center_of_mass_old.c, this iteration no longer scans all NumPart
        * entries looking for the same particle twice. We march directly over the
        * P_list[] subset, resolve each entry to its particle slot via Phash, and
        * accumulate only if it belongs to `galid` and lies inside the current
        * shrinking sphere (rsmall * rsmall). That keeps the mass-weighted sum identical
        * to the old code but removes an O(counter3 × NumPart) loop and prevents any
        * out-of-bounds writes. */

        int gidx = hash_lookup(hashtable, P[partidx].ID);
        if (gidx != galid) continue;

        float dx = P[partidx].Pos[0] - AllGal[galid].CM_Pos[0];
        float dy = P[partidx].Pos[1] - AllGal[galid].CM_Pos[1];
        float dz = P[partidx].Pos[2] - AllGal[galid].CM_Pos[2];
        float r2 = dx*dx + dy*dy + dz*dz;

        if (r2 <= rsmall * rsmall) {
            CM_Px += P[partidx].Mass * P[partidx].Pos[0];
            CM_Py += P[partidx].Mass * P[partidx].Pos[1];
            CM_Pz += P[partidx].Mass * P[partidx].Pos[2];
            CM_M  += P[partidx].Mass;
        }
    }

    /* Reduce as double for better numeric stability */
    double send_px = (double)CM_Px, send_py = (double)CM_Py, send_pz = (double)CM_Pz, send_m = (double)CM_M;
    double pcmx = 0.0, pcmy = 0.0, pcmz = 0.0, cmass = 0.0;
    MPI_Allreduce(&send_px, &pcmx, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&send_py, &pcmy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&send_pz, &pcmz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&send_m,  &cmass,1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    /* Previous (orig) CM with sensible fallbacks */
    double cmx_orig = isnan(AllGal[galid].CM_Pos[0]) ? (double)AllGal[galid].Pos[0] : (double)AllGal[galid].CM_Pos[0];
    double cmy_orig = isnan(AllGal[galid].CM_Pos[1]) ? (double)AllGal[galid].Pos[1] : (double)AllGal[galid].CM_Pos[1];
    double cmz_orig = isnan(AllGal[galid].CM_Pos[2]) ? (double)AllGal[galid].Pos[2] : (double)AllGal[galid].CM_Pos[2];

    int diff = 0;
    if (cmass > 0.0) {
        double cmx = pcmx / cmass;
        double cmy = pcmy / cmass;
        double cmz = pcmz / cmass;

        AllGal[galid].CM_Pos[0] = (float)cmx;
        AllGal[galid].CM_Pos[1] = (float)cmy;
        AllGal[galid].CM_Pos[2] = (float)cmz;

        double dx = fabs(cmx - cmx_orig);
        double dy = fabs(cmy - cmy_orig);
        double dz = fabs(cmz - cmz_orig);

        /* tolerance: 0.01 (adjust for your units) */
        diff = (dx <= 0.01 && dy <= 0.01 && dz <= 0.01) ? 0 : 1;
    } else {
        /* No mass in shrinking sphere: treat as converged (or choose to expand/abort) */
        diff = 0;
    }

    return diff;
}

