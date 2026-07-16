# Hotfix Notes: Simulation Crash and Performance Fixes
**Date:** 2026-07-15

This document details the critical fixes and optimizations applied to resolve the simulation crash at Snapshot 254 (and the subsequent snapshot 252 restart crash) and improve execution efficiency.

---

## 1. Uninitialized Variable Fix in `forcetree.c`
* **File:** [forcetree.c](file:///u/bbhattar/no_backup/SimulationOutputs/center_of_mass_c_new_test/new_code/code/forcetree.c)
* **Problem:** Local stack variables `diskscaler` and `zo` inside `force_gal()` were accessed before being computed (at line 1879 and line 1918). They contained dirty stack memory left by the deep recursive SAGE run, causing NaN errors or floating-point exceptions (SIGFPE) on certain ranks.
* **Fix:** Initialized both variables at the top of the function. Moved the calculation of `diskmass`, `diskscaler`, and `zo = diskscaler / 9.4;` to execute *before* they are accessed at line 1879.

---

## 2. Typo in Extrapolation Timestep Index in `extrap.c`
* **File:** [extrap.c](file:///u/bbhattar/no_backup/SimulationOutputs/center_of_mass_c_new_test/new_code/code/extrap.c)
* **Problem:** In `calc_Radius_slopes()`, the loop variable `i` (which equals `All.SnapshotFileCount` after loop termination) was used to index the active galaxy structures `AllGal[i].timestep` and `OldGal[i].timestep` instead of the active galaxy arguments `newgal` and `oldgal`.
* **Fix:** Replaced the index `i` with `newgal` and `oldgal` to correctly compute:
  `scalediff2 = AllGal[newgal].timestep - OldGal[oldgal].timestep;`

---

## 3. Nested Loop Performance Optimization in `extrap.c`
* **File:** [extrap.c](file:///u/bbhattar/no_backup/SimulationOutputs/center_of_mass_c_new_test/new_code/code/extrap.c)
* **Problem:** The function `Gal_search()` executed `rotation_angles(j)` inside the outer loop over `i` ($N^2$ complexity). Because `rotation_angles(j)` does not depend on `i`, this ran the calculation $243$ million times redundantly, locking all 512 ranks at 100% CPU at peak memory for 5 minutes, triggering OS Out-of-Memory (OOM) kills.
* **Fix:** Moved the `rotation_angles(j)` loop to run once at the top of `Gal_search()`. This reduces the complexity to $O(N)$ (15,590 runs), dropping execution time from 5 minutes to <0.01 seconds.

---

## 4. Uninitialized Variable `b` in `disk_orient.c`
* **File:** [disk_orient.c](file:///u/bbhattar/no_backup/SimulationOutputs/center_of_mass_c_new_test/new_code/code/disk_orient.c)
* **Problem:** In `Calculate_Disk_VerF()`, the potential scaling factor variable `b` was used uninitialized at line 416 (`x_MN = b/diskscaler;`). This copy-paste error caused vertical disk forces to be computed with arbitrary stack garbage.
* **Fix:** Assigned `b = fabs(zo);` right after `zo` is calculated, mirroring the correct radial force counterpart function `Calculate_Disk_RadF()`.

---

## 5. Zero-Diskscaler Guards for Zero-Disk Galaxies
* **Files:** [forcetree.c](file:///u/bbhattar/no_backup/SimulationOutputs/center_of_mass_c_new_test/new_code/code/forcetree.c) & [disk_orient.c](file:///u/bbhattar/no_backup/SimulationOutputs/center_of_mass_c_new_test/new_code/code/disk_orient.c)
* **Problem:** Some galaxies have no disk components (hence `DiskScaleRadius == 0.0`). When we initialized `diskscaler = 0.0;` for safety, a particle exactly at the center of such galaxies (`r == 0.0`) evaluated the check `r <= 200.0*diskscaler` as `0.0 <= 0.0` (TRUE). Inside the force block, division by `diskscaler == 0.0` caused a `NaN` timestep crash.
* **Fix:** 
  1. Added `diskscaler > 1e-10` to the check in `forcetree.c`:
     `if(diskscaler > 1e-10 && r <= 200.0*diskscaler)`
  2. Added early returns in `Calculate_Disk_RadF()` and `Calculate_Disk_VerF()` inside `disk_orient.c`:
     ```c
     if(diskscaler <= 1e-10)
     {
         return 0.0;
     }
     ```
     This bypasses disk force calculations entirely for galaxies with no disk.
