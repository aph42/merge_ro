# Merge RO

A set of tools to process and merge radio occultation profiles with other datasets
to allow for radiative calculations and to carry out various analyses.

## Existing Workflow

 0. do_proc.sh: Starts processes for dl_mission.sh and proc_mission.sh
    - dl_mission.sh: Batch downloading from CDAAC
      * Downloads tarballs by mission, year, and day
    - proc_mission.sh: Batch processing of GPS output from CDAAC
      * Expands tarballs day-by-day and calls process_mission.py for each day, 
        which in turn calls process.py process_date
      * process.py: Processes netcdf files to allow for merged input. See
        open_profile() for full details of processing, but part of this
        involves ensuring profiles are all associated with unique timestamps.
        Outputs to raw/{mission}/

 1. do_merge.sh: Starts processes running merge_all.sh for different years.
    - merge_all.sh: Loops through days of the year calling merge_rad.py merge
      * merge_rad.py merge: calls merge()
        Blends RO temperature profiles with ERA Interim, MLS ozone and water
        vapour profiles with a WACCM climatology, and includes radiative heating rates,
        albedo and skin temperature from ERA Interim to enable radiative calculations.
        Outputs to merged/mls/{year}-{month}
      * merge_rad.py merge: calls merge_EI_comp()
        Outputs ERA Interim ozone and water vapour at profile locations for comparison
        Outputs to merged/ei/{year}-{month}

 2. do_rad.sh: Starts processes running rrtm_all.sh for different years
    - rrtm_all.sh: Loops through days of the year, calling merge_rad.py rrtm
      * merge_rad.py rrtm: calls either run_rrtm() or run_pyracc() to compute profile-by-profile radiative heating rates
        Outputs to rad/mls/
 
 3. do_stats.sh: Starts processes running filter_all.sh for different years
    - filter_all.sh: Runs filter.py for each day of a given year
      * filter.py: profile-by-profile based cacluation of a variety of quantities
                   general statistics are in rad/all/year
 
 4. do_grid.sh: Starts processes running grid_all.sh for different years
    - grid_all.sh: Runs grid.py for each month
      * grid.py: outputs gridded quantities

 5. do_kw.sh: Starts processes running kw_filter_all.sh
    - kw_filter_all.sh: Runs filter.py with different actions
      * e.g. calculating wave-number frequency decompositions etc. 

## Other key files

read.py: utilities for opening pygeode datasets of the various stages of the workflow

## TODO list
 - [ ] Decide on list of missions to use
 - [ ] Check on date_hash implementation to omit existing cosmic1 data
 - [ ] Process cosmic2021 re-processing
 - [ ] Write summary function that looks for raw nc files, quotes first and last days and number of missing dates within that range
 - [ ] Test merging and radiative calculation for Jan 2010 using ERA5 humidity and surface quantities
 - [ ] Discuss with Aaron if he can generate a global gridded MLS dataset that spans 2004 through 2025
