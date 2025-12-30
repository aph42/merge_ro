# Merge RO

A set of tools to process and merge radio occultation profiles with other datasets
to allow for radiative calculations and to carry out various analyses.

## Existing Workflow

 0. proc_mission.sh: Batch downloading and processing of GPS output from CDAAC
    - Downloads tarballs by mission, year, and day, extracts them.
      This calls process_mission.py for each day, which in turn calls process.py process_date
      * process.py: Processes netcdf files to allow for merged input. See open_profile() for
        full details of processing, but part of this involves ensuring profiles are all associated
        with unique timestamps. Outputs to raw/{mission}/

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
      *

 5. do_kw.sh: Starts processes running kw_filter_all.sh
    - kw_filter_all.sh: Runs filter.py with different actions
      *e.g. 

## Other key files

read.py: utilities for opening pygeode datasets of the various stages of the workflow

## TODO list
 - [ ] Run merging of a single month/profile to ensure code still works. Choose example
       profile from figure.
 - [ ] Adapt to merge with ERA 5
