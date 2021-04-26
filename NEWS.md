# bayesmove 0.2.0 (2021-04-26)

* Add cluster_obs() function to make observation-level inference of behavioral states.
* Add Shiny app as shiny_tracks() function to dynamically explore telemetry data within and among individuals.
* Add insert_NAs() function to insert NAs for creation of regular time series.
* Export get_MAP_internal() function for use within workflow for observation-level behavior estimation using cluster_obs().
* Update assign_tseg_internal() to account for non-consecutive `time1` variable.
* Update plot_breakpoints_behav() to account for non-consecutive `time` variable.
* Update some function arguments to be compatible with {furrr} package (>= 0.2.0).
* Update API for viewing progress bar in segment_behavior() via {progressr}.
* Update prep_data() to also calculate net-squared displacement
* Update method by which round_track_time() calculates new date and add argument for units


# bayesmove 0.1.0 (2020-09-21)

* First release of bayesmove package.


# bayesmove 0.0.0.9000 (2020-06-02)

* Bayesmove package in development.
