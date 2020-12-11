# bayesmove 0.2.0 (2020-12-01)

* Add cluster_obs() function to make observation-level inference of behavioral states.
* Add Shiny app as shiny_tracks() function to dynamically explore telemetry data within and among individuals.
* Add arguments to get_behav_hist() to specify order of states and iteration of MAP estimate for mixture-model results (observation-level behavioral states).
* Add insert_NAs() function to insert NAs for creation of regular time series.
* Export get_MAP_internal() function for use within workflow for observation-level behavior estimation using cluster_obs().
* Update assign_tseg_internal() to account for non-consecutive `time1` variable.
* Update plot_breakpoints_behav() to account for non-consecutive `time` variable.
* Update some function arguments to be compatible with {furrr} package (>= 0.2.0).
* Update prep_data() to also calculate net-squared displacement


# bayesmove 0.1.0 (2020-09-21)

* First release of bayesmove package.


# bayesmove 0.0.0.9000 (2020-06-02)

* Bayesmove package in development.
