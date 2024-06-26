# polycrystalline-PUF-project
Models: 
series_model.py (note: the LUTs don't calculate the min/max voltage so they musn't be exceeded eg. keep GB voltage below 0.9V for 1e16 to 2e16 GB)
array_model.py
(some of the paths may need altering as they used to be in separate directories)

Visualisation:
To generate colour plots in griddata and show where  (sorry, it's a bit convoluted):
Use colour_plot_divider.py to generate the voltage/current map, then use griddata_plot_not_adjusted_position.py with voltage divider data set to True. This saves all the points where it exceeds a defined limit (which can be obtained from the min/max voltage of a LUT).
Then use colour_plot.py to regenerate voltage/current plot. griddata_plot_not_adjusted_position.py can then be used to generate the plot showing where the voltage limits are exceeded. (If you skip the divider step, it will still work but won't highlight if it those points are exceeded).
NOTE: the position calculation does not return reliable results, it's mostly just a good way of visualising it in 2D.

energy_band_plot.py (plot for single boundary, and also used to calculate the voltage split data and write to csv)
energy_band_plot_mult_boundaries.py (uses .raw sim file, assumes no voltage dropped across grain which could be added)

simple_graphs.py used to plot the single GB IV plots
contact_plots.py used to plot the 2D array model results, there are different plots commented eg. different temperatures, and the quantisation calculations/plots are also at the bottom

effective_mobility.py lots of different versions commented out

voltage_split_plot.py uses csv files from energy_band_plot.py
