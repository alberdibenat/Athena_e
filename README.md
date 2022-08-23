# Athena_e
Scripts to run bandwith acceptance simulations and analysis for the interferometric beam-size monitor.

First, we calculate the interference pattern for a given V and a different BW in the visibility calculator. The different visibilities are achieved by changing the geometric parameter L:
L = 0.04m -> V = 0.1008
L = 0.05m -> V = 0.2303
L = 0.07m -> V = 0.4728
L = 0.1m -> V = 0.6928
L = 0.3m -> V = 0.9600

Then one can calculate the retrieved visibilities by fitting the visibility equation to the interference patterns using the script "retrieve visibility". The fitted equation, however, does not include the BW contribution, so we have to see up to which BW the retrieved visibility corresponds to the real visibility. This is seen in the plot.

The plot_spectrum script plots the spectrum extracted from spectra.
