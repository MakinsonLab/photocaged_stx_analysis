This repository contains  analysis code for _Behavioral control through the direct, focal silencing of neuronal activity_ by Elleman _et. al_, 2024. Each folder contains code used for data analysis associated with the corresponding figure. 

****

**Figure 4: Activation of STX-pbc 13 abolishes action potentials in layer 4 cortical neurons** 

Igor Pro scripts were used for AP detection and analysis of AP waveform . R scripts were use for statistical analysis and creation of plots. 

The Igor Pro scripts contain a functions to import .abf files written by Rothman JS and Silver RA. NeuroMatic: An Integrated Open-Source Software Toolkit for Acquisition, Analysis and Simulation of Electrophysiological Data. _Front Neuroinform._ 2018 Apr 4;12:14. [10.3389/fninf.2018.00014](https://www.frontiersin.org/articles/10.3389/fninf.2018.00014/full)


`Figure_4_AP_characteristics.ipf` contains functions used to analyze the AP waveform, and export example traces. **,A - B, D - J, M**

`Figure_4_AP_phase_plane.ipf` was used to generate the phase plane plot.**C**

`Figure_4_AP_raster.ipf` was used to create an action potential raster plot.**K - M**

****
**Figure S12 Silencing of cortical network activity by STX-bpc 13**

Python scripts were used to calculate current source densities from LFP recordings.

`Figure_S12_CSD_quantification.ipynb` is a Jupyter Notebook which contains code used to calculate current sources and sinks and generate the figure.











