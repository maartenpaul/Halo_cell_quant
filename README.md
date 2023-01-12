# Halo_cell_quant

This repository contains scripts to quantify the concentration of fluorescently tagged proteins using 3D confocal microscopy.
The analysis requires several scripts in the following order:



1. stardist3D_batch.py - This script can be used to segment nuclei in 3D using Fiji (https://www.fiji.sc) with Trackmate (https://imagej.net/plugins/trackmate/) using Stardist for 3D nuclear segmentation (https://github.com/stardist/stardist) 
2. 3DsegmentationBRCA2v3.cppipe - This pipeline needs to be run in CellProfiler: https://cellprofiler.org/
3. Halo cell quant analysis script.R
