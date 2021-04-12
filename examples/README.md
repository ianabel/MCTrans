
# Examples 

This directory contains some sample inputs for MCTrans++ and reference outputs.
To compare reference outputs and the output of your current copy of MCTrans++ the expected behaviour is:
`
MCTrans/examples$ ../MCTrans++ Mirror_1keV.conf > Mirror_1kev.out
MCTrans/examples$ diff -s Mirror_1keV.report Mirror_1keV.out
Files Mirror_1keV.report and Mirror_1keV.out are identical
`

## Warm Plasma Experiment -- `Mirror_1keV.conf`

This configuration file describes a simple mirror experiment designed to achieve 1 keV ion and electron temperatures. The reference output is in `Mirror_1keV.report`. 
This configuration does not include alpha particles, nuclear physics, or neutral physics. The main variant of this configuration targets the temperature and reports the 
required voltage and Mach number. This is the simplest mode of operation of MCTrans and can be used as a template from which to make your own configuration files.

