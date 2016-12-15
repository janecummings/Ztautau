# Ztautau

The measurements of the tau polarization in W decays to a tau and neutrino and that in Z decays to tau lepton pairs in data collected by the ATLAS detector at the LHC in the years 2010 and 2012, respectively, comprise the subject of my PhD thesis. The code in this repository is a representative sample of the analysis framework that I wrote while working on the measurement of tau polarization in Z decays. This was the first analysis framework that I developed in Python. My main responsibilities on this analysis team were developing and validating the event selection, performing investigations of the comparability of the observed data and simulated distributions, developing an understanding of the generated polarization and kinematics in simulated signal samples, and, finally, developing and validating the fit model for the measurement. I helped to write and edit the paper that documents this analysis and is yet to be published. 

* base/Drawer.py - plotting class and plot instances
* base/Sample.py - Sample and Group classes to define file handling behavior; 
			      Region class to define analysis regions (defined by event selections / cuts )
* base/stream.py - stream class to join sample groups and regions
* base/Groups.py - instances of samples and groups
* base/streams.py - analysis region definitions and stream instances 
* base/trueTau.py - methods for working with tau truth event information 
 
* run/makeOut.py - main program for running analysis over data and simulation samples 
* run/run.py - definitions of analysis specs for different analysis settings

* fit/makeData.py - generate pseudo data for fit development from simulations
* fit/makeFit.py - produce fit with ROOT HistFactory framwork
* fit/toyFits.py - produce toy fits to pseudodata 
* fit/systematics.py - study the impact of systematic uncertainties in the fit model
* fit/Minuit/ - early fit development 

* draw/ - plotting code 
* tools/ - various helper code, mostly for formatting

* plots/ - some representative plots ( not official public plots, but those included in thesis ) 
 