Geant4Medipix 
-------------

 Pixel based Medipix geometry with digitiser for simulating Monte Carlo and charge transport simulation as
 well as chip electronics.

## REFERENCES

 Feel free to use and modify the code. If you intend to publish, please use these references:

 [1] A. Schübel, D. Krapohl, E. Fröjdh, C. Fröjdh, and G. Thungström, “A Geant4 based framework for pixel  
 detector simulation,” J. Instrum., vol. 9, no. 12, 2014.
 http://dx.doi.org/10.1088/1748-0221/9/12/C12018

 [2] D. Krapohl, A. Schubel, E. Fröjdh, G. Thungstrom, and C. Fröjdh, “Validation of Geant4 Pixel Detector  
 Simulation Framework by Measurements with the Medipix Family Detectors,” IEEE Trans. Nucl. Sci., vol. 63, 
 no. 3, pp. 1874–1881, 2016.
 http://dx.doi.org/10.1109/TNS.2016.2555958
 
## Authors
A. Schübel, David Kraphol, Erik Fröjdh, Lucas Roussel, Paul van Schayck

## License
MIT 

## GEOMETRY DEFINITION

standard setup 
```
╔═════════════════════════════════════════════════════════════════════════════╗
║                                                            ║                ║
║                                                         ░:█║                ║
║                                                         ░:█║<-Medipix       ║
║                                                         ░:█║  -sensor       ║
║                 ║                 ╔═╗                   ░:█║  -bumps(false) ║
║ beam            ║                 ║ ║                   ░:█║  -chip(false)  ║
║ ======>         ║                 ╚═╝                   ░:█║  -pcb (false)  ║
║                 ║                  ^                    ░:█║                ║
║                 ^                  Object               ░:█║                ║
║                 Filter(or XRF)                          ░:█║                ║
║                                                         ░:█║                ║
║                                                            ║                ║
╚═════════════════════════════════════════════════════════════════════════════╝
```
The definition of the Medipix detector is flexible in terms of pixel size
and number of pixels
 
## REQUIREMENTS
 
Geant4-10.00 or higher is required in order to profit from multi-threading
The Boost library is used to read ini-files
HDF5 is used to write output files

CAD libraries are optional

	    
## PHYSICS LISTS
 
Physics lists can be local (eg. in this example) or from G4 kernel
physics_lists subdirectory.

From geant4/source/physics_lists/builders:	 
- "emstandard_opt0" recommended standard EM physics for LHC
- "emstandard_opt1" best CPU performance standard physics for LHC
- "emstandard_opt2"     
- "emstandard_opt3" best current advanced EM options. 
- "emlivermore"  low-energy EM physics using Livermore data
- "empenelope"   low-energy EM physics implementing Penelope models
- "dna"	     low-energe DNA physics    
Physics lists and options can be (re)set with UI commands

Please, notice that options set through G4EmProcessOPtions are global, eg
for all particle types. In G4 builders, it is shown how to set options per
particle type.
 				
## VISUALIZATION
 
  The Visualization Manager is set in the main().
  The initialisation of the drawing is done via the commands :
  /vis/... in the macro vis.mac. In interactive session:
  PreInit or Idle > /control/execute vis.mac
 	
  The default view is a longitudinal view of the calorimeter.
 	
  The tracks are drawn at the end of event, and erased at the end of run.
  Optionaly one can choose to draw all particles, only the charged one, or none.
  This command is defined in EventActionMessenger class.

## HOW TO START ?
  - Install Geant4 and configure the environment (source pathto/geant4.sh)
  - Create a build folder and execute 
  ```
    % cmake ..
    % make -j <nothreads>
```
  - Execute G4Medipix in 'batch' mode from macro files
 ```
 % G4Medipix  -m gamma.mac -t <nothreads>
```
 		
  - Execute G4Medipix in 'interactive mode' with visualization
 ```
% G4Medipix
....
Idle> type your commands. For instance:
Idle> /control/execute gamma.mac
....
Idle> exit
 ```   
There are example mac files in the macro folder as well as tube spectra
  
## OUTPUT FILES
 
G4Medipix can produce HDF5 files containing a `/g4medipix` and `/trajectories` group containing the pixel output 
matrix and full trajectories of the incident particle.   

 

 