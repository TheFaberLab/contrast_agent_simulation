# IONP Simulations
Still in developement phase!

#What is needed to run the Simulations

1) C++ compiler (C++17 or newer)
2) cmake (min 3.16.)

#How to compile the Simulations

1) create build folder in the simulation folder
2) open terminal in build folder and run:
      cmake ../
3) then run:
      make
4) after that the C++ program is compiled and the executable file is in "/simulation/build/src/ionp/" and is called "IONP"

#How to set parameters for a Simulation

Write the parameters in .txt/.dat file and put the file in the "simulation/ini" folder
Important: All .txt or .dat files in this folder are read for simulations - it is possible to put multiple files in this folder to run multiple simulations one after the other!!!

Possible parameters are: Everything in SI-Units (except concentration)!

- Ionp radius (double):                                                                                                 IONP_r
- Standard deviation of IONP radius (lognormal distribution) (double):                                                  IONP_r_std
- Thickness of IONP coating (double):                                                                                   Coating
- Standard deviation of IONP coating (lognormal distribution) (double):                                                 Coating_std
- Aggregate radius (double):                                                                                            Agg_r
- Standard deviation of Aggregate radius (lognormal distribution) (double):                                             Agg_r_std
- Origin Position (double):                                                                                             Origin or Matrix_x0, Matrix_y0, Matrix_z0
- Simulation edge length (double):                                                                                      Length or Matrix_xlength, Matrix_ylength, Matrix_zlength
- Amount of Nodes (for B-Field) (int):                                                                                  Nodes
- IONP number (int):                                                                                                    IONP_n
- Aggregate number (int):                                                                                               Agg_n
- Iron concentration in simulation volume (in µmol/L) (double):                                                         Conc
- Proton number (int):                                                                                                  Proton_n
- Duration of the Simulation (double):                                                                                  echo_time
- Size of single time steps (double):                                                                                   time_step
- Echo spacing (only for multi SE) (double):                                                                            echo_spacing
- Amounts of threads used (if none given uses maximum threads available!!) (int):                                       Threads
- T2-constant of surrounding (double):                                                                                  T2
- Second IONP radius for bimodal lognormal distribution (double):                                                       IONP_r_1
- Standard deviation of second IONP radius (bimodal lognormal distribution) (double):                                   IONP_r_std_1
- For bimodal lognormal distribution of radii probability (%) of a radius to be drwan from first distribution (double): IONP_per
- Ionp magnetic moment (double):                                                                                        magnetic_moment
- Which MR sequence is used (FID, SE, MSE):                                                                             Sequence
- Is an DLS file for the radius distribution available or not (bool):                                                   DLS_file

- Protons that are bound to IONP and don't diffuse:
      -> Number of protons bound to IONP (int):                                                                         Bound_Prot_n

- Additional option for simulations of IONP movement:
      -> velocity of IONP (double):                                                                                     vIONP
      -> movement direction of IONP (as Vector double x y z):                                                           Direction

- Option to define a smaller cube inside simulation volume where signal is calculated:
      -> Origin of the cube (vector double):                                                                            Cube_Orig
      -> Size of the cube edges (double):                                                                               Cube_size


Not all parameters have to be defined! For Example:
- General: Not specifically defined parameters are 0
- No IONP number: The IONP_pos file for manual placement of IONP is used
- No aggregate number but a aggregate radius: Agg_pos file used for manual Aggregate placement
- No aggregate radius but a aggregate number: Aggregates are placed as dense Aggregates
- No nodes: B-Field will be calculated at every diffusion step of every proton (time consuming)
- No time steps: Variable time steps are used - smaller if protons are close to IONP

Under "/include/mptmacros.h" parameters are defined that are usually not changes that often (as global variables):
- Gyromagnetic ratio
- Diffusion constant
- IONP magnetization (used if no magnetic moment is given)

Examples for ini files of usual simulations in /ini!


In "simulation/ini/IONP_pos" it is possible to define the position of the IONPs manually:

- Put a .txt file in this folder with the positions of the IONP like
  x1 y1 z1
  x2 y2 z2
  ...

The same concept for "simulation/ini/Agg_pos", where you can define position of Aggregates manually

In "simulation/ini/DLS_File" it is possible to include a Histogramm file that contains a distribution of IONP radii

- Put a .txt file in this folder with columns radius, intensity like 
      r1 i1
      r2 i2
      ...


Output is in "simulation/data" folder with a log in "simulation/data/log"
and the actual files in "simulation/data/tmp"

#How to run Simulation

1) run the executable IONP file in folder "/simulation/build/src/ionp/" (depending on operating system what kind of file)

#Contact

lauritz.kluender@uni-muenster.de