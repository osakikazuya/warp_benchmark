#ReadMe file of warp simulation single particle simulations for nonlinear solenoid#

#Structure of warp_benchmark#

warp_benchmark/warp_numerical_error_converge_nonlinear_solenoid.pptx
			  /Fig_files/Change_dr/
					    /Change_dz/
                        /Change_nstep/					   
		      /warp_simulation_script/add_sol_grid_sample.sh
		   			                 /add_sol_grid_sample.py
					                 /Bgrid_data/

#Explanation of each directory and file#

There are two directories and one pptx file in warp_benchmark, "Fig_files", "warp_simulation_script", and "warp_numerical_error_converge_nonlinear_solenoid.pptx".

The "warp_numerical_error_converge_nonlinear_solenoid.pptx" is a power point file of systematic warp simulations.

We can see and make figures of "warp_numerical_error_converge_nonlinear_solenoid.pptx" in the "Fig_files" directory.
In the directory "Change_dr", I summarized simulation results of several transverse grid size of magnetic field data.
There are simulation data files 'Change_dr.tar.bz2', a mathematica file 'warp_numerical_error_nr1', and several figures in "Change_dr". 
In the directory "Change_dz", I summarized simulation results of several longitudinal grid size of magnetic field data.
There are simulation data files 'Change_dz.tar.bz2', a mathematica file 'warp_numerical_error_nz1', and several figures in "Change_dz". 
In the directory "Change_nstep", I summarized simulation results of several simulation step size.
There are simulation data files 'Change_nstep.tar.bz2', a mathematica file 'warp_numerical_error_nstep1', and several figures in "Change_nstep". 

In the "warp_simulation_script", there are a sample inputfile 'add_sol_grid_sample.py' of warp and a script file 'add_sol_grid_sample.sh'.
We should specify several simulation conditions in 'add_sol_grid_sample.sh' and we can conduct a simulation by "./add_sol_grid_sample.sh".
There are magnetic grid data files with several grid size in the directory "Bgrid_data", and when we conduct a simulation, we copy magnetic data files from here.