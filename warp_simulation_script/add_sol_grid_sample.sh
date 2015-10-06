#!/bin/sh 

for bz in 1.0 # center of magnetic field [T]
do
for ns in 500 # step number of simulation
do
for na in 1 # number of advance
do
for gnr in 100 #grid number in warp, not important this time
do
for gnz in 100 #grid number in warp, not important this time
do
for data in Lz20_nr101_nz3001 # set grid data
do
for vz in 3e7 # initial velosity of longitudinal direction
do

rm *dat
rm *cgm*

cp Bgrid_data/B_table_R0.08_L0.3_"$data".dat .

sed s/BZFIELD/"$bz"/g  add_sol_grid_sample.py > add_sol_grid1.py

sed s/NSTEP/"$ns"/g  add_sol_grid1.py > add_sol_grid2.py

sed s/NADVANCE/"$na"/g  add_sol_grid2.py > add_sol_grid3.py

sed s/GRIDZ/"$gnz"/g  add_sol_grid3.py > add_sol_grid4.py

sed s/GRIDR/"$gnr"/g  add_sol_grid4.py > add_sol_grid5.py

sed s/B_DATA/"$data"/g  add_sol_grid5.py > add_sol_grid6.py

sed s/VELOSITY/"$vz"/g  add_sol_grid6.py > add_sol_grid7.py

/usr/local/bin/python add_sol_grid7.py

tar jcf Solenoid_add_grid_"$data"_Vz"$vz"_Bz"$bz"_NStep"$ns"_NAd"$na".tar.bz2 orbit*dat *cgm add_sol_grid7.py

rm *dat
rm *cgm*

done
done
done
done
done
done
done
