scons config=gnu_release.py layers=85 flux_solver=upwind -j4

dmax=14

command="bin/samoa_darcy_notasks_upwind_gnu_l85 -sections 1 -threads 4 -dmax "$dmax" -xmloutput -tout 86.4e4 -tmax 172.8e6 -epsilon 1.0e-5 -S_ref_th 0.5 -p_ref_th 5.0 -courant 1.0"

echo $command > "output/SPE10_d"$dmax".log"
$command | tee -a "output/SPE10_d"$dmax".log"

