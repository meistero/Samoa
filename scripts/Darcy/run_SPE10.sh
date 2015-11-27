layers=32
dmax=16
epsilon=3.0e-5
perm_averaging='harmonic'

echo "Compiling (skip with Ctrl-C)..."

scons config=gnu_release.py flux_solver='upwind' layers=$layers perm_averaging=$perm_averaging -j4

virtual_cores=$(lscpu | grep "^CPU(s)" | grep -oE "[0-9]+" | tr "\n" " ")
threads_per_core=$(lscpu | grep "^Thread(s) per core" | grep -oE "[0-9]+" | tr "\n" " ")
cores=$(( $virtual_cores / $threads_per_core ))

command="bin/samoa_darcy_notasks_upwind_gnu_l"$layers" -sections 1 -threads "$cores" -dmax "$dmax" -xmloutput -tout 86.4e4 -tmax 345.6e5 -epsilon $epsilon -S_ref_th 0.5 -p_ref_th 5.0e-9 -courant 0.5"

echo "Running..."

echo $command > "output/SPE10_l"$layers"_d"$dmax"_"$perm_averaging".log"
$command | tee -a "output/SPE10_l"$layers"_d"$dmax"_"$perm_averaging".log"

