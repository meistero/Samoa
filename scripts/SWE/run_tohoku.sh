max_depths='18 20 22 24 26'
flux_solvers='llfbath aug_riemann'
fdispl="data/tohoku_static/displ.nc"
fbath="data/tohoku_static/bath_2014.nc"

echo "SWE scenario exeuction script."
echo ""
echo "bathymetry file   : "$fdispl
echo "displacement file : "$fbath
echo "max depths        : "$max_depths
echo "flux solvers      : "$flux_solvers
echo ""

echo ""
echo "Compiling (skip with Ctrl-C)..."
echo ""

for flux_solver in $flux_solvers
do
    scons config=gnu_release.py scenario=swe flux_solver=$flux_solver -j4
    
    for max_depth in $max_depths
    do
        mkdir -p "output/tohoku_"$flux_solver"_d"$max_depth
    done
done

echo ""
echo "Running..."
echo ""

for max_depth in $max_depths
do
    for flux_solver in $flux_solvers
    do
        if [ $flux_solver = "aug_riemann" ] ; then
            exe='bin/samoa_swe_notasks_gnu'
        else
            exe='bin/samoa_swe_notasks_'$flux_solver'_gnu'
        fi

        command=$exe' -sections 1 -threads 4 -tout 20.0 -dmax '$max_depth' -tmax 10.8e3 -fdispl '$fdispl' -fbath '$fbath' -stestpoints "545735.266126 62716.4740303,935356.566012 -817289.628677,1058466.21575 765077.767857"'

        echo $command > "output/tohoku_"$flux_solver"_d"$max_depth"/tohoku_"$flux_solver"_d"$max_depth".log"
        $command | tee -a "output/tohoku_"$flux_solver"_d"$max_depth"/tohoku_"$flux_solver"_d"$max_depth".log"
        mv output/swe* "output/tohoku_"$flux_solver"_d"$max_depth"/"
    done
done
