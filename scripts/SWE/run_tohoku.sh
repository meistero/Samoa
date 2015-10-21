compiler='gnu'
fdispl="data/tohoku_static/displ.nc"
fbath="data/tohoku_static/bath_2014.nc"
max_depths='18 20 22 24 26'
flux_solvers='llfbath fwave aug_riemann'

echo "SWE scenario exeuction script."
echo ""
echo "compiler          : "$compiler
echo "bathymetry file   : "$fdispl
echo "displacement file : "$fbath
echo "max depths        : "$max_depths
echo "flux solvers      : "$flux_solvers
echo ""

echo ""
echo "Compiling (skip with Ctrl-C)..."
echo ""

exe_base='bin/samoa_swe'

if [ $compiler = "intel" ] ; then
    for flux_solver in $flux_solvers
    do
        scons config=intel_release.py scenario=swe flux_solver=$flux_solver -j4
    done

    comp_suffix=''
else
    for flux_solver in $flux_solvers
    do
        scons config=gnu_release.py scenario=swe flux_solver=$flux_solver -j4
    done

    exe_base=$exe_base'_notasks'
    comp_suffix='_gnu'
fi


echo ""
echo "Running..."
echo ""

for max_depth in $max_depths
do
    for flux_solver in $flux_solvers
    do

        if [ $flux_solver = "aug_riemann" ] ; then
            exe=$exe_base
        else
            exe=$exe_base'_'$flux_solver''$comp_suffix
        fi

        command=$exe' -sections 1 -threads 4 -tout 20.0 -dmin 0 -dmax '$max_depth' -tmax 10.8e3 -fdispl '$fdispl' -fbath '$fbath' -stestpoints "545735.266126 62716.4740303,935356.566012 -817289.628677,1058466.21575 765077.767857"'

        echo $command > "output/tohoku_"$flux_solver"_d"$max_depth".log"
        $command | tee -a "output/tohoku_"$flux_solver"_d"$max_depth".log"
        
        mkdir -p "output/tohoku_"$flux_solver"_d"$max_depth
        mv output/tohoku*.log "output/tohoku_"$flux_solver"_d"$max_depth"/"
        mv output/swe* "output/tohoku_"$flux_solver"_d"$max_depth"/"
    done
done
