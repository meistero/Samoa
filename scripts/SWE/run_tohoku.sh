compiler='gnu'
max_depths='18 20 22 24 26'
fdispl="data/tohoku_static/displ.nc"
fbath="data/tohoku_static/bath_2014.nc"
flux_solvers='llfbath fwave aug_riemann'
vtk='.true.'

virtual_cores=$(lscpu | grep "^CPU(s)" | grep -oE "[0-9]+" | tr "\n" " ")
threads_per_core=$(lscpu | grep "^Thread(s) per core" | grep -oE "[0-9]+" | tr "\n" " ")
cores=$(( $virtual_cores / $threads_per_core ))

echo "Tsunami scenario execution script."
echo ""
echo "compiler          : "$compiler
echo "max depths        : "$max_depths
echo "flux solvers      : "$flux_solvers
echo "bathymetry file   : "$fdispl
echo "displacement file : "$fbath
echo ""

echo ""
echo "Compiling (skip with Ctrl-C)..."
echo ""

for flux_solver in $flux_solvers
do
    exe="samoa_swe_"$compiler"_"$flux_solver
    scons config=$compiler"_release.py" scenario=swe flux_solver=$flux_solver exe=$exe -j4
done


echo ""
echo "Running..."
echo ""

for max_depth in $max_depths
do
    for flux_solver in $flux_solvers
    do
        exe="samoa_swe_"$compiler"_"$flux_solver

        if  [ $max_depth -gt '22' ] ; then
            vtk='.false.'
        fi

        output_dir="output/tohoku_"$compiler"_"$flux_solver"_d"$max_depth
        mkdir -p $output_dir

        command='bin/'$exe' -sections 1 -threads '$cores' -courant 0.95 -tout 20.0 -dmin 0 -dmax '$max_depth' -tmax 10.8e3 -fdispl '$fdispl' -fbath '$fbath' -xmloutput '$vtk' -stestpoints "545735.266126 62716.4740303,935356.566012 -817289.628677,1058466.21575 765077.767857" -output_dir '$output_dir

        echo $command > $output_dir"/tohoku_"$compiler"_"$flux_solver"_d"$max_depth".log"
        $command | tee -a $output_dir"/tohoku_"$compiler"_"$flux_solver"_d"$max_depth".log"
    done
done
