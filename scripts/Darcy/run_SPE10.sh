compiler='intel'
perm_averagings='geometric'
fperm="data/darcy_five_spot/spe_perm_renamed.nc"
fpor="data/darcy_five_spot/spe_phi_renamed.nc"
layers=85
max_depths='14'
epsilon=1.0e-4
S_ref_th=1.0e2

if [ $layers = 0 ] ; then
    p_ref_th=3.0e-9
else
    p_ref_th=2.0e-8
fi

virtual_cores=$(lscpu | grep "^CPU(s)" | grep -oE "[0-9]+" | tr "\n" " ")
threads_per_core=$(lscpu | grep "^Thread(s) per core" | grep -oE "[0-9]+" | tr "\n" " ")
cores=$(( $virtual_cores / $threads_per_core ))

echo "SPE10 scenario execution script."
echo ""
echo "compiler          : "$compiler
echo "max depths        : "$max_depths
echo "perm. averaging   : "$perm_averagings
echo "perm. file        : "$fperm
echo "porosity file     : "$fpor
echo "layers            : "$layers
echo ""

echo ""
echo "Compiling (skip with Ctrl-C)..."
echo ""

for perm_averaging in $perm_averagings
do
    exe="samoa_darcy_"$compiler"_l"$layers"_"$perm_averaging
    scons config=$compiler"_release.py" flux_solver='upwind' layers=$layers perm_averaging=$perm_averaging exe=$exe -j4      
done

echo "Running..."

for dmax in $max_depths
do
    for perm_averaging in $perm_averagings
    do
        exe="samoa_darcy_"$compiler"_l"$layers"_"$perm_averaging

        output_dir="output/SPE10_"$compiler"_l"$layers"_"$perm_averaging"_d"$dmax
        mkdir -p $output_dir

        command="bin/$exe -fperm $fperm -fpor $fpor -sections 1 -threads $cores -dmax $dmax -xmloutput -tout 86.4e3 -tmax 172.8e6 -epsilon $epsilon -S_ref_th $S_ref_th -p_ref_th $p_ref_th -courant 1.0 -output_dir $output_dir"

        echo $command > $output_dir"/SPE10_"$compiler"_l"$layers"_"$perm_averaging"_d"$dmax.log
        $command | tee -a $output_dir"/SPE10_"$compiler"_l"$layers"_"$perm_averaging"_d"$dmax.log
    done
done
