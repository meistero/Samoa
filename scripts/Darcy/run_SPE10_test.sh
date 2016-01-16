compiler='gnu'
perm_averagings='arithmetic geometric harmonic'
fperm="data/darcy_five_spot/spe_perm_renamed.nc"
fpor="data/darcy_five_spot/spe_phi_renamed.nc"
epsilon=3.0e-5
S_ref_th=1.0e2
layers="85"

virtual_cores=$(lscpu | grep "^CPU(s)" | grep -oE "[0-9]+" | tr "\n" " ")
threads_per_core=$(lscpu | grep "^Thread(s) per core" | grep -oE "[0-9]+" | tr "\n" " ")
cores=$(( $virtual_cores / $threads_per_core ))

echo "" > list.tmp

for layers in $layers ; do
    for perm_averaging in $perm_averagings ; do
        exe="samoa_darcy_"$compiler"_l"$layers"_"$perm_averaging
        echo "scons config="$compiler"_release.py flux_solver='upwind' layers="$layers" perm_averaging="$perm_averaging" exe="$exe >> list.tmp
        #scons config=$compiler"_release.py" flux_solver='upwind' layers=$layers perm_averaging=$perm_averaging exe=$exe &
    done
done

cat list.tmp | parallel

echo "" > list.tmp

for layers in $layers ; do
    if [ $layers = 0 ] ; then
        p_ref_th=3.0e-9
    else
        p_ref_th=2.0e-8
    fi

    for dmax in {8..16..1} ; do
        for perm_averaging in $perm_averagings ; do
            exe="samoa_darcy_"$compiler"_l"$layers"_"$perm_averaging

            output_dir="output/SPE10test_"$compiler"_l"$layers"_"$perm_averaging"_d"$dmax

            command="bin/$exe -fperm $fperm -fpor $fpor -sections 2 -threads 1 -dmax $dmax -xmloutput -tout 86.4e4 -tmax 0.0 -epsilon $epsilon -S_ref_th $S_ref_th -p_ref_th $p_ref_th -courant 1.0 -output_dir $output_dir"

            echo "mkdir -p "$output_dir" ; echo "$command" > "$output_dir"/SPE10test_"$compiler"_l"$layers"_"$perm_averaging"_d"$dmax".log ; "$command" | tee -a "$output_dir"/SPE10test_"$compiler"_l"$layers"_"$perm_averaging"_d"$dmax".log" >> list.tmp
            #$command | tee -a $output_dir"/SPE10test_"$compiler"_l"$layers"_"$perm_averaging"_d"$dmax.log &
        done
    done
done

cat list.tmp | parallel

rm list.tmp
