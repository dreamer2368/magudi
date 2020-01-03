prefix="AcousticMonopole"
controlRegions=()
controlRegions+=("controlRegion")

rm previous.*.dat
rm ${prefix}.norm_${controlRegions[0]}.dat
rm ${prefix}.gradient_${controlRegions[0]}.dat
rm ${prefix}.control_forcing_${controlRegions[0]}.dat
rm ${prefix}.conjugate_gradient_${controlRegions[0]}.dat

for k in {0..2}
do
    rm ${k}/${prefix}-${k}.control_forcing_${controlRegions[0]}.dat
    rm ${k}/${prefix}-${k}.gradient_${controlRegions[0]}.dat
done

for dir in step
do
    rm ${prefix}.control_forcing_${controlRegions[0]}.dat
    for k in {0..2}
    do
        rm ${dir}/${k}/${prefix}-${k}.control_forcing_${controlRegions[0]}.dat
        rm ${dir}/${k}/${prefix}-${k}.gradient_${controlRegions[0]}.dat
    done
done
