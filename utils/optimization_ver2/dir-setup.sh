prefix="MultiblockJet"
controlRegions=()
controlRegions+=("controlRegion.E")
controlRegions+=("controlRegion.W")
controlRegions+=("controlRegion.N")
controlRegions+=("controlRegion.S")

BUILDS="~/Builds/magudi-release/bin"

for k in {0..29}
do
    cp bc.dat ${k}/
    cp decomp.map ${k}/
    cp magudi.inp ${k}/magudi-${k}.inp
    
    cp ${prefix}.xyz ${k}/ &
    cp ${prefix}.target.q ${k}/ &
    cp ${prefix}.control_mollifier.f ${k}/ &
    cp ${prefix}.target_mollifier.f ${k}/ &
    cp ${prefix}.mean_pressure.f ${k}/ &

    for j in {1..6}
    do
        cp ${prefix}-0${j}.eigenmode_real.q ${k}/${prefix}-${k}-0${j}.eigenmode_real.q &
        cp ${prefix}-0${j}.eigenmode_imag.q ${k}/${prefix}-${k}-0${j}.eigenmode_imag.q &
    done

    ln -s ${BUILDS}/forward ${k}/
    ln -s ${BUILDS}/adjoint ${k}/
done

for dir in {a,b,c,x}
do
    for k in {0..29}
    do
        cp bc.dat ${dir}/${k}/
        cp decomp.map ${dir}/${k}/
        cp magudi.inp ${dir}/${k}/magudi-${k}.inp
        
        cp ${prefix}.xyz ${dir}/${k}/ &
        cp ${prefix}.target.q ${dir}/${k}/ &
        cp ${prefix}.control_mollifier.f ${dir}/${k}/ &
        cp ${prefix}.target_mollifier.f ${dir}/${k}/ &
        cp ${prefix}.mean_pressure.f ${dir}/${k}/ &
    
        for j in {1..6}
        do
            cp ${prefix}-0${j}.eigenmode_real.q ${dir}/${k}/${prefix}-${k}-0${j}.eigenmode_real.q &
            cp ${prefix}-0${j}.eigenmode_imag.q ${dir}/${k}/${prefix}-${k}-0${j}.eigenmode_imag.q &
        done
    
        ln -s ${BUILDS}/forward ${dir}/${k}/
        ln -s ${BUILDS}/adjoint ${dir}/${k}/
    done
done
