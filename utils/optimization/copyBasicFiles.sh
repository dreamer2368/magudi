for i in {1..13}
do
cp MultiblockJet.ic.q ${i}/
cp MultiblockJet.target.q ${i}/
cp MultiblockJet.*_mollifier.f ${i}/
cp MultiblockJet*.eigenmode_*.q ${i}/
cp bc.dat ${i}/
cp magudi.inp ${i}/
cp quartz.sh ${i}/
cp forward ${i}/
cp MultiblockJet.mean_pressure.f ${i}/
cp decomp.map ${i}/
cp MultiblockJet.xyz ${i}
cp MultiblockJet.inflow.xyz ${i}
done
