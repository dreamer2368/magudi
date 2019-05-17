#!/bin/bash
FAIL=0

nodeList=()
while IFS= read -r line 
do
nodeList+=("$line")
done < <( scontrol show hostnames ${SLURM_JOB_NODELIST} )

for i in {1..10}
do
    let "nodeIndex=30*${i}-30"
    nodeListString="${nodeList[${nodeIndex}]}"
    for j in {1..29}
    do
        let "nodeIndex=30*${i}-30+${j}"
        nodeListString+=",${nodeList[${nodeIndex}]}"
    done

    cd ${i}/
    rm *.forward_run.txt
    echo "srun -N 30 -n 1056 -w ${nodeListString} ./forward &"
    srun -N 30 -n 1056 -w ${nodeListString} ./forward &>result.out &
    pids[${i}]=$!
    cd ../
done

for pid in ${pids[*]}
do
echo $pid
wait $pid || let "FAIL+=1"
done

echo $FAIL
if [ $FAIL -gt 0 ]
then
exit -1
fi
