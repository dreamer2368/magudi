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
wait $pid || let "FAIL+=1"
echo $pid
done

echo "Number of failures: $FAIL"
if [ $FAIL -eq 0 ]; then
echo "intermediate forward runs succeeded."
exit 0
else
echo "intermediate forward runs failed."
exit -1
fi
