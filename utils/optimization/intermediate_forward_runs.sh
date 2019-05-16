FAIL=0

for i in {1..10}
do
cd ${i}/
rm *.forward_run.txt
srun -N 1 -n 36 ./forward &
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
