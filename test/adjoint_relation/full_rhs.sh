echo "full_rhs test."
echo "working dir:" $(pwd)
echo "binary dir:" $1

exampleDir="AcousticMonopole"

cp -r ../../examples/${exampleDir} ./
cp ../../utils/python/plot3dnasa.py ./${exampleDir}
cd ${exampleDir}
python2 config.py
mpirun -n 2 $1/full_rhs
if [ $? -eq 0 ]; then
	cd ..
	rm -rf ${exampleDir}
else
	exit -1
fi
