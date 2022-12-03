#/usr/bin/bash
EXAMPLE="OneDWave"
export commandFile="${EXAMPLE}.command.sh"

cd build

BINARYDIR="$(pwd)/bin"

mkdir ${EXAMPLE}

echo "Setting up the python files for multi-point optimization."
cp ../utils/optimization_ver3/*.py ./${EXAMPLE}/
if [ $? -ne 0 ]; then exit -1; fi

echo "Copying input files for ${EXAMPLE}"

cp -r ../examples/${EXAMPLE}/* ./${EXAMPLE}/
if [ $? -ne 0 ]; then exit -1; fi

cp ../utils/python/plot3dnasa.py ./${EXAMPLE}/
if [ $? -ne 0 ]; then exit -1; fi

echo "Generating grid and solutions for ${EXAMPLE}"

cd ${EXAMPLE}
python2 config.py
if [ $? -ne 0 ]; then exit -1; fi

echo "Generating links to executables in ${EXAMPLE}"

ln -s ${BINARYDIR}/forward ./
ln -s ${BINARYDIR}/control_space_norm ./

echo "Setting up the baseline for ${EXAMPLE}"

./forward
if [ $? -ne 0 ]; then exit -1; fi

./control_space_norm
if [ $? -ne 0 ]; then exit -1; fi

echo "Setting up the directories for multi-point optimization."

python3 optimization.py optim.yml --mode setup
if [ $? -ne 0 ]; then exit -1; fi

bash ${commandFile}
if [ $? -ne 0 ]; then exit -1; fi

echo "Check gradient accuracy: checkGradientAccuracy.py"
python3 checkGradientAccuracy.py optim.yml
