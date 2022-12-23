#/usr/bin/bash
EXAMPLE="OneDWave"
export commandFile="${EXAMPLE}.command.sh"
export checkResultFile="check_result.py"

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
python3 config.py
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

export REF_ERROR="5.2613616409101998E-08"
cat <<EOF > ${checkResultFile}
from inputs import InputParser
from optimizer import Optimizer
config = InputParser('./optim.yml')
optim = Optimizer(config)

J, subJ = optim.base.QoI('b')
error = abs(J - ${REF_ERROR})
if (error > 1.0e-10):
    raise RuntimeError('Error: %.15E' % error)
else:
    print('Optimization test passed.')
EOF

for k in {1..15}
do
    echo "Optimization Iteration $k"

    python3 optimization.py optim.yml --mode schedule
    if [ $? -ne 0 ]; then
        echo "Scheduling is not run successfully."
        exit -1
    fi
    bash $commandFile
    export RESULT=$?
    if [ $RESULT -eq 1 ]; then
      echo "Optimization finished."
      python3 ${checkResultFile}
      if [ $? -eq 0 ]; then exit 0; else exit -1; fi
    else
      echo $RESULT
      python3 optimization.py optim.yml --mode log_result --result $RESULT
    fi
done
