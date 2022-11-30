#/usr/bin/bash
EXAMPLE="OneDWave"
export commandFile="${EXAMPLE}.command.sh"
export decisionMakerCommandFile="${EXAMPLE}.command.python.sh"
export nextDecisionMakerCommandFile="${EXAMPLE}.command.python.ready.sh"
export checkResultFile="check_result.py"

cd build

BINARYDIR="$(pwd)/bin"

mkdir ${EXAMPLE}

echo "Setting up the python files for multi-point optimization."
cp ../utils/optimization_ver2/*.py ./${EXAMPLE}/
if [ $? -ne 0 ]; then exit -1; fi

mv ./${EXAMPLE}/base_extension.py ./${EXAMPLE}/base.py
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
ln -s ${BINARYDIR}/adjoint ./
ln -s ${BINARYDIR}/control_space_norm ./
ln -s ${BINARYDIR}/qfile_zaxpy ./
ln -s ${BINARYDIR}/spatial_inner_product ./
ln -s ${BINARYDIR}/zxdoty ./
ln -s ${BINARYDIR}/zwxmwy ./
ln -s ${BINARYDIR}/zaxpy ./
ln -s ${BINARYDIR}/slice_control_forcing ./
ln -s ${BINARYDIR}/paste_control_forcing ./
ln -s ${BINARYDIR}/patchup_qfile ./

echo "Setting up the baseline for ${EXAMPLE}"

./forward
if [ $? -ne 0 ]; then exit -1; fi

./control_space_norm
if [ $? -ne 0 ]; then exit -1; fi

echo "Setting up the directories for multi-point optimization."

python3 optimization.py -h

for dir in {a,b,c,x,x0}
do
  mkdir $dir
  for k in {0..23}
  do
    mkdir -p ${dir}/${k}
    ln -s ${BINARYDIR}/forward ./${dir}/${k}/
    ln -s ${BINARYDIR}/adjoint ./${dir}/${k}/
    cp ./magudi.multi-point.inp ./${dir}/${k}/magudi-${k}.inp
  done
done
if [ $? -ne 0 ]; then exit -1; fi

echo "Setting up the magudies for multi-point optimization."

python3 magudi-setup.py
if [ $? -ne 0 ]; then exit -1; fi
bash magudi-setup.sh
if [ $? -ne 0 ]; then exit -1; fi

echo "Setting up the initial conditions for multi-point optimization."

for k in {0..23}
do
  let "idx=20*$k"
  printf -v oldfile "OneDWave-%08d.q" $idx
  printf -v newfile "OneDWave-%d.ic.q" $k
  cp $oldfile x0/${newfile}
done
if [ $? -ne 0 ]; then exit -1; fi

echo "Initial forward and adjoint run."

python3 initial.py
if [ $? -ne 0 ]; then exit -1; fi
bash initial-forward.sh
if [ $? -ne 0 ]; then exit -1; fi
bash initial-adjoint.sh
if [ $? -ne 0 ]; then exit -1; fi

cat <<EOF > ${nextDecisionMakerCommandFile}
python3 optimization.py 1 -initial_cg -zero_baseline
EOF

export REF_ERROR="5.2613616409101998E-08"
cat <<EOF > ${checkResultFile}
from base import *
J, subJ = QoI('b')
error = abs(J - ${REF_ERROR})
if (error > 1.0e-10):
    raise RuntimeError('Error: %.15E' % error)
else:
    print('Optimization test passed.')
EOF

for k in {1..15}
do
    echo "Optimization Iteration $k"

    mv $nextDecisionMakerCommandFile $decisionMakerCommandFile
    bash $decisionMakerCommandFile
    if [ $? -ne 0 ]; then
        echo "$decisionMakerCommandFile is not run successfully."
        exit -1
    fi
    bash $commandFile
    export RESULT=$?
    if [ $RESULT -eq 1 ]; then
      echo "Optimization finished."
      python3 ${checkResultFile}
    elif [ $RESULT -ne 0 ]; then
      echo $RESULT
      echo "$commandFile is not run successfully."
      exit -1
    fi
done
