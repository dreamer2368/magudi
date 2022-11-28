EXAMPLE="OneDWave"

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

mkdir out grad diff txt

for dir in {x,x0}
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

echo "Check gradient accuracy: checkGradientAccuracy.py"
python3 checkGradientAccuracy.py
