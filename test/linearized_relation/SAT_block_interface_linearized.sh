echo "SAT_block_interface script test."
echo "working dir:" $(pwd)
echo "binary dir:" $1

cd ../adjoint_relation/block_interface_inputs
$1/SAT_block_interface_linearized
