echo "control_functional script test."
echo "working dir:" $(pwd)
echo "binary dir:" $1

cd control_functional
$1/control_functional
if [ $? -eq 0 ]; then
	rm test-*.q test.*.dat
else
	exit -1
fi
