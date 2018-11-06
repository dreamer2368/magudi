import numpy as np
import subprocess

prefix = 'AcousticMonopole'
forwardFile = prefix+'.forward_run.txt'
adjointFile = prefix+'.adjoint_run.txt'
numProcs = 2

subprocess.check_call('mpirun -np '+str(numProcs)+' ./forward '+forwardFile, stdout=open('test.out', 'w'), stderr=subprocess.STDOUT, shell=True)
fQoI0 = open(forwardFile,'r')
QoI0 = float(fQoI0.read())
fQoI0.close()
print "QoI0: ", QoI0

subprocess.check_call('mpirun -np '+str(numProcs)+' ./adjoint '+adjointFile, stdout=open('adjoint_test.out', 'w'), stderr=subprocess.STDOUT, shell=True)
fGrad0 = open(adjointFile,'r')
Grad0 = float(fGrad0.read())
fGrad0.close()
print "Grad0: ", Grad0

actuation_amount = 1.0
growth = 10.0**(-0.25)
QoIk = np.zeros((32,),dtype=np.double)

grad_file = prefix+'.gradient_controlRegion.dat'
control_forcing_file = prefix+'.control_forcing_controlRegion.dat'
for k in range(32):
    command = './zaxpy '+control_forcing_file+' '+str(actuation_amount)+' '+grad_file
    subprocess.check_call(command, stdout=open('zaxpy_test.out', 'w'), stderr=subprocess.STDOUT, shell=True)
    subprocess.check_call('mpirun -np '+str(numProcs)+' ./forward '+forwardFile, stdout=open('test.out', 'w'), stderr=subprocess.STDOUT, shell=True)
    fQoI0 = open(forwardFile,'r')
    QoIk[k] = float(fQoI0.read())
    fQoI0.close()
    print actuation_amount, QoIk[k]
    actuation_amount *= growth

Ak = 10.0**(-0.25*np.array(range(32)))
Gradk = (QoIk-QoI0)/Ak
ek = abs( (Gradk-Grad0)/Grad0 )
fId = open(prefix+'.gradient_accuracy.txt','w')
for k in range(32):
    fId.write(str(Ak[k])+"\t"+str(QoIk[k])+"\t"+str(Gradk[k])+"\t"+str(ek[k])+"\n")
fId.close()
