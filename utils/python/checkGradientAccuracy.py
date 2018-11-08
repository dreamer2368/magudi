import numpy as np
import subprocess

prefix = 'AcousticMonopole'
forwardFile = prefix+'.forward_run.txt'
adjointFile = prefix+'.adjoint_run.txt'
numProcs = 32

#subprocess.check_call('mpirun -np '+str(numProcs)+' ./forward '+forwardFile, shell=True)
fQoI0 = open(forwardFile,'r')
QoI0 = float(fQoI0.read())
fQoI0.close()
print "QoI0: ", QoI0

#subprocess.check_call('mpirun -np '+str(numProcs)+' ./adjoint '+adjointFile, shell=True)
fGrad0 = open(adjointFile,'r')
Grad0 = float(fGrad0.read())
fGrad0.close()
print "Grad0: ", Grad0

Nk = 24
Ak = 10.0**(-0.25*np.array(range(Nk)))
QoIk = np.zeros((Nk,),dtype=np.double)
Gradk = np.zeros((Nk,),dtype=np.double)
ek = np.zeros((Nk,),dtype=np.double)

forwardFile = prefix+'.forward_run1.txt'
#baseline_control_forcing_file = prefix+'.control_forcing_controlRegion0.dat'
baseline_control_forcing_file = ''
grad_file = prefix+'.gradient_controlRegion.dat'
control_forcing_file = prefix+'.control_forcing_controlRegion.dat'
for k in range(Nk):
    actuation_amount = Ak[k]
    command = './zaxpy '+control_forcing_file+' '+str(actuation_amount)+' '+grad_file+' '+baseline_control_forcing_file
    subprocess.check_call(command, shell=True)
    subprocess.check_call('srun -n '+str(numProcs)+' ./forward '+forwardFile, shell=True)
    fQoI0 = open(forwardFile,'r')
    QoIk[k] = float(fQoI0.read())
    fQoI0.close()

    Gradk[k] = (QoIk[k]-QoI0)/Ak[k]
    ek[k] = abs( (Gradk[k]-Grad0)/Grad0 )
    print actuation_amount, QoIk[k], Gradk[k], ek[k]

fId = open(prefix+'.gradient_accuracy.txt','w')
for k in range(Nk):
    fId.write(str(Ak[k])+"\t"+str(QoIk[k])+"\t"+str(Gradk[k])+"\t"+str(ek[k])+"\n")
fId.close()
