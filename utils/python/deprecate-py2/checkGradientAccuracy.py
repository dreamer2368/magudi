import numpy as np
import subprocess

prefix = 'MultiblockJet'
forwardFile = prefix+'.forward_run.txt'
adjointFile = prefix+'.adjoint_run.txt'
numNodes = 30
numProcs = 1056

#subprocess.check_call('srun -N '+str(numNodes)+' -n '+str(numProcs)+' ./forward '+forwardFile, shell=True)
fQoI0 = open(forwardFile,'r')
QoI0 = float(fQoI0.read())
fQoI0.close()
print "QoI0: ", QoI0

#subprocess.check_call('srun -N '+str(numNodes)+' -n '+str(numProcs)+' ./adjoint '+adjointFile, shell=True)
fGrad0 = open(adjointFile,'r')
Grad0 = float(fGrad0.read())
fGrad0.close()
print "Grad0: ", Grad0

Nk = 10
Ak = 10.0**(-3.0-2.5-0.25*np.array(range(Nk)))
QoIk = np.zeros((Nk,),dtype=np.double)
Gradk = np.zeros((Nk,),dtype=np.double)
ek = np.zeros((Nk,),dtype=np.double)

forwardFile = prefix+'.forward_run.txt'
#baseline_control_forcing_file = [prefix+'.control_forcing_controlRegion0.E.dat',
#                                 prefix+'.control_forcing_controlRegion0.N.dat',
#                                 prefix+'.control_forcing_controlRegion0.S.dat',
#                                 prefix+'.control_forcing_controlRegion0.W.dat']
baseline_control_forcing_file = ['','','','']
grad_file = [prefix+'.gradient_controlRegion.E.dat',
             prefix+'.gradient_controlRegion.N.dat',
             prefix+'.gradient_controlRegion.S.dat',
             prefix+'.gradient_controlRegion.W.dat']
control_forcing_file = [prefix+'.control_forcing_controlRegion.E.dat',
                        prefix+'.control_forcing_controlRegion.N.dat',
                        prefix+'.control_forcing_controlRegion.S.dat',
                        prefix+'.control_forcing_controlRegion.W.dat']

fId = open(prefix+'.gradient_accuracy.txt','w')
for k in range(Nk):
    actuation_amount = Ak[k]
#    for i in range(4):
#        command = 'msub ./ZAXPY.sh '+str(k+1)+'/'+control_forcing_file[i]+' '+"{:.16E}".format(actuation_amount)+' '+grad_file[i]+' '+baseline_control_forcing_file[i]
#        subprocess.check_call(command, shell=True)
    fQoI0 = open(str(k+1)+'/'+forwardFile,'r')
    QoIk[k] = float(fQoI0.read())
    fQoI0.close()

    Gradk[k] = (QoIk[k]-QoI0)/Ak[k]
    ek[k] = abs( (Gradk[k]-Grad0)/Grad0 )
    fId.write("{:.16E}".format(actuation_amount)+'\t'
                +"{:.16E}".format(QoIk[k])+'\t'
                +"{:.16E}".format(Gradk[k])+'\t'
                +"{:.16E}".format(ek[k]))
    print "{:.16E}".format(actuation_amount), "{:.16E}".format(QoIk[k]), "{:.16E}".format(Gradk[k]), "{:.16E}".format(ek[k])

fId.close()
