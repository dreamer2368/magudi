import numpy as np
import subprocess

prefix = 'AcousticMonopole'
forwardFile = prefix+'.forward_run.txt'
deltaFile = prefix+'.linearized.txt'
forwardFile1 = prefix+'.forward_run.1.txt'
numNodes = 1
numProcs = 4
numTimesteps = 480

filename = '%s-%08d'%(prefix,numTimesteps)

baseline = False

forwardFile = prefix+'.forward_run.txt'
#baseline_control_forcing_file = [prefix+'.control_forcing_controlRegion0.E.dat',
#                                 prefix+'.control_forcing_controlRegion0.N.dat',
#                                 prefix+'.control_forcing_controlRegion0.S.dat',
#                                 prefix+'.control_forcing_controlRegion0.W.dat']
baseline_control_forcing_file = ['']
grad_file = [prefix+'.gradient_controlRegion.dat']
control_forcing_file = [prefix+'.control_forcing_controlRegion.dat']

if (baseline):
    subprocess.check_call('cp %s.ic.0.q %s.ic.q' % (prefix,prefix), shell=True)
    subprocess.check_call('./forward --output '+forwardFile, shell=True)
    subprocess.check_call('cp %s.q %s.0.q' % (filename,filename), shell=True)
    subprocess.check_call('./linearized --output '+forwardFile, shell=True)
    subprocess.check_call('./spatial_inner_product %s.linearized.q %s.0.q --output %s'                       \
                            %(filename,filename,forwardFile), shell=True)
    fQoI0 = open(forwardFile,'r')
    QoI0 = float(fQoI0.read())
    fQoI0.close()
    print ("QoI0: %.15E" % QoI0)

    subprocess.check_call('./spatial_inner_product %s.linearized.q %s.linearized.q --output %s'                       \
                            %(filename,filename,deltaFile), shell=True)
    fGrad0 = open(deltaFile,'r')
    Grad0 = float(fGrad0.read())
    fGrad0.close()
    print ("Grad0: %.15E" % Grad0)
else:
    fQoI0 = open(forwardFile,'r')
    QoI0 = float(fQoI0.read())
    fQoI0.close()
    print ("QoI0: %.15E" % QoI0)

    fGrad0 = open(deltaFile,'r')
    Grad0 = float(fGrad0.read())
    fGrad0.close()
    print ("Grad0: %.15E" % Grad0)

    Nk = 10
    Ak = 10.0**(-3.0-2.5-0.25*np.array(range(Nk)))
    QoIk = np.zeros((Nk,),dtype=np.double)
    Gradk = np.zeros((Nk,),dtype=np.double)
    ek = np.zeros((Nk,),dtype=np.double)

    fId = open(prefix+'.gradient_accuracy.txt','w')
    for k in range(Nk):
        actuation_amount = Ak[k]
        for i in range(1):
            command = './zaxpy %s %.16E %s %s' % (control_forcing_file[i], actuation_amount, grad_file[i], baseline_control_forcing_file[i])
            subprocess.check_call(command, shell=True)
        command = './qfile_zaxpy %s.ic.q %.16E %s.ic.linearized.q %s.ic.0.q' % (prefix, actuation_amount, prefix, prefix)
        subprocess.check_call(command, shell=True)
        subprocess.check_call('./forward --output '+forwardFile, shell=True)
        subprocess.check_call('./spatial_inner_product %s.linearized.q %s.q --output %s'                       \
                                %(filename,filename,forwardFile1), shell=True)
        fQoI1 = open(forwardFile1,'r')
        QoIk[k] = float(fQoI1.read())
        fQoI1.close()

        Gradk[k] = (QoIk[k]-QoI0)/Ak[k]
        ek[k] = abs( (Gradk[k]-Grad0)/Grad0 )
        fId.write("%.16E\t%.16E\t%.16E\t%.16E\n"%(actuation_amount,QoIk[k],Gradk[k],ek[k]))
        print ("%.16E %.16E %.16E %.16E"%(actuation_amount,QoIk[k],Gradk[k],ek[k]))

    fId.close()
