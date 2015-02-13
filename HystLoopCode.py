from numpy import arange, array, pi, zeros, empty,sqrt, copy, sin, cos, nonzero, ravel, append,tanh, add
from matplotlib import pyplot as plt
# from pylab import *
import math
import scipy.optimize as opt
from random import randrange,shuffle
from matplotlib.ticker import MultipleLocator, MaxNLocator

print('init')
theta=pi/4
phis = [pi/2,0.] #first entry is  the initial phi for Ni, second entry is the initial phi for Gd

folder='4'

Tlist=[5.] # Tlist=[5.,70.]
for T in Tlist:
    Tc=[631.,28.] #(28 is pretty good)

    mu0 = 4*pi*1E-7
    N_Ni = 9.14E28
    N_Gd = 3.02E28
    muB = 9.274E-24
    M = array([1.3, 7.])/2 #(1.3, 7.)
    MGd_min=0.3414
    MT = array([1.3, (7.-MGd_min)/2*(1+tanh(0.1143*(Tc[1]-T)))+MGd_min])/2
    K = array([4000./N_Ni*2, 17500./N_Gd*(MT[1]/7*2)**3/4])*3/1.5 #4000, 17500
    Js = array([1.2E-22, 2.625E-24])*4.45/2
    Jex = -5.E-22
    Thickness_Ni=2.5 #7.#       at Ni4. and Gd3. they both flip?
    Thickness_Gd=2.5 #7.9

    # print('phis:', phis)
    print('K:',K, " T: ",T)
    print('MT:',MT)
    print('E_ex:',Js[0]*M[0]**2,Js[1]*M[1]**2)

    # Define the thickness of the magnetic stack

    nGd=int(round(Thickness_Gd/0.3,0))
    nNihalf=int(round(Thickness_Ni/0.4,0))
    nNi=nNihalf*2


    
    Ni = [0]
    Gd = [1]
    layers = nNihalf*Ni + (nGd)*Gd + nNihalf*Ni

    print('gd: ', nGd, ' nh ', nNihalf, ' m2 ', nNihalf*2, ' layers ', layers)
    def create_stack(n,layers,x):
        x_vector = []
        for i in range(0,n,1):
            x_vector = append(x_vector,x[layers[i]])
        x_vector = append(append(0,x_vector),0)
        # print('xvector: ', x_vector)
        return  x_vector #a LIST of given order is returned here.

    n = len(layers)
    n_stack = n+2
    Ku = create_stack(n,layers,K)
    Ms = create_stack(n,layers,M)
    Mt = create_stack(n,layers,MT)
    phi0 = create_stack(n,layers,phis)

    print("\n")    
    print('ku ', Ku, ' ms ', Ms, ' mt ', Mt, ' phi10 ', phi0)
    print("\n")
    if n_stack%2==1:
        half_point=int((n_stack+1)/2)
    else:
        half_point=int(n_stack/2)

    k=2*pi/10/2 #the period of the variation of JGd is 10*2
    f=1
    J=zeros(n_stack)
    for i in range(half_point):
        if i==nNi/2:
            J[i]=Jex
        elif i<=nNi/2-1 and i>=nNi/2-5:
            J[i] = Js[0]
        elif i>=nNi/2+1 and i<=nNi/2+5:
            J[i] = Js[1]
        elif i>nNi/2+5:
            J[i] = Js[1]
        else:
            J[i] = Js[0]
        J[len(Ms)-2-i]=J[i]

    B = arange(-0.5,0.5,.01)
    H = B/mu0

    Mtot_Ni = empty(len(H))  #generate an array of random small number of the size len(H)
    Mtot_Gd = empty(len(H))
    Mtot_Ni_bw = empty(len(H))
    Mtot_Gd_bw = empty(len(H))
    phi_check = empty(n_stack) #generate an array of random small number of the size len(phi0), which is the same as number of layers
                               #phi0 is the array that stands for the initial phi of the layers
    phi_store = empty(n_stack)


    """In the iteration process, for a particular field Hk we find the energy and angle phi layer by layer using i.
        What we do is we start from layer 1, keeping phi of all other layers fixed(at the initial value or the value from last iteration), then we
        rotate the spin in layer 1 around and find the phi_min with minimum E for layer 1. Put phi_min as the new angle for layer 1, then we go on
        for layer 2, etc. After we found the phi_min for whole stack, we will see how much the new phi is differ from last phi, and if the difference
        is small, we know that we are closer to the stationary configuration of phi. We can do this because we know that at the equilibrium configuration
        of the stack, if we vary the angle of ANY layers, it will lead to a increase in E. Hence, if the stack is in the stationary configuration,
        and we put it through the iteration, it should remain unchange. The closer the stack is to the stationary configuration, the smaller the change
        it will undergo when we put it through the iteration."""


    def magnetic_energy(H,Ku,theta,Mt,Ms,Msup,phiup,Jup,Msdn,phidn,Jdn):

        Asin=Jup*Msup*sin(phiup)+Jdn*Msdn*sin(phidn)
        Acos=Jup*Msup*cos(phiup)+Jdn*Msdn*cos(phidn)
        Atot=sqrt(Asin**2+Acos**2)
        Btot=mu0*2*muB*H
        Ctot1=sqrt((Ms*Acos+Mt*Btot+Ku*cos(theta))**2+(Ms*Asin+Ku*sin(theta))**2)
        absphase1=math.asin((Ms*Asin+Ku*sin(theta))/Ctot1)

        #let phi=phase to have minimum energy
        if Ms*Acos+Mt*Btot+Ku*cos(theta)>=0:
            phase1=absphase1
        else:
            phase1=pi-absphase1

        #for pi/2+theta<phi<(3/2)pi+theta
        Ctot2=sqrt((Ms*Acos+Mt*Btot-Ku*cos(theta))**2+(Ms*Asin-Ku*sin(theta))**2)
        absphase2=math.asin((Ms*Asin-Ku*sin(theta))/Ctot2)

        #let phi=phase to have minimum energy
        if Ms*Acos+Mt*Btot-Ku*cos(theta)>=0:
            phase2=absphase2
        else:
            phase2=pi-absphase2

        if Ctot2>Ctot1:
            phase=phase2
        else:
            phase=phase1

        return phase

    def hysterisis(phi0,RANGE):

        for k in RANGE:
            iterations = 0
            delta = array(range(n_stack))

            while delta[2:-2].max()>0.1:
                phi_check=copy(phi0)

                for i in range(1, half_point,1):
                    # (H,Ku,theta,Mt,Ms, \
                    # Msup,phiup,Jup,Msdn,phidn,Jdn)
                    phi0[i]= magnetic_energy(H[k],Ku[i],theta,Mt[i],Ms[i],\
                        Ms[i-1],phi0[i-1],J[i-1],Ms[i+1],phi0[i+1],J[i])
                    phi0[len(Ms)-1-i]=phi0[i]

                delta = abs(phi_check - phi0)

                iterations+=1

            Mtot_Ni[k] = sum(Mt[1:nNi/2+1]*cos(phi0[1:nNi/2+1])) + \
                         sum(Mt[nNi/2+nGd+1:-0]*cos(phi0[nNi/2+nGd+1:-0]))
            Mtot_Gd[k] = sum(Mt[nNi/2+1:nNi/2+nGd+1]*cos(phi0[nNi/2+1:nNi/2+nGd+1]))


            #print ('Ground state found after ' + str(iterations) + ' iterations')

        return Mtot_Ni, Mtot_Gd, phi0

    hysterisis_fw1=hysterisis(phi0,range(len(arange(-0.5,0,0.1)),len(H),1))
    phi0=hysterisis_fw1[2]
    hysterisis_bw=hysterisis(phi0,range(len(H)-1,-1,-1))

    Mtot_Ni_bw=hysterisis_bw[0]
    Mtot_Gd_bw=hysterisis_bw[1]
    Mtot_Ni=-Mtot_Ni_bw
    Mtot_Gd=-Mtot_Gd_bw

    print ('T = ' + str(T))

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    plt.rcParams['xtick.major.pad'] = 8
    plt.rcParams['ytick.major.pad'] = 8

    fig = plt.figure()
    axes = fig.add_subplot(111)
    plt.plot(-mu0*H*1E3,array(Mtot_Ni),'g')
    plt.plot(mu0*H*1E3,array(Mtot_Ni_bw),'r')

    ax1 = gca()
    for label in ax1.get_xticklabels() + ax1.get_yticklabels():
        label.set_fontsize(30)
    ax1.set_xlabel('Magnetic Field (mT)', fontsize = 30, labelpad = 15)
    ax1.set_ylabel('Magnetisation ($\mu_B$/atom)', fontsize=30, labelpad = 3)
    plt.locator_params(axis = 'x', nbins = 7)
    plt.subplots_adjust(bottom=0.18, left=0.18) #ensures xlabel doesn't get cut off
    axes.autoscale(True,'x',True)
    T=int(T)
    #plt.savefig('/home/tdch4/Dropbox/Uni/PhD/Physics/Presentations/THiggs_Pres_02Oct14/Figures/Ni_thin_%i.pdf'%T)

    #plt.show()
    #plt.clf()

    saveNi = zip(*[-mu0*H*1E3,Mtot_Ni_bw])
#     #np.savetxt('/home/tdch4/Dropbox/Uni/PhD/Physics/PubPapers/PRLSept14/ModelData/Ni%iK.txt'%T, saveNi)
#
#     saveNi_up = zip(*[-mu0*H*1E3,Mtot_Ni])
#     #np.savetxt('/home/tdch4/Dropbox/Uni/PhD/Physics/PubPapers/PRLSept14/ModelData/Ni%iK_up.txt'%T, saveNi_up)
#
#
#     fig = plt.figure()
#     axes = fig.add_subplot(111)
#     plt.plot(-mu0*H*1E3,array(Mtot_Gd),'g')
#     plt.plot(mu0*H*1E3,array(Mtot_Gd_bw),'r')
#     plt.title('Gd')
#     #plt.xlabel('$\mu0$ H (mT)')
#     #plt.ylabel('magnetic moment ($\mu_B$/atom)')
#     #plt.title('Ni')
#     #plt.xlim([min(1E3*mu0*H), max(1E3*mu0*H)])
#     #plt.grid(True)
#     ax1 = gca()
#     for label in ax1.get_xticklabels() + ax1.get_yticklabels():
#         label.set_fontsize(30)
#     ax1.set_xlabel('Magnetic Field (mT)', fontsize = 30, labelpad = 15)
#     ax1.set_ylabel('Magnetisation ($\mu_B$/atom)', fontsize=30, labelpad = 3)
#     plt.locator_params(axis = 'x', nbins = 7)
#     plt.subplots_adjust(bottom=0.18, left=0.18) #ensures xlabel doesn't get cut off
#     axes.autoscale(True,'x',True)
#     T=int(T)
#     #plt.savefig('/home/tdch4/Dropbox/Uni/PhD/Physics/Presentations/THiggs_Pres_02Oct14/Figures/Gd_thin_%i.pdf'%T)
#
#     plt.show()
#
#     saveGd = zip(*[-mu0*H*1E3,Mtot_Gd_bw])
#     #np.savetxt('/home/tdch4/Dropbox/Uni/PhD/Physics/PubPapers/PRLSept14/ModelData/Gd%iK.txt'%T, saveGd)
#
#     saveGd_up = zip(*[-mu0*H*1E3,Mtot_Gd])
#     #np.savetxt('/home/tdch4/Dropbox/Uni/PhD/Physics/PubPapers/PRLSept14/ModelData/Gd%iK_up.txt'%T, saveGd_up)
#
#     '''
#     plt.figure(1)
#     plt.plot(-mu0*H*1E3,array(Mtot_Gd),'b')
#     plt.plot(mu0*H*1E3,array(Mtot_Gd_bw),'r')
#     plt.xlabel('$\mu0$ H (mT)')
#     plt.ylabel('magnetic moment ($\mu_B$/atom)')
#     plt.title('Gd')
#     plt.xlim([min(1E3*mu0*H), max(1E3*mu0*H)])
#     plt.grid(True)
#
#     plt.figure(2)
#     plt.plot(-mu0*H*1E3, add(Mtot_Ni,Mtot_Gd), 'b')
#     plt.plot(mu0*H*1E3, add(Mtot_Ni_bw,Mtot_Gd_bw), 'r')
#     plt.xlabel('$\mu0$ H (mT)')
#     plt.ylabel('magnetic moment ($\mu_B$/atom)')
#     plt.title('Total')
#     plt.xlim([min(1E3*mu0*H), max(1E3*mu0*H)])
#     plt.show()
#     '''
  



