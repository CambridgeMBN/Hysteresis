from numpy import arange, array, pi, zeros, empty,sqrt, copy, sin, cos, nonzero, ravel, append,tanh, add, savetxt
#from matplotlib import pyplot as plt
#from pylab import *
import math
#import scipy.optimize as opt
#from random import randrange,shuffle
#from matplotlib.ticker import MultipleLocator, MaxNLocator

theta=pi/4
phis = [pi/2,0.] #first entry is  the initial phi for Ni, second entry is the initial phi for Gd

savepath = 'Arrow_Data/'

Tlist=[6.4]#, 30., 37.5, 40., 50.]#5.,10.,15.,20.,25.,45.,55.,60.,65.,70.,75.]
for T in Tlist:
    #Tc=[631,33.05]
    Tc=[631.,28.] #(28 is pretty good)
#    kb=1.3806488E-23
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
    Thickness_Ni=7. #7.#   
    Thickness_Gd=7.9 #7.9    

    nGd=int(round(Thickness_Gd/0.3,0))
    nNihalf=int(round(Thickness_Ni/0.4,0))
    nNi=nNihalf*2
    
    Ni = [0]
    Gd = [1]
    layers = nNihalf*Ni + (nGd)*Gd + nNihalf*Ni

    def create_stack(n,layers,x):
        x_vector = []
        for i in range(0,n,1): 
            x_vector = append(x_vector,x[layers[i]])
        x_vector = append(append(0,x_vector),0)
        return  x_vector #a LIST of given order is returned here.
        
    n = len(layers)
    n_stack = n+2
    Ku = create_stack(n,layers,K)
    Ms = create_stack(n,layers,M)
    Mt = create_stack(n,layers,MT)
    phi0 = create_stack(n,layers,phis)
    
    # print 'layers ', layers, ' ku ', Ku, ' ms ', Ms, ' ', Mt, ' ', phi0
    
    
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

    # print 'J ', J
    B = arange(-0.5,0.5,.01)
    H = B/mu0

    Mtot_Ni = empty(len(H))  #generate an array of random small number of the size len(H)
    Mtot_Gd = empty(len(H))
    Mtot_Ni_bw = empty(len(H))
    Mtot_Gd_bw = empty(len(H))
    phi_check = empty(n_stack) #generate an array of random small number of the size len(phi0), which is the same as number of layers
                               #phi0 is the array that stands for the initial phi of the layers
    phi_store = empty(n_stack)

    def magnetic_energy(H,Ku,theta,Mt,Ms,Msup,phiup,Jup,Msdn,phidn,Jdn, i):
        # print locals.get('arg')
        # if i > 16:
        #     print 'i ', i, ' H ', H, ' Ku ', Ku, ' theta ', theta, ' Mt ', Mt, ' Ms ', Ms, ' MSup', Msup, ' phiup ', phiup, \
        #         ' Jup ', Jup, ' MsDn ', Msdn, ' PhiDn ', phidn, ' JDn ', Jdn
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
        
        print 'i ', i, ' phase: ', phase
        # if i > 16:
        #     print 'i ', i, 'k ', Ku, ' as ', Asin, ' acos ', Acos, ' part ', Jdn, ' abs1 ', absphase1, ' abs2 ', absphase2, ' ct1 ', Ctot2, ' ct2 ', Ctot2, ' p1 ', phase1, ' p2 ',phase2, ' p ', phase
        # print phase
        return phase
        
    def hysterisis(phi0,RANGE):
        RANGE = [0]
        for k in RANGE:
            iterations = 0
            delta = array(range(n_stack))

            # while delta[2:-2].max()>0.1: # 0.0001
            while (iterations <= 0):
                print " \n"
                phi_check=copy(phi0)
                
                for i in range(1, half_point, 1): # range(1, half_point,1): 
                    phi0[i]= magnetic_energy(H[k],Ku[i],theta,Mt[i],Ms[i],\
                        Ms[i-1],phi0[i-1],J[i-1],Ms[i+1],phi0[i+1],J[i], i)
                    a = phi0
                    # print "p01 ", phi0
                    phi0[len(Ms)-1-i]=phi0[i]
                    # b = phi0
                    # print "p02 ", phi0 , " diff ", a-b
                     
                delta = abs(phi_check - phi0) 
                
                iterations+=1

            b = int(B[k]*100)
            t = int(T)
            savetxt(savepath + 'Arrows_%i_%i' %(b, t), phi0)
            
            
            Mtot_Ni[k] = sum(Mt[1:nNi/2+1]*cos(phi0[1:nNi/2+1])) + \
                         sum(Mt[nNi/2+nGd+1:-0]*cos(phi0[nNi/2+nGd+1:-0]))
            Mtot_Gd[k] = sum(Mt[nNi/2+1:nNi/2+nGd+1]*cos(phi0[nNi/2+1:nNi/2+nGd+1]))

        return Mtot_Ni, Mtot_Gd, phi0

    hysterisis_fw1=hysterisis(phi0,range(len(arange(-0.5,0,0.01)),len(H),1))
    phi0=hysterisis_fw1[2]
    hysterisis_bw=hysterisis(phi0,range(len(H)-1,-1,-1))

    Mtot_Ni_bw=hysterisis_bw[0]
    Mtot_Gd_bw=hysterisis_bw[1]
    Mtot_Ni=-Mtot_Ni_bw
    Mtot_Gd=-Mtot_Gd_bw

    print ('T = ' + str(T))
