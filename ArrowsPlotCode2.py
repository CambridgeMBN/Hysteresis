
#def clearall():
# """clear all globals"""
# for uniquevar in [var for var in globals().copy() if var[0] != "_" and var != 'clearall']:
# del globals()[uniquevar]
#clearall()

from numpy import arange, array, pi, zeros, empty,sqrt, copy, sin, cos, nonzero, ravel, append,tanh, add, savetxt
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import math, os
import scipy.optimize as opt
from random import randrange,shuffle
from matplotlib.ticker import MultipleLocator
#import time
#theta = pi/6
theta=pi/4
phis = [1/2*pi,0] #first entry is  the initial phi for Ni, second entry is the initial phi for Gd
#phis = [phi[randrange(1000)], phi[randrange(1000)]]

T = 6.4     #T transit around 55,60K
#Tc=[631,33.05]
Tc=[631,28.]
kb=1.3806488E-23
mu0 = 4*pi*1E-7
N_Ni = 9.14E28
N_Gd = 3.02E28
muB = 9.274E-24
#M = array([1.3, 7.])/2
M = array([1.3, 7.])/2
MGd_min=0.3414
MT = array([1.3, (7.-MGd_min)/2*(1+tanh(0.1143*(Tc[1]-T)))+MGd_min])/2 
#MT = array([1.3, 2.0])/2 
#K = array([4200./N_Ni*2, 16500./N_Gd*(MT[1]/7*2)**3/4])*3/1.5
K = array([4000./N_Ni*2, 17500./N_Gd*(MT[1]/7*2)**3/4])*3/1.5
Js = array([2.4E-22, 3.5E-24])*4.45/2
Jex = -5E-22
Thickness_Ni=7.        
Thickness_Gd=7.9

#drt='/home/tdch4/Dropbox/Uni/PhD/Physics/Cai/ArrowGifs/%inmNi%inmGd%iK/'%(Thickness_Ni,Thickness_Gd,T)
drt='/home/tdch4/Dropbox/Uni/PhD/Physics/PubPapers/APLSept14/Figures/Arrows/%iNi%inmGd/'%(Thickness_Ni, Thickness_Gd)

#print('phis:', phis)
#print('K:',K)
#print('MT:',MT)
#print('E_ex:',Js[0]*M[0]**2,Js[1]*M[1]**2)

# Define the thickness of the magnetic stack

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
    return  x_vector #an LIST of given order is returned here.
    

n = len(layers)
n_stack=n+2
Ku = create_stack(n,layers,K)
Ms = create_stack(n,layers,M)
Mt = create_stack(n,layers,MT)
phi0 = create_stack(n,layers,phis)
color=['b','r']
color_stack = create_stack(n,layers,color)

#print(n_stack, len(Ms))

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
        J[i]=(Js[0]-f*Js[0])*(cos(k*(i-nNi/2+1))**2)+f*Js[0]
    elif i>=nNi/2+1 and i<=nNi/2+5:
        J[i]=(Js[1]-f*Js[1])*(cos(k*(i-nNi/2-1))**2)+f*Js[1]
    elif i>nNi/2+5:
        J[i]=f*(Js[1])
    else:
        J[i]=f*(Js[0])
    J[len(Ms)-2-i]=J[i]

#print (J)
    
B = arange(-0.5,0.5,.01)
H = B/mu0

Mtot_Ni = empty(len(H))  #generate an array of random small number of the size len(H)
Mtot_Gd = empty(len(H))
Mtot_Ni_bw = empty(len(H))
Mtot_Gd_bw = empty(len(H))

phi_plot = [0]*len(H) #generate an array of 0 of the size len(H)

"""deleted a line about " e = [0]*len(H)" here as it has no use."""
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
    #Zeeman=-Btot*S*cos(phi)
    #ex=-S*(Asin*sin(phi)+Acos*cos(phi))
    #An=Ku*sin(phi-theta)=Ku*sin(phi)cos(theta)-Ku*cos(phi)*sin(theta)
    #Etot=-{[S*Acos+S*Btot+K*sin(theta)]*cos(phi)+[S*Asin-K*cos(theta)]*sin(phi)}
    
    #for phi<pi/2+theta or phi>(3/2)pi+theta:
    Ctot1=sqrt((Ms*Acos+Mt*Btot+Ku*cos(theta))**2+(Ms*Asin+Ku*sin(theta))**2)
    absphase1=math.asin((Ms*Asin+Ku*sin(theta))/Ctot1)
    #let phi=phase to have minimum energy
    #Etot=-Ctot
    if Ms*Acos+Mt*Btot+Ku*cos(theta)>=0:
        phase1=absphase1
    else:
        phase1=pi-absphase1

    #for pi/2+theta<phi<(3/2)pi+theta
    Ctot2=sqrt((Ms*Acos+Mt*Btot-Ku*cos(theta))**2+(Ms*Asin-Ku*sin(theta))**2)
    absphase2=math.asin((Ms*Asin-Ku*sin(theta))/Ctot2)
    #let phi=phase to have minimum energy
    #Etot=-Ctot
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
        while delta[2:-2].max()>0.001:
            phi_check=copy(phi0)
            
            for i in range(1, half_point,1): 
                phi0[i]= magnetic_energy(H[k],Ku[i],theta,Mt[i],Ms[i],\
                    Ms[i-1],phi0[i-1],J[i-1],Ms[i+1],phi0[i+1],J[i])
                phi0[len(Ms)-1-i]=phi0[i]
                 
            delta = abs(phi_check - phi0)
            
            iterations+=1

        #print phi0
        plt.figure(k)
        for i in range(1,n_stack-1,1):
            plt.arrow(0,1*i,3*cos(phi0[i]), \
                  3*sin(phi0[i]),shape='full',length_includes_head=False,color=color_stack[i],head_width=1.0)

        plt.ylim([-2,n_stack+3])
        #plt.xlim([-2-n_stack/2,n_stack/2+1])
        plt.xlim([-n_stack/4,n_stack/4])
        #plt.text(1-n_stack/2,n_stack,'{:1.2f}T'.format(B[k]),size=20)
        plt.axis('off')
        #plt.savefig('%s/Fig_%0*d.png'%(drt,3,k), bbox_inches='tight')
        plt.savefig('%s/Fig_%d.pdf'%(drt,B[k]*100))
        plt.close()#'figure(k)')

        phi_plot[k] = phi0

        if b1[k] == 49:
            phi_dict[b1[k]] = phi0
        
        Mtot_Ni[k] = sum(Mt[1:nNi/2+1]*cos(phi0[1:nNi/2+1])) + \
                     sum(Mt[nNi/2+nGd+1:-0]*cos(phi0[nNi/2+nGd+1:-0]))
        Mtot_Gd[k] = sum(Mt[nNi/2+1:nNi/2+nGd+1]*cos(phi0[nNi/2+1:nNi/2+nGd+1]))
        
        
        #print ('Ground state found after ' + str(iterations) + ' iterations')

    return Mtot_Ni,Mtot_Gd,phi0

print 'Gd thickness = ' + str(Thickness_Gd) + ' nm'
print 'Ni thickness = ' + str(Thickness_Ni) + ' nm'
print '-------------------------------------------'
print 'Starting T = ' + str(T) + ' K'
print '...'
print 'Finished T = ' + str(T) + ' K'

if not os.path.exists('%s'%drt):
    os.makedirs('%s'%drt)

phi_dict = {}
b1 = []
for i in B:
    b1.append(int(i*100))

#hysterisis_fw1=hysterisis(phi0,range(len(arange(-0.5,0,0.01)),len(H),1))
#phi0=hysterisis_fw1[2]
hysterisis_bw=hysterisis(phi0,range(len(H)-1,-1,-1))

print 'Saving ...'

#path = '/home/tdch4/Dropbox/Uni/PhD/Physics/PubPapers/PRLSept14/ModelData/Arrows/'
#t = int(T)

#for i in phi_dict:
#    savetxt(path + 'Arrows%i_%iK'%(i,t),phi_dict[i])


Mtot_Ni_bw=hysterisis_bw[0]
Mtot_Gd_bw=hysterisis_bw[1]
Mtot_Ni=-Mtot_Ni_bw
Mtot_Gd=-Mtot_Gd_bw

#for k in range(len(H)):
#
#    fig = plt.figure(k+3)
#    phi_plotk=phi_plot[k]
#    for i in range(1,n_stack-1,1):
#        plt.arrow(0,1*i,3*cos(phi_plotk[i]), \
#                  3*sin(phi_plotk[i]),shape='full',length_includes_head=False,color=color_stack[i],head_width=1.0)
#    plt.ylim([-2,n_stack+1])
#    #plt.xlim([-2-n_stack/2,n_stack/2+1])
#    plt.xlim([-n_stack/4,n_stack/4])
#    plt.axis('off')
    #fig.tight_layout(pad=0)
    #plt.savefig('%s/Fig_%0*d.png'%(drt,3,k))#, bbox_inches='tight')
    #plt.savefig('%s/Fig_%d.png'%(drt,B[k]*100))
    #print B[k]
#    plt.close()

#plt.show()
#convert -delay 10 -loop 0 -deconstruct Fig_*.png anim.gif

'''
plt.clf()
plt.figure(1111)
plt.plot(-mu0*H*1E3,array(Mtot_Ni),'g')
plt.plot(mu0*H*1E3,array(Mtot_Ni_bw),'r')
plt.xlabel('$\mu0$ H (mT)')
plt.ylabel('magnetic moment ($\mu_B$/atom)')
plt.title('Ni')
plt.xlim([min(1E3*mu0*H), max(1E3*mu0*H)])
plt.grid(True)
#plt.savefig('%s/NiLoop.png'%(drt))

plt.clf()
plt.figure(1112)
plt.plot(-mu0*H*1E3,array(Mtot_Gd),'b')
plt.plot(mu0*H*1E3,array(Mtot_Gd_bw),'r')
plt.xlabel('$\mu0$ H (mT)')
plt.ylabel('magnetic moment ($\mu_B$/atom)')
plt.title('Gd')
plt.xlim([min(1E3*mu0*H), max(1E3*mu0*H)])
plt.grid(True)
#plt.savefig('%s/GdLoop.png'%(drt))

plt.clf()
plt.figure(1113)
plt.plot(-mu0*H*1E3, add(Mtot_Ni,Mtot_Gd), 'b')
plt.plot(mu0*H*1E3, add(Mtot_Ni_bw,Mtot_Gd_bw), 'r')
plt.xlim([min(1E3*mu0*H), max(1E3*mu0*H)])
#plt.savefig('%s/MTotLoop.png'%(drt))

plt.clf()
plt.figure(1114)
plt.text(2,8,'theta:     {:1.3f}'.format(theta),size=20)
plt.text(2,7,'phis: {:1.3f}    {:1.3f}'.format(phis[0],phis[1]),size=20)
plt.text(2,6,'K: {:1.2e}     {:1.2e}'.format(K[0],K[1]),size=20)
plt.text(2,5,'MT: {:1.3f}     {:1.3f}'.format(MT[0],MT[1]),size=20)
plt.text(2,4,'Js: {:1.2e}     {:1.2e}'.format(Js[0],Js[1]),size=20)
plt.text(2,3,'E_ex: {:1.2e}     {:1.2e}'.format(Js[0]*M[0]**2,Js[1]*M[1]**2),size=20)
plt.text(2,2,'f: {:1.2e}'.format(f),size=20)
plt.text(2,1,'Thickness_Ni: {:1.2f}'.format(Thickness_Ni),size=20)
plt.text(2,0,'Thickness_Gd: {:1.2f}'.format(Thickness_Gd),size=20)
plt.text(2,9,'Temp: {:1.1f}'.format(T),size=20)
plt.ylim([0,10])
plt.xlim([0,10])
plt.axis('off')
#plt.savefig('%s/data.png'%drt)

#plt.show()
'''

  



