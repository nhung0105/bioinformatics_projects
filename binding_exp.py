import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

#Oct = A; Nanog = B
def reaction_system(state,t):
    [A,B,AB] =  state #define state
    d_A =  -kp*A*B + km*AB
    d_B = -kp*A*B
    d_AB = kp*A*A - km*AB
    return [d_A, d_B, d_AB]

#Parameters of binding 
kp = 0.002 #rate of binding
km = 0.1 #rate of dissociation
K_D = km/kp #dissociation constant
print('K_D = ', K_D)

#initial conditions
A_T = 2
B_T = 1
init_state = [A_T, B_T, 0]

#formulas from the analytical solution
AB_free = A_T * B_T / K_D
AB_eq =(A_T + B_T + K_D - pow(((A_T + B_T + K_D)*(A_T + B_T + K_D) - 4*A_T*B_T),0.5))/2
# AB equilibrium => [AB] = 0 => quadratic equation

t = np.arange(0,40,0.1)
#creates an array of time points from 0 to 40 with an interval of 0.1
state = odeint(reaction_system, init_state, t)

#Plotting the results:
'''''
fig = plt.figure()
#plt.ylim(0,2)
plt.xlim(0,t[-1])
#plt.plot(t,state[:, 0],'k-',label = 'A')
#plt.plot(t,state[:, 1],'b-',label = 'B')
plt.plot(t,state[:, 2],'r-',label = '[Oct4/Nanog complex]')
plt.plot([0,t[-1]],[AB_free,AB_free],'--y',label = '[Oct4/Nanog] (free ligand)')
plt.plot([0,t[-1]],[AB_eq,AB_eq],'--g',label = '[Oct4/Nanog] in equilibrium')
#plt.plot(t,t*kp*A_T*B_T,'--r',label = 'initial rate')
plt.ylabel('concentrations')
plt.xlabel('time')
plt.legend()
plt.show()
'''
#The higher the dissociation constant the lower the complex concentration, and vice versa

#get final state for a range of initial Nanog (B) values, increasing initial Oct 4 (A)
A_T = np.arange(0.2, 10, 0.2)
B_T = np.arange(0.2, 10, 0.2)
t = np.arange(0, 10, 0.1)

j=0
i=0
abf = []

A_T = np.arange(0.2, 10, 0.2)
B_T = np.arange(0.2, 10, 0.2)
t = np.arange(0, 10, 0.1)

i=0
j=0
abf = []            
for b in B_T:
    abf.append([])
    for a in A_T:
        init_state = [a, b, 0]
        state = odeint(reaction_system, init_state, t)
        abf[i].append(state[-1,2])
    i+=1

 
fig = plt.figure()
#plt.xlim(0,t[-1])
i=0
for b in B_T:
    if (i % 10 == 0):
        lb = 'Nanog_T = ' + str (A_T[i])
        plt.plot(A_T,abf[i],'k-',label = lb, color=plt.cm.YlOrBr(A_T[i]/max(A_T)))
    i+=1
plt.ylabel('[Oct4/Nanog] (au)')
plt.xlabel('[Oct4_T] (au)')
plt.legend()
plt.show()