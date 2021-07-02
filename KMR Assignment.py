import numpy as np # Since we will define the equation variables and coefficients in array form 

from numdifftools.nd_algopy import Jacobian # For use of grad operator and finding jacobian of a matrix with ease 

import math # for using mathematical funcs 

from matplotlib.pyplot import * #for plotting the links for animation 

ct = 0 

# Given data [taken as per Student's wish as Instructed] 

r1 = 60 

o2 = math.radians(0) #varries 

r2 = 40 

r6 = 60 

# Defining Func Eqns in Python  

#r3,r4,r5,o3 = 87.17,130.75,51.96,math.radians(66.58) for theta = 30 degrees 

eqn = [] 

eq1 = lambda x : 40*np.sin(o2)-x[0]*np.sin(x[3])+60 

eqn.append(eq1) 

eq2 = lambda x : 40*np.cos(o2)-x[0]*np.cos(x[3]) 

eqn.append(eq2) 

eq3 = lambda x : 60-x[1]*np.sin(x[3])+60 

eqn.append(eq3) 

eq4 = lambda x : x[2]-x[1]*np.cos(x[3]) 

eqn.append(eq4) 

#Jacobians Defined 

jacob1 = Jacobian(eq1) 

jacob2 = Jacobian(eq2) 

jacob3 = Jacobian(eq3) 

jacob4 = Jacobian(eq4) 

#plot 

def ploter(ct): 

    #Naming file 

    if ct<10 : 

        filename = 'QuickReturn00'+str(ct)+'.jpg' 

    elif ct<100 : 

        filename = 'QuickReturn0'+str(ct)+'.jpg' 

    else : 

        filename = 'QuickReturn'+str(ct)+'.jpg'   

     

    #Plotting 

    figure() 

    xlim(-200, 200) 

    ylim(0, 125) 

    plot((s0[0],s3[0]),(s0[1],s3[1]),linewidth='4') 

    plot((s1[0],s2[0]),(s1[1],s2[1]),linewidth='4') 

    plot((-200,200),(120,120),linestyle='dashed') 

    plot((0,0),(0,120),linestyle='dashed') 

    savefig(filename) 

    ct=ct+1  

    return ct 

# Newton Raphson Method for Multivariable systems of Non-Linear equation to be implemented 

  

def root_finder(): 

    i = 0 

    er = 100 

    tol = 0.000001 

    maxiter = 100 

    M = 4 

    N = 4 

    x0 = np.array([1,1,1,1],dtype=float).reshape(N,1) 

    while np.any(abs(er)>tol) and i< maxiter: 

        func_eval = np.array([eq1(x0),eq2(x0),eq3(x0),eq4(x0)]).reshape(M,1) 

        flat_x0 = x0.flatten() 

        jac = np.array([jacob1(flat_x0),jacob2(flat_x0),jacob3(flat_x0),jacob4(flat_x0)]) 

        jac = jac.reshape(N,M) 

        x_new = x0 - np.linalg.inv(jac)@func_eval 

        er = x_new - x0 

        x0 = x_new 

        i+=1 

    return x_new 

#Itering through theta to find all solutions 

for i in range(360): 

    o2 = np.radians(i) 

    x_new = root_finder() 

    #co-ordinates 

    o3 = x_new[3][0] 

    r4 = x_new[1][0] 

    s0 = [0,0] 

    s1 = [0,r1] 

    s2 = [s1[0]+r2*math.cos(o2),s1[1]+r2*math.sin(o2)] 

    s3 = [r4*math.cos(o3),r4*math.sin(o3)] 

    if i ==30 : 

        print(x_new) 

    clf() 

    ct = ploter(ct) 