import numpy as np
import matplotlib.pyplot as plt

with open('convergence_err.txt') as f:
    lines = f.readlines()

size = len(lines)
nx =[]; err_L1 =[]; err_L2 =[]; err_Li=[]

for k in range(0,int(size/6)):
    nx.append(      int(lines[6*k+1][9:15])  )
    err_L1.append(  float(lines[6*k+2][12:25])) 
    err_L2.append(  float(lines[6*k+3][12:25])) 
    err_Li.append(  float(lines[6*k+4][12:25])) 
    
plt.loglog(nx,err_L1,'b' , label="err_L1")
plt.loglog(nx,err_L2,'r', label="err_L2")
plt.loglog(nx,err_Li,'g', label="err_Linf")
plt.legend()
plt.grid()
plt.show()
