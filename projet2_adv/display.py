import numpy as np
import matplotlib.pyplot as plt

with open('file_data.txt') as f:
    lines = f.readlines()

nt = int(lines[0][5:11])
nx = int(lines[1][5:11])
sm = int(lines[2][12:16])

T=[]
for k in range(3,3+sm):
    T.append(float(lines[k][12:22]))

with open('file_sol.txt') as f:
    lines = f.readlines()

X   = np.zeros(nx)
dec = 0

U       = np.zeros((2,sm,nx))
U_ex    = np.zeros((2,sm,nx))
err = np.zeros((sm,nx))
err_L2 = np.zeros(sm)

for k in range(0,sm):
    print(k)
    for i in range(0,(nx)):
        
        if(lines[k*(nx+1) +i][1] != '*') :
            X[i] = lines[k*(nx+1) + i][1:10]
            
        if(lines[k*(nx+1) +i][11] != '*') :
            U[0][k][i]   = lines[k*(nx+1) +i][10:22]
        else : U[0][k][i] = 0
        
        if(lines[k*(nx+1) +i][23] != '*') :
            U[1][k][i]     = lines[k*(nx+1) +i][23:34] 
        else : U[1][k][i] = 0
        
        
        if(lines[k*(nx+1) +i][36] != '*') :
            U_ex[0][k][i]     = lines[k*(nx+1) +i][36:47] 
        else : U_ex[0][k][i] = 0
        
        
        if(lines[k*(nx+1) +i][49] != '*') :
            U_ex[1][k][i]     = lines[k*(nx+1) +i][48:59] 
        else : U_ex[1][k][i] = 0
        
        
    plt.title("ondes à temps :"+str(T[k]))    
    plt.plot(X,U[0][k],'g')
    #plt.plot(X,U_ex[0][k],'b')
    plt.plot(X,U[1][k],'k')
    #plt.plot(X,U[0][k] - U_ex[0][k],'r')
    #plt.plot(X,U[1][k]/U[0][k],'b')
    plt.ylim(-2,3); plt.xlim(X[1],X[-1])
    plt.show()
    
    err_L2[k] = np.sqrt(sum((U[0][k]-U_ex[0][k])**2) + sum((U[1][k]-U_ex[1][k])**2))

print(max(err_L2))
    