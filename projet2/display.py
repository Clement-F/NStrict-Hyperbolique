import numpy as np
import matplotlib.pyplot as plt

with open('file_data.txt') as f:
    lines = f.readlines()

nt = int(lines[0][5:11])
nx = int(lines[1][5:11])
sm = int(lines[2][12:16])

T=[]
for k in range(3,3+sm):
    T.append(float(lines[k][14:22]))

with open('file_sol.txt') as f:
    lines = f.readlines()

X   = np.zeros(nx)
dec = 0

U1 = np.zeros((sm,nx))
U2 = np.zeros((sm,nx))
err = np.zeros((sm,nx))
err_L2 = np.zeros(sm)

for k in range(0,sm):
    print(k)
    for i in range(0,(nx)):
        X[i] = lines[k*(nx+1) + i][1:11]
        U1[k][i]   = lines[k*(nx+1) +i][11:23]
        U2[k][i]     = lines[k*(nx+1) +i][23:35] 
        
        
    plt.plot(X,U1[k],'g')
    plt.plot(X,U2[k]/U1[k],'b')
    plt.ylim(-2,3)
    plt.show()

print(max(err_L2))
    