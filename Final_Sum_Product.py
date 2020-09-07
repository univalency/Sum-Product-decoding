# WHAT THIS CODE DOES:
#generates a random binary message (vector of size k) using a given LDPC generator matrix G.
#modulates (changes 0 to -1)
#adds Gaussian noise to every entry 
#the target is to recover the original message using Sum-Product algorithm which is done by 
#thinking of parity check matrix H as a Tanner graph (variable nodes, check nodes, etc.)




import numpy as np
import math
import random

def check_codeword(Y, H):
    N= np.sum(Y.dot(H)) 
    return N

def de_modul(Y):
    for i in range(len(Y)) :
        if Y[i] > 0: Y[i] = 0
        else: Y[i] = 1
    return np.asarray(Y)
    
def modul(Y):
    for i in range(len(Y)):
        if Y[i] == 0:
            Y[i]=1
        else:
            Y[i]=-1
    return Y

def LLR(x):  # function for calculating log likelihood ratios
    sigma = 1
    return 2 * x / (sigma**2)

def SumProduct_iteration(L_v, L_ch, d_check, d_var,L_c):
    k = 120
    N = 520
    M = 400
    d_check_msgs = {}#messages that a check_node receives from each var_node
    for i in range (0,M):
        a = []
        for j in range(len(d_check[i])):
            a.append(L_v[d_check[i][j]]-L_c[i][j])
        d_check_msgs[i] = a
   # print("d_check_msgs:", d_check_msgs[1][:3])
    
    
    L_c.clear()
    for k in range (M):# for each check_node we calculate a message that it must send
        a=[]
        for i in range ( len(d_check[k])):
            product = 1
            for j in range(len(d_check[k])):# j is the number of an adjacent vertex. 
                if(j!=i):
                    product *= np.tanh(d_check_msgs[k][j]/2)# message which k's сheck_node received from j's var_node
            if(product>=0.99999999999999994) :
                product=0.99999999999999994
            if(product<=-0.99999999999999994) :
                product=-0.99999999999999994
            product = 2*np.arctanh(product)
            a.append(product)
        L_c[k]=a
        

    #print(L_v) 
    for i in range (0,N):# for each var_node calculate its value
        sum = L_ch[i]# add value from channel to i'th element 
        for j in d_var[i]:# loop across all adjacent check nodes 
            for k in range(len( d_check [j])):
                if d_check [j][k] == i:
                    x = L_c [j][k]# receive message from check node 
           
            sum+= x
            L_v[i] = sum
            
            
            
G = np.loadtxt("G_new.txt")
H = np.loadtxt('H_new.txt').astype(np.int)
np.sum(np.mod((H).dot(G.transpose()),2))

k = G.shape[0]
N = H.shape[1]
M = H.shape[0]

message = (np.random.randint(0, 2, size=(k, ))).astype(np.float64)
code = modul( np.mod(message.dot(G),2) ).astype(np.float64)

noise = np.random.normal(0,1,N)
noisy_code = (code + noise)


L_ch = [LLR(noisy_code[i]) for i in range(N)]
L_v = L_ch[:] 
d_check = {}
for i in range (0,M):
    a = []
    for j in range (0,N):
        if H[i,j] == 1:
            a.append(j)
    d_check[i] = a
d_var = {}
for i in range (0,N):
    b = []
    for j in range (0,M):
        if H[j,i] == 1:
            b.append(j)
    d_var[i] = b
    i=0
Check = np.array(L_v[:]).transpose()
L_c = {}
for k in range (M):# for each check node calculate messages which in must send.
    a=[]
    for i in range ( len(d_check[k])): #go along adjacent var_nodes, i -  var_node
        for j in range(len(d_check[k])):#j- adjacent vertex
            if(j!=i):
                a.append(0)
    L_c[k]=a
x=np.mod((H).dot(de_modul(Check)),2)
J = np.sum(x) 
print("initially = ",J, "\n")
i = 0
while( J>0 ):
    i += 1
    SumProduct_iteration(L_v, L_ch, d_check, d_var, L_c)
    Check = L_v[:]
    J = np.sum(np.mod((H).dot(de_modul(Check)),2)) 
    print("after",i, "iteration :","J =",J,"\n")

T = (Check)
print("initial:",de_modul(code)[:24])   
print("with noise",de_modul(noisy_code)[:24])   
print("output",T[:24],"\n")

print("difference of noisy from original", np.sum(abs(noisy_code-(code))))



