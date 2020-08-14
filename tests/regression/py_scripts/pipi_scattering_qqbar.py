import numpy as np
import h5py as h5
from gamma import g_5 as G

'''To make this code consistent with the notes, we note:
    
G = gamma5 matrix
S = quark propagator
M = G S G
D = delta matrix
Q = qqbar object'''

'''The following spin indices corresponds to the notes, since we don't want to type Greek letters in the einsum.
alpha  = i
beta   = j
sigma  = k
tau    = l
alpha' = m
beta'  = n
sigma' = o
tau'   = p '''

'''Denote primes
x' = q
y' = r
z' = s
t' = u

 '''

# Remember this program uses unsmeared propagator!
# For the current purpose, let us treat the four propagators as the same.

t_0 = 3
Nt = 8

D = np.zeros([3, 3], dtype = np.complex128)
for i in range(3):
    D[i,i] = 1.

path = '/Users/haobo/lattice_qcd/test_run/lalibe/test_propagator2.h5'
f = h5.File(path, 'r')
S = f['sh_sig2p0_n5/PS_prop'][()]
f.close()

M = np.einsum('li,tzyxikab,kj->tzyxljab', G, S, G)

correlator  = np.zeros([8,4,4,4,4,4,4], dtype = np.complex128)
correlator1 = np.zeros([8,4,4,4,4,4,4], dtype = np.complex128)
correlator2 = np.zeros([8,4,4,4,4,4,4], dtype = np.complex128)
correlator3 = np.zeros([8,4,4,4,4,4,4], dtype = np.complex128)
correlator4 = np.zeros([8,4,4,4,4,4,4], dtype = np.complex128)

quarkloop  = np.zeros([8,4,4,4,4,4,3,3], dtype = np.complex128)
for t in range(8):
    for z in range(4):
        for y in range(4):
            for x in range(4):
                quarkloop[t,z,y,x,:,:,:,:] =  np.einsum('mjba,mn,nibc->jiac', np.conj(M[t,z,y,x,:,:,:,:]), G, S[t,z,y,x,:,:,:,:], optimize='greedy')

p2max = 1
ptotmax = 2
N=10

momlist = []
for px in range(-N, N):
    for py in range(-N, N):
        for pz in range(-N, N):
            if (px ** 2 + py ** 2 + pz ** 2) <= p2max:
                momlist.append([px,py,pz])

# Construct qqbar objects
Q  = np.zeros([8,len(momlist),4,4,3,3], dtype = np.complex128)

print('qqbar starts.')
pnum = 0
for p in momlist:
    for t in range(8):
        t_relative = t - t_0
        if (t_relative < 0):
            t_relative += Nt
        for z in range(4):
            for y in range(4):
                for x in range(4):
                    Q[t_relative,pnum,:,:,:,:] += np.exp(1j*(np.pi/2)*(p[0]*x+p[1]*y+p[2]*z)) * quarkloop[t,z,y,x,:,:,:,:]
    pnum += 1

print('qqbar completed.')

correlator_FTed = {}
p1num = 0
for p1 in momlist:
    p2num = 0
    for p2 in momlist:
        if ((p1[0]+p2[0]) ** 2 + (p1[2]+p2[2]) ** 2 + (p1[2]+p2[2]) ** 2) <= ptotmax:            
            # Let me omit the origin fix at this moment
            tmpFTlist = np.zeros(8, dtype = np.complex128)
            for t in range(8):
                tmpFTlist[t] += np.einsum('jiaa,ij', Q[t,p1num,:,:,:,:], G, optimize='greedy') * np.einsum('jiaa,ij', Q[t,p2num,:,:,:,:], G, optimize='greedy')
                tmpFTlist[t] -= np.einsum('jkac,kl,lica,ij', Q[t,p1num,:,:,:,:], G, Q[t,p2num,:,:,:,:], G, optimize='greedy')
                tmpFTlist[t] -= np.einsum('jkac,kl,lica,ij', Q[t,p2num,:,:,:,:], G, Q[t,p1num,:,:,:,:], G, optimize='greedy')
                tmpFTlist[t] += np.einsum('jiaa,ij', Q[t,p2num,:,:,:,:], G, optimize='greedy') * np.einsum('jiaa,ij', Q[t,p1num,:,:,:,:], G, optimize='greedy')
            correlator_FTed[(p1[0],p1[1],p1[2],p2[0],p2[1],p2[2])] = tmpFTlist
        p2num += 1
    p1num += 1

# Python version with qqbar
PYQ = correlator_FTed[(0,0,1,0,0,0)]
