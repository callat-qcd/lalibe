import numpy as np
import h5py as h5
from gamma import g_5 as G

'''To make this code consistent with the notes, we note:
    
G = gamma5 matrix
S = quark propagator
M = G S G
D = delta matrix'''

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

Nt = 8
origin_list=[[0,0,1,0],[2,0,0,0]]
t_0 = origin_list[0][3]

D = np.zeros([3, 3], dtype = np.complex128)
for i in range(3):
    D[i,i] = 1.

path = '/Users/haobo/lattice_qcd/test_run/lalibe/pipi_propagators.h5'
f = h5.File(path, 'r')
S = f['l_prop_x0_y0_z1_t0'][()]
f.close()

path = '/Users/haobo/lattice_qcd/test_run/lalibe/pipi_propagators.h5'
f = h5.File(path, 'r')
S2 = f['l_prop_x2_y0_z0_t0'][()]
f.close()
'''
path = '/Users/haobo/lattice_qcd/test_run/lalibe/test_propagator1.h5'
f = h5.File(path, 'r')
S = f['sh_sig2p0_n5/PS_prop'][()]
f.close()

path = '/Users/haobo/lattice_qcd/test_run/lalibe/test_propagator2.h5'
f = h5.File(path, 'r')
S2 = f['sh_sig2p0_n5/PS_prop'][()]
f.close()'''

M = np.einsum('li,tzyxikab,kj->tzyxljab', G, S, G)
M2 = np.einsum('li,tzyxikab,kj->tzyxljab', G, S2, G)

correlator  = np.zeros([8,4,4,4,4,4,4], dtype = np.complex128)
correlator1 = np.zeros([8,4,4,4,4,4,4], dtype = np.complex128)
correlator2 = np.zeros([8,4,4,4,4,4,4], dtype = np.complex128)
correlator3 = np.zeros([8,4,4,4,4,4,4], dtype = np.complex128)
correlator4 = np.zeros([8,4,4,4,4,4,4], dtype = np.complex128)

cnt = 0
tot = 8
for t in range(8):
    for z in range(4):
        for y in range(4):
            for x in range(4):
                for s in range(4):
                    for r in range(4):
                        for q in range(4):
                            correlator1[t,z,y,x,s,r,q] =  np.einsum('mjba,mn,oldc,op,pkdc,kl,niba,ij', np.conj(M[t,z,y,x,:,:,:,:]), G, np.conj(M2[t,s,r,q,:,:,:,:]), G, S2[t,s,r,q,:,:,:,:], G, S[t,z,y,x,:,:,:,:], G, optimize='greedy')
                            correlator2[t,z,y,x,s,r,q] = -np.einsum('mjba,mn,nkbc,op,oldc,kl,pida,ij', np.conj(M[t,z,y,x,:,:,:,:]), G, S2[t,z,y,x,:,:,:,:], G, np.conj(M2[t,s,r,q,:,:,:,:]), G, S[t,s,r,q,:,:,:,:], G, optimize='greedy')
                            correlator3[t,z,y,x,s,r,q] = -np.einsum('mlbc,mn,niba,op,pkdc,kl,ojda,ij', np.conj(M2[t,z,y,x,:,:,:,:]), G, S[t,z,y,x,:,:,:,:], G, S2[t,s,r,q,:,:,:,:], G, np.conj(M[t,s,r,q,:,:,:,:]), G, optimize='greedy')
                            correlator4[t,z,y,x,s,r,q] =  np.einsum('mlbc,mn,nkbc,op,pida,kl,ojda,ij', np.conj(M2[t,z,y,x,:,:,:,:]), G, S2[t,z,y,x,:,:,:,:], G, S[t,s,r,q,:,:,:,:], G, np.conj(M[t,s,r,q,:,:,:,:]), G, optimize='greedy')
                            correlator[t,z,y,x,s,r,q] = correlator1[t,z,y,x,s,r,q] + correlator2[t,z,y,x,s,r,q] + correlator3[t,z,y,x,s,r,q] +correlator4[t,z,y,x,s,r,q]
    cnt += 1
    print('Contraction {:.2%} completed.'.format(cnt/tot))

p2max = 1
ptotmax = 2
N=10

momlist = []
for px in range(-N, N):
    for py in range(-N, N):
        for pz in range(-N, N):
            if (px ** 2 + py ** 2 + pz ** 2) <= p2max:
                momlist.append([px,py,pz])

#correlator_FTed = np.zeros([8,4,4,4,4,4,4], dtype = np.complex128)
correlator_FTed1 = {}
correlator_FTed2 = {}
correlator_FTed3 = {}
correlator_FTed4 = {}
correlator_FTed = {}
tot =len(momlist) ** 2

print('Fourier transform starts.')
for p1 in momlist:
    for p2 in momlist:
        if ((p1[0]+p2[0]) ** 2 + (p1[2]+p2[2]) ** 2 + (p1[2]+p2[2]) ** 2) <= ptotmax:
            tmpFTlist1 = np.zeros(8, dtype = np.complex128)
            tmpFTlist2 = np.zeros(8, dtype = np.complex128)
            tmpFTlist3 = np.zeros(8, dtype = np.complex128)
            tmpFTlist4 = np.zeros(8, dtype = np.complex128)
            tmpFTlist = np.zeros(8, dtype = np.complex128)
            for t in range(8):
                tmpFT1 = 0
                tmpFT2 = 0
                tmpFT3 = 0
                tmpFT4 = 0
                tmpFT = 0
                for z in range(4):
                    for y in range(4):
                        for x in range(4):       
                            for s in range(4):
                                for r in range(4):
                                    for q in range(4):
                                        tmpFT1 += np.exp(1j*(np.pi/2)*(p1[0]*x+p1[1]*y+p1[2]*z + p2[0]*q+p2[1]*r+p2[2]*s)) * correlator1[t,z,y,x,s,r,q]
                                        tmpFT2 += np.exp(1j*(np.pi/2)*(p1[0]*x+p1[1]*y+p1[2]*z + p2[0]*q+p2[1]*r+p2[2]*s)) * correlator2[t,z,y,x,s,r,q]
                                        tmpFT3 += np.exp(1j*(np.pi/2)*(p1[0]*x+p1[1]*y+p1[2]*z + p2[0]*q+p2[1]*r+p2[2]*s)) * correlator3[t,z,y,x,s,r,q]
                                        tmpFT4 += np.exp(1j*(np.pi/2)*(p1[0]*x+p1[1]*y+p1[2]*z + p2[0]*q+p2[1]*r+p2[2]*s)) * correlator4[t,z,y,x,s,r,q]
                                        tmpFT += np.exp(1j*(np.pi/2)*(p1[0]*x+p1[1]*y+p1[2]*z + p2[0]*q+p2[1]*r+p2[2]*s)) * correlator[t,z,y,x,s,r,q]

                t_relative = t - t_0
                if (t_relative < 0):
                    t_relative += Nt
                tmpFTlist1[t_relative] = tmpFT1
                tmpFTlist2[t_relative] = tmpFT2
                tmpFTlist3[t_relative] = tmpFT3
                tmpFTlist4[t_relative] = tmpFT4
                tmpFTlist[t_relative] = tmpFT
            
            origin_phases = 0
            for p_comp in range(3):
                origin_phases -= (p1[p_comp] * origin_list[0][p_comp] + p2[p_comp] * origin_list[1][p_comp]) * np.pi/2;
            
            print(p1,p2,origin_phases)

            correlator_FTed1[(p1[0],p1[1],p1[2],p2[0],p2[1],p2[2])] = tmpFTlist1 * np.exp(1j*origin_phases)
            correlator_FTed2[(p1[0],p1[1],p1[2],p2[0],p2[1],p2[2])] = tmpFTlist2 * np.exp(1j*origin_phases)
            correlator_FTed3[(p1[0],p1[1],p1[2],p2[0],p2[1],p2[2])] = tmpFTlist3 * np.exp(1j*origin_phases)
            correlator_FTed4[(p1[0],p1[1],p1[2],p2[0],p2[1],p2[2])] = tmpFTlist4 * np.exp(1j*origin_phases)
            correlator_FTed[(p1[0],p1[1],p1[2],p2[0],p2[1],p2[2])] = tmpFTlist * np.exp(1j*origin_phases)
print('Fourier transform completed.')

# Python version
PY = correlator_FTed[(0,0,1,0,0,0)]
# C++ version
f = h5.File('/Users/haobo/lattice_qcd/test_run/lalibe/lalibe_pipi_spectrum.h5','r')
CXX = f['PS']['pipi12']['x_0_0_1_0__y_2_0_0_0']['ptotx0_ptoty0_ptotz1']['px0_py0_pz1_qx0_qy0_qz0'][()]
f.close()
print(PY)
print(CXX)