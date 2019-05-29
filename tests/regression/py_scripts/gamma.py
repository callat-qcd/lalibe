import numpy as np
verbose=False
test=True

one = np.zeros([4,4],dtype=np.complex128)
for i in range(4):
    one[i,i] = 1.


g_1 = np.zeros([4,4],dtype=np.complex128) #this is real and imag Float32
g_1[0,3] = 0 + 1.j
g_1[1,2] = 0 + 1.j
g_1[2,1] = 0 - 1.j
g_1[3,0] = 0 - 1.j
# g_2 = -i |   0   tau_2 |
#          | -tau_2   0  |
g_2 = np.zeros([4,4],dtype=np.complex128)
g_2[0,3] = -1. + 0.j
g_2[1,2] = 1. + 0.j
g_2[2,1] = 1. + 0.j
g_2[3,0] = -1. + 0.j

# g_3 = i |   0   tau_3 |
#         | -tau_3   0  |
g_3 = np.zeros([4,4],dtype=np.complex128)
g_3[0,2] = 0 + 1.j
g_3[1,3] = 0 - 1.j
g_3[2,0] = 0 - 1.j
g_3[3,1] = 0 + 1.j

# g_4 =   | 0   1 |
#         | 1   0 |
g_4 = np.zeros([4,4],dtype=np.complex128)
g_4[0,2] = 1. + 0.j
g_4[1,3] = 1. + 0.j
g_4[2,0] = 1. + 0.j
g_4[3,1] = 1. + 0.j

# g_5 = g_1 g_2 g_3 g_4
g_5 = np.dot(g_1,np.dot(g_2,np.dot(g_3,g_4)))

''' This DR_to_DP is in the CODE of chroma (vs the notes)
    it gives an all (-) gamma_5 while the notes would give an all +
    g_mu(DP) = U+ g_mu(DR) U
    q(DP)    = U+ q(DR)
    S(DP)    = U+ S(DR) U
'''
U_DR_to_DP = np.dot(g_5,np.dot(g_4,np.dot(g_1, g_3))) - np.dot(g_1,g_3)
U_DR_to_DP = U_DR_to_DP / np.sqrt(2)
U_DR_to_DP_adj = np.real(np.conj(U_DR_to_DP.T))

if __name__ == '__main__':
    if verbose:
        print('\ng_1 = ')
        print(g_1)
        print('\ng_2 = ')
        print(g_2)
        print('\ng_3 = ')
        print(g_3)
        print('\ng_4 = ')
        print(g_4)
        print('\ng_5 = ')
        print(g_5)

    if test:
        print('\n 1')
        print(one)

        print('\nreal[ (g_5 g_4 -1) g_1 g_3 ] * sqrt(2)')
        print(np.real(U_DR_to_DP)*np.sqrt(2))

        print('\nUadj')
        print(U_DR_to_DP_adj*np.sqrt(2))

        print('\nUadj . U')
        print(np.real(np.dot(U_DR_to_DP_adj,U_DR_to_DP)))

        U = U_DR_to_DP
        Uadj = np.conj(U_DR_to_DP.T)

        print('\nU+ g_4 U')
        print(np.real(np.dot(Uadj,np.dot(g_4,U))))

        print('\nU+ g_5 U')
        print(np.real(np.dot(Uadj,np.dot(g_5,U))))

        U_DR_to_DP = -np.dot(g_4,np.dot(g_1, g_3)) - np.dot(g_5,np.dot(g_1,g_3))
        print('\nreal[ (g_4 -g_5) g_1 g_3 ]')
        print(np.real(U_DR_to_DP))
        U_DR_to_DP = U_DR_to_DP / np.sqrt(2)

        print('\nU+ g_4 U')
        U = U_DR_to_DP
        Uadj = -U_DR_to_DP
        print(np.real(np.dot(Uadj,np.dot(g_4,U))))

        print('\nU+ g_5 U')
        U = U_DR_to_DP
        Uadj = -U_DR_to_DP
        print(np.real(np.dot(Uadj,np.dot(g_5,U))))
