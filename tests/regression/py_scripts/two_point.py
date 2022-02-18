import numpy as np
import h5py as h5
import contractions
import gamma
import sys

U    = gamma.U_DR_to_DP
Uadj = gamma.U_DR_to_DP_adj

known_path = sys.argv[1]

f = h5.File(known_path+'/test_propagator.h5','r')

U_prop = f['sh_sig2p0_n5/PS_up'][()]
D_prop = f['sh_sig2p0_n5/PS_dn'][()]
f.close()

# Rotate from Degrand-Rossi to Dirac-Pauli basis
U_DP = np.einsum('ik,tzyxklab,lj->tzyxijab',Uadj,U_prop,U)
D_DP = np.einsum('ik,tzyxklab,lj->tzyxijab',Uadj,D_prop,U)
print(known_path+'/test_propagator.h5/sh_sig2p0_n5/PS_up shape')
print(U_DP.shape)
Nt = U_DP.shape[0]

print('\nPS')
print('(proton_up - known) / (proton_up + known) > 1.e-7?')

f = h5.File(known_path+'/lalibe_2pt_spectrum.h5','r')
for corr in ['proton','proton_np']:
    for spin in ['up','dn']:
        print(corr,spin)
        proton = contractions.proton_spin_contract(U_DP,U_DP,D_DP,corr,spin)
        proton_time = np.einsum('tzyx->t',proton)
        '''
        for t in range(Nt):
            print(t,proton_up_time[t])
        '''
        known_proton = f['PS/'+corr+'/spin_'+spin+'/x0_y0_z0_t0/px0_py0_pz0'][()]
        if np.any(np.real(proton_time - known_proton)/np.real(proton_time + known_proton) > 1.e-7):
            print('    FAIL')
        else:
            print('    PASS')
f.close()
