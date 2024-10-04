import numpy as np
import h5py as h5
import contractions
import baryons
import gamma
import sys
import argparse
import os

# Name: (isospin, isospin_z, strangeness, spin)
h5_baryons = {'delta_m': (3/2, -3/2, 0, 3/2),
              'delta_z': (3/2, -1/2, 0, 3/2),
              'delta_p': (3/2, 1/2, 0, 3/2),
              'delta_pp': (3/2, 3/2, 0, 3/2),
              'lambda_z': (0, 0, -1, 1/2),
              'neutron': (1/2, -1/2, 0, 1/2),
              'proton': (1/2, 1/2, 0, 1/2),
              'omega_m': (0, 0, -3, 3/2),
              'sigma_m': (1, -1, -1, 1/2),
              'sigma_z': (1, 0, -1, 1/2),
              'sigma_p': (1, 1, -1, 1/2),
              'sigma_star_m': (1, -1, -1, 3/2),
              'sigma_star_z': (1, 0, -1, 3/2),
              'sigma_star_p': (1, 1, -1, 3/2),
              'xi_m': (1/2, -1/2, -2, 1/2),
              'xi_z': (1/2, 1/2, -2, 1/2),
              'xi_star_m': (1/2, -1/2, -2, 3/2),
              'xi_star_z': (1/2, 1/2, -2, 3/2)}

spin_zs = {1/2: [('dn',-1/2), ('up',1/2)],
           3/2: [('dndn',-3/2), ('dn',-1/2), ('up',1/2), ('upup',3/2)]}

def main():
    cur_path = os.path.dirname(os.path.realpath(__file__))
    default_known_path = os.path.join(cur_path, '..', 'known_results')

    parser = argparse.ArgumentParser(description="Create baryon two point correlation functions from quark propagators and check results against an h5 file.")
    parser.add_argument("--known_file", type=str, help="h5 file containing known results", default=os.path.join(default_known_path,'lalibe_2pt_spectrum.h5'))
    parser.add_argument("--quark_prop_file", type=str, help="h5 file containing quark propagators", default=os.path.join(default_known_path,'test_propagator.h5'))
    parser.add_argument("--up_quark_node", type=str, help="Node in quark_prop_file containg the up quark propagator", default='/sh_sig2p0_n5/PS_up')
    parser.add_argument("--dn_quark_node", type=str, help="Node in quark_prop_file containg the down quark propagator", default='/sh_sig2p0_n5/PS_dn')
    parser.add_argument("--strange_quark_node", type=str, help="Node in quark_prop_file containg the strange quark propagator", default='/sh_sig2p0_n5/PS_strange')
    parser.add_argument("-t", "--tol", type=float, help="Relative tolerance to accept or reject the new values compared to the known values for the baryon propagators", default=1.e-7)
    parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
    args = parser.parse_args()

    U    = gamma.U_DR_to_DP
    Uadj = gamma.U_DR_to_DP_adj

    f = h5.File(args.quark_prop_file)

    ps_up = f[args.up_quark_node][()]
    ps_dn = f[args.dn_quark_node][()]
    ps_strange = f[args.strange_quark_node][()]
    f.close()

    # Rotate from Degrand-Rossi to Dirac-Pauli basis
    ps_up_DP = np.einsum('ik,tzyxklab,lj->tzyxijab',Uadj,ps_up,U)
    ps_dn_DP = np.einsum('ik,tzyxklab,lj->tzyxijab',Uadj,ps_dn,U)
    ps_strange_DP = np.einsum('ik,tzyxklab,lj->tzyxijab',Uadj,ps_strange,U)
    ps_quarks_DP = (ps_up_DP, ps_dn_DP, ps_strange_DP)

    if args.verbose:
        print(args.quark_prop_file+args.up_quark_node+' shape')
        print(ps_up_DP.shape)
        print(args.quark_prop_file+args.dn_quark_node+' shape')
        print(ps_dn_DP.shape)
        print(args.quark_prop_file+args.strange_quark_node+' shape')
        print(ps_strange_DP.shape)

    passed = 0
    failed = 0
    signed = 0

    known_results_file = h5.File(args.known_file)

    for baryon, qnums in h5_baryons.items():
        for parity in [1, -1]:
            for s_name, spin_z in spin_zs[qnums[3]]:

                if args.verbose: print(baryon,s_name)

                up = int(3/2 + 0.5*qnums[2] + qnums[1])
                strange = -qnums[2]
                flavors = np.ones(3, dtype=int)
                flavors[:up] = 0
                if strange != 0: flavors[-strange:] = 2

                q1 = ps_quarks_DP[flavors[0]]
                q2 = ps_quarks_DP[flavors[1]]
                q3 = ps_quarks_DP[flavors[2]]

                #proton = contractions.proton_spin_contract(ps_up_DP,ps_up_DP,ps_dn_DP,corr,spin)
                baryon_correlator = baryons.two_point_correlator(isospin=qnums[0], isospin_z=qnums[1], strangeness=qnums[2], spin=qnums[3], spin_z=spin_z, parity=parity)
                baryon_tensor = baryon_correlator(q1, q2, q3)
                baryon_time = np.einsum('tzyx->t',baryon_tensor)
                '''
                for t in range(Nt):
                    print(t,proton_up_time[t])
                '''
                p_name = ""
                if parity==-1:
                    p_name = "_np"

                known_baryon = known_results_file['PS/'+baryon+p_name+'/spin_'+s_name+'/x0_y0_z0_t0/px0_py0_pz0'][()]

                known_ratio = np.real(baryon_time/known_baryon)
                
                if np.all(np.abs(known_ratio - 1) <= args.tol):
                    if args.verbose: print('    PASS')
                    passed += 1
                elif np.all(np.abs(known_ratio + 1) <= args.tol):
                    if args.verbose: print('    SIGN ERROR')
                    signed += 1
                else:
                    if args.verbose: 
                        print('    FAIL')
                        print('Error ratio:',max(np.abs(known_ratio)))
                    failed += 1

    known_results_file.close()

    print("\n# Passed:",passed)
    print("# Sign error:",signed)
    print("# Failed",failed)

if __name__ == '__main__':
    main()

