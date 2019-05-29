import numpy as np
import sys

def two_eps_color_contract(q1,q2,q3):
    ''' take 3 quark props of definite spin and perform color contractions
        e.g. q1[:,:,:,:,sf,si,:,:]
        eps_a,b,c eps_d,e,f q1[a,d] q2[b,e] q3[c,f]
    '''
    result  = q1[:,:,:,:,0,0] * q2[:,:,:,:,1,1] * q3[:,:,:,:,2,2]
    result -= q1[:,:,:,:,0,0] * q2[:,:,:,:,1,2] * q3[:,:,:,:,2,1]
    result -= q1[:,:,:,:,0,1] * q2[:,:,:,:,1,0] * q3[:,:,:,:,2,2]
    result += q1[:,:,:,:,0,1] * q2[:,:,:,:,1,2] * q3[:,:,:,:,2,0]
    result += q1[:,:,:,:,0,2] * q2[:,:,:,:,1,0] * q3[:,:,:,:,2,1]
    result -= q1[:,:,:,:,0,2] * q2[:,:,:,:,1,1] * q3[:,:,:,:,2,0]
    result -= q1[:,:,:,:,0,0] * q2[:,:,:,:,2,1] * q3[:,:,:,:,1,2]
    result += q1[:,:,:,:,0,0] * q2[:,:,:,:,2,2] * q3[:,:,:,:,1,1]
    result += q1[:,:,:,:,0,1] * q2[:,:,:,:,2,0] * q3[:,:,:,:,1,2]
    result -= q1[:,:,:,:,0,1] * q2[:,:,:,:,2,2] * q3[:,:,:,:,1,0]
    result -= q1[:,:,:,:,0,2] * q2[:,:,:,:,2,0] * q3[:,:,:,:,1,1]
    result += q1[:,:,:,:,0,2] * q2[:,:,:,:,2,1] * q3[:,:,:,:,1,0]
    result -= q1[:,:,:,:,1,0] * q2[:,:,:,:,0,1] * q3[:,:,:,:,2,2]
    result += q1[:,:,:,:,1,0] * q2[:,:,:,:,0,2] * q3[:,:,:,:,2,1]
    result += q1[:,:,:,:,1,1] * q2[:,:,:,:,0,0] * q3[:,:,:,:,2,2]
    result -= q1[:,:,:,:,1,1] * q2[:,:,:,:,0,2] * q3[:,:,:,:,2,0]
    result -= q1[:,:,:,:,1,2] * q2[:,:,:,:,0,0] * q3[:,:,:,:,2,1]
    result += q1[:,:,:,:,1,2] * q2[:,:,:,:,0,1] * q3[:,:,:,:,2,0]
    result += q1[:,:,:,:,1,0] * q2[:,:,:,:,2,1] * q3[:,:,:,:,0,2]
    result -= q1[:,:,:,:,1,0] * q2[:,:,:,:,2,2] * q3[:,:,:,:,0,1]
    result -= q1[:,:,:,:,1,1] * q2[:,:,:,:,2,0] * q3[:,:,:,:,0,2]
    result += q1[:,:,:,:,1,1] * q2[:,:,:,:,2,2] * q3[:,:,:,:,0,0]
    result += q1[:,:,:,:,1,2] * q2[:,:,:,:,2,0] * q3[:,:,:,:,0,1]
    result -= q1[:,:,:,:,1,2] * q2[:,:,:,:,2,1] * q3[:,:,:,:,0,0]
    result += q1[:,:,:,:,2,0] * q2[:,:,:,:,0,1] * q3[:,:,:,:,1,2]
    result -= q1[:,:,:,:,2,0] * q2[:,:,:,:,0,2] * q3[:,:,:,:,1,1]
    result -= q1[:,:,:,:,2,1] * q2[:,:,:,:,0,0] * q3[:,:,:,:,1,2]
    result += q1[:,:,:,:,2,1] * q2[:,:,:,:,0,2] * q3[:,:,:,:,1,0]
    result += q1[:,:,:,:,2,2] * q2[:,:,:,:,0,0] * q3[:,:,:,:,1,1]
    result -= q1[:,:,:,:,2,2] * q2[:,:,:,:,0,1] * q3[:,:,:,:,1,0]
    result -= q1[:,:,:,:,2,0] * q2[:,:,:,:,1,1] * q3[:,:,:,:,0,2]
    result += q1[:,:,:,:,2,0] * q2[:,:,:,:,1,2] * q3[:,:,:,:,0,1]
    result += q1[:,:,:,:,2,1] * q2[:,:,:,:,1,0] * q3[:,:,:,:,0,2]
    result -= q1[:,:,:,:,2,1] * q2[:,:,:,:,1,2] * q3[:,:,:,:,0,0]
    result -= q1[:,:,:,:,2,2] * q2[:,:,:,:,1,0] * q3[:,:,:,:,0,1]
    result += q1[:,:,:,:,2,2] * q2[:,:,:,:,1,1] * q3[:,:,:,:,0,0]

    return result

def proton_spin_contract(q1,q2,q3,corr,spin):
    src_weights = np.zeros([2],dtype=np.complex128)
    src_weights[0] = 1./np.sqrt(2)
    src_weights[1] = -1./np.sqrt(2)
    snk_weights = np.zeros([4],dtype=np.complex128)
    snk_weights[0] =  1./np.sqrt(2)
    snk_weights[1] = -1./np.sqrt(2)
    snk_weights[2] =  1./np.sqrt(2)
    snk_weights[3] = -1./np.sqrt(2)

    src_spins = np.zeros([2,3],dtype=np.int)
    snk_spins = np.zeros([4,3],dtype=np.int)
    if corr == 'proton':
        if spin == 'up':
            src_spins[0,0] = 0; src_spins[0,1] = 0; src_spins[0,2] = 1;
            src_spins[1,0] = 0; src_spins[1,1] = 1; src_spins[1,2] = 0;

            snk_spins[0,0] = 0; snk_spins[0,1] = 0; snk_spins[0,2] = 1;
            snk_spins[1,0] = 0; snk_spins[1,1] = 1; snk_spins[1,2] = 0;
            snk_spins[2,0] = 0; snk_spins[2,1] = 0; snk_spins[2,2] = 1;
            snk_spins[3,0] = 1; snk_spins[3,1] = 0; snk_spins[3,2] = 0;
        elif spin == 'dn':
            src_spins[0,0] = 1; src_spins[0,1] = 0; src_spins[0,2] = 1;
            src_spins[1,0] = 1; src_spins[1,1] = 1; src_spins[1,2] = 0;

            snk_spins[0,0] = 1; snk_spins[0,1] = 0; snk_spins[0,2] = 1;
            snk_spins[1,0] = 1; snk_spins[1,1] = 1; snk_spins[1,2] = 0;
            snk_spins[2,0] = 0; snk_spins[2,1] = 1; snk_spins[2,2] = 1;
            snk_spins[3,0] = 1; snk_spins[3,1] = 1; snk_spins[3,2] = 0;
        else:
            print('unrecognized spin - aborting',spin)
            sys.exit(-1)
    elif corr == 'proton_np':
        if spin == 'up':
            src_spins[0,0] = 2; src_spins[0,1] = 2; src_spins[0,2] = 3;
            src_spins[1,0] = 2; src_spins[1,1] = 3; src_spins[1,2] = 2;

            snk_spins[0,0] = 2; snk_spins[0,1] = 2; snk_spins[0,2] = 3;
            snk_spins[1,0] = 2; snk_spins[1,1] = 3; snk_spins[1,2] = 2;
            snk_spins[2,0] = 2; snk_spins[2,1] = 2; snk_spins[2,2] = 3;
            snk_spins[3,0] = 3; snk_spins[3,1] = 2; snk_spins[3,2] = 2;
        elif spin == 'dn':
            src_spins[0,0] = 3; src_spins[0,1] = 2; src_spins[0,2] = 3;
            src_spins[1,0] = 3; src_spins[1,1] = 3; src_spins[1,2] = 2;

            snk_spins[0,0] = 3; snk_spins[0,1] = 2; snk_spins[0,2] = 3;
            snk_spins[1,0] = 3; snk_spins[1,1] = 3; snk_spins[1,2] = 2;
            snk_spins[2,0] = 2; snk_spins[2,1] = 3; snk_spins[2,2] = 3;
            snk_spins[3,0] = 3; snk_spins[3,1] = 3; snk_spins[3,2] = 2;
        else:
            print('unrecognized spin - aborting',spin)
            sys.exit(-1)
    else:
        print('unrecognized corr',corr)
        sys.exit(-1)

    nt,nz,ny,nx = q1.shape[0:4]
    result = np.zeros([nt,nz,ny,nx],dtype=np.complex128)
    for sf,wf in enumerate(snk_weights):
        for si,wi in enumerate(src_weights):
            tmp1 = q1[:,:,:,:,snk_spins[sf,0],src_spins[si,0]]
            tmp2 = q2[:,:,:,:,snk_spins[sf,1],src_spins[si,1]]
            tmp3 = q3[:,:,:,:,snk_spins[sf,2],src_spins[si,2]]
            result += two_eps_color_contract(tmp1,tmp2,tmp3) * wf * wi

    return result
