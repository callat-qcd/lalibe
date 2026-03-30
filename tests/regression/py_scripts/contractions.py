import numpy as np
import collections
import sys

'''
Function which performs an arbitrary amount of epsilon contractions between specified axes of three nd arrays.
Input:
    a: ndarray of dtype complex128 which is the first of three objects undergoing the epsilon contraction.
    b: ndarray of dtype complex128 which is the second of three objects undergoing the epsilon contraction.
    c: ndarray of dtype complex128 which is the third of three objects undergoing the epsilon contraction.
    a_axes: The axes of a to contract. The number of elements in this array determines how many epsilon contractions total are going to be performed. 
            If an int, a single contraction will be performed.
            If b_axes and c_axes are undefined, they will be set to a_axes. This will ensure a, b, and c will be contracted over identical axes.
    b_axes: The axes of b to contract. The length of this array must match a_axes, or be an int if a_axes is an int.
    c_axes: The axes of c to contract. The length of this array must match a_axes, or be an int if a_axes is an int.
Output:
    An ndarray of dtype complex128 which is the output of the multi-contraction. It's shape will be what is left over after contracting the a_axes on a.
    By necessity, the residual shapes of b and c must match a. If this is not the case an error will be thrown.
'''
def multi_eps_contract(a,b,c,a_axes,b_axes=None,c_axes=None):
    if c_axes==None:
        if b_axes!=None:
            print("b_axes defined but c_axes is not. Either only a_axes is defined or a_axes, b_axes, and c_axes must be.")
            sys.exit(-1)
        b_axes = a_axes
        c_axes = a_axes
    if type(a_axes)==int:
        if (type(b_axes)!=int) | (type(c_axes)!=int):
            print("a_axes is an int while either one or both of b_axes and c_axes is not.")
            sys.exit(-1)
        a_axes=[a_axes]
        b_axes=[b_axes]
        c_axes=[c_axes]
    elif (len(a_axes) != len(b_axes)) | (len(b_axes) != len(c_axes)):
        print("Mismatched number of axes for a,b,c")
        sys.exit(-1)
    if len(a_axes) == 0:
        print("The axes' parameters have length 0. Must have length 1 or more.")
        sys.exit(-1)
    
    n = len(a_axes)
    elms = 6**n
    a_indices = np.array([slice(None) for _ in range(max(a_axes)+1)])
    b_indices = np.array([slice(None) for _ in range(max(b_axes)+1)])
    c_indices = np.array([slice(None) for _ in range(max(c_axes)+1)])

    a_indices[a_axes] = 0
    b_indices[b_axes] = 0
    c_indices[c_axes] = 0

    out_shape_a = np.shape(a[tuple(a_indices)])
    out_shape_b = np.shape(b[tuple(b_indices)])
    out_shape_c = np.shape(c[tuple(c_indices)])

    if (out_shape_a != out_shape_b) | (out_shape_b != out_shape_c):
        print("Remaining indices of a,b,c after contraction do not match.")
        print("a residual shape: ", out_shape_a)
        print("b residual shape: ", out_shape_b)
        print("c residual shape: ", out_shape_c)
        sys.exit(-1)

    ijk = np.array([[i for _ in range(n)] for i in range(3)])
    levels = [6**i for i in range(n)]
    parity = 1
    out = np.zeros(out_shape_a, dtype=np.complex128)
    for i in range(elms):
        a_indices[a_axes] = ijk[0]
        b_indices[b_axes] = ijk[1]
        c_indices[c_axes] = ijk[2]

        out += parity * a[tuple(a_indices)] * b[tuple(b_indices)] * c[tuple(c_indices)]

        if i == elms-1:
            return out

        for j,l in enumerate(levels):
            if i%l==0:

                r = (i/l)%3
                if r==0:
                    ijk[0,j],ijk[1,j] = ijk[1,j],ijk[0,j]
                elif r==1:
                    ijk[1,j],ijk[2,j] = ijk[2,j],ijk[1,j]
                else:
                    ijk[0,j],ijk[2,j] = ijk[2,j],ijk[0,j]

                parity *= -1

            else:
                break

'''
This function solves the general problem of contracting a term of 3 pairs of quark/antiquark creation opperators over their color indices. These operators will have definite
flavor and dirac indices. These terms are a common structure among all baryon correlation functions and having a way to generally contract these terms is helpful.

The general structure of these quark terms are as follows:
eps_abc * eps_a'b'c' * (q1)^a_d1 * (q2)^b_d2 * (q3)^c_d3 * (aq1)^a'_ad1 * (aq2)^b'_ad2 * (aq3)_c'_ad3
where the eps_ijk terms are Levi-Civita symbols, and the indices a, b, c, a', b', c' represent color indices to be summed over in the Einstien convention.
q1, q2, q3 are the quark annihilation operators where the flavor of each is determined by the flavor_indices array. 
d1, d2, d3 are the dirac indices of the quark annihilation operators dermined by the dirac_indices array.
aq1, aq2, aq3 are the antiquark annihilation operators where the flavor of each determined by the anti_flavor_indices array.
ad1, ad2, ad3 are the dirac indices of the antiquark annihilation operators determined by the anti_dirac_indices array.

Inputs:
q1, q2, q3: These are the three quark two-point correlation functions. The must be passed into this function "flavor sorted". This means passing any up quarks first,
    down quarks second, and strange quarks last.
flavor_indices: This is a 3 element array of the integers 0, 1, or 2. This array encodes the order of the flavor of the quark annihilation operators in the general term showcased above.
    0 represents an up quark, 1 a down quark, and 2 a strange quark. The number of each of these quarks must match the number of each flavor passed in q1, q2, q3. There is no test for this
    to be true in the function, so if this rule is broken you will likely get an incorrect result. 
dirac_indices: This is a 3 element array of the integers 0, 1, 2, or 3. This array encodes the dirac index of each of the three quark annihilation operators in the general term showcased above.
    0 represents a spin up, postive parity quark. 1 is a spin down, positive partiy quark. 2 is a spin up, negative parity quark. 3 is a spin down, negative parity quark.
anti_flavor_indices: This is a 3 element array of the integers 0, 1, or 2. This encodes the order of the flavor of the antiquark annihilation operators showcased in the general term above.
    See flavor_indices to see the meaning of each integer.
anti_dirac_indices: This is the analogous pair to dirac_indices for the antiquark annihilation operators in the general term above.

Output:
A tensor with the same shape as the first four axes of the quark correlation functions q1, q2, q3 passed as inputs. The tensor is equal to the general quark term described in the first
section of this comment. 
'''
def contract_quark_term(q1, q2, q3, flavor_indices, dirac_indices, anti_flavor_indices, anti_dirac_indices):

    qs = (q1,q2,q3)
    q_sort = {i: flavor_indices[i] for i in range(3)}
    q_key = sorted(q_sort, key=q_sort.get)
    q_sort = {q_key[i]: i for i in range(3)}

    qa = qs[q_sort[0]]
    qb = qs[q_sort[1]]
    qc = qs[q_sort[2]]

    F = np.full((3,3), False)
    for i,f in enumerate(flavor_indices):
        for j,af in enumerate(anti_flavor_indices):
            F[i,j] = (f==af)
    
    out = np.zeros(q1.shape[:4],dtype=np.complex128)
    d = dirac_indices
    ad = anti_dirac_indices

    if F[0,0]&F[1,1]&F[2,2]:
        out += multi_eps_contract(qa[:,:,:,:,d[0],ad[0]],qb[:,:,:,:,d[1],ad[1]],qc[:,:,:,:,d[2],ad[2]],[4,5])
    if F[0,0]&F[1,2]&F[2,1]:
        out += multi_eps_contract(qa[:,:,:,:,d[0],ad[0]],qb[:,:,:,:,d[1],ad[2]],qc[:,:,:,:,d[2],ad[1]],[4,5])
    if F[0,1]&F[1,0]&F[2,2]:
        out += multi_eps_contract(qa[:,:,:,:,d[0],ad[1]],qb[:,:,:,:,d[1],ad[0]],qc[:,:,:,:,d[2],ad[2]],[4,5])
    if F[0,1]&F[1,2]&F[2,0]:
        out += multi_eps_contract(qa[:,:,:,:,d[0],ad[1]],qb[:,:,:,:,d[1],ad[2]],qc[:,:,:,:,d[2],ad[0]],[4,5])
    if F[0,2]&F[1,0]&F[2,1]:
        out += multi_eps_contract(qa[:,:,:,:,d[0],ad[2]],qb[:,:,:,:,d[1],ad[0]],qc[:,:,:,:,d[2],ad[1]],[4,5])
    if F[0,2]&F[1,1]&F[2,0]:
        out += multi_eps_contract(qa[:,:,:,:,d[0],ad[2]],qb[:,:,:,:,d[1],ad[1]],qc[:,:,:,:,d[2],ad[0]],[4,5])

    return out


# Depreciated. See contractions.multi_eps_contract for an updated, general form of this function.
def two_eps_color_contract(q1,q2,q3):
    ''' take 3 quark props of definite spin and perform color contractions
        e.g. q1[:,:,:,:,sf,si,:,:]
        eps_a,b,c eps_d,e,f q1[a,d] q2[b,e] q3[c,f]
    '''
    return multi_eps_contract(q1,q2,q3,[4,5])

# Depreciated. See baryons.two_point_correlator for an updated, general form of this function.
def proton_spin_contract(q1,q2,q3,corr,spin):
    src_weights = np.zeros([2],dtype=np.complex128)
    src_weights[0] = 1./np.sqrt(2)
    src_weights[1] = -1./np.sqrt(2)
    snk_weights = np.zeros([4],dtype=np.complex128)
    snk_weights[0] =  1./np.sqrt(2)
    snk_weights[1] = -1./np.sqrt(2)
    snk_weights[2] =  1./np.sqrt(2)
    snk_weights[3] = -1./np.sqrt(2)
    if corr == 'proton':
        if spin == 'up':

            src_spins = [[0,0,1],
                         [0,1,0]]
        
            snk_spins = [[0,0,1],
                         [0,1,0],
                         [0,0,1],
                         [1,0,0]]
            
        elif spin == 'dn':

            src_spins = [[1,0,1],
                         [1,1,0]]
        
            snk_spins = [[1,0,1],
                         [1,1,0],
                         [0,1,1],
                         [1,1,0]]

        else:
            print('unrecognized spin - aborting',spin)
            sys.exit(-1)
    elif corr == 'proton_np':
        if spin == 'up':

            src_spins = [[2,2,3],
                         [2,3,2]]

            snk_spins = [[2,2,3],
                         [2,3,2],
                         [2,2,3],
                         [3,2,2]]
            
        elif spin == 'dn':

            src_spins = [[3,2,3],
                         [3,3,2]]

            snk_spins = [[3,2,3],
                         [3,3,2],
                         [2,3,3],
                         [3,3,2]]
        else:
            print('unrecognized spin - aborting',spin)
            sys.exit(-1)
    else:
        print('unrecognized corr',corr)
        sys.exit(-1)

    snk_spins = np.array(snk_spins)
    src_spins = np.array(src_spins)

    nt,nz,ny,nx = q1.shape[0:4]
    result = np.zeros([nt,nz,ny,nx],dtype=np.complex128)
    for sf,wf in enumerate(snk_weights):
        for si,wi in enumerate(src_weights):
            tmp1 = q1[:,:,:,:,snk_spins[sf,0],src_spins[si,0]]
            tmp2 = q2[:,:,:,:,snk_spins[sf,1],src_spins[si,1]]
            tmp3 = q3[:,:,:,:,snk_spins[sf,2],src_spins[si,2]]
            result += two_eps_color_contract(tmp1,tmp2,tmp3) * wf * wi

    return result
