import numpy as np
import sys
import contractions

'''
This function takes in the quantum numbers isospin, isospin z component, and strangeness to determine
which baryon is being referenced. The allowed baryons are:
Proton:     (1/2, 1/2, 0)
Neutron:    (1/2,-1/2, 0)
Delta++:    (3/2, 3/2, 0)
Delta+:     (3/2, 1/2, 0)
Delta0:     (3/2,-1/2, 0)
Delta-:     (3/2,-3/2, 0)
Omega:      (  0,   0,-3)
Lambda:     (  0,   0,-1)
Sigma+:     (  1,   1,-1)
Sigma0:     (  1,   0,-1)
Sigma-:     (  1,  -1,-1)
Cascade0:   (1/2, 1/2,-2)
Cascade-:   (1/2,-1/2,-2)

Output: A tuple of an (n,3) ndarray of integers 0,1,2, and a (n) vector of coeffecients.
    These represent n terms of 3 quark creation operators composing the creation operator of the desired baryon.
    In the flavor array, 0 is an up quark, 1 is a down quark, and 2 is a strange quark.
    The dirac and color indices are implied and will be tracked in a later step; This tracks the 
    flavor order of the quark creaton operators. The vector of coeffecients are the factors for the terms.

Reference:
Basak, S., Edwards, R., Fleming, G. T., Heller, U. M., Morningstar, C., Richards, D., Sato, I., & Wallace, S. J. (2005). Clebsch-Gordan construction of lattice interpolated fields for excited baryons. Physical Review D, 72(7). https://doi.org/10.1103/physrevd.72.074501 
'''
def flavor_terms(isospin, isospin_z, strangeness):

    qnums = (isospin, isospin_z, strangeness)

    flavors = 0
    coeffs = 0

    if qnums==(0,0,-1):
        # Lambda
        flavors = [[0,1,2], [1,0,2]]
        coeffs = [1,-1]/np.sqrt(2)
    elif qnums==(0,0,-3):
        # Omega
        flavors = [[2,2,2]]
        coeffs = [1]
    elif qnums==(0.5,0.5,0):
        # Proton
        flavors = [[0,1,0],[1,0,0]]
        coeffs = [1,-1]/np.sqrt(2)
    elif qnums==(0.5,-0.5,0):
        # Neutron
        flavors = [[0,1,1],[1,0,1]]
        coeffs = [1,-1]/np.sqrt(2)
    elif qnums==(0.5,0.5,-2):
        # Cascade0
        flavors = [[2,2,0]]
        coeffs = [1]
    elif qnums==(0.5,-0.5,-2):
        # Cascade-
        flavors = [[2,2,1]]
        coeffs = [1]
    elif qnums==(1,1,-1):
        # Sigma+
        flavors = [[0,0,2]]
        coeffs = [1]
    elif qnums==(1,0,-1):
        # Sigma0
        flavors = [[0,1,2],[1,0,2]]
        coeffs = [1,1]/np.sqrt(2)
    elif qnums==(1,-1,-1):
        # Sigma-
        flavors = [[1,1,2]]
        coeffs = [1]
    elif qnums==(1.5,1.5,0):
        # Delta++
        flavors = [[0,0,0]]
        coeffs = [1]
    elif qnums==(1.5,0.5,0):
        # Delta+
        flavors = [[0,0,1],[0,1,0],[1,0,0]]
        coeffs = [1,1,1]/np.sqrt(3)
    elif qnums==(1.5,-0.5,0):
        # Delta0
        flavors = [[0,1,1],[1,0,1],[1,1,0]]
        coeffs = [1,1,1]/np.sqrt(3)
    elif qnums==(1.5,-1.5,0):
        # Delta-
        flavors = [[1,1,1]]
        coeffs = [1]
    else:
        print("Error: Disallowed quantum numbers.")
        print("Isospin:",isospin)
        print("Isospin_z:",isospin_z)
        print("Strangeness:",strangeness)
        print("Check baryons.flavor_terms documentation for allowed quantum numbers.")
        sys.exit(-1)

    flavors = np.array(flavors)
    coeffs = np.array(coeffs, dtype=np.complex128)

    return flavors, coeffs

'''
This function takes in the quantum numbers isospin, strangeness, spin, spin_z, and parity
to both determine the baryon and spin/parity state. The allowed baryons are:
            (isospin, strangeness, spin)
Nucleon:    (1/2, 0,1/2)
Delta:      (3/2, 0,3/2)
Omega:      (  0,-3,3/2)
Lambda:     (  0,-1,1/2)
Sigma:      (  1,-1,1/2)
Cascade:    (1/2,-2,1/2)
Sigma*:     (  1,-1,3/2)
Cascade*:   (1/2,-2,3/2)
The allowed spin_z values are 1/2, -1/2 for spin = 1/2 particles and 3/2,1/2,-1/2,-3/2 for spin=3/2 particles.
Parity is either +1 or -1.

Output: A tuple of an (n,3) ndarray of integers 0,1,2,3, and a (n) vector of coeffecients.
    These represent n terms of baryon creation operators composing the creation operator of the desired baryon
    with the desired spin properties.
    In the spin array, 0 is spin up with positive parity, 1 is spin down with positive parity, 
    2 is spin up with negative parity, and 3 is spin down with negative parity.
    The color indices are implied and will be tracked in a later step; Each vector in the spin array tracks the 
    dirac indices of the quarks composing a baryon creaton operator. The vector of coeffecients are the factors
    for the linear combination of baryon creaton operators.
Reference:
Basak, S., Edwards, R., Fleming, G. T., Heller, U. M., Morningstar, C., Richards, D., Sato, I., & Wallace, S. J. (2005). Clebsch-Gordan construction of lattice interpolated fields for excited baryons. Physical Review D, 72(7). https://doi.org/10.1103/physrevd.72.074501 
'''
def spin_terms(isospin, strangeness, spin, spin_z, parity):

    def throw_spin_error():
        print("Error: Disallowed spin/parity quantum numbers.")
        print("Spin:",spin)
        print("Spin_z:",spin_z)
        print("Parity:",parity)
        print("Check baryons.spin_terms() documentation for allowed spin/parity quantum numbers.")
        sys.exit(-1)

    bnums = (isospin, strangeness)
    snums = (spin, spin_z)

    spins = 0
    coeffs = 0

    if (bnums==(0.5,0)) | (bnums==(0,-1)):
        # Nucleon, Lambda
        if snums==(0.5,0.5):
            spins = [[0,1,0]]
            coeffs = [1]
        elif snums==(0.5,-0.5):
            spins = [[0,1,1]]
            coeffs = [1]
        else:
            throw_spin_error()
    elif (bnums==(1.5,0)) | (bnums==(0,-3)):
        # Delta, Omega
        if snums==(1.5,1.5):
            spins = [[0,0,0]]
            coeffs = [1]
        elif snums==(1.5,0.5):
            spins = [[0,0,1]]
            coeffs = [np.sqrt(3)]
        elif snums==(1.5,-0.5):
            spins = [[0,1,1]]
            coeffs = [np.sqrt(3)]
        elif snums==(1.5,-1.5):
            spins = [[1,1,1]]
            coeffs = [1]
        else:
            throw_spin_error()
    elif (bnums==(1,-1)) | (bnums==(0.5,-2)):
        # Sigma, Cascade
        if snums==(0.5,0.5):
            spins = [[0,0,1],[0,1,0]]
            coeffs = np.sqrt(2/3)*np.array([1,-1])
        elif snums==(0.5,-0.5):
            spins = [[1,1,0],[0,1,1]]
            coeffs = np.sqrt(2/3)*np.array([-1,1])
        elif snums==(1.5,1.5):
            spins = [[0,0,0]]
            coeffs = [1]
        elif snums==(1.5,0.5):
            spins = [[0,0,1],[0,1,0]]
            coeffs = np.array([1,2])/np.sqrt(3)
        elif snums==(1.5,-0.5):
            spins = [[0,1,1],[1,1,0]]
            coeffs = np.array([2,1])/np.sqrt(3)
        elif snums==(1.5,-1.5):
            spins = [[1,1,1]]
            coeffs = [1]
        else:
            throw_spin_error()
    else:
        print("Error: Disallowed baryon quantum numbers.")
        print("Isospin:",isospin)
        print("Strangeness:",strangeness)
        print("Check baryons.spin_terms documentation for allowed quantum numbers.")
        sys.exit(-1)
    
    spins = np.array(spins)
    coeffs = np.array(coeffs, dtype=np.complex128)

    if parity==-1:
        spins += 2
    elif parity!=1:
        throw_spin_error()

    return spins, coeffs

'''
This function takes in six quantum numbers which uniquely determine the baryon, its spin state, and its parity. It returns a two point correlation function for this baryon state,
which in turn takes in three quark correlation functions and outputs the correlation function of the specified baryon state.

Inputs:
isospin: The total isospin of the baryon. See the documentation of baryons.flavor_terms for allowed values.
isospin_z: The z component of the isospin. See the documentation of baryons.flavor_terms for allowed values.
strangeness: The number of strange quarks composing the baryon, times -1. See the documentation of baryons.flavor_terms for allowed values.
spin: The total spin of the baryon. See the documentation of baryons.spin_terms for allowed values.
spin_z: The z component of the spin. See the documentation of baryons.spin_terms for allowed values.
parity: The eigenvalue of the baryon under a parity transformation. Must be either +1 or -1.

Output:
A two-point correlation function for the particular octuplet or dectuplet baryon with the specified spin/parity state. This function has the following properties:
    Input:
    q1, q2, q3: These are the three quark two-point correlation functions. The must be passed into this function "flavor sorted". This means passing any up quarks first,
        down quarks second, and strange quarks last.
    
    Output: 
    A tensor with the same shape as the first four axes of the quark correlation functions q1, q2, q3 passed as inputs. This tensor is the two-point correlator for the specified baryon state.
'''
def two_point_correlator(isospin, isospin_z, strangeness, spin, spin_z, parity):

    flavors, flavor_coeffs = flavor_terms(isospin, isospin_z, strangeness)
    spins, spin_coeffs = spin_terms(isospin, strangeness, spin, spin_z, parity)

    def correlator(q1, q2, q3):

        space_shape = q1.shape[:4]
        space_shape2 = q2.shape[:4]
        space_shape3 = q3.shape[:4]

        if (space_shape!=space_shape2) | (space_shape!=space_shape3):
            print("Error: quark propagators do not have the same shape for the space-time indices. These are expected to be the first four indices of the propagator.")
            print("q1 spacetime shape:", space_shape)
            print("q2 spacetime shape:", space_shape2)
            print("q3 spacetime shape:", space_shape3)
            sys.exit(-1)
        
        out = np.zeros(space_shape, dtype=np.complex128)

        for ai,af in enumerate(flavors):
            for i,f in enumerate(flavors):
                f = f[::-1]
                for aj,al in enumerate(spins):
                    for j,l in enumerate(spins):
                        l = l[::-1]

                        out += flavor_coeffs[ai] * flavor_coeffs[i] * spin_coeffs[aj] * spin_coeffs[j] * contractions.contract_quark_term(q1,q2,q3,f,l,af,al)

        return out

    return correlator