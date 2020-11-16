import numpy as np

# NOTE - props are in T-Z-Y-X order in the file

# NOTE - we assume all quarks have the same spacetime source, x

''' NOTE: The sign convenction used in numpy FFT is opposite of Chroma
    numpy: FFT A = sum_n exp{ -i 2pi/L n } A(n)
    So, we can use the inverse FFT provided we also normalize properly
    L**3 * IFFT A = sum_n exp{ +i 2pi/L n} A(n)
'''

def pipi(q1,aq1,q2,aq2):
    return d1(q1,aq1,q2,aq2) + d2(q1,aq1,q2,aq2) + d3(q1,aq1,q2,aq2) + d4(q1,aq1,q2,aq2)

def piplus_kplus(q1,aq1,q2,aq2):
    ''' d3 exchanges the quarks '''
    return d1(q1,aq1,q2,aq2) + d3(q1,aq1,q2,aq2)

def piminus_kminus(q1,aq1,q2,aq2):
    ''' d2 exchanges the anti-quarks '''
    return d1(q1,aq1,q2,aq2) + d2(q1,aq1,q2,aq2)


def d1(q1,aq1,q2,aq2):
    ''' Diagram 1
        w ---<--- q1
          --->--- aq1

        z --->--- aq2
          ---<--- q2
    '''
    L = q1.shape[1]

    w  = np.einsum('tzyxrsab,tzyxrsab->tzyx',np.conjugate(aq1), q1)
    z  = np.einsum('tzyxrsab,tzyxrsab->tzyx',np.conjugate(aq2), q2)
    # single sink time, two different spatial locations
    wz = np.einsum('twvu,tzyx->twvuzyx', w, z)
    # make w->p and z->q momentum blocks
    pz = np.fft.ifftn(wz, axes=(1,2,3)) * L**3
    pq = np.fft.ifftn(pz, axes=(4,5,6)) * L**3

    return pq

def d2(q1,aq1,q2,aq2):
    ''' Diagram 2, this one has a relative (-) sign from fermion exchange
        it exchanges the anti-quarks
        w ---<--- q1
          \    / aq1
           \  /
            \/
            /\
           /  \
          /    \ aq2
        z ---<--- q2
    '''
    L = q1.shape[1]

    w  = np.einsum('tzyxirca,tzyxiscb->tzyxrsab',np.conjugate(aq2), q1)
    z  = np.einsum('tzyxirca,tzyxiscb->tzyxrsab',np.conjugate(aq1), q2)
    # single sink time, two different spatial locations
    wz = -np.einsum('tzyxrsab,twvusrba->tzyxwvu', w, z)
    # FFT
    pz = np.fft.ifftn(wz,axes=(1,2,3)) * L**3
    pq = np.fft.ifftn(pz,axes=(4,5,6)) * L**3

    return pq

def d3(q1,aq1,q2,aq2):
    ''' Diagram 3
        This diagram is the same as 2 with p <--> q
        it exchanges the quarks
        w\      / q1
          --->--- aq1
           \  /
            \/
            /\
           /  \
        z --->--- aq2
          /    \  q2
    '''
    L = q1.shape[1]

    w  = np.einsum('tzyxirca,tzyxiscb->tzyxrsab',np.conjugate(aq1), q2)
    z  = np.einsum('tzyxirca,tzyxiscb->tzyxrsab',np.conjugate(aq2), q1)
    # single sink time, two different spatial locations
    wz = -np.einsum('tzyxrsab,twvusrba->tzyxwvu', w, z)
    # FFT
    pz = np.fft.ifftn(wz,axes=(1,2,3)) * L**3
    pq = np.fft.ifftn(pz,axes=(4,5,6)) * L**3

    return pq

def d4(q1,aq1,q2,aq2):
    ''' Diagram 4
        This diagram is the same as 1 with p <--> q
        w\\     //q1
          \\   // aq1
           \\ //
            \\/
            //\
           // \\
          //   \\ aq2
        z//     \\ q2
    '''
    L = q1.shape[1]

    w  = np.einsum('tzyxrsab,tzyxrsab->tzyx',np.conjugate(aq2), q2)
    z  = np.einsum('tzyxrsab,tzyxrsab->tzyx',np.conjugate(aq1), q1)
    # single sink time, two different spatial locations
    wz = np.einsum('twvu,tzyx->twvuzyx', w, z)
    # make w->p and z->q momentum blocks
    pz = np.fft.ifftn(wz, axes=(1,2,3)) * L**3
    pq = np.fft.ifftn(pz, axes=(4,5,6)) * L**3

    return pq
