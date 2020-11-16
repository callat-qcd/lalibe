import os, sys
import h5py as h5
import numpy as np
np.set_printoptions(linewidth=120)
import diagrams

iso_symmetric=False

src = 'x0_y1_z0_t0'

''' NOTE - the ordering of the propagator array is
           T, Z, Y, X, sf, si, cf, ci
'''
x_origin = {'x0_y1_z0_t0':[0,1,0]}
x0 = np.array(x_origin[src])

#with h5.File('known_results/test_propagator.h5') as p5:
with h5.File('pipi_propagator.h5') as p5:
    prop_u = p5['sh_sig2p0_n5/PS_prop_u'][()]
    if iso_symmetric:
        prop_d = prop_u
    else:
        prop_d = p5['sh_sig2p0_n5/PS_prop_d'][()]
    prop_s = p5['sh_sig2p0_n5/PS_prop_s'][()]
L = prop_u.shape[1]
T = prop_u.shape[0]

''' NOTE: The sign convenction used in numpy FFT is opposite of Chroma
    numpy: FFT A = sum_n exp{ -i 2pi/L n } A(n)
    So, we can use the inverse FFT provided we also normalize properly
    L**3 * IFFT A = sum_n exp{ +i 2pi/L n} A(n)
'''
phase = np.zeros([L,L,L],dtype=np.complex128)
for iz,pz in enumerate([0,1,2,-1]):
    for iy,py in enumerate([0,1,2,-1]):
        for ix,px in enumerate([0,1,2,-1]):
            p = np.array([pz,py,px])
            phase[iz,iy,ix] = np.exp(-2j*np.pi/L * np.sum(p*x0))

''' pipi '''
pipi_pp = diagrams.pipi(prop_u,prop_d,prop_u,prop_d)
# add phase for shifted origin
pipi_pp = np.einsum('tzyxwvu,zyx->tzyxwvu', pipi_pp, phase)
pipi_pp = np.einsum('tzyxwvu,wvu->tzyxwvu', pipi_pp, phase)

pipi_mm = diagrams.pipi(prop_d,prop_u,prop_d,prop_u)
pipi_mm = np.einsum('tzyxwvu,zyx->tzyxwvu', pipi_mm, phase)
pipi_mm = np.einsum('tzyxwvu,wvu->tzyxwvu', pipi_mm, phase)

''' pi k '''
pik_pp = diagrams.piplus_kplus(prop_u,prop_d,prop_u,prop_s)
pik_pp = np.einsum('tzyxwvu,zyx->tzyxwvu', pik_pp, phase)
pik_pp = np.einsum('tzyxwvu,wvu->tzyxwvu', pik_pp, phase)

pik_p0 = diagrams.piminus_kminus(prop_u,prop_d,prop_s,prop_d)
pik_p0 = np.einsum('tzyxwvu,zyx->tzyxwvu', pik_p0, phase)
pik_p0 = np.einsum('tzyxwvu,wvu->tzyxwvu', pik_p0, phase)

pik_mm = diagrams.piminus_kminus(prop_d,prop_u,prop_s,prop_u)
pik_mm = np.einsum('tzyxwvu,zyx->tzyxwvu', pik_mm, phase)
pik_mm = np.einsum('tzyxwvu,wvu->tzyxwvu', pik_mm, phase)

''' k k '''
kk_pp = diagrams.pipi(prop_u,prop_s,prop_u,prop_s)
kk_pp = np.einsum('tzyxwvu,zyx->tzyxwvu', kk_pp, phase)
kk_pp = np.einsum('tzyxwvu,wvu->tzyxwvu', kk_pp, phase)

kk_00 = diagrams.pipi(prop_d,prop_s,prop_d,prop_s)
kk_00 = np.einsum('tzyxwvu,zyx->tzyxwvu', kk_00, phase)
kk_00 = np.einsum('tzyxwvu,wvu->tzyxwvu', kk_00, phase)

kk_mm = diagrams.pipi(prop_s,prop_u,prop_s,prop_u)
kk_mm = np.einsum('tzyxwvu,zyx->tzyxwvu', kk_mm, phase)
kk_mm = np.einsum('tzyxwvu,wvu->tzyxwvu', kk_mm, phase)

def get_mom(pmom):
    pq = pmom.split('_')
    px = int(pq[0].split('x')[1])
    py = int(pq[1].split('y')[1])
    pz = int(pq[2].split('z')[1])
    qx = int(pq[3].split('x')[1])
    qy = int(pq[4].split('y')[1])
    qz = int(pq[5].split('z')[1])

    return px,py,pz,qx,qy,qz

l_file = 'lalibe_pipi.h5'
d1 = h5.File(l_file,'r')

corrs = {'pip_pip':pipi_pp, 'pim_pim':pipi_mm, 'pip_kp':pik_pp, 'pim_km':pik_mm, 'kp_kp':kk_pp, 'km_km':kk_mm, 'k0_k0':kk_00, 'pip_k0b':pik_p0}
with h5.File('py_pipi.h5','w') as f5:
    for corr in ['pip_pip', 'pip_kp', 'kp_kp', 'k0_k0', 'pip_k0b']:
        print(corr)
        for ptot in d1['PS/'+corr+'/x0_y1_z0_t0']:
            print(ptot)
            for prel in d1['PS/'+corr+'/x0_y1_z0_t0/'+ptot]:
                px,py,pz,qx,qy,qz = get_mom(prel)
                py_data = corrs[corr][:,pz,py,px,qz,qy,qx]
                f5.create_dataset('PS/'+corr+'/x0_y1_z0_t0/'+ptot+'/'+prel,data=py_data)



if False:
    for prel in d1['PS/pipi/x0_y1_z0_t0/ptotx0_ptoty0_ptotz0']:
        px,py,pz,qx,qy,qz = get_mom(prel)
        py = pipi_mm[:,pz,py,px,qz,qy,qx]
        ll = d1['PS/pipi/x0_y1_z0_t0/ptotx0_ptoty0_ptotz0/'+prel][()]
        print("%30s" %prel, (py/ll).real)

    print('pipi_++')
    print('andre')
    print(pipi_pp[:,0,0,0,0,0,0].real)
    print('lalibe')
    print(d1['PS/pipi/x0_y1_z0_t0/ptotx0_ptoty0_ptotz0/px0_py0_pz0_qx0_qy0_qz0'][()].real)

    print('\npik_++')
    print('andre')
    print(pik_pp[:,0,0,0,0,0,0].real)
    print('lalibe')
    print(d1['PS/pik/x0_y1_z0_t0/ptotx0_ptoty0_ptotz0/px0_py0_pz0_qx0_qy0_qz0'][()].real)

    print('\nkk_++')
    print('andre')
    print(kk_pp[:,0,0,0,0,0,0].real)
    print('lalibe')
    print(d1['PS/kk/x0_y1_z0_t0/ptotx0_ptoty0_ptotz0/px0_py0_pz0_qx0_qy0_qz0'][()].real)
