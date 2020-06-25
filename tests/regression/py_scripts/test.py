import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

f = h5.File('/Users/haobo/lattice_qcd/test_run/lalibe/lalibe_pipi_spectrum.h5','r')
g = h5.File('/Users/haobo/lattice_qcd/test_run/lalibe/lalibe_pipi_spectrum_0625.h5','r')

A = f['PS']['kk']['diagram0']['x_0_0_0_0__y_2_2_2_0']['ptotx0_ptoty0_ptotz1']['px0_py0_pz1_qx0_qy0_qz0'][()]
B = g['PS']['kk']['diagram0']['x_0_0_0_0__y_2_2_2_0']['ptotx0_ptoty0_ptotz1']['px0_py0_pz1_qx0_qy0_qz0'][()]
print(A)
print(B)

f.close()
g.close()
