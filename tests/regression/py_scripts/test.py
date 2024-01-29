import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

f = h5.File('/Users/haobo/lattice_qcd/test_run/lalibe/lalibe_pipi_spectrum.h5','r')
g = h5.File('/Users/haobo/lattice_qcd/test_run/lalibe/lalibe_pipi_spectrum_7.15saving.h5','r')
h = h5.File('/Users/haobo/lattice_qcd/test_run/lalibe/lalibe_pipi_spectrum_7.15saving_not_cleaned.h5','r')

A = f['PS']['pik']['diagram0']['0_1_0_0']['ptotx-1_ptoty0_ptotz1']['px0_py0_pz1_qx-1_qy0_qz0'][()]
B = g['PS']['pik']['diagram0']['x_0_1_0_0__y_0_1_0_0']['ptotx-1_ptoty0_ptotz1']['px0_py0_pz1_qx-1_qy0_qz0'][()]
C = h['PS']['pik']['diagram0']['0_1_0_0']['ptotx-1_ptoty0_ptotz1']['px0_py0_pz1_qx-1_qy0_qz0'][()]

print(A)
print(B)
print(C)

f.close()
g.close()
h.close()