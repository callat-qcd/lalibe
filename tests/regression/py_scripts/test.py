import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

f = h5.File('/Users/haobo/lattice_qcd/test_run/lalibe/lalibe_pipi_spectrum.h5','r')
'''g = h5.File('/Users/haobo/lattice_qcd/test_run/lalibe/lalibe_pipi_spectrum.h5','r')
C1 = f['PS']['pipi12']['x_0_0_0_0__y_2_2_2_0']['ptotx1_ptoty-1_ptotz0']['px1_py0_pz0_qx0_qy-1_qz0'][()]
C2 = f['PS']['kpi12']['x_0_0_0_0__y_2_2_2_0']['ptotx1_ptoty-1_ptotz0']['px0_py-1_pz0_qx1_qy0_qz0'][()]
C3 = f['PS']['kk11']['x_0_0_0_0__y_0_0_0_0']['ptotx1_ptoty-1_ptotz0']['px0_py-1_pz0_qx1_qy0_qz0'][()]
D1 = g['PS']['pipi12']['x_0_0_0_0__y_2_2_2_0']['ptotx1_ptoty-1_ptotz0']['px1_py0_pz0_qx0_qy-1_qz0'][()]
D2 = g['PS']['kpi12']['x_0_0_0_0__y_2_2_2_0']['ptotx1_ptoty-1_ptotz0']['px0_py-1_pz0_qx1_qy0_qz0'][()]
D3 = g['PS']['kk11']['x_0_0_0_0__y_0_0_0_0']['ptotx1_ptoty-1_ptotz0']['px0_py-1_pz0_qx1_qy0_qz0'][()]
'''
X = f['PS']['kk']['diagram0']['x_0_0_0_0__y_2_2_2_0']['ptotx0_ptoty0_ptotz1']['px0_py0_pz1_qx0_qy0_qz0'][()]
print(X)

f.close()
#g.close()

'''print(C1-D1)
print(C2-D2)
print(C3-D3)'''
