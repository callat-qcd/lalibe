import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

f = h5.File('/Users/haobo/lattice_qcd/test_run/lalibe/phaes_good_with_kk_and_kpi.h5','r')
C1 = f['PS']['kk12']['x_0_0_0_0__y_2_2_2_0']['ptotx1_ptoty-1_ptotz0']['px1_py0_pz0_qx0_qy-1_qz0'][()]
C2 = f['PS']['kpi12']['x_0_0_0_0__y_2_2_2_0']['ptotx1_ptoty-1_ptotz0']['px0_py-1_pz0_qx1_qy0_qz0'][()]
C3 = f['PS']['kpi21']['x_2_2_2_0__y_0_0_0_0']['ptotx1_ptoty-1_ptotz0']['px0_py-1_pz0_qx1_qy0_qz0'][()]

C4 = f['PS']['kk12']['x_0_0_0_0__y_2_2_2_0']['ptotx0_ptoty0_ptotz1']['px0_py0_pz1_qx0_qy0_qz0'][()]
f.close()

print(C1)
print(C2)
print(C3)
print(C4)