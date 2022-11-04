import numpy as np
from collections import OrderedDict

"""
mass bins for rebinning mT and fits
"""

mass_bins_w = np.array([20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 120.0])

mass_bins_z = np.array([60.0, 68.0, 76.0, 84.0, 88.0, 92.0, 96.0, 104.0, 112.0, 120.0])

mass_bins_test = OrderedDict()
mass_bins_test[0] = np.array([0., 10.0, 20, 30, 40, 50, 60, 70, 80, 90, 120])
mass_bins_test[1] = np.array([5., 15.0, 25, 35, 45, 55, 65, 75, 85, 95, 120])
mass_bins_test[2] = np.array([10., 20, 30, 40, 50, 60, 70, 80, 90, 120])
mass_bins_test[3] = np.array([15., 25, 35, 45, 55, 65, 75, 85, 95, 120])
mass_bins_test[4] = np.array([20., 30, 40, 50, 60, 70, 80, 90, 120])
mass_bins_test[5] = np.array([25., 35, 45, 55, 65, 75, 85, 95, 120])
mass_bins_test[6] = np.array([30., 40, 50, 60, 70, 80, 90, 100, 120])
mass_bins_test[7] = np.array([35., 45, 55, 65, 75, 85, 95, 120])
mass_bins_test[8] = np.array([40., 50, 60, 70, 80, 90, 120])
mass_bins_test[9] = np.array([45., 55, 65, 75, 85, 95, 120])
mass_bins_test[10] = np.array([50., 60, 70, 80, 90, 120])

Vptbins = np.array([0, 2.0, 4, 6, 8.0, 10.0, 14.0, 16.0, 20., 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 90.0, 110.0, 15.0, 200.0, 13000.0])
