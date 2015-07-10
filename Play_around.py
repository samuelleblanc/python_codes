# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%config InlineBackend.rc = {}
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
%matplotlib inline
import numpy as np
import scipy.io as sio
from mpl_toolkits.basemap import Basemap

# <codecell>

import os

# <codecell>

os.path.exists

# <codecell>

os.path.exists('C:\\Users\\sleblan2\\Research\\python_codes\\')

# <codecell>

fu = np.array([[550.000 , 1.286162e-06 , 3.552874e+02 , 2.316014e+01 , 1.298847e-07 , 4.876075e+01 , 6.735294e+00], 
     [550.000 , 8.816101e+02 , 1.048253e+02 , 5.573513e+02 , 8.903053e+01 , 2.529635e+01 , 9.411757e+01], 
     [550.000 , 1.114990e+03 , 3.609851e-05 , 5.377008e+02 , 1.125986e+02 , 1.569951e-05 , 8.701024e+01]])

sbdart = np.array([[250.000 , 1.296362e-03 , 3.430699e+05 , 1.931734e+04 , 1.309148e-04 , 4.699350e+04 , 6.008736e+03],
          [250.000 , 8.569704e+05 , 1.053688e+05 , 5.385019e+05 , 8.654226e+04 , 2.511391e+04 , 9.108232e+04],
          [250.000 , 1.110150e+06 , 3.641137e-02 , 5.153732e+05 , 1.121100e+05 , 1.576419e-02 , 8.285177e+04]])

# <codecell>

fu.shape

# <codecell>

plt.plot(fu,'+-')
plt.plot(sbdart/1000.0,'x-')

