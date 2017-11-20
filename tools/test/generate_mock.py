'''

171120 - Gabriele Girelli
Project: GPSeq validation with FISH

Aim:
	Generate mock figure to check dotter2gpseq.

'''

# DEPENDENCIES =================================================================

import numpy as np
from skimage.io import imread, imsave
from scipy.ndimage.measurements import center_of_mass
from skimage.morphology import ball

# PARAMETERS ===================================================================

# FUNCTIONS ====================================================================

# Generate image
i = np.zeros((1, 1000, 1000, 1000))

# Add sphere
i[0, 49:250, 49:250, 49:250] = ball(100)
i[i == 1] = 100

# Check center
#center_of_mass(i)

# Export figure
imsave("dapi_001.tif", i[0:1, 0:300, 0:300, 0:300].astype('u4'))

# Import figure
f = open("dapi_001.tif", 'rb')
i2 = imread(f)
f.close()

i2.shape

# RUN ==========================================================================

# END ==========================================================================

################################################################################
