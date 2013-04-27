#! /usr/bin/env python

import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt

from myplot import getfield
from netCDF4 import Dataset
from scipy import ndimage

def disk_structure(n):
    struct = np.zeros((2 * n + 1, 2 * n + 1))
    x, y = np.indices((2 * n + 1, 2 * n + 1))
    mask = (x - n)**2 + (y - n)**2 <= n**2
    struct[mask] = 1
    return struct.astype(np.bool)


def granulometry(data, sizes=None):
    s = max(data.shape)
    if sizes == None:
        sizes = range(1, s/2, 2)
    granulo = [ndimage.binary_opening(data, \
            structure=disk_structure(n)).sum() for n in sizes]
    return granulo


#np.random.seed(1)
#n = 10
#l = 256
#im = np.zeros((l, l))
#points = l*np.random.random((2, n**2))
#im[(points[0]).astype(np.int), (points[1]).astype(np.int)] = 1
#im = ndimage.gaussian_filter(im, sigma=l/(4.*n))
#mask = im > im.mean()

namefile = "/home/aymeric/Big_Data/psfc_f18.nc"
nc  = Dataset(namefile)
psfc_full=getfield(nc,'PSFC')

#mask = psfc_full[250,100:200,100:200]
mask = psfc_full[250,:,:]
mask = ndimage.laplace(mask)
#mask = mask > 0.3

print mask.std()
print mask.mean()

lim = mask.std() * 12.

ind1 = np.where(mask > lim)
ind2 = np.where(mask < lim)
mask[ind1] = 1
mask[ind2] = 0

print mask


labeled_array, num_features = ndimage.label(mask)

print num_features

plt.pcolor(labeled_array)
plt.show()
exit()




mask = mask < mask.mean() - 0.3

#mask = mask > 0.3


#plt.pcolor(mask)
#plt.show()

granulo = granulometry(mask, sizes=[2,5,20])

print granulo

plt.figure(figsize=(6, 2.2))

plt.subplot(121)
plt.imshow(mask, cmap=plt.cm.gray)
opened = ndimage.binary_opening(mask, structure=disk_structure(10))
opened_more = ndimage.binary_opening(mask, structure=disk_structure(14))
plt.contour(opened, [0.5], colors='b', linewidths=2)
plt.contour(opened_more, [0.5], colors='r', linewidths=2)
plt.axis('off')
plt.subplot(122)
#plt.plot(np.arange(2, 19, 4), granulo, 'ok', ms=8)


plt.subplots_adjust(wspace=0.02, hspace=0.15, top=0.95, bottom=0.15, left=0, right=0.95)
plt.show()
