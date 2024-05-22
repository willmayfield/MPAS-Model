from netCDF4 import Dataset
import sys
import numpy as np

filename=(sys.argv[1])
f=Dataset(filename,'r+')

f.DX=3000.0
print(f.DX)

print(type(f.SF_SURFACE_PHYSICS))

f.SF_SURFACE_PHYSICS=np.int32(3)
print(f.SF_SURFACE_PHYSICS)


print(type(f.SF_SURFACE_PHYSICS))

f.close()

