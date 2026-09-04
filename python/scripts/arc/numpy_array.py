#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

a = np.array([1,2]) 
A = np.array([[1,2],[3,4]])
start = np.array([0,0])
end = np.array([3,4])



n = 5

t=np.linspace(0.0,1.0,n,endpoint=False)

seg = start[None,:] + t[:, None]*(end - start)[None,:]

seg = start[None, :] + t[:, None] @ (end-start)[None, :]

print( t[:, None]*(end - start)[None,:])

print( t[:, None] @ (end-start)[None, :] )

# print(start[None,:].shape)
# print(start[None,:])