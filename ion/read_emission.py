import numpy as np
import pdb

'''
Read an athena array (flattened 3D array) into a 3D numpy array. Output the
array to a text file with `printf("%d %d %d %g\n", k, j, i, emission(k, j, i))`
and remove header and footer output beforehand.
'''

#A(n,k,j,i) = A[i + N1*(j + N2*(k + N3*n))]

_, _, _, emis = np.loadtxt('out.txt').T
n1, n2, n3 = 128, 64, 64
emission = np.zeros([n3, n2, n1])
for k in range(0, n3):
    for j in range(0, n2):
        for i in range(0, n1):
            ind = i + n1*(j + n2*(k))
            emission[k, j, i] = emis[ind]

np.save("emission.npy", emission, allow_pickle=True)
