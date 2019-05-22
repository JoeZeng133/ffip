import h5py

f = h5py.File('output.h5')
print(len(f.keys()))