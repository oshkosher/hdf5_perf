# hdf5_perf
Performance comparison HDF5 vs direct MPI-IO


h5_collective <filename> <rows> <cols>

Write a 2-d dataset using HDF5 and MPI-IO directly. Compare the performance. 
Each process writes a vertical stripe of data, so collective IO is needed
to write the datafile.

extract_dxt_file.py and sequence_extract.py are used to parse the Darshan DxT
output produced by this program.

