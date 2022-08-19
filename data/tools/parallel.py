from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
    idata = 1
    comm.send(idata, dest=1)
    print ('This is process {}'.format(rank), '\nData send to process 1 successfully!')
elif rank == 1:
    idata = comm.recv(source=0)
    print ('This is process {}, data is '.format(rank),idata)
