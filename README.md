# parallel-CCL

The project includes three versions of connected component labelling algorithm. One is serial and the other two are parallel CCL implemented with OpenMPI and CUDA.

For detailed description about the algorithm, please refer to the PDF file.

To run the program, you need to install CUDA and OpenMPI to your machine.

Firstly use image2input.py to convert a gray image to binary image stored with plain text:

./image2input 1.jpg source

Then compile all files with command: Make

To run the algorithm:

Serial: ./ccl_le_cpu source > result

MPI: mpirun -np 32 ./ccl_le_mpi source > result

CUDA: ./ccl_le_gpu source > result

