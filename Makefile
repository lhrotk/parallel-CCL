# Makfile for CCL

FLAGS += -O3 -w

all:
	nvcc $(FLAGS) ccl_le.cu -o ccl_le_gpu
	g++ $(FLAGS) ccl_le.cpp -o ccl_le_cpu
	mpic++ ${FLAGS} -std=c++11 ccl_le_mpi.cpp -o ccl_le_mpi

clean:
	rm -f *~ \#* ccl_le_gpu ccl_le_cpu ccl_le_mpi
#	rm -f *~ \#* $(TARGET)
