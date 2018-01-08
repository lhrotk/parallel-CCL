#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <list>
#include <algorithm>
#include <utility>
#include <cmath>
#include <functional>
#include <cstring>
#include <cmath>
#include <limits>
#include "mpi.h"

#define NOMINMAX

#ifdef _MSC_VER		//for MS C
#include <ctime>
inline double get_time()
{
        return static_cast<double>(std::clock()) / CLOCKS_PER_SEC;
}
#else
#include <sys/time.h>	// for linux C
inline double get_time()
{
        timeval tv;
        gettimeofday(&tv, 0);
        return tv.tv_sec + 1e-6 * tv.tv_usec;
}
#endif

using namespace std;

void init_CCL(int Label[], int Reference[], int N)
{
	for (int id = 0; id < N; id++) Label[id] = Reference[id] = 0;
}

int scanning(int D[], int L[], int R[], int N, int W, int th){

	for (int id = 0; id < N; id++) {
		if(D[id]==255){
			if (id - W >= 0 && D[id-W] ==255 && id%W != 0 && D[id-1]==255){
				if(L[id-W]==L[id-1])
					L[id] = L[id-1];
				else if(L[id-W]<L[id-1]){
					int ref;
					int label1 = L[id-1];
					int label2 = L[id-W];
					do { label1 = R[ref = label1]; } while (ref ^ label1);
					do { label2 = R[ref = label2]; } while (ref ^ label2);
					R[max(label1, label2)] = min(label1, label2);
					L[id] = min(label1, label2);
				}else{	
					int ref;
					int label1 = L[id-1];
					int label2 = L[id-W];
					do { label1 = R[ref = label1]; } while (ref ^ label1);
					do { label2 = R[ref = label2]; } while (ref ^ label2);
					R[max(label1, label2)] = min(label1, label2);
					L[id] = min(label1, label2);
				}
			}else if(id - W >= 0 && D[id-W] ==255){
				L[id] = L[id-W];
			}else if(id%W != 0 && D[id-1]==255){
				L[id] = L[id-1];
			}else{
				L[id] = ++th;
				R[th] = th;
			}
		}
	}

	return th;
}


void labeling(int L[], int R[], int N)
{
	for (int id = 0; id < N; id++) {
		int label = L[id];
		int ref;
		do { label = R[ref = label]; } while (ref != label);
		L[id] = label;
	}
}

void labelingF(int L[], int globalR[], int N){
	for(int id = 0; id< N; id++){
		int label = L[id];
		int ref;
		do { label = globalR[ref = label]; } while( label != 0);
		L[id] = ref;
	}
}

class CCL {
public:
	vector<int> ccl(vector<int>& image, int W, int degree_of_connectivity, int th, int rank, int size);
};

vector<int> CCL::ccl(vector<int>& image, int W, int degree_of_connectivity, int th, int rank, int size)
{
	int N = image.size();
	int* globalR;
	int startIndex = 0;
	int* buf = new int[1]{0};
	for(int i=0; i<rank; i++){
		startIndex += (N+size-i-1)/size;
	}
	N = (N+size-rank-1)/size;
	int* D = static_cast<int*>(&image[startIndex]);
	th = 0;
	int* Label = new int[N];
	int* Reference = new int[N];
	int* boundary = new int[W];
	int* refLabel = new int[W];
	int* eqLabel = new int[W+1]{0};
	int  acc;

	MPI_Barrier(MPI_COMM_WORLD);
	double start = get_time();
	init_CCL(Label, Reference, N);
	
	//local serial ccl
	th = scanning(D, Label, Reference, N, W, th);
	labeling(Label, Reference, N);
	
	if(rank==0){
		acc = th;
		for(int i=1; i<size; i++){
			MPI_Send(&acc, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Recv(buf, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			acc += *buf;
		}
		acc++;
		globalR = new int[acc]{0};
	}else{
		
		MPI_Send(&th, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Recv(buf, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		for(int i=0; i<N; i++){
			if(Label[i]!=0)
				Label[i]+= buf[0];
		}
	}
	// the boundary to its previous process
	if(rank!=0){
		MPI_Send(D, W, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
		MPI_Send(Label, W, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
	}
	//receive the boundary
	if(rank!=size-1){
		MPI_Recv(boundary, W, MPI_INT, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(refLabel, W, MPI_INT, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		int p = 0;
		for(int i=0; i<W; i++){
			if(D[N-W+i]!=0&&boundary[i]!=0){
				eqLabel[p++] = refLabel[i];
				eqLabel[p++] = Label[N-W+i];
				i++;//two neighboring column can not generate different equal label couple
			}
		}

		//tell process 0 equal label couples
		if(rank!=0){
			MPI_Send(eqLabel, W+1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}else{
			//process generate global equal label array and tell each process
			for(int j=0; j<W; j+=2){
                                        int label = eqLabel[j];
                                        int ref1, ref2;
                                        do{
                                                label = globalR[ref1 = label];
                                        }while(label!=0);
                                        label = eqLabel[j+1];
                                        do{
                                                label = globalR[ref2 = label];
                                        }while(label!=0);
                                        if(ref1!=ref2){
                                                globalR[max(ref1, ref2)] = min(ref1, ref2);
                                        }

                        }
			for(int i=1; i<size-1; i++){
				MPI_Recv(eqLabel, W+1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				for(int j=0; j<W; j+=2){
					int label = eqLabel[j];
					int ref1, ref2;
					do{
						label = globalR[ref1 = label];
					}while(label!=0);
					label = eqLabel[j+1];
					do{
						label = globalR[ref2 = label];
					}while(label!=0);
					if(ref1!=ref2){
						globalR[max(ref1, ref2)] = min(ref1, ref2);
					}
					
				}
			}
		}
	}

	//process 0 broadcast equal labels
	MPI_Bcast(&acc, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(rank!=0){
		globalR = new int[acc]{0};
	}
	MPI_Bcast(globalR, acc, MPI_INT, 0, MPI_COMM_WORLD);

	//find connect component across partition lines
	labelingF(Label, globalR, N);

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
		double end = get_time();
		cerr << "Time: " << end - start << endl;
	}

	//check the result
        bool judge = true;

        if(rank!=0){
                MPI_Send(D, W, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
                MPI_Send(Label, W, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
        }
        //receive the boundary
        if(rank!=size-1){
                MPI_Recv(boundary, W, MPI_INT, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(refLabel, W, MPI_INT, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                int p = 0;
                for(int i=0; i<W; i++){
                        if(D[N-W+i]!=0&&boundary[i]!=0&&refLabel[i]!=Label[N-W+i]){
                                judge = false;
                        }
                }
        }

        for (int i = 0; i< N; i++) {
                        if(i%W!=0&&D[i-1]==255&&D[i]==255&&Label[i]!=Label[i-1]){
                                judge = false;
                        }
                        if(i-W>=0&&D[i]==255&&D[i-W]==255&&Label[i]!=Label[i-W]){
                                judge = false;
                        }
        }

        if(judge){
                cerr << rank << ": result is correct" <<endl; 
        }


	vector<int> result(Label, Label + N);

	
	delete [] Label;
	delete [] Reference;

	return result;
}

void read_data(const string filename, vector<int>& image, int& W, int& degree_of_connectivity, int& threshold)
{
	fstream fs(filename.c_str(), ios_base::in);
	string line;
	stringstream ss;
	int data;

	getline(fs, line);
	ss.str(line);
	ss >> W >> degree_of_connectivity >> threshold;
	getline(fs, line);
	ss.str("");  ss.clear();
	for (ss.str(line); ss >> data; image.push_back(data));
}

int main(int argc, char* argv[])
{
	ios_base::sync_with_stdio(false);
	
	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " input_file" << endl;
		exit(1);
	}
	
	int rank, size;	
	vector<int> image;
	int W, degree_of_connectivity, th;//W: width
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	read_data(argv[1], image, W, degree_of_connectivity, th);
	CCL ccl;

	
	vector<int> result(ccl.ccl(image, W, degree_of_connectivity, th, rank, size));

	MPI_Finalize();
	return 0;
}
