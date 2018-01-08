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



#define NOMINMAX

#ifdef _MSC_VER
#include <ctime>
inline double get_time()
{
        return static_cast<double>(std::clock()) / CLOCKS_PER_SEC;
}
#else
#include <sys/time.h>
inline double get_time()
{
        timeval tv;
        gettimeofday(&tv, 0);
        return tv.tv_sec + 1e-6 * tv.tv_usec;
}
#endif

using namespace std;

//const int BLOCK = 128;
const int BLOCK = 1024;

__global__ void init_CCL(int L[], int N)
{
	for(int id = blockIdx.x * blockDim.x + threadIdx.x; id<N; id+=blockDim.x*gridDim.x){
		L[id] = id;
	}
}

__global__ void scanning(int D[], int L[], int N, int W, int th)
{
	for(int id = blockIdx.x * blockDim.x + threadIdx.x; id<N; id+=blockDim.x*gridDim.x){
		if(D[id]==0)
			L[id]=0;
		else{
			if(id-W>=0 && D[id-W]!=0)
				L[id] = min(L[id], L[id-W]);
			if(id%W!=N && D[id-1]!=0)
				L[id] = min(L[id], L[id-1]);
		}
	}
}

__global__ void analysis(int D[], int L[], int W, int N)
{
	for(int id = blockIdx.x * blockDim.x + threadIdx.x; id<N; id+=blockDim.x*gridDim.x){
		int ref;
		if(D[id]!=0&&id-W>=0 && D[id-W]!=0&&id%W!=0 && D[id-1]!=0){
			int label1 = L[id-1];
			int label2 = L[id-W];
			do { label1 = L[ref = label1]; } while (ref ^ label1);
			do { label2 = L[ref = label2]; } while (ref ^ label2);
			if(label1>label2){
				atomicAdd(L+label1, label2-label1);
			}else if(label1<label2){
				atomicAdd(L+label2, label1-label2);
			}else{
				L[id] = label1;
			}
		}else{
			int label1 = L[id];
			do { label1 = L[ref = label1]; } while (ref ^ label1);
			L[id] = label1;
		}
	}
}

__global__ void labeling(int D[], int L[], int N)
{
	for (int id = blockIdx.x * blockDim.x + threadIdx.x; id<N; id+=blockDim.x*gridDim.x) {
		int label = L[id];
		int ref;
		do { label = L[ref = label]; } while (ref != label);
		L[id] = label;
	}
}

class CCL {
private:
	int* Dd;
	int* Ld;
	int* Rd;

public:
	vector<int> cuda_ccl(vector<int>& image, int W, int degree_of_connectivity, int threshold);
};

vector<int> CCL::cuda_ccl(vector<int>& image, int W, int degree_of_connectivity, int threshold)
{
	int* D = static_cast<int*>(&image[0]);
	int N = image.size();
	int* Label = new int[N];
	cudaMalloc((void**)&Ld, sizeof(int) * N);
	//cudaMalloc((void**)&Rd, sizeof(int) * N);
	cudaMalloc((void**)&Dd, sizeof(int) * N);

	
	dim3 grid(6, 1, 1);
	dim3 threads(1024, 1, 1);
	double start = get_time();
	
	cudaMemcpy(Dd, D, sizeof(int) * N, cudaMemcpyHostToDevice);

	
	init_CCL<<<grid, threads>>>(Ld, N);
	scanning<<<grid, threads>>>(Dd, Ld, N, W, threshold);
	analysis<<<grid, threads>>>(Dd, Ld, W, N);
	labeling<<<grid, threads>>>(Dd, Ld, N);
	cudaDeviceSynchronize();
	cudaMemcpy(Label, Ld, sizeof(int) * N, cudaMemcpyDeviceToHost);
	
	double end = get_time();
	cerr << "Time: " << end - start << endl;


	cudaFree(Dd);
	cudaFree(Ld);
	cudaFree(Rd);
	vector<int> result(Label, Label + N);

	delete [] Label;
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

	vector<int> image;
	int W, degree_of_connectivity, threshold;
	read_data(argv[1], image, W, degree_of_connectivity, threshold);

	CCL ccl;

	//double start = get_time();
	vector<int> result(ccl.cuda_ccl(image, W, degree_of_connectivity, threshold));
	//double end = get_time();
	//cerr << "Time: " << end - start << endl;

	cout << result.size() << endl; /// number of pixels
	cout << W << endl; /// width
	bool judge = true;
	for (int i = 0; i < image.size()/W; i++) {
		for (int j = 0; j < W; j++) cout << result[i*W+j] << " ";
		cout << endl;
	}
	for (int i = 0; i < image.size()/W-1; i++) {
		for (int j = 0; j < W-1; j++){
			if(image[i*W+j]==255&&image[i*W+W+j]==255&&result[i*W+W+j]!=result[i*W+j]){
				judge = false;
				cerr << i <<" "<< j << " " << i+1 <<" "<< j << endl;
			}
			if(image[i*W+j]==255&&image[i*W+1+j]==255&&result[i*W+1+j]!=result[i*W+j]){
				judge = false;
				cerr << i <<" "<< j << " " << i <<" "<< j+1 << endl;
			}
		}
	}
	if(judge)
		cerr << "result is correct" << endl;
	return 0;
}
