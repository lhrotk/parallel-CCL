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

class CCL {
public:
	vector<int> ccl(vector<int>& image, int W, int degree_of_connectivity, int th);
};

vector<int> CCL::ccl(vector<int>& image, int W, int degree_of_connectivity, int th)
{
	int* D = static_cast<int*>(&image[0]);
	int N = image.size();
	int* Label = new int[N];
	int* Reference = new int[N];
	th = 0;

	double start = get_time();
	init_CCL(Label, Reference, N);

	th = scanning(D, Label, Reference, N, W, th);
	labeling(Label, Reference, N);
	
	double end = get_time();
	cerr << "Time: " << end - start << endl;

	vector<int> result(Label, Label + N);

	cerr << th <<endl;

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

	vector<int> image;
	int W, degree_of_connectivity, threshold;//W: width
	read_data(argv[1], image, W, degree_of_connectivity, threshold);

	CCL ccl;

	
	vector<int> result(ccl.ccl(image, W, degree_of_connectivity, threshold));

	cout << image.size() << endl; /// number of pixels
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
