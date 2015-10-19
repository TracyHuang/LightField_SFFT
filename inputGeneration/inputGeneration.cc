#include <typeinfo>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <errno.h>
#include <cstring>
#include <cmath>
#include <complex>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fftw3.h>
//#include <Magick++.h>
#include <sstream>
#include <cv.h>
#include <climits>
#include <highgui.h>
#include "opencv2/core/core_c.h"
#include "opencv2/core/core.hpp"
#include "opencv2/flann/miniflann.hpp"
#include "opencv2/imgproc/imgproc_c.h"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/video/video.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/objdetect/objdetect.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/ml/ml.hpp"
#include "opencv2/highgui/highgui_c.h"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/contrib/contrib.hpp"

using namespace std;
//using namespace cv;

typedef std::complex<double> Complex;
typedef Complex * ComplexPtr;


#define MAX_DIRECTORY_NAME_LENGTH 100

namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}



int isFile(const char *path)
{
    struct stat buf;
    stat(path, &buf);
    return S_ISREG(buf.st_mode);
}


int getFileList(string dir, vector<string> &files)
{
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }

	string filePath;
    while ((dirp = readdir(dp)) != NULL) {
		filePath = dir + string(dirp->d_name);
		if (isFile(filePath.c_str())) {
			files.push_back(filePath);
		}
    }
	sort( files.begin(), files.end() );
    closedir(dp);
    return 0;
}

template <typename T, typename M>
void oneDimensionToTwoDimension(T * ptr_1d, M ** ptr_2d, int n1, int n2) {
	for (int i = 0; i < n1; i ++) {
		for (int j = 0; j < n2; j ++) {
			ptr_2d[i][j] = ptr_1d[i * n2 + j];
		}
	}
}

template <typename T, typename M>
void twoDimensionToOneDimension(T ** ptr_2d, M * ptr_1d, int n1, int n2) {
	for (int i = 0; i < n1; i ++) {
		for (int j = 0; j < n2; j ++) {
			ptr_1d[i * n2 + j] = ptr_2d[i][j];
		}
	}
}

void scaleUp(ComplexPtr x, int n, Complex a) {
	ComplexPtr x_ptr = x;
	//cout << "scale up by " << real(a) << "     " << imag(a) << endl;
	for (int i = 0; i < n; i ++, x_ptr ++) {
		*x_ptr = a * (*x_ptr);
	} 
}

void getIndex(int start_index, int * row_index, int * column_index, int n3, int n4) {
	*column_index = (int) floor(((double) start_index) / ((double) n3));
	*row_index = start_index % n3;
	if (*row_index >= n3 || *column_index >= n4) {
		cerr << "Fatal error, the row_index/column_index exceeds the limit: n3 " << n3 << "n4 " << n4 << "row_index " << *row_index << "column_index " << *column_index << endl;
		exit(-1);  
	}
}

void getNextIndex(int * row_index, int * column_index, int n3) {
	(*row_index) ++;
	if (*row_index >= n3) {
		(*column_index) ++;
		*row_index = 0;
	}
}

void writeData(FILE * file, Complex **** input, int num_slices_per_split, int index, int n1, int n2, int n3, int n4) {
	
	int start_index = num_slices_per_split * index;
	int row_index, column_index;
	getIndex(start_index, &row_index, &column_index, n3, n4);
	//cout << row_index << '\t' << column_index << endl;	
	double temp = 0.0;
	for (int i = 0; i < num_slices_per_split; i ++) {
		//cout << "row " << row_index << " column " << column_index << endl;
		// write real part into file
		for (int j = 0; j < n1; j ++) {
			for (int k = 0; k < n2; k ++) {
				temp = real(input[j][k][row_index][column_index]);
				fwrite(&temp, sizeof(double), 1, file);	
			}
		}

		// write imaginary part into file
		for (int j = 0; j < n1; j ++) {
			for (int k = 0; k < n2; k ++) {
				temp = imag(input[j][k][row_index][column_index]);
				fwrite(&temp, sizeof(double) , 1, file);
			}
		}
		
		getNextIndex(&row_index, &column_index, n3);
		//cout << row_index << '\t' << column_index << endl;
		//cin >> a;
	}
	//cout << num_slices_per_split << endl;
}





int main(int argc, char* argv[]) {


//InitializeMagick(*argv);


if (argc != 9) {
	cerr << "Error usage, the correct way should be ./getInput <directory_name> <n1> <n2> <n3> <n4> <num_splits> <if_copy> <channel>" << endl;
	exit(-1);
}

//set up parammeters
char directory[MAX_DIRECTORY_NAME_LENGTH];
int n1, n2, n3, n4, num_splits;
// if_copy means whether we will generate full input or we will only use the first half of data and copy their results to the other half
// for more information regarding copying, please check document
// channel means which channel of data will be generated, for current setting the channel 0 is Y, channel 1 is U, and channel 2 is V
// and by defalut channel U and V is rescaled by half to speed up calculation
bool if_copy;
int channel;

strcpy(directory, argv[1]);
int length = strlen(directory);
if (directory[length - 1] != '/') {
	directory[length + 1] = directory[length]; // copy \0
	directory[length] = '/';
}
string dir(directory);

n1 = atoi(argv[2]);
n2 = atoi(argv[3]);
n3 = atoi(argv[4]);
n4 = atoi(argv[5]);
num_splits = atoi(argv[6]);
if_copy = atoi(argv[7]);
channel = atoi(argv[8]);  


if (channel != 0) {
	if ((n3 % 2 == 1) || (n4 % 2 == 1)) {
		cerr << "Error cannot rescale the image to half as n3 or n4 is odd" << endl;
		exit(-1);
	}
	n3 /= 2;
	n4 /= 2;
	cout << "n3 rescaled to " << n3 << " and n4 rescaled to " << n4 << " to speed up calculation" << endl;
}

vector<string> files = vector<string>();
getFileList(dir, files);



// variables declaration
double ** grey_image;
grey_image = (double **) malloc(sizeof(double*) * n3);
for (int i = 0; i < n3; i ++) {
	grey_image[i] = (double *) malloc(sizeof(double) * n4);
}

double **** data;
data = (double ****) malloc(sizeof(double ***) * n1);
for (int i = 0; i < n1; i ++) {
	data[i] = (double ***) malloc(sizeof(double **) * n2);
	for (int j = 0; j < n2; j ++) {
		data[i][j] = (double **) malloc(sizeof(double *) * n3);
		for (int k = 0; k < n3; k ++) {
			data[i][j][k] = (double *) malloc(sizeof(double) * n4);
		}
	}
}

int row, column;
fftw_plan plan;
ComplexPtr fftw_input;
fftw_input = (ComplexPtr)fftw_malloc(sizeof(Complex) * n3 * n4);
ComplexPtr fftw_output; // used to store the output of the fft
fftw_output = (ComplexPtr)fftw_malloc(sizeof(Complex) * n3 * n4);
bool init_plan = false;


// generate grey scale image and store data
//Image image;
//string resize_pattern = patch::to_string(n4);
//resize_pattern += "x";
//resize_pattern += patch::to_string(n3);
//Geometry newSize = Geometry(resize_pattern);
//PixelPacket * pixels;
//Color color;
//double max = 0.0;
cv::Mat image, resized_image;
cv::Vec3b color;
for (int32_t i = 0; i < files.size(); i ++) {
	image = cv::imread(files[i].c_str(), CV_LOAD_IMAGE_COLOR);
	cv::resize(image, resized_image, cv::Size(n3, n4)); // resize the image according to the channel
	//image.resize(newSize);
	//string out_name = patch::to_string(i) + ".png";
	//image.write(out_name);
	//cout << resize_pattern << endl;
	//cout << "width " << image.columns() << "  height   " << image.rows() << endl;
	//pixels = image.getPixels(0, 0, n4, n3);   // columns rows order 

	//double red, green, blue;	

	//cout << " cureent channel is " << channel << endl;
	for (int j = 0; j < n3; j ++) {
		for (int k = 0; k < n4; k ++) {
			//cout << "j " << j << " k " << k << endl;
			//color = pixels[j * n4 + k];
			color = resized_image.at<cv::Vec3b>(j, k);

			
			// note that opencv uses a BGR channel order
			if (channel == 0) {
				grey_image[j][k] = round(0.257 * (double) color.val[2]  + 0.504 * (double) color.val[1] + 0.098 * (double) color.val[0]  +  16.0);
			}
			else if (channel == 1) {
				grey_image[j][k] = round(-0.148 * (double) color.val[2] - 0.291 * (double) color.val[1] + 0.439 * (double) color.val[0] + 128.0);
			}
			else if (channel == 2) {
				grey_image[j][k] = round(0.439 * (double) color.val[2] - 0.368 * (double) color.val[1] - 0.071 * (double) color.val[0] + 128.0);
			}
		}
	}
	row = floor(i / n2);
	column = i - row * n2; 
    //cout << "row column " << row << "   " << column << endl;	
	for (int j = 0; j < n3; j ++) {
		for (int k = 0; k < n4; k ++) {
			data[row][column][j][k] = grey_image[j][k];
		}
	}
}

//cout << " max is " << max;

//cout << "finished reading files\n";






Complex **** input;
input = (Complex ****) malloc(sizeof(Complex ***) * n1);
for (int i = 0; i < n1; i ++) {
	input[i] = (Complex ***) malloc(sizeof(Complex **) * n2);
	for (int j = 0; j < n2; j ++) {
		input[i][j] = (Complex **) malloc(sizeof(Complex *) * n3);
		for (int k = 0; k < n3; k ++) {
			input[i][j][k] = (Complex *) malloc(sizeof(Complex) * n4);
		}
	}
}

// do 2d fft on data and store it in input
for (int i = 0; i < n1; i ++) {
	for (int j = 0; j < n2; j ++) {
		//2Dto1D(data[i][j], fftw_input, n1, n2);
		if (!init_plan) {
			plan = fftw_plan_dft_2d(n3 , n4, ((fftw_complex *)(void *) fftw_input), ((fftw_complex *)(void *)fftw_output), FFTW_FORWARD, FFTW_MEASURE);
		}
		// fill in input
		twoDimensionToOneDimension(data[i][j], fftw_input, n3, n4);	
		fftw_execute(plan);
		// copy result
		scaleUp(fftw_output, n3 * n4, 1.0 / sqrt(n3 * n4));
		oneDimensionToTwoDimension(fftw_output, input[i][j], n3, n4);
	}
}

/*
out.open("fft_result.txt");
for (int i = 0; i < n3; i ++) {
	for (int j = 0; j < n4; j ++){
		out << "(" << real(input[0][0][i][j]) << "," << imag(input[0][0][i][j]) << ") "; 
	}
	out << '\n';
}
out.close();
*/

// split the input and store them into different files
string output_file_prefix = "x_00";
int num_slices = n3 * n4;
cout << "Total number of slices is " << num_slices << endl;
if (if_copy) {
	num_slices = n3 * ( ceil((double) n4 / 2) + 2 );
}
cout << "Number of slices to be calculated is " << num_slices << endl;
int num_slices_per_split = ceil( (double) num_slices / (double) num_splits);
//cout << "num_slices is " << num_slices << "  num_slices_per_split is " << num_slices_per_split << endl;
string file_name = output_file_prefix;
FILE * file;
int index_length = 0;
int dims = 3;

for (int i = 0; i < num_splits - 1; i ++) {  // the last one is a special case as it needs to deal with remainder problem

	if (i < 10) {
		file_name += "0" + patch::to_string(i) + ".dat";
	}
	else {
		file_name += patch::to_string(i) + ".dat";
	}
	file = fopen(file_name.c_str(), "w");
	
	// write header
	fwrite(&num_slices_per_split, sizeof(int), 1, file);
	fwrite(&index_length, sizeof(int), 1, file);
	fwrite(&dims, sizeof(int), 1, file);
	fwrite(&n1, sizeof(int), 1, file);
	fwrite(&n2, sizeof(int), 1, file);
	fwrite(&num_slices_per_split, sizeof(int), 1, file);	


	writeData(file, input, num_slices_per_split, i, n1, n2, n3, n4);	
	fclose(file);
	file_name = output_file_prefix;

}

// deal with the special last one

int last_index = num_splits - 1;

if (last_index < 10) {
	file_name += "0" + patch::to_string(last_index) + ".dat";
}
else {
	file_name += patch::to_string(last_index) + ".dat";
}

file = fopen(file_name.c_str(), "w");

int num_slices_left = num_slices - num_slices_per_split * (num_splits - 1);


fwrite(&num_slices_left, sizeof(int), 1, file);
fwrite(&index_length, sizeof(int), 1, file);
fwrite(&dims, sizeof(int), 1, file);
fwrite(&n1, sizeof(int), 1, file);
fwrite(&n2, sizeof(int), 1, file);
fwrite(&num_slices_left, sizeof(int), 1, file);

int start_index = num_slices_per_split * last_index;
int row_index, column_index;
getIndex(start_index, &row_index, &column_index, n3, n4);
//cout << row_index << '\t' << column_index << endl;
	
double temp = 0.0;
for (int i = 0; i < num_slices_left; i ++) {
	//cout << "row " << row_index << " column " << column_index << endl;	

	// write real part into file
	for (int j = 0; j < n1; j ++) {
		for (int k = 0; k < n2; k ++) {
			temp = real(input[j][k][row_index][column_index]);
			fwrite(&temp, sizeof(double), 1, file);	
		}
	}

	// write imaginary part into file
	for (int j = 0; j < n1; j ++) {
		for (int k = 0; k < n2; k ++) {
			temp = imag(input[j][k][row_index][column_index]);
			fwrite(&temp, sizeof(double) , 1, file);
		}
	}
	
	getNextIndex(&row_index, &column_index, n3);
	//cout << row_index << '\t' << column_index << endl;
	//cin >> a;
}


	
fclose(file);

//cout << "finished writing to files\n";


return 0;

}

 
