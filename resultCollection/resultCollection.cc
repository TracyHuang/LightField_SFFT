#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string>
#include <cstring>
#include <errno.h>
#include <vector>
#include <stdlib.h>
#include <complex>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fftw3.h>
#include <sstream>
#include <cv.h>
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
#include <climits>

using namespace std;

typedef std::complex<double> Complex;
typedef Complex * ComplexPtr;
const Complex I(0.0, 1.0);

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


void readData(FILE * file, Complex *** data, int num_slices) {    // each output file has a data structure of 3D

	if (file == NULL) {
		cerr << "The output file does not exist" << endl;
		exit(-1);
	}

	int index_length, num_data_dims;
	//fread(&num_slices, sizeof(int), 1, file);
	fread(&index_length, sizeof(int), 1, file);
	fread(&num_data_dims, sizeof(int), 1, file);
	//cout << "num_slices " << num_slices << " index length " << index_length << " num_data_dims " << num_data_dims << endl;
	if (num_data_dims != 3) {
		cerr << "The data dims from the output file is not 3" << endl;
		exit(-1);
	}
	
	int * data_dims = (int *) malloc(sizeof(int) * num_data_dims);
	for (int i = 0; i < num_data_dims; i ++) {
		fread(data_dims + i, sizeof(int), 1, file);
	}

	//cout << "data_dims is " << data_dims[0] << "," << data_dims[1] << "," << data_dims[2] << endl;
    /*
	data = (Complex ***) malloc(sizeof(Complex **) * data_dims[0]);
	for (int i = 0; i < data_dims[0]; i ++) {
		data[i] = (Complex **) malloc(sizeof(ComplexPtr) * data_dims[1]);
		for (int j = 0; j < data_dims[1]; j ++) {
			data[i][j] = (ComplexPtr) malloc(sizeof(Complex) * data_dims[2]);
		}
	}
	*/
	double ** real_d = (double **) malloc(sizeof(double *) * data_dims[0]);
	for (int i = 0; i < data_dims[0]; i ++) {
		real_d[i] = (double *) malloc(sizeof(double) * data_dims[1]);
	}

	double ** imag_d = (double **) malloc(sizeof(double *) * data_dims[0]);
	for (int i = 0; i < data_dims[0]; i ++) {
		imag_d[i] = (double *) malloc(sizeof(double) * data_dims[1]);
	}

	// sanity check
	if (num_slices != data_dims[2]) {
		cerr << "number of slices in current output file does not accord with data_dims recorded in it" << endl;
		exit(-1);
	}	

	for (int i = 0; i < num_slices; i ++) {

		// original we deal with something here as index_length != 0 
		// but now as we assume that index_length is always 0
		// so we just ignore all things here
		// for reference, please check projection_directory/matlab/ReadData.m  line 38

		// read in real values
		for (int j = 0; j < data_dims[0]; j ++) {
			for (int k = 0; k < data_dims[1]; k ++) {
				fread( &(real_d[j][k]), sizeof(double), 1, file);
			}	
		}

		// read in imaginary values
		for (int j = 0; j < data_dims[0]; j ++) {
			for (int k = 0; k < data_dims[1]; k ++) {
				fread( &(imag_d[j][k]), sizeof(double), 1, file);
			}
		}	

		// fill in values into current slices
		for (int j = 0; j < data_dims[0]; j ++) {
			for (int k = 0; k < data_dims[1]; k ++) {
				data[j][k][i] = real_d[j][k] + I * imag_d[j][k];
			}
		}

		for (int j = 0; j < data_dims[0]; j ++) {
			for (int k = 0; k < data_dims[1]; k ++) {
				data[j][k][i] = data[j][k][i] + 1.0 - 1.0;
			}
		}

	}



}

void findIndex(int slice_index, int n3, int * row_index, int * column_index) {

	*column_index = slice_index / n3;   // yes we want to use interger division here
	*row_index = slice_index % n3;
}

void scaleUp(ComplexPtr x, int n, Complex a) {
	ComplexPtr x_ptr = x;
	//cout << "scale up by " << real(a) << "     " << imag(a) << endl;
	for (int i = 0; i < n; i ++, x_ptr ++) {
		*x_ptr = a * (*x_ptr);
	} 
}


int32_t getBlockNum(int x, int y, int n) {
	if( n % 2 == 0) {
		if (y >= n / 2) { // block 2 or 4
			if(x >= n/2) {
			return 4;
			}
			else {
			return 2;
			}
		}
		else {   // block 1 or 3
			if (x >= n/2) {
			return 3;
			}	
			else {
				return 1;
			}
		}
	} 
	else if ( n % 2 == 1) {
		if ( y >= (n + 1) / 2) {
			if ( x >= (n + 1) / 2) { // block 4 or 2
			return 4;
			}
			else {
			return 2;
			}
		}
		else {
			if( x >= (n + 1) / 2) { // block 3 or 1
			return 3;
			}
			else {
			return 1;
			}
		}
	}

	return -1;
}    



void getCorrespondPos(int * old_x, int * old_y, int * new_x, int * new_y, int n) {

	//first decide which block it is
	int block_num = getBlockNum(*old_x, *old_y, n);
	if(block_num == -1) {
		fprintf(stderr, "Fatal error, block_num is -1\n");
	}


	if (n % 2 == 0) {
		int off_x, off_y;
		switch(block_num) {
			case 1:
				off_x = *old_x;
				off_y = *old_y;
				if (off_x == 0 && off_y == 0) {    // must calc, no corresponding point
				*new_x = *old_x;
				*new_y = *old_y;
				}
				else if (off_x == 0) {             // block 2
				*new_x = 0;
				*new_y = n - off_y;
				}
				else if (off_y == 0) {             // block 3
				*new_x = n - off_x;
				*new_y = 0;
					}
				else {                             // block 4
				*new_x = n - off_x;
				*new_y = n - off_y;
				}
			break;
				case 2:
				off_x = *old_x;
				off_y = *old_y - n / 2;
					if (off_x == 0 && off_y == 0) {     // must calc, no corresponding point
				*new_x = *old_x;    
				*new_y = *old_y;
				}
				else if (off_x == 0) {             // block 1
				*new_x = 0;
				*new_y = n / 2 - off_y;
				}
				else if (off_y == 0) {             // block 4
				*new_x = n / 2;
				*new_y = n - off_x;
				}
				else {                             // block 3
				*new_x = n - off_x;
				*new_y = n / 2 - off_y;
				}
			break;
			case 3:
				off_x = *old_x - n / 2;
				off_y = *old_y;
				if (off_x == 0 && off_y == 0) {     // must calc, no corresponding point
				*new_x = *old_x;
				*new_y = *old_y;
				}
				else if (off_x == 0) {              // block 4
				*new_x = n / 2;
				*new_y = n - off_y;
				}
				else if (off_y == 0) {              // block 1
				*new_x = n / 2 - off_x;
				*new_y = 0;
				}
				else {                              // block 2
				*new_x = n / 2 - off_x;
				*new_y = n - off_y;
				}
				break;
			case 4:
				off_x = *old_x - n / 2;
				off_y = *old_y - n / 2;
				if (off_x == 0 && off_y == 0) {     // must calc, no corresponding point
				*new_x = *old_x;
				*new_y = *old_y;
				}
				else if (off_x == 0) {              // block 3
				*new_x = n / 2;
				*new_y = n / 2 - off_y;
				}
				else if (off_y == 0) {              // block 2
				*new_x = n / 2 - off_x;
				*new_y = n / 2;
				}
				else {
				*new_x = n / 2 - off_x;         // block 1
				*new_y = n / 2 - off_y;
				}
			break;
				default: 
			break;
		}
	}
	else if (n % 2 == 1) {
		int off_x, off_y;
		switch (block_num) {
			case 1:
			off_x = *old_x;
			off_y = *old_y;
			if (off_x == 0 && off_y == 0) {     // must calc, no corresponding point
			*new_x = *old_x;
			*new_y = *old_y;
			}
			else if (off_x == 0) {              // block 2
			*new_x = 0;
			*new_y = n - off_y;
			}
			else if (off_y == 0) {              // block 3
			*new_x = n - off_x;
			*new_y = 0;
			}
			else {                              // block 4
			*new_x = n - off_x;
			*new_y = n - off_y;
			}
			break;
			case 2:
			off_x = *old_x;
			off_y = *old_y - (n + 1) / 2;
			if(off_x == 0) {                    // block 1
			*new_x = 0;
			*new_y = (n - 1) / 2 - off_y;
			}
			else {                              // block 3
			*new_x = n - off_x;
			*new_y = (n - 1) / 2 - off_y;
			}
			break;
			case 3:
			off_x = *old_x - (n + 1) / 2;
			off_y = *old_y;
			if (off_y == 0) {                   // block 1
			*new_x = (n - 1) / 2 - off_x;
			*new_y = 0;
			}
			else {                              // block 2
			*new_x = (n - 1) / 2 - off_x;
			*new_y = n - off_y;
			}
			break;
			case 4:                             // block 1
			off_x = *old_x - (n + 1) / 2;
			off_y = *old_y - (n + 1) / 2;
			*new_x = (n - 1) / 2 - off_x;
			*new_y = (n - 1) / 2 - off_y;
			break;
			default:
			break;
		}
	}
	return;

}



void copyData(bool **** used, Complex **** data, int n1, int n2, int n3, int n4) {

	// go through each slice 
	// if one of the slice is not calculated
	// we first find the corresponding slice it corresponds to 
	// and then do corresponding again in each slice
	// luckily they are using the same matching pattern, thank god 

	int new_x, new_y, old_x, old_y;
	int new_pos_x, new_pos_y, old_pos_x, old_pos_y;
	for (int i = 0 ; i < n3; i ++) {
		for (int j = 0; j < n4; j ++) {
			if (!used[0][0][i][j]) {
				old_x = i;
				old_y = j;
				getCorrespondPos(&old_x, &old_y, &new_x, &new_y, n3);
				//cout << "old : " <<	old_x << ", " << old_y << "  new : " << new_x << " , " << new_y << endl;
				//cout << "used" << used[0][0][new_x][new_y] << endl;
				if (!used[0][0][new_x][new_y]) {
					cerr << "Error the slice copied from is not calculated" << endl;
					exit(-1);
				}	
				// now we need to fill in the correspond positions in two different slices
				for (int k = 0; k < n1; k ++) {
					for (int l = 0; l < n2; l ++) {
						old_pos_x = k;
						old_pos_y = l;
						getCorrespondPos(&old_pos_x, &old_pos_y, &new_pos_x, &new_pos_y, n1);
						//cout << "old pos " << old_pos_x << ", " << old_pos_y << " new_pos " << new_pos_x << ", " << new_pos_y << endl;
						data[k][l][i][j] = conj(data[new_pos_x][new_pos_y][new_x][new_y]);
					}
				}
			}
		}
	}

	cout << "copy finished\n";
	


}



void uvPosCorrespond(int old_x, int old_y, int * new_x, int * new_y) {
	int temp_x = old_x % 2 == 1 ? old_x - 1 : old_x;
	int temp_y = old_y % 2 == 1 ? old_y - 1 : old_y;
	*new_x = temp_x / 2;
	*new_y = temp_y / 2;
}








int main(int argc, char *argv[]) {

	if (argc != 7) {
		cerr << "Error usage, the correct way should be ./collectResult <directory_name> <n1> <n2> <n3> <n4> <channel>" << endl;
		exit(-1);
	}

	//set up parammeters
	char directory[MAX_DIRECTORY_NAME_LENGTH];
	int n1, n2, n3, n4, channel;

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
	channel = atoi(argv[6]);
		
	if (channel != 0) {       // for u/v channel we operate on lower resolution
		n3 /= 2;
		n4 /= 2;
		cout << "n3 is rescaled to " << n3 << " n4 is rescaled to " << n4 << endl;
	}
	
	int total_slices = n3 * n4;
	int copy_slices = n3 * n4;

	// this variable is to show whether we can use copy or not
	bool copy_possible = false;
	if ((n1 == n2) && (n3 == n4)) {
		copy_possible = true;
	}

	vector<string> files = vector<string>();
	getFileList(dir, files);

	Complex **** data = (Complex ****) malloc(sizeof(Complex ***) * n1);
	for (int i = 0; i < n1; i ++) {
		data[i] = (Complex ***) malloc(sizeof(Complex **) * n2);
		for (int j = 0; j < n2; j ++) {
			data[i][j] = (Complex **) malloc(sizeof(ComplexPtr) * n3);
			for (int k = 0; k < n3; k ++) {
				data[i][j][k] = (ComplexPtr) malloc(sizeof(Complex) * n4);
			}
		}
	}

	bool **** used = (bool ****) malloc(sizeof(bool ***) * n1);
	for (int i = 0; i < n1; i ++) {
		used[i] = (bool ***) malloc(sizeof(bool **) * n2);
		for (int j = 0; j < n2; j ++) {
			used[i][j] = (bool **) malloc(sizeof(bool *) * n3);
			for (int k = 0; k < n3; k ++) {
				used[i][j][k] = (bool *) malloc(sizeof(bool) * n4);
				for (int l = 0; l < n4; l ++) {
					used[i][j][k][l] = 0;
				}
			}
		}
	}	


	Complex *** buffer = NULL;
	FILE * file;
	int num_slices;
	int slice_index = 0;
	int row_index = 0;
	int column_index = 0;

	for (int i = 0; i < files.size(); i ++) {
		
		file = fopen(files[i].c_str(), "rb");
		fread(&num_slices, sizeof(int), 1, file);
		copy_slices -= num_slices;		

		// allocate space for buffer
		buffer = (Complex ***) malloc(sizeof(Complex **) * n1);
		for (int j = 0; j < n1; j ++) {
		buffer[j] = (Complex **) malloc(sizeof(ComplexPtr) * n2);
			for (int k = 0; k < n2; k ++) {
				buffer[j][k] = (ComplexPtr) malloc(sizeof(Complex) * num_slices);
			}
		}

		
		// read data from output file
		readData(file, buffer, num_slices);
			
		for (int j = 0; j < num_slices; j ++) {
			findIndex(slice_index, n3, &row_index, &column_index);
			//cout << "row_index " << row_index << " column_index " << column_index <<  endl;
			for (int k = 0; k < n1; k ++) {
				for (int l = 0; l < n2; l ++) {
					data[k][l][row_index][column_index] = buffer[k][l][j];
					used[k][l][row_index][column_index] = 1;
				}
			}
			slice_index ++;
		}

		//cout << "file index is " << i << ", slice index is " << slice_index << endl;
		
		fclose(file);

		// free buffer
		for (int j = 0; j < n1; j ++) {
			for (int k = 0; k < n2; k ++) {
				if (buffer[j][k]) {
					free(buffer[j][k]);
					buffer[j][k] = NULL;
				}
			}
			if (buffer[j]) {
				free(buffer[j]);
				buffer[j] = NULL;
			}
		}
		if(buffer) {
			free(buffer);
			buffer = NULL;
		}
	
	}

	// now we need to decide whether we need to copy or not
	bool need_copy = false;
	for (int i = 0; i < n1; i ++) {
		for (int j = 0; j < n2; j ++) {
			for (int k = 0; k < n3; k ++) {
				for (int l = 0; l < n4; l ++) {
					if (!used[i][j][k][l]) {
						need_copy = true;
					}
				}
			}
		}
	}


	if (need_copy) {

		if (!copy_possible) {
			cerr << "Error we need to copy but dimension does not match\n";
			exit(-1);
		}

		cout << " we need to copy-------" << endl;
		copyData(used, data, n1, n2, n3, n4);

	}	





	// after collect output and copy we can process data to finish the task
	fftw_plan plan;
	ComplexPtr fftw_input;
	fftw_input = (ComplexPtr)fftw_malloc(sizeof(Complex) * n1 * n2 * n3 * n4);
	ComplexPtr fftw_output; // used to store the output of the fft
	fftw_output = (ComplexPtr)fftw_malloc(sizeof(Complex) *n1 * n2 * n3 * n4);
	
	int n[4] = {n1, n2, n3, n4};
	
	plan = fftw_plan_dft(4, n, ((fftw_complex *) (void *) fftw_input), ((fftw_complex *) (void *) fftw_output), FFTW_BACKWARD, FFTW_MEASURE); 
	
	// fill in input
	for (int i = 0; i < n1; i ++) {
		for (int j = 0; j < n2; j ++) {
			for (int k = 0; k < n3; k ++) {
				for (int l = 0; l < n4; l ++) {
					fftw_input[i * n2 * n3 * n4 + j * n3 * n4 + k * n4 + l] = data[i][j][k][l];
				}
			}
		}
	}

	fftw_execute(plan);	
	scaleUp(fftw_output, n1 * n2 * n3 * n4, 1 / sqrt((double)n1 * n2 * n3 * n4));

	// create space for output
	// we are using double here because we don't want to face the problem of 
	double **** image = (double ****) malloc(sizeof(double ***) * n1);
	for (int i = 0; i < n1; i ++) {
		image[i] = (double ***) malloc(sizeof(double **) * n2);
		for (int j = 0; j < n2; j ++) {
			image[i][j] = (double **) malloc(sizeof(double *) * n3);
			for (int k = 0; k < n3; k ++) {
				image[i][j][k] = (double *) malloc(sizeof(double) * n4);
			}
		}
	}

	
	for (int i = 0; i < n1; i ++) {
		for (int j = 0; j < n2; j ++) {
			for (int k = 0; k < n3; k ++) {
				for (int l = 0; l < n4; l ++) {
					//image[i][j][k][l] = (int) ( real(fftw_output[i * n2 * n3 *n4 + j * n3 * n4 + k * n4 + l]) > INT_MAX * 1.0 ? INT_MAX : real(fftw_output[i * n2 * n3 * n4 + j * n3 * n4 + k * n4 + l])) ;
					image[i][j][k][l] = real(fftw_output[i * n2 * n3 * n4 + j * n3 * n4 + k * n4 + l]);
				}
			}
		}
	}

	// write output into binary files
	string outFilePrefix;
	string outFileName;
	FILE * out_file;
	for (int i = 0; i < n1; i ++) {
		for (int j = 0; j < n2; j ++) {
			
			// first get to know the file name
			if (channel == 0) {
				outFilePrefix = "y_";
			}
			else if (channel == 1) {
				outFilePrefix = "u_";
			}
			else if (channel == 2) {
				outFilePrefix = "v_";
			}
			if (i < 10) {
				if (j < 10) {
					outFileName = outFilePrefix + "0" + patch::to_string(i) + "_0" + patch::to_string(j) + ".dat";
				}
				else {
					outFileName = outFilePrefix + "0" + patch::to_string(i) + "_" + patch::to_string(j) + ".dat";
				}
			}
			else {
				if (j < 10) {
					outFileName = outFilePrefix + patch::to_string(i) + "_0" + patch::to_string(j) + ".dat";
				}
				else {
					outFileName = outFilePrefix + patch::to_string(i) + "_" + patch::to_string(j) + ".dat";
				}
			}

			// create output file
			out_file = fopen(outFileName.c_str(), "wb");

			if (channel == 0) {
				for (int k = 0; k < n3; k ++) {
					for (int l = 0; l < n4; l ++) {
						fwrite(&(image[i][j][k][l]), sizeof(double), 1, out_file);
					}
				}
			}
			
			else {
				// here because we only have data for the lower resolution one
				// so we need to use one entry several times to fill in data file
				int actual_k, actual_l;
				for (int k = 0; k < n3 * 2; k ++) {
					for (int l = 0; l < n4 * 2; l ++) {
						uvPosCorrespond(k, l, &actual_k, &actual_l);
						fwrite(&(image[i][j][actual_k][actual_l]), sizeof(double), 1, out_file);
					}
				}
			}
			
			fclose(out_file);

		}
	}




	/*
	// transfer image array into pictures
	IplImage * out_image = NULL;
	IplImage * resized_image = NULL;
	double * image_data;
	//image_data = (unsigned char *) malloc(sizeof(unsigned char) * n3 * n4);
	image_data = (double * ) malloc(sizeof(double) * n3 * n4);
	string outFilePrefix = "out_";
	string outFileName;
	CvSize image_size;
	image_size.width = n4;
	image_size.height = n3;
	for (int i = 0; i < n1; i ++) {
		for (int j = 0; j < n2; j ++) {
			
			out_image = cvCreateImageHeader( image_size, IPL_DEPTH_64F, 1); // channel is one because we only want a one dimensional graph
			//image_data = out_image->imageData;
			// copy image data into img
			for (int k = 0; k < n3; k ++) {
				for (int l = 0; l < n4; l ++) {
					
					//if (image[i][j][k][l] > UCHAR_MAX) {
					//	image[i][j][k][l] = UCHAR_MAX;
					//}
					
					image_data[k * n4 + l] = image[i][j][k][l]; //(unsigned char) (image[i][j][k][l] < 0 ? 0 : image[i][j][k][l])
				}
			}
			cvSetData(out_image, image_data, sizeof(double) * n4); 
			if (i < 10) {
				if (j < 10) {
					outFileName = outFilePrefix + "0" + patch::to_string(i) + "_0" + patch::to_string(j) + ".png";
				}
				else {
					outFileName = outFilePrefix + "0" + patch::to_string(i) + "_" + patch::to_string(j) + ".png";
				}
			}
			else {
				if (j < 10) {
					outFileName = outFilePrefix + patch::to_string(i) + "_0" + patch::to_string(j) + ".png";
				}
				else {
					outFileName = outFilePrefix + patch::to_string(i) + "_" + patch::to_string(j) + ".png";
				}
			}
			if (channel == 0) {
				if (!cvSaveImage(outFileName.c_str(), out_image)) cerr << "the image " << i << "_" << j << "is not created successfully" << endl;
				cvReleaseImageHeader(&out_image);
			}
			else {                     // if not y channel, we need to resize the image
				resized_image = cvCreateImageHeader(cvSize(n4 * 2, n3 * 2), out_image->depth, out_image->nChannels);			
				cvResize(out_image, resized_image, CV_INTER_LINEAR);
				if (!cvSaveImage(outFileName.c_str(), resized_image)) cerr << "the image " << i << "_" << j << "is not created successfully" << endl;
				cvReleaseImageHeader(&out_image);
				cvReleaseImageHeader(&resized_image);
			}
		}
	}
	*/

	cout << "total number of slices is " << total_slices << "   number of slices that are copied is " << copy_slices << endl;

	return 0;


}
      

