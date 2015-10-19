#include <iostream>
#include <stdio.h>
#include <string>
#include <cstring>
#include <errno.h>
#include <vector>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
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
using namespace cv;

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

unsigned char double2uchar(double a) {

	//cout << a << '\t';
	if (a < 0) {
		//cout << 0 << endl;
		return 0;
	}
	if (a > UCHAR_MAX * 1.0) {
		//cout << UCHAR_MAX;
		return UCHAR_MAX;
	}
	//cout << (int) ((unsigned char) a) << endl;
	return (unsigned char) a;
}


int main(int argc, char *argv[]) {

	if (argc != 8) {
		cerr << "Error usage, the correct way should be ./rgbReconstruct <Y_directory> <U_directory> <V_directory> <n1> >n2> <n3> <n4>" << endl;
		exit(-1);
	}
	
	char y_dir[MAX_DIRECTORY_NAME_LENGTH];
	char u_dir[MAX_DIRECTORY_NAME_LENGTH];
	char v_dir[MAX_DIRECTORY_NAME_LENGTH];
	int n1, n2, n3, n4;

	strcpy(y_dir, argv[1]);
	strcpy(u_dir, argv[2]);
	strcpy(v_dir, argv[3]);

	n1 = atoi(argv[4]);
	n2 = atoi(argv[5]);
	n3 = atoi(argv[6]);
	n4 = atoi(argv[7]);

	vector<string> y_files = vector<string>();
	vector<string> u_files = vector<string>();
	vector<string> v_files = vector<string>();

	getFileList(y_dir, y_files);
	getFileList(u_dir, u_files);
	getFileList(v_dir, v_files);

	

	// note that the desired output storage order is BGR
	// The conversion formula is as below
	// R = 1.164 * ( Y - 16 ) + 1.596 * ( V - 128 )
	// G = 1.164 * ( Y - 16 ) - 0.813 * ( V - 128 ) - 0.391 * ( U - 128 )
	// B = 1.164 * ( Y - 16 ) +                       2.018 * ( U - 128 )
	
	//Mat y_image, u_image, v_image;
	// we need 3 2D arrays to store the information 
	double ** y = (double **) malloc(sizeof(double *) * n3);
	for (int i = 0; i < n3; i ++) {
		y[i] = (double *) malloc(sizeof(double) * n4);
		for (int j = 0; j < n4; j ++) {
			y[i][j] = 0.0;
		}
	} 

	double ** u = (double **) malloc(sizeof(double *) * n3);
	for (int i = 0; i < n3; i ++) {
		u[i] = (double *) malloc(sizeof(double) * n4);
		for (int j = 0; j < n4; j ++) {
			u[i][j] = 0.0;
		}
	}

	double ** v = (double **) malloc(sizeof(double *) * n3);
	for (int i = 0; i < n3; i ++) {
		v[i] = (double *) malloc(sizeof(double) * n4);
		for (int j = 0; j < n4; j ++) {
			v[i][j] = 0.0;
		}
	}
	Mat img = Mat::zeros(n3, n4, CV_8UC3);
	Vec3b color;
	string outFilePrefix = "result_";
	string outFileName;
	FILE * input_file;
	int x_index, y_index;
	for (int i = 0; i < y_files.size(); i ++) {

		
		//y_image = imread(y_files[i].c_str(), 0); // 0 means read it into grayscale
		//u_image = imread(u_files[i].c_str(), 0);
		//v_image = imread(v_files[i].c_str(), 0);
	
		// read files
		// read y
		input_file = fopen(y_files[i].c_str(), "rb");	
		if (input_file == NULL) {
			cerr << "File error " << y_files[i] << " is empty" << endl;	
			exit(-1);
		}
		for (int j = 0;  j < n3; j ++) {
			for (int k = 0; k < n4; k ++) {
				fread(&(y[j][k]), sizeof(double), 1, input_file);
			}
		}
		fclose(input_file);

		// read u
		input_file = fopen(u_files[i].c_str(), "rb");
		if (input_file == NULL) {
			cerr << "File error " << u_files[i] << " is empty" << endl;	
			exit(-1);
		}
		for (int j = 0;  j < n3; j ++) {
			for (int k = 0; k < n4; k ++) {
				fread(&(u[j][k]), sizeof(double), 1, input_file);
			}
		}
		fclose(input_file);

		// read v
		input_file = fopen(v_files[i].c_str(), "rb");
		if (input_file == NULL) {
			cerr << "File error " << v_files[i] << " is empty" << endl;	
			exit(-1);
		}
		for (int j = 0;  j < n3; j ++) {
			for (int k = 0; k < n4; k ++) {
				fread(&(v[j][k]), sizeof(double), 1, input_file);
			}
		}
		fclose(input_file);

		//int a;
		//cin >> a;


		for (int j = 0; j < n3; j ++) {       // rows
			for (int k = 0; k < n4; k ++) {	  // columns
				//y = y_image.at<double>(j, k);
				//u = u_image.at<double>(j, k);
				//v = v_image.at<double>(j, k);				

				color = img.at<Vec3b>(j, k);
				// storage is in order of BGR 
				//color.val[0] = double2uchar(y[j][k] + 2.03211 * u[j][k]);
				//color.val[1] = double2uchar(y[j][k] - 0.39465 * u[j][k] - 0.58060 * v[j][k]);
				//color.val[2] = double2uchar(y[j][k] + 1.13983 * v[j][k]);
				
				color.val[0] = double2uchar(1.164 * (y[j][k] - 16.0)                       + 2.018 * (u[j][k] - 128.0));
				color.val[1] = double2uchar(1.164 * (y[j][k] - 16.0) - 0.813 * (v[j][k] - 128.0) - 0.391 * (u[j][k] - 128.0));
				color.val[2] = double2uchar(1.164 * (y[j][k] - 16.0) + 1.596 * (v[j][k] - 128.0));
				
				// store the color back
				img.at<Vec3b>(j, k) = color;
                        
			}
		}
	
		x_index = i / n2;
		y_index = i % n2;
		
		if (x_index < 10) {
			if (y_index < 10) {
                outFileName = outFilePrefix + "0" + patch::to_string(x_index) + "_0" + patch::to_string(y_index) + ".png";
            }
            else {
                outFileName = outFilePrefix + "0" + patch::to_string(x_index) + "_" + patch::to_string(y_index) + ".png";
            }
        }
        else {
            if (y_index < 10) {
                outFileName = outFilePrefix + patch::to_string(x_index) + "_0" + patch::to_string(y_index) + ".png";
            }
            else {
                outFileName = outFilePrefix + patch::to_string(x_index) + "_" + patch::to_string(y_index) + ".png";
            }
        }
		
	
		imwrite(outFileName, img);

	
	}

	cout << " finish generating rgb images" << endl;
	return 0;

}

 
