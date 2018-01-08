#include<iostream>
#include<cstdio>
#include<string>
#include<vector>
#include<opencv2/core/core.hpp>
#include<opencv2/highgui/highgui.hpp>
#include<opencv2/imgproc/imgproc.hpp>
#include<opencv2/objdetect/objdetect.hpp>

using namespace cv;
using namespace std;

const char* _windowname = "JPEG Compressed Image";
const char* _dctwindow = "DCT Image";
const char* _filename = "src.jpg";

//Data for quantization matrix
int data[8][8] = {
	{ 16, 11, 10, 16, 24, 40, 51, 61 },
	{ 12, 12, 14, 19, 26, 58, 60, 55 },
	{ 14, 13, 16, 24, 40, 57, 69, 56 },
	{ 14, 17, 22, 29, 51, 87, 80, 62 },
	{ 18, 22, 37, 56, 68, 109, 103, 77 },
	{ 24, 35, 55, 64, 81, 104, 113, 92 },
	{ 49, 64, 78, 87, 103, 121, 120, 101 },
	{ 72, 92, 95, 98, 112, 100, 103, 99 }
};

int main()
{
	Mat image = imread(_filename, CV_LOAD_IMAGE_UNCHANGED);
	if (image.empty()) {
		cout << "Can't load the image from " << _filename << endl;
		return -1;
	}

	int x = 0;
	int y = 0;

	namedWindow(_windowname, CV_WINDOW_AUTOSIZE);
	moveWindow(_windowname, x, y);
	imshow(_windowname, image);

	x += 300;
	y += 300;

	int height = image.size().height;
	int width = image.size().width;
	Mat dctImage = image.clone();

	//Convert the 2D array data to image matrix
	Mat quant = Mat(8, 8, CV_8UC1, &data);

	//Discrete Cosine Tranform
	for (int i = 0; i < height; i += 8) {
		for (int j = 0; j < width; j += 8) {
			Mat block = dctImage(Rect(j, i, 8, 8));
			vector<Mat> planes;
			split(block, planes);
			vector<Mat> outplanes(planes.size());
			for (int k = 0; k < planes.size(); k++) {
				planes[k].convertTo(planes[k], CV_32FC1);
				dct(planes[k], outplanes[k]);
				outplanes[k].convertTo(outplanes[k], CV_8UC1);
			}

			//Quantization
			cvtColor(block, block, CV_BGR2GRAY);
			divide(block, quant, block);
			multiply(block, quant, block);
			merge(outplanes, block);
		}
	}

	//Inverse Discrete Cosine Transform
	/*for(int i=0; i < height; i+=8) {
	for(int j=0; j < width; j+=8) {
	Mat block = dctImage(Rect(j,i,8,8));
	vector<Mat> planes;
	split(block,planes);
	vector<Mat> outplanes(planes.size());
	for(int k=0; k < planes.size(); k++) {
	planes[k].convertTo(planes[k],CV_32FC1);
	idct(planes[k],outplanes[k]);
	outplanes[k].convertTo(outplanes[k],CV_8UC1);
	}
	merge(outplanes,block);
	}
	}*/

	namedWindow(_dctwindow, CV_WINDOW_AUTOSIZE);
	moveWindow(_dctwindow, x, y);
	imshow(_dctwindow, dctImage);
	imwrite("dst.jpg", dctImage);

	waitKey(0);

	destroyWindow(_windowname);
	destroyWindow(_dctwindow);
}