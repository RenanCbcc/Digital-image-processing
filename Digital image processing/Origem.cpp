
#include<opencv2/opencv.hpp>

#include<iostream>
#include<vector>
#include<conio.h>  

void rgbTOhsv(const cv::Mat & );
void rgbTOblue(const cv::Mat &);
void rgbTOgray(const cv::Mat &);
void getPixel(const cv::Mat &,int , int );


int main() {
	cv::Mat imgOriginal;
	imgOriginal = cv::imread("belem.jpg");   

	if (imgOriginal.empty()) {                                 
		std::cout << "error: image not read from file\n\n";    
		_getch();                                              
		return(0);                                             
	}
	//getPixel(imgOriginal, 100, 100);
	//rgbTOgray(imgOriginal);
	rgbTOblue(imgOriginal);
	return 0;
}

void getPixel(const cv::Mat &input, int x, int y) {
	if (!(x <= 0 || x>=800) && !(y <= 0 || y>=533)) {
		cv::Vec3b intensity = input.at< cv::Vec3b >(y, x);
		// ( Note: We could write img_rgb.at< cv::Vec3b >(x,y)[0] )
		//
		uchar blue = intensity[0];
		uchar green = intensity[1];
		uchar red = intensity[2];
		std::cout << "At (x,y) = (" << x << ", " << y <<
			"): (blue, green, red) = (" <<
			(unsigned int)blue <<
			", " << (unsigned int)green << ", " <<
			(unsigned int)red << ")" << std::endl;
		
	}
	else {
		std::cout << "Error: coordenates out of range \n";
	}

	cv::namedWindow("img_in", cv::WINDOW_AUTOSIZE);
	cv::imshow("img_in", input);

	cv::waitKey(0);

}

void rgbTOblue(const cv::Mat &input) {
	
	cv::Mat output;
	std::vector< cv::Mat> channel;
	cv::split(input, channel);

	cv::Mat b = channel[0], g = channel[1], r = channel[2];

	channel[1] = cv::Mat::zeros(input.rows, input.cols, CV_8UC1); // green channel is set to 0
	channel[2] = cv::Mat::zeros(input.rows, input.cols, CV_8UC1);// red channel is set to 0
	
	cv::merge(channel, output);

	cv::namedWindow("img_in", cv::WINDOW_AUTOSIZE);
	cv::namedWindow("img_out", cv::WINDOW_AUTOSIZE);

	cv::imshow("img_in", input);
	cv::imshow("img_out", output);
	// Wait for the user to hit a key, windows will self destruct
	//
	cv::waitKey(0);

}
void rgbTOgray(const cv::Mat & input) {
	// Create some windows to show the input
	// and output images in.
	cv::namedWindow("img_in", cv::WINDOW_AUTOSIZE);
	cv::namedWindow("img_out", cv::WINDOW_AUTOSIZE);
	// Create a window to show our input image
	cv::imshow("img_in", input);
	// Create an image to hold the smoothed output

	/*cv::Mat makes life simpler for us; we just instantiate an output matrix, out, and it will automatically
	resize/reallocate and deallocate itself as necessary as it is used.*/
	cv::Mat output;

	cv::cvtColor(input, output, cv::COLOR_BGR2GRAY);

	cv::imshow("img_out", output);
	// Wait for the user to hit a key, windows will self destruct
	//
	cv::waitKey(0);
}

inline float min(float a, float b) { return (a < b) ? a : b; }
inline float max(float a, float b) { return (a > b) ? a : b; }

void rgbTOhsv(const cv::Mat & input) {
	// Create some windows to show the input
	// and output images in.
	cv::namedWindow("img_in", cv::WINDOW_AUTOSIZE);
	cv::namedWindow("img_out", cv::WINDOW_AUTOSIZE);
	// Create a window to show our input image
	//
	cv::imshow("img_in", input);
	// Create an image to hold the smoothed output
	
	/*cv::Mat makes life simpler for us; we just instantiate an output matrix, out, and it will automatically
	resize/reallocate and deallocate itself as necessary as it is used.*/
	cv::Mat output; 

	/* RGB and HSV values are in the range from 0 to 1.0 */
	const float noHue = -1.0;

	float minRGB = min(r, min(g, b)), maxRGB = max(r, max(g, b));
	float deltaRGB = maxRGB - minRGB;
	v = maxRGB;
	if (maxRGB != 0.0)
		s = deltaRGB / maxRGB;
	else
		s = 0.0;
	Color Models and Color Applications
		582if (s <= 0.0)
		h = noHue;
		else {
			if (r == maxRGB)
				h = (g - b) / deltaRGB;
			else
				if (g == maxRGB)
					h = 2.0 + (b - r) / deltaRGB;
				else
					if (b == maxRGB)
						h = 4.0 + (r - g) / deltaRGB;
			h *= 60.0;
			if (h < 0.0)
				h += 360.0;
			h /= 360.0;
		}
	cv::imshow("img_out", output);
	// Wait for the user to hit a key, windows will self destruct
	//
	cv::waitKey(0);
}