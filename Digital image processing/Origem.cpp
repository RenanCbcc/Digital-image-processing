
#include<opencv2/opencv.hpp>

#include<iostream>
#include<vector>
#include<conio.h>  

void rgbTOhsv(const cv::Mat & );
void rgbTOblue(const cv::Mat &);
void rgbTOgray(const cv::Mat &);
void getPixel(const cv::Mat &,int , int );
void negative(const cv::Mat &);
void contrast(const cv::Mat &);
void histogram(const cv::Mat &);
void equalization(const cv::Mat &);
int computeOutput(int, int, int, int, int);

int nain() {
	cv::Mat origin;
	origin = cv::imread("belem.jpg");

	if (origin.empty()) {
		std::cout << "error: image not read from file\n\n";    
		_getch();                                              
		return(0);                                             
	}
	//getPixel(origin, 100, 100);
	//rgbTOgray(origin);
	//rgbTOblue(origin);
	//rgbTOhsv(origin);
	//histogram(origin);
	//negative(origin);
	//contrast(origin);//r1 = 70 s1 = 0 r2 = 140 s2 = 255
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

void negative(const cv::Mat &input)
{
	cv::Mat output = cv::Mat::zeros(input.size(), input.type());
	// create a matrix with all elements equal to 255 for subtraction
	cv::Mat sub_mat = cv::Mat::ones(input.size(), input.type()) * 64;

	cv::subtract(
		sub_mat, // First input array
		input, // Second input array
		output // Result array
	);

	cv::namedWindow("img_in", cv::WINDOW_AUTOSIZE);
	cv::namedWindow("img_out", cv::WINDOW_AUTOSIZE);

	cv::imshow("img_in", input);
	cv::imshow("img_out", output);
	// Wait for the user to hit a key, windows will self destruct
	//
	cv::waitKey(0);

}

void contrast(const cv::Mat & input)
{
	cv::Mat output = input.clone();

	int r1, s1, r2, s2;
	std::cout << "Enter r1: " << std::endl; std::cin >> r1;
	std::cout << "Enter s1: " << std::endl; std::cin >> s1;
	std::cout << "Enter r2: " << std::endl; std::cin >> r2;
	std::cout << "Enter s2: " << std::endl; std::cin >> s2;

	for (int y = 0; y < input.rows; y++) {
		for (int x = 0; x < input.cols; x++) {
			for (int c = 0; c < 3; c++) {
				int result = computeOutput(input.at<cv::Vec3b>(y, x)[c], r1, s1, r2, s2);
				output.at<cv::Vec3b>(y, x)[c] = cv::saturate_cast<uchar>(result);
			}
		}
	}

	cv::namedWindow("Original Image", cv::WINDOW_AUTOSIZE);
	cv::imshow("Original Image", input);

	cv::namedWindow("New Image", cv::WINDOW_AUTOSIZE);
	cv::imshow("New Image", output);

	cv::waitKey();
}


void histogram(const cv::Mat &input)
{
	// allcoate memory for no of pixels for each intensity value
	int histogram[256];

	// initialize all intensity values to 0
	for (int i = 0; i < 255; i++)
	{
		histogram[i] = 0;
	}

	// calculate the no of pixels for each intensity values
	for (int y = 0; y < input.rows; y++)
		for (int x = 0; x < input.cols; x++)
			histogram[(int)input.at<uchar>(y, x)]++;

	
	// draw the histograms
	int hist_w = 512; int hist_h = 400;
	int bin_w = cvRound((double)hist_w / 256);

	cv::Mat histImage(hist_h, hist_w, CV_8UC1, cv::Scalar(255, 255, 255));

	// find the maximum intensity element from histogram
	int max = histogram[0];
	for (int i = 1; i < 256; i++) {
		if (max < histogram[i]) {
			max = histogram[i];
		}
	}

	// normalize the histogram between 0 and histImage.rows

	for (int i = 0; i < 255; i++) {
		histogram[i] = ((double)histogram[i] / max)*histImage.rows;
	}


	// draw the intensity line for histogram
	for (int i = 0; i < 255; i++)
	{
		line(histImage, cv::Point(bin_w*(i), hist_h),
			cv::Point(bin_w*(i), hist_h - histogram[i]),
			cv::Scalar(0, 0, 0), 1, 8, 0);
	}

	// display histogram
	cv::namedWindow("Intensity Histogram", CV_WINDOW_AUTOSIZE);
	cv::imshow("Intensity Histogram", histImage);

	cv::namedWindow("img_in", CV_WINDOW_AUTOSIZE);
	cv::imshow("img_in", input);
	cv::waitKey();

}

void 
equalization(const cv::Mat &)
{

}

int computeOutput(int x, int r1, int s1, int r2, int s2)
{
	float result;
	if (0 <= x && x <= r1) {
		result = s1 / r1 * x;
	}
	else if (r1 < x && x <= r2) {
		result = ((s2 - s1) / (r2 - r1)) * (x - r1) + s1;
	}
	else if (r2 < x && x <= 255) {
		result = ((255 - s2) / (255 - r2)) * (x - r2) + s2;
	}
	return (int)result;
}

void rgbTOblue(const cv::Mat &input) {
	
	cv::Mat output;
	std::vector< cv::Mat> channel;
	cv::split(input, channel);

	cv::Mat b = channel[0], g = channel[1], r = channel[2];

	channel[1] = cv::Mat::zeros(input.rows, input.cols, CV_8UC1); // green channel is set to 0
	channel[2] = cv::Mat::zeros(input.rows, input.cols, CV_8UC1);// red channel is set to 0
	
	cv::merge(channel, output);

	
	cv::namedWindow("Original Image", cv::WINDOW_AUTOSIZE);
	cv::imshow("Original Image", input);

	cv::namedWindow("New Image", cv::WINDOW_AUTOSIZE);
	cv::imshow("New Image", output);

	cv::waitKey();

}
void rgbTOgray(const cv::Mat & input) {
	
	/*cv::Mat makes life simpler for us; we just instantiate an output matrix, out, and it will automatically
	resize/reallocate and deallocate itself as necessary as it is used.*/
	cv::Mat output;

	cv::cvtColor(input, output, cv::COLOR_BGR2GRAY);

	cv::namedWindow("Original Image", cv::WINDOW_AUTOSIZE);
	cv::imshow("Original Image", input);

	cv::namedWindow("New Image", cv::WINDOW_AUTOSIZE);
	cv::imshow("New Image", output);

	cv::waitKey();

}

inline float min(float a, float b) { return (a < b) ? a : b; }
inline float max(float a, float b) { return (a > b) ? a : b; }


inline float con(float a, float b) { 
/*
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

		*/
}


void rgbTOhsv(const cv::Mat & input) {

	/*cv::Mat makes life simpler for us; we just instantiate an output matrix, out, and it will automatically
	resize/reallocate and deallocate itself as necessary as it is used.*/
	cv::Mat output;
	std::vector< cv::Mat> channel;
	cv::split(input, channel);

	/* RGB and HSV values are in the range from 0 to 1.0 */
	/*
	const float noHue = -1.0;
	cv::Mat r = channel[0], g = channel[1], b = channel[2];
	cv::merge(channel, output);

	for (int i = 0; i < channel.size(); i++) {
		
	
	}
	
	*/
	cvtColor(input, output, CV_RGB2HSV);
	cv::namedWindow("Original Image", cv::WINDOW_AUTOSIZE);
	cv::imshow("Original Image", input);

	cv::namedWindow("New Image", cv::WINDOW_AUTOSIZE);
	cv::imshow("New Image", output);

	cv::waitKey();
}