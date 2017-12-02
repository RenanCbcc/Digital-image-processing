#pragma once
#include<opencv2/opencv.hpp>

#include<iostream>
#include<vector>  
#include<string>
class Filter
{
public:
	Filter(std::string );
	void rgbTOhsv(const cv::Mat &);
	void rgbTOblue(const cv::Mat &);
	void rgbTOgray(const cv::Mat &);
	void getPixel(const cv::Mat &, int, int);
	void negative(const cv::Mat &);
	void contrast(const cv::Mat &);
	void histogram(const cv::Mat &);
	int computeOutput(int x, int r1, int s1, int r2, int s2);
	float min(float, float);
	float max(float, float);


private:
		cv::Mat image;
};

