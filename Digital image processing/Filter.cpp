#include "Filter.h"



Filter::Filter(std::string path)
{
	 image =cv::imread(path);

	if (image.empty()) {
		std::cout << "Error: image not read from file\n\n";
		
	}
}


