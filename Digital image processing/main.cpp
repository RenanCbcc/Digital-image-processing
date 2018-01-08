
//Boost https://www.youtube.com/watch?v=5AmwIwedTCM

#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include "Image.h"


int main(int argc, char *argv[]) {

	if (argc < 2) {
		std::cout << "No filename was written" << std::endl;
		return 0;
	}

	std::string ppmFilename = argv[1];
	auto img = loadPPM(ppmFilename);

	std::string jpgFilename;

	if (argc < 3)
	{
		jpgFilename = "noname.jpg";
	}
	else
	{
		jpgFilename = argv[2];
	}

	img.writeJPEG(jpgFilename);

	return 0;
}