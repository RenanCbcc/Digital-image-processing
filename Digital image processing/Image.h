#pragma once

#include <cassert>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <string>
#include <exception>
#include <fstream>
#include <cctype>
#include <iostream>
#include <array>

#include <boost/numeric/ublas/matrix.hpp>

#include "Coding.hpp"
#include "Huffman.hpp"

typedef unsigned int uint;
typedef uint8_t Byte;
typedef double PixelDataType;

using boost::numeric::ublas::matrix;

class Image;

// load a ppm file (P3 or P6 version)
Image loadPPM(std::string path);

// fast version of atoi. No error checking, nothing.
int fast_atoi(const char * str);

// image class handling three matrix<PixelDataType>s (RGB, YUV, whatever) with one byte pixels
class Image
{
	// NESTED ENUMS
public:
	enum ColorSpace
	{
		RGB,
		YCbCr
	};

	enum SubsamplingMode
	{
		S444,       // full sampling
		S422,       // every second pixel in a row
		S411,       // every fourth pixel in a row
		S420,       // every second pixel in every second row
		S420_m,     // between vertical and horizontal pixels
		S420_lm,    // between vertical pixels
	};

	enum DCTMode
	{
		Simple,
		Matrix,
		Arai
	};

	// INTERFACE
public:
	// CTORS
	explicit Image(uint w, uint h, ColorSpace color);   // ctor
	Image(const Image& other);                          // copy ctor
	Image(Image&& other);                               // move ctor

	~Image(); // dtor

			  // ASSIGNMENTS
	Image& operator=(const Image &other);   // copy assignment
	Image& operator=(Image &&other);        // move assignment

											// METHODS
											// returns a new image object, this object won't be modified
	Image convertToColorSpace(ColorSpace target_space) const;

	// apply subsampling to the color matrix<PixelDataType>s (Cb, Cr)
	void applySubsampling(SubsamplingMode mode);

	void applyDCT(DCTMode mode);
	void applyQuantization(const matrix<Byte>& q_table_y, const matrix<Byte>& q_table_c);
	void applyDCdifferenceCoding();
	void doZigZagSorting();
	void doRLEandCategoryCoding();
	void doHuffmanEncoding(SymbolCodeMap &Y_DC,
		SymbolCodeMap &Y_AC,
		SymbolCodeMap &C_DC,
		SymbolCodeMap &C_AC);

	// JPEG SEGMENTS
	void writeJPEG(std::string file);

	// HELPER
private:
	struct Mask;
	void subsample(matrix<PixelDataType>&, int, int, Mask&, bool, SubsamplingMode);

	// ACCESSORS
public:
	uint width, height;
	uint real_width, real_height;
	uint subsample_width, subsample_height;
	matrix<PixelDataType> &R, &G, &B;
	matrix<PixelDataType> &Y, &Cb, &Cr;

	// HIDDEN MEMBERS
private:
	ColorSpace color_space_type;
	matrix<PixelDataType> one, two, three;
	matrix<PixelDataType> DctY, DctCb, DctCr;
	matrix<int> QY, QCb, QCr;
	matrix<std::vector<Category_Code>> CategoryCodeY, CategoryCodeCb, CategoryCodeCr;
	matrix<Bitstream> BitstreamY, BitstreamCb, BitstreamCr;
};