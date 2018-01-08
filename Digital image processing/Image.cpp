#include "Image.h"
#include <chrono>
#include <array>
#include <streambuf>
#include <sstream>
#include <future>
#include <omp.h>
#include <future>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "JpegSegments.hpp"
#include "Dct.hpp"

using boost::numeric::ublas::matrix_range;
using boost::numeric::ublas::range;
using boost::numeric::ublas::subrange;
using boost::numeric::ublas::zero_matrix;

using namespace std::chrono;

static const auto debug = false;

//
// CONSTRUCTORS
//
// ctor
Image::Image(uint w, uint h, ColorSpace color)
	: color_space_type(color),
	width(w), height(h),
	real_width(w), real_height(h),
	subsample_width(w), subsample_height(h),
	one(h, w), two(h, w), three(h, w),
	Y(one), Cb(two), Cr(three),
	R(one), G(two), B(three)
{
	if (debug) std::cout << "Image::Image()\n" << std::endl;
}

// copy ctor
Image::Image(const Image& other)
	: color_space_type(other.color_space_type),
	width(other.width), height(other.height),
	real_width(other.real_width), real_height(other.real_height),
	subsample_width(other.subsample_width), subsample_height(other.subsample_height),
	one(other.one), two(other.two), three(other.three),
	Y(one), Cb(two), Cr(three),
	R(one), G(two), B(three)
{
	if (debug) std::cout << "Image::Image(Image&)\n" << std::endl;
}

// move ctor
Image::Image(Image&& other)
	: color_space_type(other.color_space_type),
	width(other.width), height(other.height),
	real_width(other.real_width), real_height(other.real_height),
	subsample_width(other.subsample_width), subsample_height(other.subsample_height),
	one(std::move(other.one)), two(std::move(other.two)), three(std::move(other.three)),
	Y(one), Cb(two), Cr(three),
	R(one), G(two), B(three)
{
	if (debug) std::cout << "Image::Image(Image&&)\n" << std::endl;
}

// dtor
Image::~Image()
{}

//
// ASSIGNMENTS
//
// copy assignment
Image& Image::operator = (const Image &other) {
	if (this != &other) {
		one = other.one;
		two = other.two;
		three = other.three;
		width = other.width;
		height = other.height;
		real_width = other.real_width;
		real_height = other.real_height;
		subsample_width = other.subsample_width;
		subsample_height = other.subsample_height;
		color_space_type = other.color_space_type;
		if (debug) std::cout << "Image::operator=(Image&)\n" << std::endl;
	}
	return *this;
}

// move assignment
Image& Image::operator=(Image &&other) {
	if (this != &other) {
		one = std::move(other.one);
		two = std::move(other.two);
		three = std::move(other.three);
		width = other.width;
		height = other.height;
		real_width = other.real_width;
		real_height = other.real_height;
		subsample_width = other.subsample_width;
		subsample_height = other.subsample_height;
		color_space_type = other.color_space_type;
		if (debug) std::cout << "Image::operator=(Image&&)\n" << std::endl;
	}
	return *this;
}

Image Image::convertToColorSpace(ColorSpace target_color_space) const {
	// no converting if already in target color space
	if (color_space_type == target_color_space)
		return *this;

	Image converted(*this);
	const auto num_pixel = converted.width * converted.height;

	switch (target_color_space) {
	case ColorSpace::YCbCr:
	{
		assert(color_space_type == ColorSpace::RGB);

		// matrix factors
		/*
		0.299  0.587  0.114
		-0.169 -0.331  0.500
		0.500 -0.419 -0.081
		*/
		static const float Flat[]{ .0f, 256 / 2.f, 256 / 2.f };
		static const float Yv[]{ .299f,   .587f,   .114f };
		static const float Cb[]{ -.1687f, -.3312f,  .5f };
		static const float Cr[]{ .5f,    -.4186f, -.0813f };

		for (uint x = 0; x < num_pixel; ++x) {
			auto& r = R.data()[x];
			auto& g = G.data()[x];
			auto& b = B.data()[x];

			converted.Y.data()[x] = Flat[0] + (Yv[0] * r + Yv[1] * g + Yv[2] * b) - 128;
			converted.Cb.data()[x] = Flat[1] + (Cb[0] * r + Cb[1] * g + Cb[2] * b) - 128;
			converted.Cr.data()[x] = Flat[2] + (Cr[0] * r + Cr[1] * g + Cr[2] * b) - 128;
		}

		converted.color_space_type = ColorSpace::YCbCr;
	}
	break;
	case ColorSpace::RGB:
	{
		assert(color_space_type == ColorSpace::YCbCr);

		// matrix factors
		/*
		1.000  0.000  1.402
		1.000 -0.344 -0.714
		1.000  1.772  0.000
		*/
		static const float Flat[]{ .0f, 256 / 2.f, 256 / 2.f };
		static const float r[]{ 1.f,  .0f,    1.402f };
		static const float g[]{ 1.f, -.344f,  -.714f };
		static const float b[]{ 1.f, 1.772f,   .0f };

		for (uint x = 0; x < num_pixel; ++x) {
			auto y = Y.data()[x] + 128;
			auto cb = Cb.data()[x] + 128;
			auto cr = Cr.data()[x] + 128;

			converted.R.data()[x] = (r[0] * y + r[1] * cb + r[2] * cr);
			converted.G.data()[x] = (g[0] * y + g[1] * cb + g[2] * cr);
			converted.B.data()[x] = (b[0] * y + b[1] * cb + b[2] * cr);
		}

		converted.color_space_type = ColorSpace::RGB;
	}
	break;
	}
	return converted;
}

// simple "Matrix" used by the subsampling method
struct Image::Mask
{
	std::vector<Byte> row;
	bool scanline_jump = false;

	Mask() {}

	Mask(std::vector<Byte>&& r1, bool _scanline_jump)
		: row(std::move(r1))
	{
		scanline_jump = _scanline_jump;
	}

	uint rowsize() const { return static_cast<uint>(row.size()); }
};

void Image::subsample(matrix<PixelDataType>& chan, int hor_res_div, int vert_res_div, Mask& mat, bool averaging, SubsamplingMode mode)
{
	matrix<PixelDataType> new_chan(chan.size1() / vert_res_div, chan.size2() / hor_res_div);

	// subsample Cr matrix<PixelDataType>
	// size1() = rows
	// size2() = columns
	auto pixidx = 0U;
	auto pixidx2 = 0U;
	for (auto y = 0U; y < chan.size1(); y += 2) {
		for (auto x = 0U; x < chan.size2(); x += mat.rowsize()) {
			PixelDataType pix_val = 0;
			for (auto m = 0U; m < mat.rowsize(); ++m) {
				pix_val += mat.row[m] * chan(y, x + m);
			}
			new_chan.data()[pixidx++] = pix_val;
		}

		if (!mat.scanline_jump) {
			// go through next scanline and average the pixels
			if (averaging) {
				for (auto x = 0U; x < chan.size2(); x += mat.rowsize()) {
					PixelDataType pix_val = 0;
					for (auto m = 0U; m < mat.rowsize(); ++m) {
						pix_val += mat.row[m] * chan(y + 1, x + m);
					}
					(new_chan.data()[pixidx2++] += (pix_val)) /= ((mode == S420_m) ? 4 : 2);
				}
			}
			// go through next scanline
			else
				--y;
		}

	}

	std::swap(chan, new_chan);
}

void Image::applySubsampling(SubsamplingMode mode)
{
	//enum SubsamplingMode
	//{
	//    S444,       // full sampling
	//    S422,       // every second pixel in a row
	//    S411,       // every fourth pixel in a row
	//    S420,       // every second pixel in every second row
	//    S420_lm,    // between vertical pixels
	//    S420_m      // between vertical and horizontal pixels
	//};

	bool averaging = false;

	int vert_res_div = 1;
	int hor_res_div = 1;

	Mask mat;

	switch (mode) {
	case Image::S444:
		// no subsampling
		// xx xx
		// xx xx
		return;
		break;
	case Image::S422:
		// x- x-
		// x- x-
		mat = Mask{ { 1, 0 }, false };

		vert_res_div = 1;
		hor_res_div = 2;
		break;
	case Image::S411:
		// x- --
		// x- --
		mat = Mask{ { 1, 0, 0, 0 }, false };

		vert_res_div = 1;
		hor_res_div = 4;
		break;
	case Image::S420:
		// x- x-
		// -- --
		mat = Mask{ { 1, 0 }, true };

		vert_res_div = 2;
		hor_res_div = 2;
		break;
	case Image::S420_m:
		// like S420, but taking the mean in vertical and horizontal direction
		// ++ ++
		// ++ ++
		mat = Mask{ { 1, 1 }, false };

		vert_res_div = 2;
		hor_res_div = 2;
		averaging = true;
		break;
	case Image::S420_lm:
		// like S420, but taking the mean only in vertical direction
		// +- +-
		// +- +-
		mat = Mask{ { 1, 0 }, false };

		vert_res_div = 2;
		hor_res_div = 2;
		averaging = true;
		break;
	default:
		break;
	}

	subsample_width = width / hor_res_div;
	subsample_height = height / vert_res_div;

	// subsampling Cr matrix<PixelDataType>
	subsample(Cr, hor_res_div, vert_res_div, mat, averaging, mode);

	// subsample Cb matrix<PixelDataType>
	subsample(Cb, hor_res_div, vert_res_div, mat, averaging, mode);
}

//
// PPM loading
//

// fast version of atoi. No error checking, nothing.
int fast_atoi(const char * str)
{
	int val = 0;
	while (*str) {
		val = val * 10 + (*str++ - '0');
	}
	return val;
};
struct PPMFileBuffer
{
	std::string& file;
	std::string::size_type file_pos;
	std::string::size_type eof;

	PPMFileBuffer(std::string& file)
		: file{ file },
		file_pos{ 0 },
		eof{ file.size() } {}

	Byte read_byte() {
		return file[file_pos++];
	}

	void read_word(std::string& buf) {
		buf.clear();

		auto c = read_byte();

		if (std::isspace(c)) // discard any leading whitespace
			c = discard_whitespace();

		auto first = file_pos - 1;
		while (1) {
			if (c == '#') {
				discard_current_line();
				first = file_pos;
			}
			else if (isspace(c)) {
				buf.insert(0, &file[first], file_pos - 1 - first);
				break;
			}
			else if (is_eof()) {
				buf.insert(0, &file[first], file_pos - first);
				break;
			}

			c = read_byte();
		}
	}

private:
	bool is_eof() {
		return file_pos >= eof;
	}

	Byte discard_whitespace() {
		while (isspace(file[file_pos]) && !is_eof())
			++file_pos;
		return read_byte();
	}

	void discard_current_line() {
		while ((file[file_pos++] != '\n') && !is_eof());
	}
};

// load a ppm file into an RGBImage
void loadP3PPM(PPMFileBuffer& file, Image& img, double scale_factor) {
	std::string buf;
	buf.reserve(32);

	const auto max = img.width * img.height;
	for (auto x = 0U; x < max; ++x) {
		file.read_word(buf);
		img.R.data()[x] = fast_atoi(buf.c_str()) * scale_factor;

		file.read_word(buf);
		img.G.data()[x] = fast_atoi(buf.c_str()) * scale_factor;

		file.read_word(buf);
		img.B.data()[x] = fast_atoi(buf.c_str()) * scale_factor;
	}
}

// Needs the ifstream opened as binary!
void loadP6PPM(PPMFileBuffer& ppm, Image& img, double scale_factor) {
	const auto max = img.width * img.height;
	for (auto x = 0U; x < max; ++x) {
		img.R.data()[x] = ppm.read_byte() * scale_factor;
		img.G.data()[x] = ppm.read_byte() * scale_factor;
		img.B.data()[x] = ppm.read_byte() * scale_factor;
	}
}

// top level load function with a path to a ppm file
Image loadPPM(std::string path) {
	auto start = high_resolution_clock::now();

	std::ifstream::sync_with_stdio(false);
	std::ifstream ppm_file{ path, std::fstream::binary };

	if (!ppm_file.is_open())
		throw std::runtime_error("Failed to open \"" + std::string(path) + "\"");

	ppm_file.seekg(0, std::ios::end);
	auto file_size = 0 + ppm_file.tellg();
	ppm_file.seekg(0, std::ios::beg);

	// using stringstream
	std::stringstream strstrbuf;
	strstrbuf << ppm_file.rdbuf();
	auto raw_buf = strstrbuf.str();

	ppm_file.close();

	PPMFileBuffer ppm{ raw_buf };

	std::string buf;
	buf.reserve(32);

	// magic number
	ppm.read_word(buf);
	std::string magic = buf;
	if (magic != "P3" && magic != "P6")
		throw std::runtime_error("Only P3 and P6 format is supported!");

	// width and height
	ppm.read_word(buf);
	uint width = std::stoi(buf);

	ppm.read_word(buf);
	uint height = std::stoi(buf);

	ppm.read_word(buf);
	auto max_color = std::stoi(buf);

	assert((max_color < 256) && "Only 1 byte colors supported for now!");

	// scale according to max_color (if max_color is 15 -> white is 15! that means we have to scale that up to 255)
	auto scale_factor = 255. / max_color;

	Image img(width, height, Image::RGB);

	if (magic == "P3")
		loadP3PPM(ppm, img, scale_factor);
	else if (magic == "P6") {
		loadP6PPM(ppm, img, scale_factor);
	}


	img.real_height = height;
	img.real_width = width;

	//adjust size of our image for using 16x16 blocks
	if (width % 16 != 0 || height % 16 != 0)
	{
		if (img.width % 16 != 0) {
			img.width += 16;
			img.width -= img.width % 16;
		}
		if (img.height % 16 != 0) {
			img.height += 16;
			img.height -= img.height % 16;
		}

		img.subsample_height = img.height;
		img.subsample_width = img.width;

		img.R.resize(img.height, img.width, true);
		img.G.resize(img.height, img.width, true);
		img.B.resize(img.height, img.width, true);

		// Fill the new Pixel with data from the border. 
		// Right
		for (auto y = 0U; y < height; ++y)
		{
			for (auto x = width; x < img.width; ++x)
			{
				img.R(y, x) = img.R(y, width - 1);
				img.G(y, x) = img.G(y, width - 1);
				img.B(y, x) = img.B(y, width - 1);
			}
		}

		// Bottom
		for (auto y = height; y < img.height; ++y)
		{
			for (auto x = 0U; x < width; ++x)
			{
				img.R(y, x) = img.R(height - 1, x);
				img.G(y, x) = img.G(height - 1, x);
				img.B(y, x) = img.B(height - 1, x);
			}
		}

		// Corner
		for (auto y = height; y < img.height; ++y)
		{
			for (auto x = width; x < img.width; ++x)
			{
				img.R(y, x) = img.R(height - 1, width - 1);
				img.G(y, x) = img.G(height - 1, width - 1);
				img.B(y, x) = img.B(height - 1, width - 1);
			}
		}

	}

	auto end = high_resolution_clock::now();
	std::cout << "PPM loading took " << duration_cast<milliseconds>(end - start).count() << " ms\n";

	return img;
}

void Image::applyDCT(DCTMode mode)
{
	std::function<void(const matrix_range<matrix<PixelDataType>>&, matrix_range<matrix<PixelDataType>>&)> dctFn;

	switch (mode) {
	case Simple:
		dctFn = dctDirect;
		break;
	case Matrix:
		dctFn = dctMat;
		break;
	case Arai:
		dctFn = dctArai;
		break;
	default:
		assert(!"This DCT mode isn't supported!");
	}

	const auto blocksize = 8;

	//omp_set_num_threads(4); // can also be set via an environment variable

	DctY = zero_matrix<PixelDataType>(Y.size1(), Y.size2());
	DctCb = zero_matrix<PixelDataType>(Cb.size1(), Cb.size2());
	DctCr = zero_matrix<PixelDataType>(Cr.size1(), Cr.size2());

#pragma omp parallel for
	for (int h = 0; h < height; h += blocksize) {
		for (int w = 0; w < width; w += blocksize) {
			// generate slices for data source and the destination of the dct result
			const matrix_range<matrix<PixelDataType>> slice_src(Y, range(h, h + blocksize), range(w, w + blocksize));
			matrix_range<matrix<PixelDataType>> slice_dst(DctY, range(h, h + blocksize), range(w, w + blocksize));
			dctFn(slice_src, slice_dst);
		}
	}

#pragma omp parallel for
	for (int h = 0; h < subsample_height; h += blocksize) {
		for (int w = 0; w < subsample_width; w += blocksize) {
			// generate slices for data source and the destination of the dct result
			const matrix_range<matrix<PixelDataType>> slice_src(Cb, range(h, h + blocksize), range(w, w + blocksize));
			matrix_range<matrix<PixelDataType>> slice_dst(DctCb, range(h, h + blocksize), range(w, w + blocksize));
			dctFn(slice_src, slice_dst);
		}
	}

#pragma omp parallel for
	for (int h = 0; h < subsample_height; h += blocksize) {
		for (int w = 0; w < subsample_width; w += blocksize) {
			// generate slices for data source and the destination of the dct result
			const matrix_range<matrix<PixelDataType>> slice_src(Cr, range(h, h + blocksize), range(w, w + blocksize));
			matrix_range<matrix<PixelDataType>> slice_dst(DctCr, range(h, h + blocksize), range(w, w + blocksize));
			dctFn(slice_src, slice_dst);
		}
	}
}

void Image::applyQuantization(const matrix<Byte>& qtable_y, const matrix<Byte>& qtable_c) {
	assert(qtable_y.size1() == 8);
	assert(qtable_y.size2() == 8);
	assert(qtable_c.size1() == 8);
	assert(qtable_c.size2() == 8);

	QY = zero_matrix<int>(DctY.size1(), DctY.size2());
	QCb = zero_matrix<int>(DctCb.size1(), DctCb.size2());
	QCr = zero_matrix<int>(DctCr.size1(), DctCr.size2());

#pragma omp parallel for
	for (int h = 0; h < height; h += blocksize) {
		for (int w = 0; w < width; w += blocksize) {
			// generate slices for data source and the destination of the dct result
			const matrix_range<matrix<PixelDataType>> slice_src(DctY, range(h, h + blocksize), range(w, w + blocksize));
			matrix_range<matrix<int>> slice_dst(QY, range(h, h + blocksize), range(w, w + blocksize));
			slice_dst.assign(quantize(slice_src, qtable_y));
		}
	}

#pragma omp parallel for
	for (int h = 0; h < subsample_height; h += blocksize) {
		for (int w = 0; w < subsample_width; w += blocksize) {
			// generate slices for data source and the destination of the dct result
			const matrix_range<matrix<PixelDataType>> slice_src(DctCb, range(h, h + blocksize), range(w, w + blocksize));
			matrix_range<matrix<int>> slice_dst(QCb, range(h, h + blocksize), range(w, w + blocksize));
			slice_dst.assign(quantize(slice_src, qtable_c));
		}
	}

#pragma omp parallel for
	for (int h = 0; h < subsample_height; h += blocksize) {
		for (int w = 0; w < subsample_width; w += blocksize) {
			// generate slices for data source and the destination of the dct result
			const matrix_range<matrix<PixelDataType>> slice_src(DctCr, range(h, h + blocksize), range(w, w + blocksize));
			matrix_range<matrix<int>> slice_dst(QCr, range(h, h + blocksize), range(w, w + blocksize));
			slice_dst.assign(quantize(slice_src, qtable_c));
		}
	}
}

void Image::applyDCdifferenceCoding() {
	int b = 0;
	for (int h = 0; h < height; h += 2 * blocksize) {
		for (int w = 0; w < width; w += 2 * blocksize) {
			// when subsampling is on, dc differences of y are NOT ordered left-right top-bottom
			auto tmp = QY(h, w);
			QY(h, w) = tmp - b;
			b = tmp;

			tmp = QY(h, w + 8);
			QY(h, w + 8) = tmp - b;
			b = tmp;

			tmp = QY(h + 8, w);
			QY(h + 8, w) = tmp - b;
			b = tmp;

			tmp = QY(h + 8, w + 8);
			QY(h + 8, w + 8) = tmp - b;
			b = tmp;
		}
	}

	b = 0;
	for (int h = 0; h < subsample_height; h += blocksize) {
		for (int w = 0; w < subsample_width; w += blocksize) {
			auto tmp = QCb(h, w);
			QCb(h, w) = tmp - b;
			b = tmp;
		}
	}

	b = 0;
	for (int h = 0; h < subsample_height; h += blocksize) {
		for (int w = 0; w < subsample_width; w += blocksize) {
			auto tmp = QCr(h, w);
			QCr(h, w) = tmp - b;
			b = tmp;
		}
	}
}

void Image::doRLEandCategoryCoding() {
	CategoryCodeY.clear();
	CategoryCodeCb.clear();
	CategoryCodeCr.clear();

	CategoryCodeY.resize(QY.size1() / 8, QY.size2() / 8);
	CategoryCodeCb.resize(QCb.size1() / 8, QCb.size2() / 8);
	CategoryCodeCr.resize(QCr.size1() / 8, QCr.size2() / 8);

	// the data in the CategoryCodeXX vectors must be sequential correct (left to right, then top to bottom),
	// so no parallel execution of the loops possible
	// but we can run the rle parallel on the different channels
	auto f1 = std::async([&]() {
		for (int h = 0; h < height; h += blocksize) {
			for (int w = 0; w < width; w += blocksize) {
				// generate slices for the data source
				const auto slice_src = subrange(QY, h, h + blocksize, w, w + blocksize);

				auto rle_data = RLE_AC(slice_src);
				auto encoded_coeffs = encode_category(rle_data);
				CategoryCodeY(h / blocksize, w / blocksize) = encoded_coeffs;
			}
		}
	});

	auto f2 = std::async([&]() {
		for (int h = 0; h < subsample_height; h += blocksize) {
			for (int w = 0; w < subsample_width; w += blocksize) {
				// generate slices for the data source
				const auto slice_src = subrange(QCb, h, h + blocksize, w, w + blocksize);

				auto rle_data = RLE_AC(slice_src);
				auto encoded_coeffs = encode_category(rle_data);
				CategoryCodeCb(h / blocksize, w / blocksize) = encoded_coeffs;
			}
		}
	});

	auto f3 = std::async([&]() {
		for (int h = 0; h < subsample_height; h += blocksize) {
			for (int w = 0; w < subsample_width; w += blocksize) {
				// generate slices for the data source
				const auto slice_src = subrange(QCr, h, h + blocksize, w, w + blocksize);

				auto rle_data = RLE_AC(slice_src);
				auto encoded_coeffs = encode_category(rle_data);
				CategoryCodeCr(h / blocksize, w / blocksize) = encoded_coeffs;
			}
		}
	});

	// wait for results
	f1.get();
	f2.get();
	f3.get();
}

void Image::doHuffmanEncoding(SymbolCodeMap &Y_DC,
	SymbolCodeMap &Y_AC,
	SymbolCodeMap &C_DC,
	SymbolCodeMap &C_AC)
{

	// the data in the BitstreamXX vectors must be sequential correct (left to right, then top to bottom),
	// so no parallel execution of the loops possible
	// but we can run the rle parallel on the different channels

	auto f1 = std::async([&]() {
		BitstreamY.clear();
		BitstreamY.resize(CategoryCodeY.size1(), CategoryCodeY.size2());

		for (int i = 0; i < CategoryCodeY.size1(); i++) {
			for (int j = 0; j < CategoryCodeY.size2(); j++) {
				Bitstream stream;
				auto& data = CategoryCodeY(i, j);
				auto& code = data[0];

				auto encoded_DC = Y_DC[code.symbol];
				stream.push_back(encoded_DC.code, encoded_DC.length);
				stream << code.code;

				for (auto it = begin(data) + 1; it != end(data); ++it) {
					auto& code = *it;
					auto encoded_AC = Y_AC[code.symbol];
					stream.push_back(encoded_AC.code, encoded_AC.length);
					stream << code.code;
				}

				BitstreamY(i, j) = stream;
			}
		}
	});

	auto f2 = std::async([&]() {
		BitstreamCb.clear();
		BitstreamCb.resize(CategoryCodeCb.size1(), CategoryCodeCb.size2());

		for (int i = 0; i < CategoryCodeCb.size1(); ++i) {
			for (int j = 0; j < CategoryCodeCb.size2(); ++j) {
				Bitstream stream;
				auto& data = CategoryCodeCb(i, j);
				auto& code = data[0];

				auto encoded_DC = C_DC[code.symbol];
				stream.push_back(encoded_DC.code, encoded_DC.length);
				stream << code.code;

				for (auto it = begin(data) + 1; it != end(data); ++it) {
					auto& code = *it;
					auto encoded_AC = C_AC[code.symbol];
					stream.push_back(encoded_AC.code, encoded_AC.length);
					stream << code.code;
				}

				BitstreamCb(i, j) = stream;
			}
		}
	});

	auto f3 = std::async([&]() {
		BitstreamCr.clear();
		BitstreamCr.resize(CategoryCodeCr.size1(), CategoryCodeCr.size2());

		for (int i = 0; i < CategoryCodeCr.size1(); ++i) {
			for (int j = 0; j < CategoryCodeCr.size2(); ++j) {
				Bitstream stream;
				auto& data = CategoryCodeCr(i, j);
				auto& code = data[0];

				auto encoded_DC = C_DC[code.symbol];
				stream.push_back(encoded_DC.code, encoded_DC.length);
				stream << code.code;

				for (auto it = begin(data) + 1; it != end(data); ++it) {
					auto& code = *it;
					auto encoded_AC = C_AC[code.symbol];
					stream.push_back(encoded_AC.code, encoded_AC.length);
					stream << code.code;
				}

				BitstreamCr(i, j) = stream;
			}
		}
	});

	// wait for results
	f1.get();
	f2.get();
	f3.get();
}

void Image::writeJPEG(std::string file)
{
	auto start = high_resolution_clock::now();

	// printing some info
	std::cout << "Processing image size: " << real_width << "x" << real_height << std::endl;

	// color conversion to YCbCr
	*this = convertToColorSpace(YCbCr);

	// Cb/Cr subsampling
	applySubsampling(SubsamplingMode::S420_m);

	applyDCT(DCTMode::Arai);
	Y = zero_matrix<PixelDataType>(0, 0);
	Cb = zero_matrix<PixelDataType>(0, 0);
	Cr = zero_matrix<PixelDataType>(0, 0);

	// quantization
	const auto qtable_y = from_vector<int>({
		16, 11, 10, 16, 24, 40, 51, 61,
		12, 12, 14, 19, 26, 58, 60, 55,
		14, 13, 16, 24, 40, 57, 69, 56,
		14, 17, 22, 29, 51, 87, 80, 62,
		18, 22, 37, 56, 68, 109, 103, 77,
		24, 35, 55, 64, 81, 104, 113, 92,
		49, 64, 78, 87, 103, 121, 120, 101,
		72, 92, 95, 98, 112, 100, 103, 99
	});
	const auto qtable_c = from_vector<int>({
		17, 18, 24, 47, 99, 99, 99, 99,
		18, 21, 26, 66, 99, 99, 99, 99,
		24, 26, 56, 99, 99, 99, 99, 99,
		47, 66, 99, 99, 99, 99, 99, 99,
		99, 99, 99, 99, 99, 99, 99, 99,
		99, 99, 99, 99, 99, 99, 99, 99,
		99, 99, 99, 99, 99, 99, 99, 99,
		99, 99, 99, 99, 99, 99, 99, 99
	});

	applyQuantization(qtable_y, qtable_c);

	DctY = zero_matrix<PixelDataType>(0, 0);
	DctCb = zero_matrix<PixelDataType>(0, 0);
	DctCr = zero_matrix<PixelDataType>(0, 0);

	// DC difference coding
	applyDCdifferenceCoding();

	// zigzag, RLE and category encoding
	doRLEandCategoryCoding();

	QY = zero_matrix<int>(0, 0);
	QCb = zero_matrix<int>(0, 0);
	QCr = zero_matrix<int>(0, 0);

	// generate Huffman tables
	std::vector<int> Y_DC_symbols, Y_AC_symbols;
	std::vector<int> C_DC_symbols, C_AC_symbols;

	// count symbols
	for (const auto& data : CategoryCodeY.data()) {
		Y_DC_symbols.push_back(data[0].symbol);
		for (auto it = begin(data) + 1; it != end(data); ++it)
			Y_AC_symbols.push_back((*it).symbol);
	}
	for (const auto& data : CategoryCodeCb.data()) {
		C_DC_symbols.push_back(data[0].symbol);
		for (auto it = begin(data) + 1; it != end(data); ++it)
			C_AC_symbols.push_back((*it).symbol);
	}
	for (const auto& data : CategoryCodeCr.data()) {
		C_DC_symbols.push_back(data[0].symbol);
		for (auto it = begin(data) + 1; it != end(data); ++it)
			C_AC_symbols.push_back((*it).symbol);
	}

	auto Y_DC_huff = generateHuffmanCode(Y_DC_symbols);
	auto Y_AC_huff = generateHuffmanCode(Y_AC_symbols);

	auto C_DC_huff = generateHuffmanCode(C_DC_symbols);
	auto C_AC_huff = generateHuffmanCode(C_AC_symbols);

	auto& Y_DC_encoder = Y_DC_huff.first;
	auto& Y_DC_Huffman_Table = Y_DC_huff.second;

	auto& Y_AC_encoder = Y_AC_huff.first;
	auto& Y_AC_Huffman_Table = Y_AC_huff.second;

	auto& C_DC_encoder = C_DC_huff.first;
	auto& C_DC_Huffman_Table = C_DC_huff.second;

	auto& C_AC_encoder = C_AC_huff.first;
	auto& C_AC_Huffman_Table = C_AC_huff.second;

	// Huffman encode symbols
	doHuffmanEncoding(Y_DC_encoder, Y_AC_encoder, C_DC_encoder, C_AC_encoder);

	// jpeg needs zigzag sorted quantization table
	auto zigzag_qtable_y = zigzag<Byte>(qtable_y);
	auto zigzag_qtable_c = zigzag<Byte>(qtable_c);

	std::ofstream jpeg(file, std::ios::binary);

	using namespace Segment;
	jpeg << sSOI()
		<< sAPP0()
		<< sDQT().pushQuantizationTable(zigzag_qtable_y, ComponentSetup::QuantizationTableID::Zero)
		<< sDQT().pushQuantizationTable(zigzag_qtable_c, ComponentSetup::QuantizationTableID::One)
		<< sSOF0()
		.setImageSizeX(real_width)
		.setImageSizeY(real_height)
		.setupY(ComponentSetup::NoSubSampling, ComponentSetup::QuantizationTableID::Zero)
		.setupCb(ComponentSetup::Half, ComponentSetup::QuantizationTableID::One)
		.setupCr(ComponentSetup::Half, ComponentSetup::QuantizationTableID::One)
		<< sDHT().pushCodeData(Y_DC_Huffman_Table, sDHT::DC, sDHT::First)
		<< sDHT().pushCodeData(Y_AC_Huffman_Table, sDHT::AC, sDHT::First)
		<< sDHT().pushCodeData(C_DC_Huffman_Table, sDHT::DC, sDHT::Second)
		<< sDHT().pushCodeData(C_AC_Huffman_Table, sDHT::AC, sDHT::Second)
		<< sSOS()
		.setupY(sDHT::First, sDHT::First) // DC, AC
		.setupCb(sDHT::Second, sDHT::Second)
		.setupCr(sDHT::Second, sDHT::Second)
		;


	Bitstream stream;

	for (auto i = 0U; i < BitstreamCb.size1(); ++i) {
		for (auto j = 0U; j < BitstreamCb.size2(); ++j) {
			stream << BitstreamY(2 * i, 2 * j);
			stream << BitstreamY(2 * i, 2 * j + 1);
			stream << BitstreamY(2 * i + 1, 2 * j);
			stream << BitstreamY(2 * i + 1, 2 * j + 1);

			stream << BitstreamCb(i, j) << BitstreamCr(i, j);
		}
	}

	stream.fill();
	jpeg << stream;
	jpeg << sEOI();

	auto end = high_resolution_clock::now();
	std::cout << "Encoding duration: " << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
}
