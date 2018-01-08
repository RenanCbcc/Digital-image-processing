#pragma once

#include <vector>
#include <cassert>
#include <utility>

#include <boost/numeric/ublas/matrix.hpp>

#include "BitstreamGeneric.hpp"

using boost::numeric::ublas::matrix;
typedef double PixelDataType;
using mat = matrix<PixelDataType>;
typedef unsigned int uint;
typedef uint8_t Byte;

template <typename T>
matrix<T> from_vector(const std::vector<T>& v) {
	assert(v.size() == 64);

	matrix<T> m(8, 8);
	for (size_t i = 0; i < m.size1(); i++) {
		for (size_t j = 0; j < m.size2(); j++) {
			m(i, j) = v[i*m.size2() + j];
		}
	}
	return m;
}

template <typename T>
std::vector<T> zigzag(matrix<T> m) {
	const auto lookup = from_vector<uint>({
		0, 1, 5, 6, 14, 15, 27, 28,
		2, 4, 7, 13, 16, 26, 29, 42,
		3, 8, 12, 17, 25, 30, 41, 43,
		9, 11, 18, 24, 31, 40, 44, 53,
		10, 19, 23, 32, 39, 45, 52, 54,
		20, 22, 33, 38, 46, 51, 55, 60,
		21, 34, 37, 47, 50, 56, 59, 61,
		35, 36, 48, 49, 57, 58, 62, 63
	});

	assert(m.size1() == 8);
	assert(m.size2() == 8);
	std::vector<T> r;
	r.resize(64);

	for (uint i = 0; i < 8; ++i) {
		for (uint j = 0; j < 8; ++j)
			r[lookup(i, j)] = m(i, j);
	}

	return r;
};

// takes an 8x8 block and the index and returns the matrix index in zigzag order
inline int zigzag(int i) {
	assert(i < 64);
	//static const auto lookup = std::vector<uint>({
	//    0, 1, 5, 6, 14, 15, 27, 28,
	//    2, 4, 7, 13, 16, 26, 29, 42,
	//    3, 8, 12, 17, 25, 30, 41, 43,
	//    9, 11, 18, 24, 31, 40, 44, 53,
	//    10, 19, 23, 32, 39, 45, 52, 54,
	//    20, 22, 33, 38, 46, 51, 55, 60,
	//    21, 34, 37, 47, 50, 56, 59, 61,
	//    35, 36, 48, 49, 57, 58, 62, 63
	//});
	const auto lookup = std::vector<uint>({
		0, 1, 1 * 8, 2 * 8, // 0 1 2 3
		1 * 8 + 1, 2, 3, 1 * 8 + 2, 2 * 8 + 1, 3 * 8, 4 * 8, // 4 - 10
		3 * 8 + 1, 2 * 8 + 2, 1 * 8 + 3, 4, 5, 1 * 8 + 4, 2 * 8 + 3, 3 * 8 + 2, 4 * 8 + 1, 5 * 8, 6 * 8, // 11 - 21
		5 * 8 + 1, 4 * 8 + 2, 3 * 8 + 3, 2 * 8 + 4, 1 * 8 + 5, 6, 7, 1 * 8 + 6, 2 * 8 + 5, 3 * 8 + 4, 4 * 8 + 3, 5 * 8 + 2, 6 * 8 + 1, 7 * 8, 7 * 8 + 1, // 22 - 36
		6 * 8 + 2, 5 * 8 + 3, 4 * 8 + 4, 3 * 8 + 5, 2 * 8 + 6, 1 * 8 + 7, 2 * 8 + 7, 3 * 8 + 6, 4 * 8 + 5, 5 * 8 + 4, 6 * 8 + 3, 7 * 8 + 2, 7 * 8 + 3, // 37 - 42 - 49
		6 * 8 + 4, 5 * 8 + 5, 4 * 8 + 6, 3 * 8 + 7, 4 * 8 + 7, 5 * 8 + 6, 6 * 8 + 5, 7 * 8 + 4, 7 * 8 + 5, // 50 - 53 - 58
		6 * 8 + 6, 5 * 8 + 7, 6 * 8 + 7, 7 * 8 + 6, 7 * 8 + 7 // 59 - 60 - 63
	});
	assert(lookup.size() == 64);

	return lookup[i];
};


inline matrix<int> quantize(const mat& m, const mat& table) {
	assert(m.size1() == 8);
	assert(m.size2() == 8);
	assert(table.size1() == 8);
	assert(table.size2() == 8);

	matrix<int> result(8, 8);

	for (uint i = 0; i < 64; ++i) {
		result.data()[i] = static_cast<int>(std::round(m.data()[i] / table.data()[i]));
	}

	return result;
}

struct RLE_PAIR {
	unsigned short num_zeros_before : 4;
	int value;

	RLE_PAIR() = default;
	RLE_PAIR(short zeros, int _value) : num_zeros_before(zeros), value(_value) { assert(zeros < 16); }
};

inline bool operator==(const RLE_PAIR &left, const RLE_PAIR &right) {
	return (left.num_zeros_before == right.num_zeros_before) && (left.value == right.value);
}

// takes the full, zigzag sorted value list with DC and AC data, first entry is DC component and will be ignored
inline std::vector<RLE_PAIR> RLE_AC(const std::vector<int> &data) {
	assert(data.size() > 1);

	std::vector<RLE_PAIR> AC_rle;

	AC_rle.push_back(RLE_PAIR(0, data[0]));

	unsigned int zero_counter = 0;
	for (auto it = begin(data) + 1; it != end(data); ++it) {
		const auto& value = *it;

		if (value == 0) {
			++zero_counter;
			continue;
		}
		else {
			if (zero_counter > 15) {
				do {
					AC_rle.push_back(RLE_PAIR(15, 0));
					zero_counter -= 16;
				} while (zero_counter > 15);
			}
			AC_rle.push_back(RLE_PAIR(zero_counter, value));
			zero_counter = 0;
		}
	}

	// EOB
	if (zero_counter > 0)
		AC_rle.push_back(RLE_PAIR(0, 0));

	return AC_rle;
}

// takes an 8x8 quantized DCT block
// and does an RLE on the zigzag sorted values
inline std::vector<RLE_PAIR> RLE_AC(const matrix<int> &data) {
	assert(data.size1() == data.size2());
	assert(data.size1() == 8);

	std::vector<RLE_PAIR> AC_rle;
	// dc part, run length is always 0
	AC_rle.push_back(RLE_PAIR(0, data(0, 0)));


	unsigned int zero_counter = 0;
	for (int i = 1; i < 64; ++i) {
		auto zigzag_idx = zigzag(i);
		const auto& value = data.data()[zigzag_idx];

		if (value == 0) {
			++zero_counter;
			continue;
		}
		else {
			if (zero_counter > 15) {
				do {
					AC_rle.push_back(RLE_PAIR(15, 0));
					zero_counter -= 16;
				} while (zero_counter > 15);
			}
			AC_rle.push_back(RLE_PAIR(zero_counter, value));
			zero_counter = 0;
		}
	}

	// EOB
	if (zero_counter > 0)
		AC_rle.push_back(RLE_PAIR(0, 0));

	return AC_rle;
}

struct Category_Code {
	uint8_t symbol;
	Bitstream code;

	Category_Code(uint8_t p, Bitstream b) : symbol(p), code(b) {}
	~Category_Code() {}
};

inline bool operator==(const Category_Code &left, const Category_Code &right) {
	return (left.symbol == right.symbol) && (left.code == right.code);
}

inline void getCategoryAndCode(int value, short &_category, Bitstream &_code) {
	if (value == 0) {
		_category = 0;
		_code = Bitstream();
		return;
	}

	unsigned short category = 1;
	auto bound = 2L;
	while (category < 16) {
		auto upper_bound = bound - 1;
		auto lower_bound = (bound >> 1);

		auto abs_val = abs(value);
		if (abs_val >= lower_bound && abs_val <= upper_bound) {
			// got category, generate code

			auto offset = 0L;
			if (value < 0)
				offset = upper_bound - abs_val;
			else
				offset = value;

			_category = category;
			_code = Bitstream(offset, category);
			return;
		}

		bound <<= 1;
		++category;
	}

	assert(!"Shouldn't happen!");
}

inline std::pair<short, Bitstream> getCategoryAndCode(int value) {
	if (value == 0)
		return std::make_pair(0, Bitstream());

	unsigned short category = 1;
	auto bound = 2l;
	while (category < 16) {
		auto upper_bound = bound - 1;
		auto lower_bound = (bound >> 1);
		auto bound_diff = upper_bound - lower_bound;

		auto abs_val = abs(value);
		if (abs_val >= lower_bound && abs_val <= upper_bound) {
			// got category, generate code

			auto offset = 0l;
			if (value < 0)
				offset = upper_bound - abs_val;
			else
				offset = value;

			return std::make_pair(category, Bitstream(offset, category));
		}

		bound <<= 1;
		++category;
	}

	assert(!"Shouldn't happen!");
	return std::make_pair(0, Bitstream());
}

// takes the encoded list of RLE_PAIRS and generates the symbol for huffman coding and a code from category encoding 
inline std::vector<Category_Code> encode_category(const std::vector<RLE_PAIR> &data) {
	std::vector<Category_Code> category_list;
	category_list.reserve(data.size());

	for (const auto& rle_pair : data) {
		short category = 0;
		Bitstream code;
		getCategoryAndCode(rle_pair.value, category, code);

		assert(rle_pair.num_zeros_before < 16);
		assert(category < 16);
		// rle_pair.num_zeros_before == 0 for the DC value
		auto symbol = (rle_pair.num_zeros_before << 4) | category;

		category_list.emplace_back(symbol, code);
	}

	return category_list;
}