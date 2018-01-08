#pragma once

#include <cassert>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "Image.h"

using boost::numeric::ublas::matrix;
using boost::numeric::ublas::matrix_range;
using boost::numeric::ublas::range;
using boost::numeric::ublas::trans;
using boost::numeric::ublas::prod;
using boost::math::constants::pi;
using boost::math::constants::root_two;
using std::cos;


const PixelDataType _pi = pi<PixelDataType>();
const PixelDataType c1 = cos(1 * _pi / 16);
const PixelDataType c2 = cos(2 * _pi / 16);
const PixelDataType c3 = cos(3 * _pi / 16);
const PixelDataType c4 = cos(4 * _pi / 16);
const PixelDataType c5 = cos(5 * _pi / 16);
const PixelDataType c6 = cos(6 * _pi / 16);
const PixelDataType c7 = cos(7 * _pi / 16);

const PixelDataType a1 = c4;
const PixelDataType a2 = c2 - c6;
const PixelDataType a3 = c4;
const PixelDataType a4 = c6 + c2;
const PixelDataType a5 = c6;

const PixelDataType s0 = 1 / (2 * root_two<PixelDataType>());
const PixelDataType s1 = 1 / (4 * c1);
const PixelDataType s2 = 1 / (4 * c2);
const PixelDataType s3 = 1 / (4 * c3);
const PixelDataType s4 = 1 / (4 * c4);
const PixelDataType s5 = 1 / (4 * c5);
const PixelDataType s6 = 1 / (4 * c6);
const PixelDataType s7 = 1 / (4 * c7);



inline void dctArai(const matrix_range<matrix<PixelDataType>>& x, matrix_range<matrix<PixelDataType>>& y) {
	assert(x.size1() == 8 && x.size2() == 8);

	matrix<PixelDataType> temp_mat(8, 8);

	for (uint j = 0; j < 8; j++) {
		auto x0 = x(0, j);
		auto x1 = x(1, j);
		auto x2 = x(2, j);
		auto x3 = x(3, j);
		auto x4 = x(4, j);
		auto x5 = x(5, j);
		auto x6 = x(6, j);
		auto x7 = x(7, j);

		auto z0 = x0 + x7;
		auto z1 = x1 + x6;
		auto z2 = x2 + x5;
		auto z3 = x3 + x4;
		auto z4 = -x4 + x3;
		auto z5 = -x5 + x2;
		auto z6 = -x6 + x1;
		auto z7 = -x7 + x0;

		auto r0 = z0 + z3;
		auto r1 = z1 + z2;
		auto r2 = z1 - z2;
		auto r3 = z0 - z3;
		auto r4 = -z4 - z5;
		auto r5 = z5 + z6;
		auto r6 = z6 + z7;
		auto r7 = z7;

		auto t0 = r0 + r1;
		auto t1 = r0 - r1;
		auto t2 = r2 + r3;
		auto t3 = r3;
		auto t4 = r4;
		auto t5 = r5;
		auto t6 = r6;
		auto t7 = r7;

		auto tmp = (t4 + t6) * a5;

		t2 *= a1;
		t4 *= a2;
		t5 *= a3;
		t6 *= a4;

		auto u0 = t0;
		auto u1 = t1;
		auto u2 = t2;
		auto u3 = t3;
		auto u4 = -t4 - tmp;
		auto u5 = t5;
		auto u6 = t6 - tmp;
		auto u7 = t7;

		auto v0 = u0;
		auto v1 = u1;
		auto v2 = u2 + u3;
		auto v3 = u3 - u2;
		auto v4 = u4;
		auto v5 = u5 + u7;
		auto v6 = u6;
		auto v7 = u7 - u5;

		auto w0 = v0;
		auto w1 = v1;
		auto w2 = v2;
		auto w3 = v3;
		auto w4 = v4 + v7;
		auto w5 = v5 + v6;
		auto w6 = -v6 + v5;
		auto w7 = v7 - v4;

		// also transposes
		temp_mat(j, 0) = w0 * s0;
		temp_mat(j, 4) = w1 * s4;
		temp_mat(j, 2) = w2 * s2;
		temp_mat(j, 6) = w3 * s6;
		temp_mat(j, 5) = w4 * s5;
		temp_mat(j, 1) = w5 * s1;
		temp_mat(j, 7) = w6 * s7;
		temp_mat(j, 3) = w7 * s3;
	}

	for (uint j = 0; j < 8; j++) {
		auto x0 = temp_mat(0, j);
		auto x1 = temp_mat(1, j);
		auto x2 = temp_mat(2, j);
		auto x3 = temp_mat(3, j);
		auto x4 = temp_mat(4, j);
		auto x5 = temp_mat(5, j);
		auto x6 = temp_mat(6, j);
		auto x7 = temp_mat(7, j);

		auto z0 = x0 + x7;
		auto z1 = x1 + x6;
		auto z2 = x2 + x5;
		auto z3 = x3 + x4;
		auto z4 = -x4 + x3;
		auto z5 = -x5 + x2;
		auto z6 = -x6 + x1;
		auto z7 = -x7 + x0;

		auto r0 = z0 + z3;
		auto r1 = z1 + z2;
		auto r2 = z1 - z2;
		auto r3 = z0 - z3;
		auto r4 = -z4 - z5;
		auto r5 = z5 + z6;
		auto r6 = z6 + z7;
		auto r7 = z7;

		auto t0 = r0 + r1;
		auto t1 = r0 - r1;
		auto t2 = r2 + r3;
		auto t3 = r3;
		auto t4 = r4;
		auto t5 = r5;
		auto t6 = r6;
		auto t7 = r7;

		auto tmp = (t4 + t6) * a5;

		t2 *= a1;
		t4 *= a2;
		t5 *= a3;
		t6 *= a4;

		auto u0 = t0;
		auto u1 = t1;
		auto u2 = t2;
		auto u3 = t3;
		auto u4 = -t4 - tmp;
		auto u5 = t5;
		auto u6 = t6 - tmp;
		auto u7 = t7;

		auto v0 = u0;
		auto v1 = u1;
		auto v2 = u2 + u3;
		auto v3 = u3 - u2;
		auto v4 = u4;
		auto v5 = u5 + u7;
		auto v6 = u6;
		auto v7 = u7 - u5;

		auto w0 = v0;
		auto w1 = v1;
		auto w2 = v2;
		auto w3 = v3;
		auto w4 = v4 + v7;
		auto w5 = v5 + v6;
		auto w6 = -v6 + v5;
		auto w7 = v7 - v4;

		// also transposes
		y(j, 0) = w0 * s0;
		y(j, 4) = w1 * s4;
		y(j, 2) = w2 * s2;
		y(j, 6) = w3 * s6;
		y(j, 5) = w4 * s5;
		y(j, 1) = w5 * s1;
		y(j, 7) = w6 * s7;
		y(j, 3) = w7 * s3;
	}
}


const auto blocksize = 8U;
// Matrix A
static const matrix<PixelDataType> A = []() {
	const auto PI = pi<PixelDataType>();
	const auto N = blocksize;
	const auto Ntimes2 = (2.*N);
	const auto scale = sqrt(2. / N);
	const auto C0 = [](unsigned int row) -> PixelDataType { return row == 0 ? 1. / root_two<PixelDataType>() : 1.; };

	matrix<PixelDataType> A(blocksize, blocksize);
	for (auto k = 0U; k < blocksize; ++k) {
		for (auto n = 0U; n < blocksize; ++n) {
			auto cos_term = (2.*n + 1.) * ((k * PI) / Ntimes2);
			A(k, n) = C0(k) * scale * cos(cos_term);
		}
	};
	return A;
}();
static const auto A_transpose = trans(A);

inline void dctDirect(const matrix_range<matrix<PixelDataType>>& X, matrix_range<matrix<PixelDataType>>& Y)
{
	assert(X.size1() == 8 && X.size2() == 8);

	// Size of blocks 8x8
	const uint N = 8;

	for (uint i = 0; i < N; ++i)
	{
		for (uint j = 0; j < N; ++j)
		{
			PixelDataType Sum = 0.0;

			// Calc the Sum
			for (uint x = 0; x < N; ++x)
			{
				for (uint y = 0; y < N; ++y)
				{
					Sum += X(y, x) * A(i, x) * A(j, y);
				}
			}
			Y(j, i) = Sum;
		}
	}
}

inline void dctMat(const matrix_range<matrix<PixelDataType>>& X, matrix_range<matrix<PixelDataType>>& Y)
{
	assert(X.size1() == blocksize && X.size2() == blocksize);
	assert(Y.size1() == blocksize && Y.size2() == blocksize);
	assert(A.size1() == blocksize && A.size2() == blocksize);
	assert(A_transpose.size1() == blocksize && A_transpose.size2() == blocksize);

	using mat = matrix<PixelDataType>;

	// Y = A * (X * A_transpose)
	mat first = prod(X, A_transpose);
	noalias(Y) = mat(prod(A, first));
}

inline matrix<PixelDataType> inverseDctMat(matrix<PixelDataType> X)
{
	const auto blocksize = 8U;
	assert(X.size1() == blocksize && X.size2() == blocksize);

	auto C = [](unsigned int row) -> PixelDataType { return row == 0 ? 1. / root_two<PixelDataType>() : 1.; };
	const auto N = blocksize;

	auto scale = sqrt(2. / N);

	// Output Matrix
	matrix<PixelDataType> A(blocksize, blocksize);
	for (auto k = 0U; k < blocksize; ++k) {
		for (auto n = 0U; n < blocksize; ++n) {
			auto co = C(k);
			auto cos_term = (2.*n + 1.) * ((k * pi<PixelDataType>()) / (2.*N));
			auto cos_val = cos(cos_term);
			A(k, n) = co * scale * cos_val;
		}
	}

	auto A_transpose = trans(A);

	// axpy_prod(A, B, C, false); // C += A * B
	matrix<PixelDataType> first = prod(A_transpose, X);
	matrix<PixelDataType> res = prod(first, A);

	return res;
}