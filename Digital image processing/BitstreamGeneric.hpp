#pragma once

#include <cstdint>
#include <vector>
#include <initializer_list>
#include <ostream>
#include <istream>
#include <cassert>

// Bitstream_Generic[0] is the least significant bit (LSB)
// Bitsream[Bitstream_Generic.size()-1] is the most significant bit (MSB)
template <typename BlockType>
class Bitstream_Generic
{
public:
	friend class BitView;
	static const auto block_size = sizeof(BlockType) * 8;

	typedef std::vector<BlockType> ContainerType;

	class BitView
	{

		friend class Bitstream_Generic<BlockType>;
		BitView(ContainerType& bstream, unsigned int block_idx, uint8_t bit_idx)
			: blocks(bstream),
			block_index(block_idx),
			bit_index(bit_idx)
		{}

		ContainerType& blocks;
		const unsigned int block_index;
		const uint8_t bit_index;

	public:
		operator bool() const
		{
			return (blocks[block_index] & ((BlockType)1 << bit_index)) != 0;
		}

		void operator=(bool val)
		{
			auto& block = blocks[block_index];
			if (val) {
				// set bit at bit_idx
				auto mask = (BlockType)1 << bit_index;
				block |= mask;
			}
			else {
				// unset bit at bit_idx
				auto mask = ~((BlockType)1 << bit_index);
				block &= mask;
			}
		}
	};

private:
	ContainerType blocks;
	uint8_t bit_idx;
	unsigned int sz; // up to 4 billion bits should be enough

private:
	// appending Blocks, don't let the user see this as it needs proper alignment,
	// only used for reading in a file 
	Bitstream_Generic& operator<<(BlockType val);

public:
	// ctors
	Bitstream_Generic();
	Bitstream_Generic(std::initializer_list<bool> args);
	Bitstream_Generic(uint32_t data, int number_of_bits);

	// appending bits / bitstreams
	Bitstream_Generic& operator<<(bool val);
	Bitstream_Generic& operator<<(std::initializer_list<bool> args);
	Bitstream_Generic& operator<<(Bitstream_Generic<BlockType>& stream);
	Bitstream_Generic& push_back(bool val);
	Bitstream_Generic& push_back(uint32_t data, int number_of_bits);
	Bitstream_Generic& push_back_LSB_mode(uint32_t data, int number_of_bits);

	// streaming
	template<typename BlockType>
	friend std::ostream& operator<<(std::ostream& out, const Bitstream_Generic<BlockType>& bitstream);
	template<typename BlockType>
	friend std::istream& operator>>(std::istream& in, Bitstream_Generic<BlockType>& bitstream);

	// accessing
	BitView operator[](unsigned int pos);

	// bit at from_position will be the MSB
	template<typename T>
	T extractT(uint8_t number_of_bits, size_t from_position);
	uint32_t extract(uint8_t number_of_bits, size_t from_position) { return extractT<uint32_t>(number_of_bits, from_position); }

	// others
	unsigned int size() const;
	void fill(); // fill remaining bits in the last block with 1s
	template<typename BlockType>
	friend bool operator==(const Bitstream_Generic<BlockType>& lhs, const Bitstream_Generic<BlockType>& rhs);
};

//
// Implementation
//
template<typename BlockType>
Bitstream_Generic<BlockType>::Bitstream_Generic()
	: bit_idx(block_size),
	blocks(),
	sz(0)
{}

template<typename BlockType>
Bitstream_Generic<BlockType>::Bitstream_Generic(std::initializer_list<bool> args)
	: Bitstream_Generic()
{
	*this << args;
}

template<typename BlockType>
Bitstream_Generic<BlockType>::Bitstream_Generic(uint32_t data, int number_of_bits)
	: Bitstream_Generic()
{
	push_back_LSB_mode(data, number_of_bits);
}

template<typename BlockType>
Bitstream_Generic<BlockType>& Bitstream_Generic<BlockType>::operator<<(bool val)
{
	// make bit stream longer if we run out of space
	if (bit_idx >= block_size) {
		blocks.push_back(0);
		bit_idx = block_size - 1;
	}

	auto& block = blocks.back();
	if (val) {
		// set bit at bit_idx
		auto mask = (BlockType)1 << bit_idx;
		block |= mask;
	}

	--bit_idx; // wraps around to 255 if it was 0
	++sz;

	return *this;
}

template<typename BlockType>
Bitstream_Generic<BlockType>& Bitstream_Generic<BlockType>::operator<<(BlockType val)
{
	blocks.push_back(val);
	sz += block_size;
	return *this;
}

template<typename BlockType>
Bitstream_Generic<BlockType>& Bitstream_Generic<BlockType>::operator<<(std::initializer_list<bool> args)
{
	// todo: perf improve: use operator<<(bool) as long as the block is not aligned
	//          then use operator<<(BlockType) if enough bits are available
	for (auto& bit : args)
		*this << bit;
	return *this;
}

template<typename BlockType>
Bitstream_Generic<BlockType>& Bitstream_Generic<BlockType>::operator<<(Bitstream_Generic<BlockType>& stream)
{
	for (auto i = 0U; i < stream.size(); ++i) {
		*this << stream[i];
	}
	return *this;
}

template<typename BlockType>
Bitstream_Generic<BlockType>& Bitstream_Generic<BlockType>::push_back(bool val)
{
	*this << val;
	return *this;
}

template<typename BlockType>
Bitstream_Generic<BlockType>& Bitstream_Generic<BlockType>::push_back(uint32_t data, int number_of_bits)
{
	// append number_of_bits from MSB
	unsigned int mask = 1 << 31;
	while (number_of_bits > 0) {
		bool bit = (data & mask) != 0;
		*this << bit;

		--number_of_bits;
		mask >>= 1;
	}
	return *this;
}

template<typename BlockType>
Bitstream_Generic<BlockType>& Bitstream_Generic<BlockType>::push_back_LSB_mode(uint32_t data, int number_of_bits)
{
	// append number_of_bits from LSB
	unsigned int mask = 1 << (number_of_bits - 1);
	while (number_of_bits > 0) {
		bool bit = (data & mask) != 0;
		*this << bit;

		--number_of_bits;
		mask >>= 1;
	}
	return *this;
}


template<typename BlockType>
std::ostream& operator<<(std::ostream& out, const Bitstream_Generic<BlockType>& bitstream)
{
	if (bitstream.sz > 0) {
		for (const auto& block : bitstream.blocks) {
			out.write((const char*)&block, sizeof(block));
			if (block == 0xFF)
				out.put(0x00);
		}
	}
	return out;
}

template<typename BlockType>
std::istream& operator>>(std::istream& in, Bitstream_Generic<BlockType>& bitstream)
{
	BlockType b = 0;
	// only reads sizeof(b) bytes, if the istream was not aligned to that, the rest is ignored
	while (in.read((char*)&b, sizeof(b)))
		bitstream << b;
	return in;
}

template<typename BlockType>
unsigned int Bitstream_Generic<BlockType>::size() const
{
	return sz;
}

template<typename BlockType>
void Bitstream_Generic<BlockType>::fill()
{
	// set all remaining bits in the last block to 1
	while (bit_idx <= block_size)
		*this << true;
}

template<typename BlockType>
bool operator==(const Bitstream_Generic<BlockType>& lhs, const Bitstream_Generic<BlockType>& rhs) {
	return lhs.blocks == rhs.blocks && lhs.size() == rhs.size();
}

template<typename BlockType>
typename Bitstream_Generic<BlockType>::BitView Bitstream_Generic<BlockType>::operator[](unsigned int pos) {
	assert(pos < size());
	auto block_idx = static_cast<unsigned int>(pos / block_size);
	auto bit_idx = static_cast<uint8_t>(block_size - (pos - block_idx * block_size) - 1);

	return BitView(blocks, block_idx, bit_idx);
}

template<typename BlockType>
template<typename T>
T Bitstream_Generic<BlockType>::extractT(uint8_t number_of_bits, size_t from_position)
{
	assert(number_of_bits <= sizeof(T) * 8);
	assert(from_position + number_of_bits - 1 < size());

	auto num_bits = number_of_bits;
	T result = 0;
	auto block_idx = from_position / block_size;

	while (number_of_bits > 0) {
		auto& block = blocks[block_idx];

		uint8_t bit_idx_from = block_size - (from_position % block_size) - 1;

		uint8_t bit_idx_to;
		if (number_of_bits > bit_idx_from)
			bit_idx_to = 0;
		else
			bit_idx_to = bit_idx_from + 1 - number_of_bits;

		// build mask, make sure that calculations aren't narrowed to int or whatever (we can have up to 8byte Blocktypes!)
		BlockType mask = ((BlockType)1 << (bit_idx_from + 1)) - 1;
		mask -= ((BlockType)1 << (bit_idx_to)) - 1;

		// extract bits
		result |= ((block & mask) >> (bit_idx_to));

		number_of_bits -= bit_idx_from + 1 - bit_idx_to;
		from_position = 0; // always start at the beginning of the next block
		++block_idx;
		if (number_of_bits > 0)
			result <<= (number_of_bits > block_size ? block_size : number_of_bits);
	}

	// make the first extracted bit the MSB
	auto shift = sizeof(T) * 8 - num_bits;
	result <<= shift;

	return result;
}

// Convenience definitions (new type alias format)
using Bitstream8 = Bitstream_Generic<uint8_t>;
using Bitstream16 = Bitstream_Generic<uint16_t>;
using Bitstream32 = Bitstream_Generic<uint32_t>;
using Bitstream64 = Bitstream_Generic<uint64_t>;
using Bitstream = Bitstream8;

typedef std::initializer_list<bool> Bits;