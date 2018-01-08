#pragma once

#include "BitstreamGeneric.hpp"

#include <vector>
#include <memory>
#include <unordered_map>
#include <numeric>
#include <queue>
#include <iostream>
#include <assert.h>
#include <iterator>
#include <utility>

using std::unordered_map;
using std::priority_queue;
using std::vector;
using std::pair;

// a huffman code for a symbol
struct Code {
	using CodeType = uint32_t;
	static const auto max_code_length = sizeof(CodeType) * 8;

	Code() : code(0), length(0) {}
	Code(CodeType code, uint8_t length)
		: code(code), length(length)
	{
		assert(length <= max_code_length);

		// make the leftmost bit of the code the MSB
		auto shift = max_code_length - length;
		this->code <<= shift;
	}

	Code(Bitstream bitstream) {
		length = bitstream.size();
		assert(length <= max_code_length);

		code = bitstream.extract(length, 0);
	}

	// like bistreams, the code begins at the msb
	CodeType code;
	uint8_t length;
};

using SymbolCodeMap = std::unordered_map<int, Code>;
using SymbolsPerLength = vector<vector<int>>;


// takes a text, caluclates the probability of every symbol and returns a map from symbols to huffman codes
pair<SymbolCodeMap, SymbolsPerLength> generateHuffmanCode(std::vector<int> text);

void preventOnlyOnesCode(SymbolsPerLength& symbols);
SymbolCodeMap generateCodes(const SymbolsPerLength& symbols);



// en- and decoding
Bitstream huffmanEncode(vector<int> text, SymbolCodeMap code_map);
vector<int> huffmanDecode(Bitstream bitstream, SymbolCodeMap code_map);

struct DecodeEntry {
	DecodeEntry(uint32_t code, uint8_t code_length, int symbol);

	// code starting at msb, rest filled with 1s
	uint32_t code;
	uint8_t code_length;
	int symbol;

	bool operator<(const DecodeEntry& other) {
		return code < other.code;
	}
};

struct Symbol {
	int symbol, frequency;

	Symbol()
		: symbol{},
		frequency{}
	{}
	Symbol(int _symbol, int _freq)
		: symbol(_symbol),
		frequency(_freq)
	{}
};

inline bool operator<(const Symbol& lhs, const Symbol& rhs) {
	return lhs.symbol < rhs.symbol;
}

struct Package {
	int weight;
	vector<Symbol> symbols;

	Package(const Symbol& _symbol) {
		symbols.push_back(_symbol);
		weight = _symbol.frequency;
	}

	Package(const Package& p1, const Package& p2) {
		weight = p1.weight + p2.weight;
		std::merge(begin(p1.symbols), end(p1.symbols),
			begin(p2.symbols), end(p2.symbols),
			back_inserter(symbols));
	}
};

// input: list of symbols (with its frequency)
//        maximum code length
// output: a code length for each symbol
inline SymbolsPerLength package_merge(vector<Symbol> symbols, int length_limit) {
	assert(symbols.size() <= pow(2, length_limit));

	auto comp = [](const Package& lhs, const Package& rhs) {
		return lhs.weight > rhs.weight;
	};
	using Level = priority_queue<Package, vector<Package>, decltype(comp)>;

	// blueprint level
	Level level;
	for (const auto& sym : symbols)
		level.push(sym);

	// allocate number of levels 2^-1 up to 2^-length_limit
	// levels[0] should be 2^-length_limit and levels[length_limit] should be 2^0
	vector<Level> levels;
	levels.reserve(length_limit);
	for (auto i = 0; i < length_limit; ++i)
		levels.push_back(level);
	levels.push_back(Level()); // last list is empty (level 2^0 or 1)

	for (auto i = 0; i < length_limit; ++i) {
		auto& level = levels[i];
		auto& next_level = levels[i + 1];

		// package the least 2 packages if there are more than 2
		while (level.size() > 1) {
			auto p1 = level.top();
			level.pop();
			auto p2 = level.top();
			level.pop();

			// package & merge
			next_level.push(Package{ p1, p2 });
		}
	}

	// count the occurences of the symbols in the remaining packages
	// the count is the code length for that symbol
	unordered_map<int, int> code_lengths;
	auto& final_level = levels[length_limit];
	while (final_level.size()) {
		const auto package = final_level.top();
		final_level.pop();

		for (const auto& sym : package.symbols)
			code_lengths[sym.symbol]++;
	}

	// put it in a vector for easy usage
	// +2 to allow for easy insertion of the 111..1 code one level deeper
	// (guess we should refactor that stupid data structure)
	SymbolsPerLength symbolsByCodeLength(length_limit + 2);
	for (auto it = code_lengths.begin(); it != code_lengths.end(); ++it) {
		// second: code length
		// first: symbol
		symbolsByCodeLength[it->second].push_back(it->first);
	}

	return symbolsByCodeLength;
}