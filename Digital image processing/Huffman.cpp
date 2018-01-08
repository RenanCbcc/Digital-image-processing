#include "Huffman.h"


Bitstream Huffman::huffmanEncode(vector<int> text, SymbolCodeMap code_map)
{
	Bitstream result;
	for (auto symbol : text) {
		const Code& code = code_map[symbol];
		result.push_back(code.code, code.length);
	}
	return result;
}

vector<int> Huffman::huffmanDecode(Bitstream bitstream, SymbolCodeMap code_map)
{
	vector<DecodeEntry> code_table;
	code_table.reserve(code_map.size());

	uint8_t max_code_length = 0;

	// prepare code_table and calc max_code_length
	for (const std::pair<int, Code>& pair : code_map) {
		uint8_t code_length = pair.second.length;
		if (code_length > max_code_length)
			max_code_length = code_length;

		code_table.emplace_back(pair.second.code, code_length, pair.first);
	}
	assert(max_code_length <= 32);

	// sort accoring to code
	std::sort(code_table.begin(), code_table.end());


	vector<int> decoded_text;
	uint8_t to_extract = max_code_length;

	// current pos in bitstream
	unsigned int pos = 0;

	while (pos < bitstream.size()) {
		// when we are at the end, we have to extract less
		if (pos + max_code_length >= bitstream.size())
			to_extract = bitstream.size() - pos;

		uint32_t bits = bitstream.extract(to_extract, pos);
		bits = fillRestWithOnes(bits, max_code_length);

		// find matching symbol for extracted code
		// todo: binary search
		// linear search for now...
		int index = -1;
		for (int i = 0; i < code_table.size(); i++) {
			// first entry that is greater than the extracted bits
			if (code_table[i].code >= bits) {
				index = i;
				break;
			}
		}
		// no entry found
		// some better error handling?
		assert(index != -1);

		decoded_text.push_back(code_table[index].symbol);

		// next position to extract from
		pos += code_table[index].code_length;
	}

	return decoded_text;
}

pair<SymbolCodeMap, SymbolsPerLength> Huffman::generateHuffmanCode(std::vector<int> text)
{
	assert(text.size() > 0);

	unordered_map<int, int> symbol_counts;
	for (auto& symbol : text) {
		++symbol_counts[symbol];
	}

	vector<Symbol> symbol_frequency;
	for (auto it = symbol_counts.begin(); it != symbol_counts.end(); ++it) {
		symbol_frequency.push_back(Symbol(it->first, it->second));
	}

	// special case when we only have one type of symbol
	if (symbol_counts.size() == 1) {
		SymbolCodeMap code_map;
		code_map.emplace(text[0], Code(Bitstream({ 0 })));

		SymbolsPerLength symbols(17);
		symbols[1] = { text[0] };

		return std::make_pair(code_map, symbols);
	}

	// list of symbols grouped by code length
	// symbols_for_length = symbols[code_length]
	SymbolsPerLength symbols = package_merge(symbol_frequency, 15);
	preventOnlyOnesCode(symbols);

	SymbolCodeMap code_map = generateCodes(symbols);

	return std::make_pair(code_map, symbols);
}
DecodeEntry::DecodeEntry(uint32_t code, uint8_t code_length, int symbol) {
	this->symbol = symbol;
	this->code_length = code_length;
	this->code = Huffman::fillRestWithOnes(code, code_length);
}

void Huffman::preventOnlyOnesCode(SymbolsPerLength & symbols)
{
	assert(symbols.back().empty());
	

	// prevents a code consiting of only ones by putting that one one level deeper
	auto it = std::find_if(symbols.rbegin(), symbols.rend(), [](const vector<int>& x) {return !x.empty(); });
	assert(it != symbols.rbegin());

	int only_ones_symbol = it->back();
	it->pop_back();
	--it;
	it->push_back({ only_ones_symbol });
}

SymbolCodeMap Huffman::generateCodes(const SymbolsPerLength& symbols)
{
	// based on the algorithm in
	// Reza Hashemian: Memory Efficient and High-speed Search Huffman Coding, 1995
	
	SymbolCodeMap code_map;

	uint32_t code = 0;
	for (int length = 1; length < symbols.size(); length++) {
		for (int symbol : symbols[length]) {
			code_map[symbol] = Code(code, length);
			++code;
		}
		code <<= 1;
	}

	return code_map;
}

uint32_t Huffman::fillRestWithOnes(uint32_t  code, uint8_t code_length)
{
	return code | ((1 << (32 - code_length)) - 1);
}

