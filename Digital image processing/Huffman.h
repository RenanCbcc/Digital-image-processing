#pragma once
#include "Huffman.hpp"

class Huffman
{
public:
	Bitstream huffmanEncode(vector<int>, SymbolCodeMap);
	vector<int> huffmanDecode(Bitstream, SymbolCodeMap);

	pair<SymbolCodeMap, SymbolsPerLength> generateHuffmanCode(std::vector<int>);
	void preventOnlyOnesCode(SymbolsPerLength&);
	SymbolCodeMap generateCodes(const SymbolsPerLength&);
	static uint32_t fillRestWithOnes(uint32_t, uint8_t);
	
};





