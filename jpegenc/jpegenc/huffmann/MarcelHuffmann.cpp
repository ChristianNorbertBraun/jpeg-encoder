//
//  MarcelHuffmann.cpp
//  jpegenc
//
//  Created by Marcel Groß on 22.11.16.
//  Copyright © 2016 FHWS. All rights reserved.
//

#include "MarcelHuffmann.hpp"

void InputWord::increase() {
	++amount;
}




void MarcelHuffmann::addToWords(std::vector<int> input) {
	for (int i = 0; i < input.size(); ++i) {
		bool notFount = true;
		for (int j = 0; j < words.size(); ++j) {
			if (input[i] == words[j].symbol) {
				words[j].increase();
				notFount = false;
				break;
			}
		}
		if (notFount) {
			words.push_back(InputWord(input[i]));
		}
		
	}
	
	sortByAppearance();
}

Node* MarcelHuffmann::generateTree() {
	Node* node = new Node();
	node->left = new Node(&words[0]);
	node->right = new Node(&words[1]);
	node->calculateValue();
	
	InputWord* currentWord;
	Node* newNode;
	for (int i = 2; i < words.size(); ++i) {
		currentWord = &words[i];
		newNode = new Node();
		newNode->depth = node->depth + 1;
		if (currentWord->amount > node->value->amount) {
			newNode->right = new Node(currentWord);
			newNode->left = node;
		} else{
			newNode->left = new Node(currentWord);
			newNode->right = node;
		}
		newNode->calculateValue();
		node = newNode;
	}
	
	return node;
}

std::map<Symbol, SymbolBits>* MarcelHuffmann::generateEncodingTable(Node* node) {
	SymbolBits bitsForSymbol;
	std::map<Symbol, SymbolBits>* map = new std::map<Symbol, SymbolBits>();
	climbTree(bitsForSymbol, node, map);
	
	return map;
}

void MarcelHuffmann::climbTree(SymbolBits bitsForSymbol, Node* node, std::map<Symbol, SymbolBits>* map) {
	Node* left = node->left;
	Node* right = node->right;
	
	if (left == nullptr && right == nullptr) {
		map->insert(std::make_pair(node->value->symbol, bitsForSymbol));
	} else{
		++bitsForSymbol.numberOfBits;
		bitsForSymbol.bits <<= 1;
		climbTree(bitsForSymbol, left, map);
		bitsForSymbol.bits |= 1;
		climbTree(bitsForSymbol, right, map);
	}
}

std::vector<Symbol> MarcelHuffmann::decode(Bitstream* bitstream, Node* rootNode) {
	std::vector<Symbol> symbols;
	size_t numberOfBits = bitstream->numberOfBits();
	Node* node = rootNode;
	bool bit;
	for (int i = 0; i < numberOfBits; ++i) {
	
		if (node->left == nullptr && node->right == nullptr) {
			symbols.push_back(node->value->symbol);
			node = rootNode;
		}
		
		bit = bitstream->read(i);
		if (bit) {
			node = node->right;
		} else {
			node = node->left;
		}
	}
	
	return symbols;
}


void MarcelHuffmann::sortByAppearance() {
	std::sort(words.begin(), words.end());
}




void Node::calculateValue() {
	Node* left = this->left;
	Node* right = this->right;
	
	if ((left != nullptr) && (right != nullptr)) {
		this->value = new InputWord(left->value->amount + right->value->amount, NULL);
	}
	
}

void Node::print() {
	
}