//
//  JPEGImage.cpp
//  jpegenc
//
//  Created by Christian Braun on 16/11/16.
//  Copyright © 2016 FHWS. All rights reserved.
//

#include "JPEGSegments.hpp"
#include <iostream>
#include "../dct/Arai.hpp"
#include "ImageDataEncoding.hpp"

using namespace JPEGSegments;

void APP0::addToStream(Bitstream &stream) {
	stream.add(type, 16);
	stream.add(length, 16);
	stream.add(JFIF[0], 8);
	stream.add(JFIF[1], 8);
	stream.add(JFIF[2], 8);
	stream.add(JFIF[3], 8);
	stream.add(JFIF[4], 8);
	stream.add(MAJOR_REVISION_NUMBER, 8);
	stream.add(MINOR_REVISION_NUMBER, 8);
	stream.add(PIXEL_UNIT, 8);
	stream.add(X_DENSITY, 16);
	stream.add(Y_DENSITY, 16);
	stream.add(PREVIEW_WIDTH, 8);
	stream.add(PREVIEW_HEIGHT, 8);
}

void StartOfFrame0::addToStream(Bitstream &stream) {
	stream.add(type, 16);
	stream.add(length, 16);
	stream.add(precision, 8);
	stream.add(height, 16);
	stream.add(width, 16);
	stream.add(numberOfComponents, 8);
	while (numberOfComponents--) {
		stream.add(0x01, 8);    // ID (Y)
		stream.add(0x22, 8);    // Subsampling
		stream.add(0x00, 8);    // Quantisierungstabelle
	}
}

void StartOfImage::addToStream(Bitstream &stream) {
    stream.add(type, 16);
}

void EndOfImage::addToStream(Bitstream &stream) {
    stream.add(type, 16);
}

void DefineHuffmanTable::addToStream(Bitstream &stream) {
    stream.add(type, 16);
    stream.add(length, 16);
    
    addTableData(0, 0, Y_DC, stream);
    addTableData(1, 1, Y_AC, stream);
    addTableData(2, 0, CbCr_DC, stream);
    addTableData(3, 1, CbCr_AC, stream);
}

void DefineHuffmanTable::addTableData(uint8_t htNumber, uint8_t htType, EncodingTable table, Bitstream &stream) {
    // HT Information
    stream.add(0, 3); // rest
    stream.add(htType, 1);
    stream.add(htNumber, 4);

    // number of symbols per level
    unsigned char symbolsPerLevel[16] = {0};

    for (const std::pair<Symbol, Encoding> &encoding : table) {
        unsigned short numberOfBits = encoding.second.numberOfBits;
        symbolsPerLevel[numberOfBits - 1] += 1;
    }

    for (int i = 0; i < 16; ++i) {
        stream.add(symbolsPerLevel[i], 8);
    }

    for (const std::pair<Symbol, Encoding> &encoding : table) {
        stream.add(encoding.first, 8);
    }
}

uint8_t* DefineQuantizationTable::sortZickZack(const uint8_t* table) {
    uint8_t *sortedTable = new uint8_t[64];
    
    for (int i = 0; i < 64; ++i) {
        sortedTable[i] = table[ZICK_ZACK_INDEXES[i]];
    }
    return sortedTable;
}

void DefineQuantizationTable::addToStream(Bitstream &stream) {
    stream.add(type, 16);
    stream.add(length, 16);
    
    stream.add(precision, 4);
    stream.add(qt_number, 4);
    
    uint8_t* sortedTable = sortZickZack(table);
    
    for (int i = 0; i < 64; ++i) {
        stream.add(sortedTable[i], 8);
    }
}

void StartOfScan::addToStream(Bitstream &stream) {
	stream.add(type, 16);
	stream.add(length, 16);
	stream.add(numberOfComponents, 8);
    
    stream.add(1, 8); // Y
    stream.add(0, 4); // Y_DC
    stream.add(1, 4); // Y_AC
    
    stream.add(2, 8); // Cb
    stream.add(2, 4); // CbCr_DC
    stream.add(3, 4); // CbCr_AC
    
    stream.add(3, 8); // Cr
    stream.add(2, 4); // CbCr_DC
    stream.add(3, 4); // CbCr_AC
    
    stream.add(0x00, 8);
    stream.add(0x3f, 8);
    stream.add(0x00, 8);
    
    // TODO: Entropy coded segment data
}

void JPEGWriter::writeJPEGImage(const char *pathToFile) {    
    StartOfImage* soi = new StartOfImage();
    segments.push_back(soi);
    
    APP0* app0 = new APP0();
    segments.push_back(app0);

    DefineQuantizationTable* Y_dqt = new DefineQuantizationTable(0, 0, luminanceQT);
    segments.push_back(Y_dqt);
    
    DefineQuantizationTable* CbCr_dqt = new DefineQuantizationTable(1, 0, chrominanceQT);
    segments.push_back(CbCr_dqt);
    
    StartOfFrame0* sof0 = new StartOfFrame0(1, image);  // 1 = numberOfComponents
    segments.push_back(sof0);
    
    ChannelData* channelData= new ChannelData(image);
    EncodedImageData encodedImageData(channelData);
    encodedImageData.initialize();
    DefineHuffmanTable* dht = new DefineHuffmanTable(encodedImageData.Y_DC, encodedImageData.Y_AC, encodedImageData.CbCr_DC, encodedImageData.CbCr_AC);
    segments.push_back(dht);
    
    // sos
    
    EndOfImage* eoi = new EndOfImage();
    segments.push_back(eoi);
    
    for (size_t i = 0; i < segments.size(); ++i) {
        segments[i]->addToStream(stream);
    }
    stream.saveToFile(pathToFile);
}

void EncodedImageData::generateYDataAndHT()
{
    ImageDataEncoding channel1Encoder(channelData->channel1->values, (unsigned int) Y_width, (unsigned int) Y_height);
    channel1Encoder.init();

    Y_DC = channel1Encoder.generateDCEncodingTable(Y_DC_encoding);
    Y_AC = channel1Encoder.generateACEncodingTable(Y_AC_byteReps, Y_AC_encoding);
}

void EncodedImageData::generateCbCrDataAndHT()
{
    ImageDataEncoding channel2Encoder(channelData->channel2->values, (unsigned int) CbCr_width, (unsigned int) CbCr_height);
    ImageDataEncoding channel3Encoder(channelData->channel3->values, (unsigned int) CbCr_width, (unsigned int) CbCr_height);
    channel2Encoder.init();
    channel3Encoder.init();
    
    Cb_DC_encoding = channel2Encoder.differenceEncoding();
    Cr_DC_encoding = channel3Encoder.differenceEncoding();
    channel2Encoder.runLengthEncoding(Cb_AC_byteReps, Cb_AC_encoding);
    channel3Encoder.runLengthEncoding(Cr_AC_byteReps, Cr_AC_encoding);
    
    generateCbCrHT();
}

void EncodedImageData::initialize() {
    Arai::transform(channelData->channel1->values, Y_width, Y_height);
    Arai::transform(channelData->channel2->values, CbCr_width, CbCr_height);
    Arai::transform(channelData->channel3->values, CbCr_width, CbCr_height);

    Quantization::run(channelData->channel1->values, Y_width, Y_height, luminanceQT);
    Quantization::run(channelData->channel2->values, CbCr_width, CbCr_height, chrominanceQT);
    Quantization::run(channelData->channel3->values, CbCr_width, CbCr_height, chrominanceQT);
}

void EncodedImageData::generateCbCrHT() {
    generateCbCrHT_DC();
    generateCbCrHT_AC();
}

void EncodedImageData::generateCbCrHT_DC() {
    Huffman DC_huffman;
    
    for (auto &current : Cb_DC_encoding) {
        DC_huffman.addSymbol(current.numberOfBits);
    }
    
    for (auto &current : Cr_DC_encoding) {
        DC_huffman.addSymbol(current.numberOfBits);
    }
    
    DC_huffman.generateNodeList();
    
    CbCr_DC = DC_huffman.canonicalEncoding(16);
}

void EncodedImageData::generateCbCrHT_AC() {
    Huffman AC_huffman;
    
    for (auto &current : Cb_AC_byteReps)
    {
        AC_huffman.addSymbol(current);
    }
    
    for (auto &current : Cr_AC_byteReps)
    {
        AC_huffman.addSymbol(current);
    }
    AC_huffman.generateNodeList();
    CbCr_AC = AC_huffman.canonicalEncoding(16);    
}
