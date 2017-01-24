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

void StartOfImage::addToStream(Bitstream &stream) {
    stream.add(type, 16);
}

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
    
    // Y
    stream.add(1, 8);    // ID
    stream.add(useSubSampling ? 0x22 : 0x11, 8); // Subsampling
    stream.add(0, 8);    // QT
    
    // Cb
    stream.add(2, 8);    // ID
    stream.add(0x11, 8); // Subsampling
    stream.add(1, 8);    // QT
    
    // Cr
    stream.add(3, 8);    // ID
    stream.add(0x11, 8); // Subsampling
    stream.add(1, 8);    // QT
}

void EndOfImage::addToStream(Bitstream &stream) {
    stream.add(type, 16);
}


void StartOfScan::addToStreamNoFF(Bitstream &stream, Encoding enc) {
	short n = enc.numberOfBits - 1;

	for (; n >= 0; --n) {
		uint8_t firstBit = (enc.code >> n) & 1;
		
		if (bufferIndex < 8) {
			buffer <<= 1;
			buffer |= firstBit;
			++bufferIndex;
		} else {
			if (buffer == 0xff) {
				stream.add(0x00, 8);
			}
			
			bufferIndex = 1;
			buffer = firstBit;
		}
		
		stream.add(firstBit);
	}
}
bool comparePairs(const std::pair<Symbol, Encoding> &pair1, const std::pair<Symbol, Encoding> &pair2) {
    return pair1.second.code < pair2.second.code;
}

void DefineHuffmanTable::addToStream(Bitstream &stream) {
    stream.add(type, 16);
    stream.add(length, 16);

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
    
    std::vector<std::pair<Symbol, Encoding>> sortedEncodings;

    for (const std::pair<Symbol, Encoding> &encoding : table) {
        sortedEncodings.push_back(encoding);
    }
    std::sort(sortedEncodings.begin(), sortedEncodings.end(), comparePairs);
    
    for (auto encoding : sortedEncodings) {
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
    
    stream.add(0, 4); // precision
    stream.add(qt_number, 4); // number
    uint8_t* qt_sorted = sortZickZack(qt);
    
    for (int i = 0; i < 64; ++i) {
        stream.add(qt_sorted[i], 8);
    }
}

void StartOfScan::addToStream(Bitstream &stream) {
    stream.add(type, 16);
	stream.add(length, 16);
	stream.add(numberOfComponents, 8);
    stream.add(1, 8); // Y
    stream.add(0, 4); // Y_DC
    stream.add(0, 4); // Y_AC
    
    stream.add(2, 8); // Cb
    stream.add(1, 4); // CbCr_DC
    stream.add(1, 4); // CbCr_AC
    
    stream.add(3, 8); // Cr
    stream.add(1, 4); // CbCr_DC
    stream.add(1, 4); // CbCr_AC
    
    stream.add(0x00, 8);
    stream.add(0x3f, 8);
    stream.add(0x00, 8);
    
    if ( useSubSampling ) {

        // Prepare Y data bitstreams
        size_t numberOf_Y_Blocks = encodedImageData->Y_DC_encoding.size();
        Bitstream* y_data = new Bitstream[numberOf_Y_Blocks];
        
        int k_y = 0;
        
        for (size_t i = 0; i < numberOf_Y_Blocks; ++i) {
            // Y_DC
            auto index = encodedImageData->Y_DC_encoding[i].numberOfBits;
            addToStreamNoFF(y_data[i], encodedImageData->Y_DC.at(index));
            if (index != 0) {
                addToStreamNoFF(y_data[i], encodedImageData->Y_DC_encoding.at(i));
            }
            
            // Y_AC
            int written_AC = 0;
            for (; encodedImageData->Y_AC_byteReps[k_y] != 0; ++k_y)
            {
                index = encodedImageData->Y_AC_byteReps[k_y];
                int leadingZeros = (index & 0xF0) >> 4;
                written_AC += leadingZeros + 1;
                addToStreamNoFF(y_data[i], encodedImageData->Y_AC.at(index));
                addToStreamNoFF(y_data[i], encodedImageData->Y_AC_encoding.at(k_y));
            }
            if (written_AC < 63) {
                // Add EOB
                addToStreamNoFF(y_data[i],encodedImageData->Y_AC.at(0));
            }
            ++k_y;
        }
        
        // Write blocks in right order
        size_t numberOf_CbCr_Blocks = encodedImageData->Cb_DC_encoding.size();
        size_t numberOf_Y_BlocksPerLine = encodedImageData->Y_width / 8;
        
        for (size_t i = 0; i < numberOf_CbCr_Blocks; ++i) {
            
            // Calculate correct 4 Y indeces
            size_t outer_offset = (i / numberOf_Y_BlocksPerLine) * numberOf_Y_BlocksPerLine * numberOf_Y_BlocksPerLine;
            size_t inner_offset = (i % numberOf_Y_BlocksPerLine) * 2;
            size_t y_1 = outer_offset + inner_offset;
            size_t y_2 = y_1 + 1;
            
            size_t y_3 = y_1 + numberOf_Y_BlocksPerLine * 2;
            size_t y_4 = y_3 + 1;

            // Copy y data to final bitstream
            for (size_t m = 0; m < y_data[y_1].numberOfBits(); ++m) {
                stream.add(y_data[y_1].read(m));
            }
            for (size_t m = 0; m < y_data[y_2].numberOfBits(); ++m) {
                stream.add(y_data[y_2].read(m));
            }
            for (size_t m = 0; m < y_data[y_3].numberOfBits(); ++m) {
                stream.add(y_data[y_3].read(m));
            }
            for (size_t m = 0; m < y_data[y_4].numberOfBits(); ++m) {
                stream.add(y_data[y_4].read(m));
            }
            
            // Add cb data
            int k_cb = 0;
            
            // Cb_DC
            auto index = encodedImageData->Cb_DC_encoding[i].numberOfBits;
            addToStreamNoFF(stream, encodedImageData->CbCr_DC.at(index));
            
            if (index != 0) {
                addToStreamNoFF(stream, encodedImageData->Cb_DC_encoding.at(i));
            }
            
            // Cb_AC
            int written_AC = 0;
            
            for (; encodedImageData->Cb_AC_byteReps[k_cb] != 0; ++k_cb)
            {
                index = encodedImageData->Cb_AC_byteReps[k_cb];
                int leadingZeros = (index & 0xF0) >> 4;
                written_AC += leadingZeros + 1;
                addToStreamNoFF(stream, encodedImageData->CbCr_AC.at(index));
                addToStreamNoFF(stream, encodedImageData->Cb_AC_encoding.at(k_cb));
            }
            if (written_AC < 63) {
                addToStreamNoFF(stream, encodedImageData->CbCr_AC.at(0));
            }
            ++k_cb;

            // Add cr data
            int k_cr = 0;
            
            // Cr_DC
            index = encodedImageData->Cr_DC_encoding[i].numberOfBits;
            addToStreamNoFF(stream, encodedImageData->CbCr_DC.at(index));
            
            if (index != 0) {
                addToStreamNoFF(stream, encodedImageData->Cr_DC_encoding.at(i));
            }
            
            // Cr_AC
            written_AC = 0;
            
            for (; encodedImageData->Cr_AC_byteReps[k_cr] != 0; ++k_cr)
            {
                index = encodedImageData->Cr_AC_byteReps[k_cr];
                int leadingZeros = (index & 0xF0) >> 4;
                written_AC += leadingZeros + 1;
                addToStreamNoFF(stream, encodedImageData->CbCr_AC.at(index));
                addToStreamNoFF(stream, encodedImageData->Cr_AC_encoding.at(k_cr));
            }
            if (written_AC < 63) {
                addToStreamNoFF(stream, encodedImageData->CbCr_AC.at(0));
            }
            ++k_cr;
        }
        delete[] y_data;
    }
    else {
        size_t numberOfBlocks = encodedImageData->Y_DC_encoding.size();
        
        int k_y  = 0;
        int k_cb = 0;
        int k_cr = 0;
        
        for (int i = 0; i < numberOfBlocks; ++i) {
            
            // Y_DC
            auto index = encodedImageData->Y_DC_encoding[i].numberOfBits;
            addToStreamNoFF(stream, encodedImageData->Y_DC.at(index));
            if (index != 0) {
                addToStreamNoFF(stream, encodedImageData->Y_DC_encoding.at(i));
            }
            
            // Y_AC
            int written_AC = 0;
            for (; encodedImageData->Y_AC_byteReps[k_y] != 0; ++k_y)
            {
                index = encodedImageData->Y_AC_byteReps[k_y];
                int leadingZeros = (index & 0xF0) >> 4;
                written_AC += leadingZeros + 1;
                addToStreamNoFF(stream, encodedImageData->Y_AC.at(index));
                addToStreamNoFF(stream, encodedImageData->Y_AC_encoding.at(k_y));
            }
            if (written_AC < 63) {
                // Add EOB
                addToStreamNoFF(stream,encodedImageData->Y_AC.at(0));
            }
            ++k_y;
            
            // Cb_DC
            index = encodedImageData->Cb_DC_encoding[i].numberOfBits;
            addToStreamNoFF(stream, encodedImageData->CbCr_DC.at(index));
            if (index != 0) {
                addToStreamNoFF(stream, encodedImageData->Cb_DC_encoding.at(i));
            }
            
            // Cb_AC
            written_AC = 0;
            for (; encodedImageData->Cb_AC_byteReps[k_cb] != 0; ++k_cb)
            {
                index = encodedImageData->Cb_AC_byteReps[k_cb];
                int leadingZeros = (index & 0xF0) >> 4;
                written_AC += leadingZeros + 1;
                addToStreamNoFF(stream, encodedImageData->CbCr_AC.at(index));
                addToStreamNoFF(stream, encodedImageData->Cb_AC_encoding.at(k_cb));
            }
            if (written_AC < 63) {
                addToStreamNoFF(stream, encodedImageData->CbCr_AC.at(0));
            }
            ++k_cb;
            
            // Cr_DC
            index = encodedImageData->Cr_DC_encoding[i].numberOfBits;
            addToStreamNoFF(stream, encodedImageData->CbCr_DC.at(index));
            if (index != 0) {
                addToStreamNoFF(stream, encodedImageData->Cr_DC_encoding.at(i));
            }
            
            // Cr_AC
            written_AC = 0;
            for (; encodedImageData->Cr_AC_byteReps[k_cr] != 0; ++k_cr)
            {
                index = encodedImageData->Cr_AC_byteReps[k_cr];
                int leadingZeros = (index & 0xF0) >> 4;
                written_AC += leadingZeros + 1;
                addToStreamNoFF(stream, encodedImageData->CbCr_AC.at(index));
                addToStreamNoFF(stream, encodedImageData->Cr_AC_encoding.at(k_cr));
            }
            if (written_AC < 63) {
                addToStreamNoFF(stream, encodedImageData->CbCr_AC.at(0));
            }
            ++k_cr;
        }
    }
	auto bitsToFill = stream.numberOfBits() % 8 == 0 ? 0 : 8 - (stream.numberOfBits() % 8);
    stream.add(0xFFF, bitsToFill);
}

void JPEGWriter::writeJPEGImage(const char *pathToFile, bool useSubSampling) {
    StartOfImage* soi = new StartOfImage();
    soi->addToStream(stream);
    
    APP0* app0 = new APP0();
    app0->addToStream(stream);

    DefineQuantizationTable* y_dqt = new DefineQuantizationTable(0, luminanceQT);
    y_dqt->addToStream(stream);
    
    DefineQuantizationTable* cbcr_dqt = new DefineQuantizationTable(1, chrominanceQT);
    cbcr_dqt->addToStream(stream);
    
    StartOfFrame0* sof0 = new StartOfFrame0(3, image, useSubSampling);
    sof0->addToStream(stream);
    
    if ( useSubSampling )
    {
        image->channel2->reduceBySubSampling(2, 2);
        image->channel3->reduceBySubSampling(2, 2);
    }
    
    ChannelData* channelData= new ChannelData(image);
    
    EncodedImageData *encodedImageData = new EncodedImageData(channelData);
    encodedImageData->initialize();
    
    DefineHuffmanTable* Y_DC_dht = new DefineHuffmanTable(0, 0, encodedImageData->Y_DC);
    Y_DC_dht->addToStream(stream);
    
    DefineHuffmanTable* Y_AC_dht = new DefineHuffmanTable(0, 1, encodedImageData->Y_AC);
    Y_AC_dht->addToStream(stream);
    
    DefineHuffmanTable* CbCr_DC_dht = new DefineHuffmanTable(1, 0, encodedImageData->CbCr_DC);
    CbCr_DC_dht->addToStream(stream);

    DefineHuffmanTable* CbCr_AC_dht = new DefineHuffmanTable(1, 1, encodedImageData->CbCr_AC);
    CbCr_AC_dht->addToStream(stream);
	
    StartOfScan* sos = new StartOfScan(3, encodedImageData, useSubSampling);
    sos->addToStream(stream);
    
    EndOfImage* eoi = new EndOfImage();
    eoi->addToStream(stream);
	
    stream.saveToFile(pathToFile);
	
    delete channelData;
    delete encodedImageData;
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
    channelData->unnormalize(255);
    
    Arai::transform(channelData->channel1->values, Y_width, Y_height);
    Arai::transform(channelData->channel2->values, CbCr_width, CbCr_height);
    Arai::transform(channelData->channel3->values, CbCr_width, CbCr_height);

    Quantization::run(channelData->channel1->values, Y_width, Y_height, luminanceQT);
    Quantization::run(channelData->channel2->values, CbCr_width, CbCr_height, chrominanceQT);
    Quantization::run(channelData->channel3->values, CbCr_width, CbCr_height, chrominanceQT);
    
    generateYDataAndHT();
    generateCbCrDataAndHT();
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
	DC_huffman.preventAllOnesPath();
    
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
	AC_huffman.preventAllOnesPath();
    CbCr_AC = AC_huffman.canonicalEncoding(16);
}
