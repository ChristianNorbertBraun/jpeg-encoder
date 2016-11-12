#include <iostream>
#include <stdlib.h>
#include "io/PPMLoader.hpp"
#include "converter/RGBToYCbCrConverter.hpp"
#include "converter/YCbCrToRGBConverter.hpp"
#include "helper/Test.hpp"

#include "bitstream/BitstreamChris.hpp"
#include "bitstream/BitstreamMarcel.hpp"
#include "bitstream/BitstreamMarv.hpp"
#include "bitstream/BitstreamOleg.hpp"


//  ---------------------------------------------------------------
// |
// |  PPM Image Processing
// |
//  ---------------------------------------------------------------

void testImage() {
	std::cout << "Loading image ..." << std::endl;
	Test::performance([]{
		PPMLoader loader;
		auto image = loader.load("data/gigantic.test.ppm");
		
//		RGBToYCbCrConverter converter1;
//		converter1.convert(image);
//
//		image->print();
//
//		image->channel2->reduceBySubSampling( image->imageSize.width, image->imageSize.height );
//		image->channel3->reduceBySubSampling( image->imageSize.width, image->imageSize.height );
//
//		YCbCrToRGBConverter converter2;
//		converter2.convert(image);
//
//		image->reduceBySubSample(2, 2);
//		image->reduceByAverage(2, 2);
//		image->print();

//		Test::performance([&loader, &image]{
//			loader.write("data/output.test.ppm", image);
//		});
	});
}

//  ---------------------------------------------------------------
// |
// |  Bitstream
// |
//  ---------------------------------------------------------------

void testChris() {
	std::cout << "Testing, Chris" << std::endl;
	std::cout << "Write single bit: ";
	Test::performance(10000000, 10, [](size_t numberOfElements){
		BitStream bitstream;
		while (numberOfElements--) {
			bitstream.add(true);
		}
	});
}

void testMarcel() {
	BitstreamMarcel bitstream;
	
	for (int i = 0; i < 16; ++i) {
		bitstream.add(true);
	}
	
	bitstream.print();
	
	cout << bitstream.read(30) << endl;
}

void testMarv() {
	
	BitStreamMarv bitStreamMarv;
	
	// Test adding bits
	bitStreamMarv.add(true);
	bitStreamMarv.add(true);
	bitStreamMarv.add(false);
	bitStreamMarv.add(true);
	bitStreamMarv.add(false);
	bitStreamMarv.add(false);
	bitStreamMarv.add(true);
	bitStreamMarv.add(false);
	
	// Test adding bytes
	char byte = (char) 0xd2;   // 11010010
	bitStreamMarv.add( byte );
	
	// Test printing
	bitStreamMarv.print();
	
	// Test reading bits
	for (size_t i = 0; i < bitStreamMarv.size(); ++i) {
		std::cout << bitStreamMarv.read(i);
	}
	std::cout << std::endl;
	
	// Test reading multiple bits
	size_t firstIndex = 0;
	size_t lastIndex = bitStreamMarv.size() - 1;
	bool *bits = bitStreamMarv.read(firstIndex, lastIndex);
	for (size_t i = 0; i < bitStreamMarv.size(); ++i) {
	std:cout << bits[i];
	}
	std::cout << std::endl;
	delete[] bits;
	
	// Test saving
	bitStreamMarv.saveToFile("/home/marv/Projects/jpeg-encoder/bitstream.txt");
	
	
	std::cout << "Testing, Marv" << std::endl;
	// Test performance adding bits (0.022s)
	std::cout << "Write single bit: "; // (0.074520 seconds)
	Test::performance(10000000, 10, [](size_t numberOfElements){
		BitStreamMarv bitStream(numberOfElements);
		while (numberOfElements--) {
			bitStream.add(true);
		}
	});
	
	// Test performance adding bytes (0.154s)
	std::cout << "Write byte bits: "; // (0.126076 seconds)
	Test::performance(10000000, 10, [](size_t numberOfElements){
		BitStreamMarv bitStream(numberOfElements);
		while (numberOfElements--) {
			bitStream.add((char) 0xd2); // 11010010
		}
	});
}

void testOleg() {
	BitstreamOleg bitstream;
//	bitstream.add(true);
//	bitstream.add(false);
//	bitstream.add(false);
//	bitstream.add(true);
	bitstream.add('A', 8); // 0100 0001
	bitstream.add('x', 8); // 0111 1000
	bitstream.add('l', 8); // 0110 1100
	bitstream.add('D', 8); // 0100 0100
	bitstream.add(' ', 8); // 0010 0000
//	for (size_t i = 0; i < 524288*8-2; ++i) {
//		bitstream.add(1);
//	}
	bitstream.add(255, 6);
//	bitstream.add(0);
//	bitstream.saveToFile("data/out.txt");
	bitstream.add(0);
	bitstream.add(1);
	bitstream.add(0);
//	bitstream.saveToFile("data/out2.txt");
//	bitstream.print();
	
	
	std::cout << "Testing, Oleg" << std::endl;
	
	std::cout << "Write single bit: ";
	Test::performance(10000000, 10, [](size_t numberOfElements){
		BitstreamOleg bitstream;
		while (numberOfElements--) {
			bitstream.add(1);
		}
	});
	
	std::cout << "Write byte bits: ";
	Test::performance(10000000, 10, [](size_t numberOfElements){
		BitstreamOleg bitstream;
		//bitstream.add(1);
		while (numberOfElements--) {
			bitstream.add(0xd2, 8);
		}
	});
	
	// create random bitstream for reading
	BitstreamOleg testStream;
	size_t fillRandom = 100000;
	while (fillRandom--)
		bitstream.add( arc4random() % 2 );
	
	std::cout << "Read single bit: ";
	Test::performance(10000000, 10, [&testStream](size_t numberOfElements){
		size_t maxRead = testStream.numberOfBits() - 2;
		size_t idx = 0;
		while (numberOfElements--) {
			testStream.read(idx++);
			if (idx > maxRead)
				idx = 0;
		}
	});
}

// ################################################################
// #
// #  Main
// #
// ################################################################

int main(int argc, const char *argv[]) {
//	testImage();
//	testChris();
//	testMarcel();
//	testMarv();
//	testOleg();
	
	return 0;
}
