#include <iostream>
#include "io/PPMLoader.hpp"
#include "converter/RGBToYCbCrConverter.hpp"
#include "converter/YCbCrToRGBConverter.hpp"
#include <stdlib.h>
#include "bitstream/BitstreamMarcel.hpp"
#include "bitstream/BitstreamMarv.hpp"
#include "bitstream/BitstreamOleg.hpp"

void testImage() {
	PPMLoader loader;
	clock_t t;
	t = clock();

	auto image = loader.load("data/gigantic.test.ppm");

	RGBToYCbCrConverter converter1;
	image = converter1.convert(image);

	image->channel2->reduceBySubSampling( image->imageSize.width, image->imageSize.height );
	image->channel3->reduceBySubSampling( image->imageSize.width, image->imageSize.height );

	YCbCrToRGBConverter converter2;
	image = converter2.convert(image);

	loader.write("data/output.test.ppm", image);

	t = clock() - t;
	printf ("It took me %lu clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}

void testChris() {
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

	// Test performance adding bits (0.022s)
	size_t numberOfRounds = 100;
	size_t numberOfElements = 10000000;
	clock_t timeStamp;
	clock_t totalTime = 0;

	for (size_t i = 0; i < numberOfRounds; ++i)
	{
		timeStamp = clock();

		BitStreamMarv bitStream(numberOfElements);
		for (size_t k = 0; k < numberOfElements; ++k)
		{
			bitStream.add(true);
		}
		timeStamp = clock() - timeStamp;
		totalTime = totalTime + timeStamp;
	}
	clock_t averageTime = totalTime / numberOfRounds;

	printf("Adding %lu single bits took %lu clicks (%f seconds) on average (%lu times).\n",numberOfElements, averageTime,((float)averageTime)/CLOCKS_PER_SEC, numberOfRounds);

	// Test performance adding bytes (0.154s)
	numberOfRounds = 100;
	numberOfElements = 10000000;
	totalTime = 0;

	for (size_t i = 0; i < numberOfRounds; ++i)
	{
		timeStamp = clock();

		BitStreamMarv bitStream(numberOfElements);
		for (size_t k = 0; k < numberOfElements; ++k)
		{
			bitStream.add((char) 0xd2); // 11010010
		}
		timeStamp = clock() - timeStamp;
		totalTime = totalTime + timeStamp;
	}
	averageTime = totalTime / numberOfRounds;

	printf("Adding %lu bytes took %lu clicks (%f seconds) on average (%lu times).\n",numberOfElements, averageTime,((float)averageTime)/CLOCKS_PER_SEC, numberOfRounds);
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
//	for (size_t i = 0; i < 131072*8-2; ++i) {
//		bitstream.add(1);
//	}
	bitstream.add(255, 6);
//	bitstream.add(0);
	bitstream.saveToFile("data/out.txt");
	bitstream.add(0);
	bitstream.add(1);
	bitstream.add(0);
	bitstream.saveToFile("data/out2.txt");
//	bitstream.print();
	
	
	// Test performance adding bits
	size_t numberOfRounds = 100;
	size_t numberOfElements = 10000000;
	clock_t timeStamp;
	clock_t totalTime = 0;
	
	for (size_t i = 0; i < numberOfRounds; ++i)
	{
		timeStamp = clock();
		
		BitstreamOleg bitStream;
		for (size_t k = 0; k < numberOfElements; ++k)
		{
			bitStream.add(true);
		}
		timeStamp = clock() - timeStamp;
		totalTime = totalTime + timeStamp;
	}
	clock_t averageTime = totalTime / numberOfRounds;
	
	printf("Adding %lu single bits took %lu clicks (%f seconds) on average (%lu times).\n",numberOfElements, averageTime,((float)averageTime)/CLOCKS_PER_SEC, numberOfRounds);
	
	// Test performance adding bytes
	numberOfRounds = 100;
	numberOfElements = 10000000;
	totalTime = 0;
	
	for (size_t i = 0; i < numberOfRounds; ++i)
	{
		timeStamp = clock();
		
		BitstreamOleg bitStream;
//		bitStream.add(1);
		for (size_t k = 0; k < numberOfElements; ++k)
		{
			bitStream.add((char) 0xd2, 8); // 11010010
		}
		timeStamp = clock() - timeStamp;
		totalTime = totalTime + timeStamp;
	}
	averageTime = totalTime / numberOfRounds;
	
	printf("Adding %lu bytes took %lu clicks (%f seconds) on average (%lu times).\n",numberOfElements, averageTime,((float)averageTime)/CLOCKS_PER_SEC, numberOfRounds);
}

int main(int argc, const char *argv[]) {
//	testImage();
//	testChris();
//	testMarcel();
//	testMarv();
//	testOleg();

	return 0;
}
