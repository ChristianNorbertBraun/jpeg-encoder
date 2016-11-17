//
//  JPEGImage.cpp
//  jpegenc
//
//  Created by Christian Braun on 16/11/16.
//  Copyright © 2016 FHWS. All rights reserved.
//

#include "JPEGSegments.hpp"
#include <iostream>
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
	
	addChannel1ToStream(stream);
}

void StartOfFrame0::addChannel1ToStream(Bitstream &stream) {
	auto numberOfPixel = image->imageSize.pixelCount;
	for(int i = 0; i < numberOfPixel; ++i) {
		stream.add(1, 8);
		stream.add(0x22, 8);
		stream.add(image->channel1->getValue(i, image->imageSize), 8);
	}
}