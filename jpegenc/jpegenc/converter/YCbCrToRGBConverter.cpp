#include "YCbCrToRGBConverter.hpp"
#include <math.h>

std::shared_ptr<Image> YCbCrToRGBConverter::convert(std::shared_ptr<Image> originalImage) {
	
	Dimension size = originalImage->imageSize;
	auto convertedImage = std::make_shared<Image>(size);
	
	size_t index = size.pixelCount;
	while (index) {
		--index;
		color y = (color) originalImage->channel1->getValue(index, size);
		color cb = (color) originalImage->channel2->getValue(index, size);
		color cr = (color) originalImage->channel3->getValue(index, size);
		color r = normalize(round(y + 1.4021 * (cr - 128)));
		color g = normalize(round(y - 0.3441 * (cb - 128) - 0.7142 * (cr - 128)));
		color b = normalize(round(y + 1.772 * (cb - 128)));
		convertedImage->channel1->setValue(index, r);
		convertedImage->channel2->setValue(index, g);
		convertedImage->channel3->setValue(index, b);
	}
	convertedImage->colorSpace = ColorSpaceRGB;
	return convertedImage;
}

color YCbCrToRGBConverter::normalize(int value) {
	if (value < 0)
		value = 0;
	if (value > 255) // TODO: only valid for 8bit
		value = 255;
	return (color) value;
}
