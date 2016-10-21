#ifndef Image_hpp
#define Image_hpp

#include <stdlib.h>
#include <stdio.h>

struct Channel {
	size_t channelSize;
	size_t *values;

	Channel(size_t size = 0) : channelSize(size) {
		values = new size_t[size];
	}

	~Channel() {
		delete[] values;
	}
};

struct Image {
	size_t width, height;
	size_t numberOfPixels;
	std::string colorSpace;
	Channel *channel1;
	Channel *channel2;
	Channel *channel3;

	Image(size_t width, size_t height) : width(width), height(height), numberOfPixels(width * height) {
		channel1 = new Channel(numberOfPixels);
		channel2 = new Channel(numberOfPixels);
		channel3 = new Channel(numberOfPixels);
	}

	~Image() {
		delete channel1;
		delete channel2;
		delete channel3;
	}

	size_t getIndex(size_t x, size_t y) const;

	size_t getValueFromChannel1(size_t x, size_t y);

	size_t getValueFromChannel2(size_t x, size_t y);

	size_t getValueFromChannel3(size_t x, size_t y);

	size_t getValueFromChannel1(size_t index);

	size_t getValueFromChannel2(size_t index);

	size_t getValueFromChannel3(size_t index);

	void setValueOnChannel1(size_t x, size_t y, size_t value);

	void setValueOnChannel2(size_t x, size_t y, size_t value);

	void setValueOnChannel3(size_t x, size_t y, size_t value);

	void setValueOnChannel1(size_t index, size_t value);

	void setValueOnChannel2(size_t index, size_t value);

	void setValueOnChannel3(size_t index, size_t value);

	void print();

	void setColorSpace(std::string colorSpace) {
		this->colorSpace = colorSpace;
	}

	void reduceBySubSamplingChannel1(size_t stepWidth);
	void reduceBySubSamplingChannel2(size_t stepWidth);
	void reduceBySubSamplingChannel3(size_t stepWidth);
};

#endif /* Image_hpp */