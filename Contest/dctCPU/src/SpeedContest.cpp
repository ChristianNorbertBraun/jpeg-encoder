#include "SpeedContest.hpp"
#include <stdlib.h>
#include <thread>
#include <math.h>
#include <mutex>
#include "Arai.hpp"
#include "DCT.hpp"

#ifdef _WIN32
#include <Windows.h>
#endif

#include <functional>
#include <chrono>

#define PerformancePrintOperationsPerSecond(__desc, __time, __count) \
printf("<%s> took %lfsec with %lu iterations (%lfms per operation)\n", __desc, __time, (unsigned long)__count, (__time / __count) * 1000);


class Timer {
public:
	Timer() : beg_(clock_::now()) {}
	void reset() { beg_ = clock_::now(); }
	double elapsed() const {
		return std::chrono::duration_cast<second_>
		(clock_::now() - beg_).count(); }
	
private:
	typedef std::chrono::high_resolution_clock clock_;
	typedef std::chrono::duration<double, std::ratio<1> > second_;
	std::chrono::time_point<clock_> beg_;
};

// 1, 7, 3, 4, 5, 4, 3, 2
// one way transform gets:
// 10.253, 0.797218, -2.19761, -0.0377379, -1.76777, -2.75264, -2.53387, -1.13403

//  ---------------------------------------------------------------
// |
// |  Helper
// |
//  ---------------------------------------------------------------

float* createTestMatrix(size_t width, size_t height) {
	float *data = new float[width * height];
	for (size_t y = 0; y < height; y++) {
		for (size_t x = 0; x < width; x++) {
			data[y * width + x]= (float)((x+y*8) % 256);
		}
	}
	return data; // remember to delete[]
}

float* createOurTestMatrix(size_t w, size_t h) {
	float data[8] = {1, 7, 3, 4, 5, 4, 3, 2}; // generate our well known test matrix
	float *vls = new float[w * h];
	for (size_t y = 0; y < h; y++) {
		for (size_t x = 0; x < w; x++) {
			float &u = vls[y * w + x];
			if (y % 8 == 0) {
				u = data[x % 8];
			} else {
				if (x % 8 == 0) {
					u = data[y % 8];
				} else {
					u = 0;
				}
			}
		}
	}
	return vls;
}

void printFloatMatrix(float* &mat, size_t w, size_t h) {
	for (size_t i = 0; i < w * h; ++i) {
		if (i % (w*8) == 0) printf("\n");
		if (mat[i] < 0.0005F && mat[i] > -0.0005F) printf("%8d ", 0);
		else printf("%8.3f ", mat[i]);
		if (i % 8 == 7)   printf("   ");
		if (i % w == w-1) printf("\n");
	}
}

inline void copyArray(float* dst, float* src, size_t size) {
	memcpy(dst, src, size * sizeof(float));
//	while (size--) {
//		*(dst++) = *(src++);
//	}
}

void verify(const char* desc, float* originalMatrix, size_t width, size_t height, std::function<float*(float*)> func) {
	size_t size = width * height;
	float *vls = new float[size];
	copyArray(vls, originalMatrix, size);
	float *result = func(vls);
	
	printf("\n%s\n", desc);
	printFloatMatrix(result, width, height);
	printf("------------------------------------------------------------> ");
	
	// Check if result is correct
	bool inverseIsCorrect = true;
	float *out = new float[size];
	DCT::inverse(result, out, width, height);
	while (size--) {
		if (fabsf(originalMatrix[size] - out[size]) > 0.0005F )
			inverseIsCorrect = false;
	}
	delete [] vls;
	
	// Break on Error
	if (inverseIsCorrect) {
		printf("Inverse: CORRECT\n");
	} else {
		printf("Inverse: FAILED\nInverse:\n");
		printFloatMatrix(out, width, height);
		exit(EXIT_FAILURE);
	}
	delete [] out;
}

// ################################################################
// #
// #  CPU - Single Core
// #
// ################################################################

void testSingleCore(const double seconds, const char* description, std::function<void()> func)
{
	unsigned long numberOfIterations = 0;
	Timer t;
	while (t.elapsed() < seconds) {
		func();
		++numberOfIterations;
	}
	double time = t.elapsed();
	
	PerformancePrintOperationsPerSecond(description, time, numberOfIterations);
}

void runCPUSingleCore(float* matrix, size_t width, size_t height, double seconds) {
	size_t size = width * height;
	float* vls = new float[size];
	float* out = new float[size];
	
	copyArray(vls, matrix, size);
	testSingleCore(seconds, "Normal DCT", [&]{
		DCT::transform(vls, out, width, height);
	});
	delete [] out;
	
	copyArray(vls, matrix, size);
	testSingleCore(seconds, "Separated DCT", [&]{
		DCT::transform2(vls, width, height);
	});
	
	copyArray(vls, matrix, size);
	testSingleCore(seconds, "Arai inline transpose", [&]{
		Arai::transformInlineTranspose(vls, width, height);
	});
	
	copyArray(vls, matrix, size);
	testSingleCore(seconds, "Arai DCT", [&]{
		Arai::transform(vls, width, height);
	});
	
	delete [] vls;
}

// ################################################################
// #
// #  CPU - Multi Core
// #
// ################################################################

void testMultiThread(const char* desc, float* matrix, size_t width, size_t height, double seconds, std::function<void(float*,float*)> func) {
	static const unsigned int hwThreadCount = std::thread::hardware_concurrency();
	unsigned long* counters = new unsigned long[hwThreadCount];
	std::thread* threads = new std::thread[hwThreadCount];
	
	// Prepare an array of input images to have coalesced memory access
	// This wouldn't be necessary in a real life situation, since input
	// should be in general different images with different memory locations
	float** multiSrc = new float*[hwThreadCount];
	float** multiDst = new float*[hwThreadCount];
	unsigned int u = hwThreadCount;
	while (u--) {
		multiSrc[u] = new float[width * height];
		multiDst[u] = new float[width * height];
		copyArray(multiSrc[u], matrix, width * height);
	}
	
	// Start as many threads as supported
	unsigned int i = hwThreadCount;
	Timer t;
	while (i--) {
		threads[i] = std::thread([multiSrc, multiDst, seconds, i, func, &counters, width, height]{
			Timer inner;
			unsigned long innerCounter = 0;
			while (inner.elapsed() < seconds) {
				func(multiSrc[i], multiDst[i]);
				++innerCounter;
			}
			counters[i] = innerCounter;
		});
	}
	
	unsigned long totalIterations = 0;
	i = hwThreadCount;
	while (i--) {
		threads[i].join();
		totalIterations += counters[i];
	}
	
	double time = t.elapsed();
	PerformancePrintOperationsPerSecond(desc, time, totalIterations);
	
	
	// Cleanup
	u = hwThreadCount;
	while (u--) {
		delete [] multiSrc[u];
		delete [] multiDst[u];
	}
	delete [] counters;
	delete [] threads;
}

#ifdef _WIN32
CRITICAL_SECTION g_critSec;
#else
std::mutex lock;
#endif

unsigned long threadSyncedCounter = 0;
inline bool multiThreadShouldProcess(unsigned long val) {
#ifdef _WIN32
	EnterCriticalSection( &g_critSec );
#else
	lock.lock();
#endif
	if (val > threadSyncedCounter) {
#ifdef _WIN32
		LeaveCriticalSection( &g_critSec );
#else
		lock.unlock();
#endif
		return false;
	} else {
		threadSyncedCounter = val + 1;
#ifdef _WIN32
		LeaveCriticalSection( &g_critSec );
#else
		lock.unlock();
#endif
		return true;
	}
}

void testMultiThread2(const char* desc, float* matrix, size_t width, size_t height, double seconds, std::function<void(float*,float*,size_t)> func) {
	static const unsigned int hwThreadCount = std::thread::hardware_concurrency();
	unsigned long* counters = new unsigned long[hwThreadCount];
	std::thread* threads = new std::thread[hwThreadCount];
	
	float* src = new float[width * height];
	float* dst = new float[width * height];
	copyArray(src, matrix, width * height);
	
	threadSyncedCounter = 0;
	
	// calc size and offset
	size_t heightPerThread = height / hwThreadCount;
	if (heightPerThread & 0b111) {
		heightPerThread -= heightPerThread & 0b111; // floor to full 8-float block
	}
	size_t threadHeight = height - (heightPerThread * (hwThreadCount - 1)); // height of last thread
	
	// Start as many threads as supported
	unsigned int i = hwThreadCount;
	Timer t;
	while (i--) {
		threads[i] = std::thread([heightPerThread, threadHeight, src, dst, seconds, i, func, &counters, width]{
			float* ptrSrc = &src[ i * width * heightPerThread ];
			float* ptrDst = &dst[ i * width * heightPerThread ];
			unsigned long innerCounter = 0;
			Timer inner;
			while (inner.elapsed() < seconds) {
				if (multiThreadShouldProcess(innerCounter)) {
					func(ptrSrc, ptrDst, threadHeight);
					++innerCounter;
				}
			}
			counters[i] = innerCounter;
		});
		threadHeight = heightPerThread;
	}
	
	i = hwThreadCount;
	while (i--) {
		threads[i].join();
	}
	
	double time = t.elapsed();
	PerformancePrintOperationsPerSecond(desc, time, threadSyncedCounter);
	
	// Cleanup
	delete [] counters;
	delete [] threads;
}

void runCPUMultiCore(float* matrix, size_t width, size_t height, double seconds) {
#ifdef _WIN32
	InitializeCriticalSection( &g_critSec );
#endif
	printf("Threading on single image:\n");
	testMultiThread2("Normal DCT", matrix, width, height, seconds, [width, height](float *in, float* out, size_t h){
		DCT::transform(in, out, width, h);
	});
	testMultiThread2("Separated DCT", matrix, width, height, seconds, [width, height](float *in, float*, size_t h){
		DCT::transform2(in, width, h);
	});
	testMultiThread2("Arai inline transpose", matrix, width, height, seconds, [width, height](float *in, float*, size_t h){
		Arai::transformInlineTranspose(in, width, h);
	});
	testMultiThread2("Arai DCT", matrix, width, height, seconds, [width, height](float *in, float*, size_t h){
		Arai::transform(in, width, h);
	});
	
	printf("\nThreading on image queue:\n");
	testMultiThread("Normal DCT", matrix, width, height, seconds, [width, height](float *in, float* out){
		DCT::transform(in, out, width, height);
	});
	testMultiThread("Separated DCT", matrix, width, height, seconds, [width, height](float *in, float*){
		DCT::transform2(in, width, height);
	});
	testMultiThread("Arai inline transpose", matrix, width, height, seconds, [width, height](float *in, float*){
		Arai::transformInlineTranspose(in, width, height);
	});
	testMultiThread("Arai DCT", matrix, width, height, seconds, [width, height](float *in, float*){
		Arai::transform(in, width, height);
	});
	
#ifdef _WIN32
	DeleteCriticalSection( &g_critSec );
#endif
}

// ################################################################
// #
// #  Main
// #
// ################################################################

void SpeedContest::run(double seconds, bool skipSingeCore) {
	size_t width = 256, height = 256;
	float *matrix = createTestMatrix(width, height); // (x+y*8) % 256;
	printf("Size: %lux%lu\n", width, height);
	if (skipSingeCore == false) {
		printf("\n== Single-Threaded ==\n");
		runCPUSingleCore(matrix, width, height, seconds);
	}
	
	printf("\n== Multi-Threading ==\n");
	printf("Threads: %d\n", std::thread::hardware_concurrency());
	runCPUMultiCore(matrix, width, height, seconds);
	
	printf("\n");
	
	delete [] matrix;
}

// ################################################################
// #
// #  Correctness Test
// #
// ################################################################

void SpeedContest::testForCorrectness(bool ourTestMatrix, bool use16x16, bool modifyData) {
	size_t width = 8, height = 8;
	if (use16x16) {
		width = 16; height = 16;
	}
	
	float* matrix = ( ourTestMatrix
					 ? createOurTestMatrix(width, height) // 1, 7, 3, 4, 5, 4, 3, 2
					 : createTestMatrix(width, height) );
	
	if (use16x16 && modifyData) {
		// modify some values to get different results, especially for 16 x 16
		matrix[8] = 0;
		matrix[8+1] = 4;
		matrix[8+width] = 4;
		matrix[2*width + 1] = 8;
		matrix[12*width + 1] = 1;
	}
	
	printf("\nInput:\n");
	printFloatMatrix(matrix, width, height);
	printf("------------------------------------------------------------------------\n");
	
	float *out = new float[width * height];
	verify("DCT Normal", matrix, width, height, [&](float* mat){
		DCT::transform(mat, out, width, height); return out;
	});
	verify("DCT Separated", matrix, width, height, [&](float* mat){
		DCT::transform2(mat, width, height); return mat;
	});
	verify("Arai", matrix, width, height, [&](float* mat){
		Arai::transform(mat, width, height); return mat;
	});
	verify("Arai Inline Transpose", matrix, width, height, [&](float* mat){
		Arai::transformInlineTranspose(mat, width, height); return mat;
	});
	
	printf("\nInverse:\n");
	DCT::transform2(matrix, width, height);
	DCT::inverse(matrix, out, width, height);
	printFloatMatrix(out, width, height);
	printf("-----------------------------------------------------------------> ALL CORRECT\n");
	delete [] out;
	delete [] matrix;
}
