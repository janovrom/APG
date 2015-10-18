#pragma once

#include <cfloat>

#define MIN_CONTEXTS 32

/**
It provides information if sglBegin was called. It is set to false
and after sglBegin call changed to true (likewise after sglEnd
it is set to false).
*/
bool hasBegun;
int offsetX, offsetY, windowWidth, windowHeight;

class SglContext {
private:
	int width;
	int height;
	/**
	Color buffer saves integers from 0 to 1.
	*/
	float *colorBuffer;
	float *depthBuffer;
	/**
	Color to which is color buffer saved.
	*/
	float r, g, b, alpha;
public:
	SglContext(int width, int height) : width{ width }, height{ height } {
		colorBuffer = new float[width*height * 3];
		depthBuffer = new float[width*height];
	};
	~SglContext() {
		delete[] colorBuffer;
		delete[] depthBuffer;
	}

	float* getColorBuffer() {
		return colorBuffer;
	}

	float* getDepthBuffer() {
		return depthBuffer;
	}

	void setClearColor(float r, float g, float b, float alpha) {
		this->alpha = alpha;
		this->r = r;
		this->g = g;
		this->b = b;
	}

	void clearColor() {
		for (int i = 0; i < width; i += 3) {
			for (int j = 0; j < height; j += 3) {
				*(colorBuffer + i*j + j) = r;
				*(colorBuffer + i*j + j + 1) = g;
				*(colorBuffer + i*j + j + 2) = b;
			}
		}
	}

	void clearDepth() {
		for (int i = 0; i < width; ++i) {
			for (int j = 0; j < height; ++j) {
				*(depthBuffer + i*j + j) = FLT_MAX;
			}
		}
	}
};

struct ContextWrapper {
	int activeContext = -1;
	SglContext* contexts[MIN_CONTEXTS];
	int size = MIN_CONTEXTS;
	int count = 0;

	SglContext* operator[] (int id) {
		if (id < size)
			return contexts[id];
		else
			return nullptr;
	}

	void clear() {
		for (int i = 0; i < size; ++i) {
			delete contexts[i];
			contexts[i] = nullptr;
		}
		count = 0;
	}

	void clear(int id) {
		if (id < size) {
			delete contexts[id];
			--count;
		}
	}

	bool empty() {
		return count == 0;
	}

	int size() {
		return size;
	}

	void add(SglContext* c);
private:
	ContextWrapper();
	int findFirstEmpty() {
		for (int i = 0; i < size; ++i) {
			if (contexts[i] == nullptr)
				return i;
		}
		return -1;
	}
} contexts;