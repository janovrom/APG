#pragma once

#include <cfloat>
#include <queue>
#include <stack>
#define MIN_CONTEXTS 32

#define _CRTDBG_MAP_ALLOC
#define _CRTDBG_MAPALLOC
#include <stdlib.h>
#include <crtdbg.h>
#ifdef _DEBUG
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW
#endif

std::stack<float *> modelViewStack;
std::stack<float *> projectionStack;

float identity[16] ={	1.0f, 0.0f, 0.0f, 0.0f,
						0.0f, 1.0f, 0.0f, 0.0f, 
						0.0f, 0.0f, 1.0f, 0.0f, 
						0.0f, 0.0f, 0.0f, 1.0f};

struct inputPoint4f {
	float x, y, z, w;
	float r, g, b, a;
};

/**
Compute dot product of line and column of left and right matrix.
*/
float dotVectors(const float* left, const float* right, int i, int j) {
	float tmp = 0;
	tmp += left[i] * right[j];
	tmp += left[i + 4] * right[j + 1];
	tmp += left[i + 8] * right[j + 2];
	tmp += left[i + 12] * right[j + 3];

	return tmp;
}
/**
Compute dot product of line of left matrix and vector.
*/
float dotVectors(const float* left, inputPoint4f& vector, int i) {
	float tmp = 0;
	tmp += left[i] * vector.x;
	tmp += left[i + 4] * vector.y;
	tmp += left[i + 8] * vector.z;
	tmp += left[i + 12] * vector.w;

	return tmp;
}
void multiplyMatrixVector(const float* matrix, inputPoint4f& vector, inputPoint4f& output) {
	output.x = dotVectors(matrix, vector, 0);
	output.y = dotVectors(matrix, vector, 1);
	output.z = dotVectors(matrix, vector, 2);
	output.w = dotVectors(matrix, vector, 3);
}

void copyMatrix(float* output, const float* input) {
	for (int i = 0; i < 16; ++i) {
		output[i] = input[i];
	}
}

void multiplyMatrix(float* left, const float* right) {
	float output[16];
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			output[i*j + j] = dotVectors(left, right, i, j);
		}
	}
	copyMatrix(left, output);
}

/**
It provides information if sglBegin was called. It is set to false
and after sglBegin call changed to true (likewise after sglEnd
it is set to false).
*/
bool hasBegun;
int offsetX, offsetY, windowWidth, windowHeight;
/**
Do not do depth test as default.
*/
bool testDepth = false;
sglEMatrixMode matrixMode = SGL_MODELVIEW;

short pointSize = 0;
float colorVertexR = 0, colorVertexG = 0, colorVertexB = 0;
float colorClearR = 0, colorClearG = 0, colorClearB = 0;

bool depthEnabled = false;
sglEElementType drawingMethod = sglEElementType::SGL_POINTS;

std::queue<inputPoint4f> queue4f;

class SglContext {
private:
	int width;
	int height;
	/**
	Color buffer saves integers from 0 to 1.
	*/
	float *colorBuffer;
	float *depthBuffer;
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

	int getWidth() {
		return width;
	}

	int getHeight() {
		return height;
	}

	void clearColor(float r, float g, float b) {
		int w = width * 3;
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < w; j += 3) {
				*(colorBuffer + i*w + j) = r;
				*(colorBuffer + i*w + j + 1) = g;
				*(colorBuffer + i*w + j + 2) = b;
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
public:
	int activeContext = -1;
	SglContext* contexts[MIN_CONTEXTS];
	int length = MIN_CONTEXTS;
	int count = 0;

	~ContextWrapper() {
		clear();
	}

	SglContext* operator[] (int id) {
		if (id < length)
			return contexts[id];
		else
			return nullptr;
	}

	void clear() {
		for (int i = 0; i < length; ++i) {
			delete contexts[i];
			contexts[i] = nullptr;
		}
		count = 0;
	}

	void clear(int id) {
		if (id < length) {
			delete contexts[id];
			--count;
		}
	}

	bool empty() {
		return count == 0;
	}

	int size() {
		return length;
	}

	int add(SglContext* c);

	int findFirstEmpty() {
		for (int i = 0; i < length; ++i) {
			if (contexts[i] == nullptr)
				return i;
		}
		return -1;
	}
} contextWrapper;