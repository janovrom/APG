#pragma once

#include <cfloat>
#include <queue>
#include <stack>
#include <cstdlib>
#include <cmath>
#define MIN_CONTEXTS 32
/*
#define _CRTDBG_MAP_ALLOC
#define _CRTDBG_MAPALLOC
#include <stdlib.h>
#include <crtdbg.h>
#ifdef _DEBUG
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW
#endif
*/
std::stack<float *> modelViewStack;
std::stack<float *> projectionStack;

// identity matrix for easy copying
float identityMatrix[16] ={	1.0f, 0.0f, 0.0f, 0.0f,
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
	tmp += left[i] * right[j * 4];
	tmp += left[i + 4] * right[j * 4 + 1];
	tmp += left[i + 8] * right[j * 4 + 2];
	tmp += left[i + 12] * right[j * 4 + 3];

	return tmp;
}
/**
Compute dot product of line of left matrix and vector.
*/
float dotVectors(const float* left, inputPoint4f* vector, int i) {
	float tmp = 0;
	tmp += left[i] * vector->x;
	tmp += left[i + 4] * vector->y;
	tmp += left[i + 8] * vector->z;
	tmp += left[i + 12] * vector->w;

	return tmp;
}

void multiplyMatrixVector(const float* matrix, inputPoint4f* vector, inputPoint4f& output) {
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

// HERE UNUSED CODE STARTS
bool invertMatrix(const float m[16], float invOut[16])
{
	float inv[16], det;
	int i;

	inv[0] = m[5] * m[10] * m[15] -
		m[5] * m[11] * m[14] -
		m[9] * m[6] * m[15] +
		m[9] * m[7] * m[14] +
		m[13] * m[6] * m[11] -
		m[13] * m[7] * m[10];

	inv[4] = -m[4] * m[10] * m[15] +
		m[4] * m[11] * m[14] +
		m[8] * m[6] * m[15] -
		m[8] * m[7] * m[14] -
		m[12] * m[6] * m[11] +
		m[12] * m[7] * m[10];

	inv[8] = m[4] * m[9] * m[15] -
		m[4] * m[11] * m[13] -
		m[8] * m[5] * m[15] +
		m[8] * m[7] * m[13] +
		m[12] * m[5] * m[11] -
		m[12] * m[7] * m[9];

	inv[12] = -m[4] * m[9] * m[14] +
		m[4] * m[10] * m[13] +
		m[8] * m[5] * m[14] -
		m[8] * m[6] * m[13] -
		m[12] * m[5] * m[10] +
		m[12] * m[6] * m[9];

	inv[1] = -m[1] * m[10] * m[15] +
		m[1] * m[11] * m[14] +
		m[9] * m[2] * m[15] -
		m[9] * m[3] * m[14] -
		m[13] * m[2] * m[11] +
		m[13] * m[3] * m[10];

	inv[5] = m[0] * m[10] * m[15] -
		m[0] * m[11] * m[14] -
		m[8] * m[2] * m[15] +
		m[8] * m[3] * m[14] +
		m[12] * m[2] * m[11] -
		m[12] * m[3] * m[10];

	inv[9] = -m[0] * m[9] * m[15] +
		m[0] * m[11] * m[13] +
		m[8] * m[1] * m[15] -
		m[8] * m[3] * m[13] -
		m[12] * m[1] * m[11] +
		m[12] * m[3] * m[9];

	inv[13] = m[0] * m[9] * m[14] -
		m[0] * m[10] * m[13] -
		m[8] * m[1] * m[14] +
		m[8] * m[2] * m[13] +
		m[12] * m[1] * m[10] -
		m[12] * m[2] * m[9];

	inv[2] = m[1] * m[6] * m[15] -
		m[1] * m[7] * m[14] -
		m[5] * m[2] * m[15] +
		m[5] * m[3] * m[14] +
		m[13] * m[2] * m[7] -
		m[13] * m[3] * m[6];

	inv[6] = -m[0] * m[6] * m[15] +
		m[0] * m[7] * m[14] +
		m[4] * m[2] * m[15] -
		m[4] * m[3] * m[14] -
		m[12] * m[2] * m[7] +
		m[12] * m[3] * m[6];

	inv[10] = m[0] * m[5] * m[15] -
		m[0] * m[7] * m[13] -
		m[4] * m[1] * m[15] +
		m[4] * m[3] * m[13] +
		m[12] * m[1] * m[7] -
		m[12] * m[3] * m[5];

	inv[14] = -m[0] * m[5] * m[14] +
		m[0] * m[6] * m[13] +
		m[4] * m[1] * m[14] -
		m[4] * m[2] * m[13] -
		m[12] * m[1] * m[6] +
		m[12] * m[2] * m[5];

	inv[3] = -m[1] * m[6] * m[11] +
		m[1] * m[7] * m[10] +
		m[5] * m[2] * m[11] -
		m[5] * m[3] * m[10] -
		m[9] * m[2] * m[7] +
		m[9] * m[3] * m[6];

	inv[7] = m[0] * m[6] * m[11] -
		m[0] * m[7] * m[10] -
		m[4] * m[2] * m[11] +
		m[4] * m[3] * m[10] +
		m[8] * m[2] * m[7] -
		m[8] * m[3] * m[6];

	inv[11] = -m[0] * m[5] * m[11] +
		m[0] * m[7] * m[9] +
		m[4] * m[1] * m[11] -
		m[4] * m[3] * m[9] -
		m[8] * m[1] * m[7] +
		m[8] * m[3] * m[5];

	inv[15] = m[0] * m[5] * m[10] -
		m[0] * m[6] * m[9] -
		m[4] * m[1] * m[10] +
		m[4] * m[2] * m[9] +
		m[8] * m[1] * m[6] -
		m[8] * m[2] * m[5];

	det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

	if (det == 0)
		return false;

	det = 1.0 / det;

	for (i = 0; i < 16; i++)
		invOut[i] = inv[i] * det;

	return true;
}
// HERE UNUSED CODE ENDS

void multiplyMatrix(float* left, const float* right) {
	float output[16];
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			output[i*4 + j] = dotVectors(left, right, j, i);
		}
	}
	copyMatrix(left, output);
	/*for (int i = 0; i < 4; i++) {
		printf("\n");
		for (int j = i; j < 16; j += 4) {
			printf("%f ", left[j]);
		}
	}
	*/
}
// HERE UNUSED CODE STARTS
bool invertedForObject = false;
float inversedProjectionMatrix[16];
float inversedViewportMatrix[16];
// HERE UNUSED CODE ENDS

float viewportMatrix[16];
float multipliedMatrix[16];

/**
It provides information if sglBegin was called. It is set to false
and after sglBegin call changed to true (likewise after sglEnd
it is set to false).
*/
bool hasBegun;
int viewportOffsetX, viewportOffsetY, viewportWidth, viewportHeight;
/**
Do not do depth test as default.
*/
bool testDepth = false;
sglEMatrixMode matrixMode = SGL_MODELVIEW;

short pointSize = 0;
float colorVertexR = 0, colorVertexG = 0, colorVertexB = 0;
float colorClearR = 0, colorClearG = 0, colorClearB = 0;

bool depthEnabled = false;
sglEElementType drawingMethod = SGL_POINTS;

std::queue<inputPoint4f*> queue4f;

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
	int activeContext;
	SglContext* contexts[MIN_CONTEXTS];
	int length;
	int count;

	ContextWrapper() {
		activeContext = -1;
		length = MIN_CONTEXTS;
		count = 0;
	}

	~ContextWrapper() {
		clear();
	}

	SglContext* operator[] (int id) {
		if (id < length)
			return contexts[id];
		else
			return NULL;
	}

	void clear() {
		for (int i = 0; i < length; ++i) {
			delete contexts[i];
			contexts[i] = NULL;
		}
		count = 0;
	}

	bool clear(int id) {
		if (id < length && id >= 0) {
			delete contexts[id];
			--count;
			return true;
		}else {
			return false;
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
			if (contexts[i] == NULL)
				return i;
		}
		return -1;
	}
} contextWrapper;