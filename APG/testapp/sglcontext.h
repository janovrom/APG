#pragma once

#include <cfloat>
#include <queue>
#include <deque>
#include <stack>
#include <cstdlib>
#include <cmath>
#define MIN_CONTEXTS 32


std::stack<float *> modelViewStack;
std::stack<float *> projectionStack;

// identity matrix for easy copying
float identityMatrix[16] ={	1.0f, 0.0f, 0.0f, 0.0f,
						0.0f, 1.0f, 0.0f, 0.0f, 
						0.0f, 0.0f, 1.0f, 0.0f, 
						0.0f, 0.0f, 0.0f, 1.0f};

struct inputPoint4f
{
	float x, y, z, w;
	float r, g, b, a;
};

struct Texture
{
private:
	float *texels;

	float WrapAndResizeValues(float *u, float *v)
	{
		*u *= width - 1;
		*v *= height - 1;
		switch (wrap)
		{
			// sgl clamp mode is default, so if someone wants to add wrap, add the case above
			case SGL_CLAMP:
			default:
				*u = *u > width - 1 ? width - *u : *u;
				*v = *v > height - 1 ? height - *v : *v;
				*u = *u < 0 ? width + *u - 1 : *u;
				*v = *v < 0 ? height + *v - 1 : *v;
				break;
		}
	}
public:
	enum WrapMode
	{
		SGL_CLAMP
	};

	enum Filtering
	{
		SGL_LINEAR,
		SGL_NEAREST
	};

	Filtering filtering;
	WrapMode wrap;
	int width, height;

	Texture(Filtering filtering, int width, int height, float *texels, WrapMode wrap)
	{
		this->filtering = filtering;
		this->height = height;
		this->width = width;
		this->texels = texels;
		this->wrap = wrap;
	}

	~Texture()
	{
		delete[] texels;
	}

	/**
	* Returns color on given position. Row and Col are normalized
	* in range [0,1].
	*/
	void Get(float u, float v, float *r, float *g, float *b)
	{
		WrapAndResizeValues(&u, &v);
		switch (filtering)
		{
			case SGL_LINEAR:
			{
				u = u - 0.5f;
				v = v - 0.5f;
				int x = (int)floor(u);
				int y = (int)floor(v);
				float u_ratio = u - x;
				float v_ratio = v - y;
				float u_opposite = 1 - u_ratio;
				float v_opposite = 1 - v_ratio;
				*r = (texels[x *width * 3 + y * 3] * u_opposite + texels[(x + 1) *width * 3 + (y)* 3] * u_ratio) * v_opposite +
					(texels[x *width * 3 + (y + 1) * 3] * u_opposite + texels[(x + 1) *width * 3 + (y + 1) *width * 3] * u_ratio) * v_ratio;
				*g = (texels[x *width * 3 + y * 3 + 1] * u_opposite + texels[(x + 1) *width * 3 + (y)* 3 + 1] * u_ratio) * v_opposite +
					(texels[x *width * 3 + (y + 1) * 3 + 1] * u_opposite + texels[(x + 1) *width * 3 + (y + 1) *width * 3 + 1] * u_ratio) * v_ratio;
				*b = (texels[x *width * 3 + y * 3 + 2] * u_opposite + texels[(x + 1) *width * 3 + (y)* 3 + 2] * u_ratio) * v_opposite +
					(texels[x *width * 3 + (y + 1) * 3 + 2] * u_opposite + texels[(x + 1) *width * 3 + (y + 1) *width * 3 + 2] * u_ratio) * v_ratio;
			} break;
			case SGL_NEAREST:
			default:
			{
				int thisRow = (int)roundf(v)  * width * 3;
				int thisCol = (int)roundf(u);
				*r = texels[thisRow + thisCol];
				*g = texels[thisRow + thisCol + 1];
				*b = texels[thisRow + thisCol + 2];
			} break;
		}
	}
};


enum MaterialType
{
	EMISSIVE,
	PHONG
};

struct Material
{
MaterialType type;
protected:


	Material()
	{
	}
};

struct PhongMaterial : public Material
{
	float r, g, b;
	float kd, ks, shine;
	float transmitance;
	float refractIndex;

	PhongMaterial()
	{
		type = MaterialType::PHONG;
	}
};

struct EmissiveMaterial : public Material
{
	float r, g, b;
	// attenuation
	float a0, a1, a2;

	EmissiveMaterial()
	{
		type = MaterialType::EMISSIVE;
	}
};

/**
* Stack for reference keeping of materials.
* The last one is always the currently used one.
* This stack should be purged during sglBeginScene.
*/
std::vector<Material*> materialStack;
/**
* Stack for reference keeping of textures since multiple
* objects and materials can have the same texture.
* The last one is always the currently used one.
* This stack should be purged during sglBeginScene.
*/
std::vector<Texture*> textureStack;

struct PointLight
{
	float x, y, z;
	float r, g, b;
};

/**
* Contains all lights in the scene.
* All lights will be deleted during sglBeginScene call.
*/
std::vector<PointLight*> lightStack;

struct Primitive
{
public:
	enum PrimitiveType
	{
		POLYGON,
		SPHERE
	};

	Material *mat;
	Texture *tex;
	PrimitiveType type;
protected:

	Primitive() 
	{
	}
};

struct Polygon : public Primitive
{
	std::deque<inputPoint4f*> points;

	~Polygon()
	{
		while (!points.empty())
		{
			delete points.front();
			points.pop_front();
		}
	}
};

struct Sphere : public Primitive
{
	float x, y, z, radius;
};

/**
* Stack, where spheres for raytracing are stored.
* They should be purged after ray tracing or during
* sglBeginScene call.
*/
std::vector<Sphere*> sphereStack;

bool gluInvertMatrix(const float m[16], float invOut[16])
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

	det = 1.0f / det;

	for (i = 0; i < 16; i++)
		invOut[i] = inv[i] * det;

	return true;
}

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

void multiplyMatrix(float* left, const float* right) {
	float output[16];
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			output[i*4 + j] = dotVectors(left, right, j, i);
		}
	}
	copyMatrix(left, output);
}

float viewportMatrix[16];
float matrixMVP[16];
float multipliedMatrix[16];
float zNear, zFar;

/**
* It provides information, whether scene is being constructed. That
* means, that during sglEnd call, no drawing is done. It will be done
* during call to sglSceneEnd.
*/
bool constructingScene = false;

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
sglEAreaMode areaMode = SGL_FILL;

short pointSize = 0;
float colorVertexR = 0, colorVertexG = 0, colorVertexB = 0;
float colorClearR = 0, colorClearG = 0, colorClearB = 0;

bool depthEnabled = false;
sglEElementType drawingMethod = SGL_POINTS;

//std::queue<inputPoint4f*> queue4f;
std::deque<Polygon*> polygonQueue;
std::vector<Polygon*> emissivePolygonStack;


class SglContext {
private:
	int width;
	int height;
	/**
	Color buffer saves integers from 0 to 1.
	*/
	float *colorBuffer;
	float *depthBuffer;
	Texture *environmentMap;
public:
	SglContext(int width, int height) : width{ width }, height{ height } {
		colorBuffer = new float[width*height * 3];
		depthBuffer = new float[width*height];
		environmentMap = NULL;
	};
	~SglContext() {
		delete[] colorBuffer;
		delete[] depthBuffer;
	}

	void SetEnvironmentMap(Texture *map)
	{
		environmentMap = map;
	}

	Texture* GetEnvironmentMap()
	{
		return environmentMap;
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
		for (int i = 0; i < height; ++i) {
			for (int j = 0; j < width; ++j) {
				*(depthBuffer + i*width + j) = FLT_MAX;
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

struct Ray
{
	float *start;
	float *dir;
	float length;
};