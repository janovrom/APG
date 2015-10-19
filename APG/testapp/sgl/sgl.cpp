//---------------------------------------------------------------------------
// sgl.cpp
// Empty implementation of the SGL (Simple Graphics Library)
// Date:  2011/11/1
// Author: Jaroslav Krivanek, Jiri Bittner CTU Prague
//---------------------------------------------------------------------------

#include "sgl.h"
#include "sglcontext.h"

/// Current error code.
static sglEErrorCode _libStatus = SGL_NO_ERROR;

static inline void setErrCode(sglEErrorCode c) 
{
  if(_libStatus==SGL_NO_ERROR)
    _libStatus = c;
}

//---------------------------------------------------------------------------
// sglGetError()
//---------------------------------------------------------------------------
sglEErrorCode sglGetError(void) 
{
  sglEErrorCode ret = _libStatus;
  _libStatus = SGL_NO_ERROR;
  return ret;
}

//---------------------------------------------------------------------------
// sglGetErrorString()
//---------------------------------------------------------------------------
const char* sglGetErrorString(sglEErrorCode error)
{
  static const char *errStrigTable[] = 
  {
      "Operation succeeded",
      "Invalid argument(s) to a call",
      "Invalid enumeration argument(s) to a call",
      "Invalid call",
      "Quota of internal resources exceeded",
      "Internal library error",
      "Matrix stack overflow",
      "Matrix stack underflow",
      "Insufficient memory to finish the requested operation"
  };

  if((int)error<(int)SGL_NO_ERROR || (int)error>(int)SGL_OUT_OF_MEMORY ) {
    return "Invalid value passed to sglGetErrorString()"; 
  }

  return errStrigTable[(int)error];
}

//---------------------------------------------------------------------------
// Initialization functions
//---------------------------------------------------------------------------

int ContextWrapper::add(SglContext* c) {
	int i = findFirstEmpty();
	if (i == -1) {
		setErrCode(sglEErrorCode::SGL_OUT_OF_RESOURCES);
		return -1;
	}
	else {
		contexts[i] = c;
		++count;
		return i;
	}
}

void sglInit(void) {
	hasBegun = false;
	offsetX = offsetY = windowWidth = windowHeight = 0;
}

void sglFinish(void) {
	contextWrapper.clear();
	contextWrapper.activeContext = -1;
	hasBegun = false;
	offsetX = offsetY = windowWidth = windowHeight = 0;
}

int sglCreateContext(int width, int height) {
	SglContext* c = new SglContext(width, height);
	if (!c) {
		setErrCode(sglEErrorCode::SGL_OUT_OF_MEMORY);
		return -2;
	}
	
	
	return contextWrapper.add(c);
}

void sglDestroyContext(int id) {
	contextWrapper.clear(id);
}

void sglSetContext(int id) {
	if (id < contextWrapper.size())
		contextWrapper.activeContext = id;
	else
		setErrCode(sglEErrorCode::SGL_INVALID_VALUE);
}

int sglGetContext(void) {
	return contextWrapper.activeContext;
}

float *sglGetColorBufferPointer(void) {
	if (contextWrapper.activeContext == -1)
		return 0;

	return contextWrapper[contextWrapper.activeContext]->getColorBuffer();
}

//---------------------------------------------------------------------------
// Drawing functions
//---------------------------------------------------------------------------

void sglClearColor (float r, float g, float b, float alpha) {
	if (contextWrapper.empty() || hasBegun) {
		setErrCode(sglEErrorCode::SGL_INVALID_OPERATION);
	}
	else {
		colorClearR = r;
		colorClearG = g;
		colorClearB = b;
	}
}

void sglClear(unsigned what) {
	if (contextWrapper.empty() || hasBegun) {
		setErrCode(sglEErrorCode::SGL_INVALID_OPERATION);
		return;
	}

	if ((what & SGL_COLOR_BUFFER_BIT) == SGL_COLOR_BUFFER_BIT) {
			contextWrapper[contextWrapper.activeContext]->clearColor(colorClearR, colorClearG, colorClearB );
	}
	else if ((what & SGL_DEPTH_BUFFER_BIT) == SGL_DEPTH_BUFFER_BIT) {
			contextWrapper[contextWrapper.activeContext]->clearDepth();
	}
	else {
			setErrCode(sglEErrorCode::SGL_INVALID_VALUE);
	}
}

void sglBegin(sglEElementType mode) 
{
	if (hasBegun) { setErrCode(sglEErrorCode::SGL_INVALID_OPERATION); return; }
	if (mode <= 0 || mode >= sglEElementType::SGL_LAST_ELEMENT_TYPE) { setErrCode(sglEErrorCode::SGL_INVALID_ENUM); return; }
	hasBegun = true;
	drawingMethod = mode;
}

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

/*
Transformations of points will be applied here in future, now it just returns input.
*/
inputPoint4f* transformThePoint(inputPoint4f& point)
{
	//printf("IMPLEMENT ME: sgl.cpp -> transformThePoint \n");
	inputPoint4f *ret = new inputPoint4f;
	ret->r = point.r;
	ret->g = point.g;
	ret->b = point.b;
	ret->a = point.a;
	switch (matrixMode) {
	case SGL_MODELVIEW:
		multiplyMatrixVector(modelViewStack.top(), point, *ret);
		break;
	case SGL_PROJECTION:
		multiplyMatrixVector(projectionStack.top(), point, *ret);
		break;
	default:
		break;
	}
	return ret;
}

void drawMeAPoint(inputPoint4f& point) 
{
	inputPoint4f* transformed = transformThePoint(point);

	int W, H, x, y;

	SglContext *cont = contextWrapper.contexts[contextWrapper.activeContext];



	W = cont->getWidth();
	H = cont->getHeight();
	x = (int)transformed->x;
	y = (int)transformed->y;

	float *colorBuffer = cont->getColorBuffer();
	int offset;

	int size = (int)((pointSize-1) / 2);

	for (int i = x - size; i <= x + size; i++)
	{
		for (int j = y - size; j <= y + size; j++)
		{
			if (i >= 0 && i < W && j >= 0 && j < H)
			{
				offset = j*W * 3 + i;
				*(colorBuffer + offset) = transformed->r;
				*(colorBuffer + offset + 1) = transformed->g;
				*(colorBuffer + offset + 2) = transformed->b;
			}
		}
	}

	delete transformed;
}

void drawMeALine(inputPoint4f& start, inputPoint4f& end)
{
	//printf("IMPLEMENT ME: sgl.cpp -> drawMeALine \n");
}

void drawPoints() 
{
	inputPoint4f tempPoint;

	while (!queue4f.empty())
	{
		tempPoint = queue4f.front();
		drawMeAPoint(tempPoint);
		queue4f.pop();
	}
}

void drawLines()
{
	inputPoint4f tempPoint1;
	inputPoint4f tempPoint2;

	while (!queue4f.empty())
	{
		tempPoint1 = queue4f.front();
		queue4f.pop();

		if (queue4f.empty()) { break; }

		tempPoint2 = queue4f.front();
		queue4f.pop();

		drawMeALine(tempPoint1, tempPoint2);
	}
}

void drawLineStrip()
{
	inputPoint4f tempPoint1;
	inputPoint4f tempPoint2;

	if (queue4f.empty()) { return; }
	tempPoint2 = queue4f.front();
	queue4f.pop();

	while (!queue4f.empty())
	{
		tempPoint1 = tempPoint2;
		tempPoint2 = queue4f.front();
		queue4f.pop();

		drawMeALine(tempPoint1, tempPoint2);
	}
}

void drawLineLoop()
{
	inputPoint4f origin;
	inputPoint4f tempPoint1;
	inputPoint4f tempPoint2;

	if (queue4f.empty()) { return; }
	tempPoint2 = origin = queue4f.front();
	queue4f.pop();
	int counter = 1;

	while (!queue4f.empty())
	{
		tempPoint1 = tempPoint2;
		tempPoint2 = queue4f.front();
		queue4f.pop();
		counter++;

		drawMeALine(tempPoint1, tempPoint2);
	}

	if (counter > 1)
	{
		drawMeALine(tempPoint2, origin);
	}
}

void sglEnd(void) 
{
	if (!hasBegun) { setErrCode(sglEErrorCode::SGL_INVALID_OPERATION); return; }

	switch (drawingMethod)
	{
	case sglEElementType::SGL_POINTS:
		drawPoints();
		break;
	default:
		break;
	}


	hasBegun = false;
}

void sglVertex4f(float x, float y, float z, float w) 
{
	inputPoint4f *point = new inputPoint4f;
	(*point).x = x;
	(*point).y = y;
	(*point).z = z;
	(*point).z = w;

	(*point).r = colorVertexR;
	(*point).g = colorVertexG;
	(*point).b = colorVertexB;

	queue4f.push(*point);
}

void sglVertex3f(float x, float y, float z) 
{
	inputPoint4f *point = new inputPoint4f;
	(*point).x = x;
	(*point).y = y;
	(*point).z = z;
	(*point).z = 1;

	(*point).r = colorVertexR;
	(*point).g = colorVertexG;
	(*point).b = colorVertexB;

	queue4f.push(*point);
}

void sglVertex2f(float x, float y) 
{
	inputPoint4f *point = new inputPoint4f;
	(*point).x = x;
	(*point).y = y;
	(*point).z = 0;
	(*point).z = 1;

	(*point).r = colorVertexR;
	(*point).g = colorVertexG;
	(*point).b = colorVertexB;

	queue4f.push(*point);
}

void sglCircle(float x, float y, float z, float radius) {}

void sglEllipse(float x, float y, float z, float a, float b) {}

void sglArc(float x, float y, float z, float radius, float from, float to) {}

//---------------------------------------------------------------------------
// Transform functions
//---------------------------------------------------------------------------

void sglMatrixMode( sglEMatrixMode mode ) {
	if (hasBegun || contextWrapper.empty()) {
		setErrCode(sglEErrorCode::SGL_INVALID_OPERATION);
		return;
	}
	if (mode != SGL_MODELVIEW || mode != SGL_PROJECTION) {
		setErrCode(SGL_INVALID_ENUM);
		return;
	}
	matrixMode = mode;
}

float* duplicateMatrix(const float* matrix) {
	float* ret = new float[16];
	for (int i = 0; i < 16; ++i) {
		ret[i] = matrix[i];
	}
	return ret;
}

void sglPushMatrix(void) {
	if (hasBegun || contextWrapper.empty()) {
		setErrCode(sglEErrorCode::SGL_INVALID_OPERATION);
		return;
	}

	switch (matrixMode) {
	case SGL_MODELVIEW:
		modelViewStack.push(duplicateMatrix(modelViewStack.top()));
		break;
	case SGL_PROJECTION:
		projectionStack.push(duplicateMatrix(projectionStack.top()));
		break;
	default:
		break;
	}
}

void sglPopMatrix(void) {
	if (hasBegun || contextWrapper.empty()) {
		setErrCode(sglEErrorCode::SGL_INVALID_OPERATION);
		return;
	}

	switch (matrixMode) {
	case SGL_MODELVIEW:
		if (modelViewStack.size() == 1) {
			setErrCode(SGL_STACK_UNDERFLOW);
		}
		else {
			delete[] modelViewStack.top();
			modelViewStack.pop();
		}
		break;
	case SGL_PROJECTION:
		if (projectionStack.size() == 1) {
			setErrCode(SGL_STACK_UNDERFLOW);
		}
		else {
			delete[] projectionStack.top();
			projectionStack.pop();
		}
		break;
	default:
		break;
	}
}

void copyMatrix(float* output, const float* input) {
	for (int i = 0; i < 16; ++i) {
		output[i] = input[i];
	}
}

void sglLoadIdentity(void) {
	if (hasBegun || contextWrapper.empty()) {
		setErrCode(sglEErrorCode::SGL_INVALID_OPERATION);
		return;
	}

	switch (matrixMode) {
	case SGL_MODELVIEW:
		copyMatrix(modelViewStack.top(), identity);
		break;
	case SGL_PROJECTION:
		copyMatrix(projectionStack.top(), identity);
		break;
	default:
		break;
	}
}

void sglLoadMatrix(const float *matrix) {
	if (hasBegun || contextWrapper.empty()) {
		setErrCode(sglEErrorCode::SGL_INVALID_OPERATION);
		return;
	}

	switch (matrixMode) {
	case SGL_MODELVIEW:
		copyMatrix(modelViewStack.top(), matrix);
		break;
	case SGL_PROJECTION:
		copyMatrix(projectionStack.top(), matrix);
		break;
	default:
		break;
	}
}

float dotVectors(const float* left, const float* right, int i, int j) {
	float tmp = 0;
	tmp += left[i] * right[j];
	tmp += left[i+4] * right[j+1];
	tmp += left[i+8] * right[j+2];
	tmp += left[i+12] * right[j+3];

	return tmp;
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

void sglMultMatrix(const float *matrix) {
	if (hasBegun || contextWrapper.empty()) {
		setErrCode(sglEErrorCode::SGL_INVALID_OPERATION);
		return;
	}

	switch (matrixMode) {
	case SGL_MODELVIEW:
		multiplyMatrix(modelViewStack.top(), matrix);
		break;
	case SGL_PROJECTION:
		multiplyMatrix(projectionStack.top(), matrix);
		break;
	default:
		break;
	}
}

void sglTranslate(float x, float y, float z) {
	if (hasBegun || contextWrapper.empty()) {
		setErrCode(sglEErrorCode::SGL_INVALID_OPERATION);
		return;
	}
	float translate[16];
	copyMatrix(translate, identity);
	translate[12] = x;
	translate[13] = y;
	translate[14] = z;
	switch (matrixMode) {
	case SGL_MODELVIEW:
		multiplyMatrix(modelViewStack.top(), translate);
		break;
	case SGL_PROJECTION:
		multiplyMatrix(projectionStack.top(), translate);
		break;
	default:
		break;
	}
}

void sglScale(float scalex, float scaley, float scalez) {
	if (hasBegun || contextWrapper.empty()) {
		setErrCode(sglEErrorCode::SGL_INVALID_OPERATION);
		return;
	}
	float scale[16];
	copyMatrix(scale, identity);
	scale[0] = scalex;
	scale[5] = scaley;
	scale[10] = scalez;
	switch (matrixMode) {
	case SGL_MODELVIEW:
		multiplyMatrix(modelViewStack.top(), scale);
		break;
	case SGL_PROJECTION:
		multiplyMatrix(projectionStack.top(), scale);
		break;
	default:
		break;
	}
}

void sglRotate2D(float angle, float centerx, float centery) {
	if (hasBegun || contextWrapper.empty()) {
		setErrCode(sglEErrorCode::SGL_INVALID_OPERATION);
		return;
	}

	float tmp[16];
	float rotate[16];
	copyMatrix(tmp, identity);
	copyMatrix(rotate, identity);
	tmp[12] = centerx;
	tmp[13] = centery;
	rotate[0] = cos(angle); rotate[4] = -sin(angle);
	rotate[1] = sin(angle); rotate[5] = cos(angle);

	multiplyMatrix(rotate, tmp);

	tmp[12] = -centerx;
	tmp[13] = -centery;

	multiplyMatrix(tmp, rotate);

	switch (matrixMode) {
	case SGL_MODELVIEW:
		multiplyMatrix(modelViewStack.top(), tmp);
		break;
	case SGL_PROJECTION:
		multiplyMatrix(projectionStack.top(), tmp);
		break;
	default:
		break;
	}
}

void sglRotateY(float angle) {
	if (hasBegun || contextWrapper.empty()) {
		setErrCode(sglEErrorCode::SGL_INVALID_OPERATION);
		return;
	}

	float rotate[16];
	copyMatrix(rotate, identity);
	rotate[0] = cos(angle); rotate[8] = -sin(angle);
	rotate[2] = sin(angle); rotate[10] = cos(angle);

	switch (matrixMode) {
	case SGL_MODELVIEW:
		multiplyMatrix(modelViewStack.top(), rotate);
		break;
	case SGL_PROJECTION:
		multiplyMatrix(projectionStack.top(), rotate);
		break;
	default:
		break;
	}
}

void sglOrtho(float left, float right, float bottom, float top, float near, float far) {
	if (hasBegun || contextWrapper.empty()) {
		setErrCode(sglEErrorCode::SGL_INVALID_OPERATION);
		return;
	}

	float ortho[16];
	copyMatrix(ortho, identity);
	ortho[0] = 2.0 / (right - left);
	ortho[5] = 2.0 / (top - bottom);
	ortho[10] = - 2.0 / (far - near);
	ortho[12] = - (right + left) / (right - left);
	ortho[13] = - (top + bottom) / (top - bottom);
	ortho[14] = -(far + near) / (far - near);

	switch (matrixMode) {
	case SGL_MODELVIEW:
		multiplyMatrix(modelViewStack.top(), ortho);
		break;
	case SGL_PROJECTION:
		multiplyMatrix(projectionStack.top(), ortho);
		break;
	default:
		break;
	}
}

void sglFrustum(float left, float right, float bottom, float top, float near, float far) {
	if (near < 0 || far < 0) {
		setErrCode(sglEErrorCode::SGL_INVALID_VALUE);
	}
	else if (hasBegun || contextWrapper.empty()) {
		setErrCode(sglEErrorCode::SGL_INVALID_OPERATION);
	}

	float A = (right + left) / (right - left);
	float B = (top + bottom) / (top - bottom);
	float C = (far + near) / (far - near);
	float D = - 2 * far * near / (far - near);

	float persp[16];
	copyMatrix(persp, identity);
	persp[0] = 2.0 * near/ (right - left);
	persp[5] = 2.0 * near / (top - bottom);
	persp[8] = A;
	persp[9] = B;
	persp[10] = C;
	persp[11] = -1.0;
	persp[14] = D;
	persp[15] = 0.0;

	switch (matrixMode) {
	case SGL_MODELVIEW:
		multiplyMatrix(modelViewStack.top(), persp);
		break;
	case SGL_PROJECTION:
		multiplyMatrix(projectionStack.top(), persp);
		break;
	default:
		break;
	}
}

void sglViewport(int x, int y, int width, int height) {
	if (width < 0 || height < 0) {
		setErrCode(sglEErrorCode::SGL_INVALID_VALUE);
	}
	else if (hasBegun || contextWrapper.empty()) {
		setErrCode(sglEErrorCode::SGL_INVALID_OPERATION);
	}
	else {
		offsetX = x;
		offsetY = y;
		windowWidth = width;
		windowHeight = height;
	}
}

//---------------------------------------------------------------------------
// Attribute functions
//---------------------------------------------------------------------------

void sglColor3f(float r, float g, float b) 
{
	colorVertexR = r;
	colorVertexG = g;
	colorVertexB = b;
}

void sglAreaMode(sglEAreaMode mode) {}

void sglPointSize(float size) 
{
	if (size <= 0) { setErrCode(sglEErrorCode::SGL_INVALID_VALUE); return; }
	if (hasBegun || contextWrapper.empty() ) { setErrCode(sglEErrorCode::SGL_INVALID_OPERATION); return; }
	pointSize = size;
}

void sglEnable(sglEEnableFlags cap) {
	if (hasBegun || contextWrapper.empty()) {
		setErrCode(sglEErrorCode::SGL_INVALID_OPERATION);
		return;
	}
	switch (cap)
	{
	case SGL_DEPTH_TEST:
		testDepth = true;
		break;
	default:
		setErrCode(sglEErrorCode::SGL_INVALID_ENUM);
		break;
	}
}

void sglDisable(sglEEnableFlags cap) {
	if (hasBegun || contextWrapper.empty()) {
		setErrCode(sglEErrorCode::SGL_INVALID_OPERATION);
		return;
	}
	switch (cap)
	{
	case SGL_DEPTH_TEST:
		testDepth = false;
		break;
	default:
		setErrCode(sglEErrorCode::SGL_INVALID_ENUM);
		break;
	}
}

//---------------------------------------------------------------------------
// RayTracing oriented functions
//---------------------------------------------------------------------------

void sglBeginScene() {}

void sglEndScene() {}

void sglSphere(const float x,
			   const float y,
			   const float z,
			   const float radius) {}

void sglMaterial(const float r,
				 const float g,
				 const float b,
				 const float kd,
				 const float ks,
				 const float shine,
				 const float T,
				 const float ior) {}

void sglPointLight(const float x,
				   const float y,
				   const float z,
				   const float r,
				   const float g,
				   const float b) {}

void sglRayTraceScene() {}

void sglRasterizeScene() {}

void sglEnvironmentMap(const int width,
					   const int height,
					   float *texels)
{}

void sglEmissiveMaterial(
						 const float r,
						 const float g,
						 const float b,
						 const float c0,
						 const float c1,
						 const float c2
						 )
{}
