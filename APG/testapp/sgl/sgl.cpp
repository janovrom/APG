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

/*
Transformations of points will be applied here in future, now it just returns input.
*/
inputPoint4f transformThePoint(inputPoint4f& point)
{
	//printf("IMPLEMENT ME: sgl.cpp -> transformThePoint \n");
	return point;
}

void drawMeAPoint(inputPoint4f& point) 
{
	inputPoint4f transformed = transformThePoint(point);

	int W, H, x, y;

	SglContext *cont = contextWrapper.contexts[contextWrapper.activeContext];
	W = cont->getWidth();
	H = cont->getHeight();
	x = (int)point.x;
	y = (int)point.y;
	float *colorBuffer = cont->getColorBuffer();

	int offset;

	int size = (int)((pointSize-1) / 2);
	int sizeCorrection = 1 - (int) pointSize % 2;

	for (int i = x - size; i <= x + (size + sizeCorrection); i++)
	{
		for (int j = y - size; j <= y + (size + sizeCorrection); j++)
		{
			if (i >= 0 && i < W && j >= 0 && j < H)
			{
				offset = j*W * 3 + i*3;
				
				*(colorBuffer + offset) = point.r;
				*(colorBuffer + offset + 1) = point.g;
				*(colorBuffer + offset + 2) = point.b;
				
			}
		}
	}
}

void drawMeALineNaive(inputPoint4f& start, inputPoint4f& end)
{
	int W, H, tempX, tempY, offset;

	SglContext *cont = contextWrapper.contexts[contextWrapper.activeContext];
	W = cont->getWidth();
	H = cont->getHeight();
	float *colorBuffer = cont->getColorBuffer();

	inputPoint4f startT = transformThePoint(start);
	inputPoint4f endT = transformThePoint(end);

	/*if (startT.x > endT.x)
	{
		inputPoint4f swapT = startT;
		startT = endT;
		endT = swapT;
	}*/

	int x0, x1, y0, y1;
	float k;

	x0 = startT.x;
	y0 = startT.y;
	x1 = endT.x;
	y1 = endT.y;


	int endValue;
	float lerpValue;
	if (std::abs((float)(y1 - y0) / (x1 - x0)) < 1){
		//x sampling
		if (x0 > x1)
		{
			int swap;
			swap = x0;
			x0 = x1;
			x1 = swap;

			swap = y0;
			y0 = y1;
			y1 = swap;
		}
		k = (float)(y1 - y0) / (x1 - x0);
		endValue = x1 - x0;
		for (int i = 0; i <= endValue; i++){
			//setPixel(x0 + i, y0 + i*k);
			lerpValue = (float)i / endValue;
			tempX = x0 + i;
			tempY = y0 + i*k;

			if (tempX >= 0 && tempX < W && tempY >= 0 && tempY < H)
			{
				offset = tempY*W * 3 + tempX*3;

				*(colorBuffer + offset) = (lerpValue)*startT.r + (1-lerpValue)*endT.r;
				*(colorBuffer + offset + 1) = (lerpValue)*startT.g + (1 - lerpValue)*endT.g;
				*(colorBuffer + offset + 2) = (lerpValue)*startT.b + (1 - lerpValue)*endT.b;
			}
		}
	}else {
		//y sampling
		if (y0 > y1)
		{
			int swap;
			swap = x0;
			x0 = x1;
			x1 = swap;

			swap = y0;
			y0 = y1;
			y1 = swap;
		}
		k = (float)(x1 - x0) / (y1 - y0);
		endValue = y1 - y0;
		for (int i = 0; i <= endValue; i++) {
			//setPixel(x0 + i*k, y0 + i);
			lerpValue = (float)i / endValue;
			tempX = x0 + i*k;
			tempY = y0 + i;

			if (tempX >= 0 && tempX < W && tempY >= 0 && tempY < H)
			{
				offset = tempY*W * 3 + tempX * 3;

				*(colorBuffer + offset) = (lerpValue)*startT.r + (1 - lerpValue)*endT.r;
				*(colorBuffer + offset + 1) = (lerpValue)*startT.g + (1 - lerpValue)*endT.g;
				*(colorBuffer + offset + 2) = (lerpValue)*startT.b + (1 - lerpValue)*endT.b;
			}
		}
	}
	//printf("IMPLEMENT ME: sgl.cpp -> drawMeALine \n");
}

void drawMeALineBresenham(inputPoint4f& start, inputPoint4f& end)

{
	int W, H, offset;
	int x0, x1, y0, y1, swap;

	SglContext *cont = contextWrapper.contexts[contextWrapper.activeContext];
	W = cont->getWidth();
	H = cont->getHeight();
	float *colorBuffer = cont->getColorBuffer();

	inputPoint4f *startT = &transformThePoint(start);
	inputPoint4f *endT = &transformThePoint(end);

	x0 = (*startT).x;
	y0 = (*startT).y;
	x1 = (*endT).x;
	y1 = (*endT).y;

	int dx, dy;

	bool drivingX = std::abs(x1 - x0) > std::abs(y1 - y0);

	if (!drivingX)
	{
		swap = x0;
		x0 = y0;
		y0 = swap;

		swap = x1;
		x1 = y1;
		y1 = swap;
	}

	if (x0 > x1)
	{
		swap = x0;
		x0 = x1;
		x1 = swap;

		swap = y0;
		y0 = y1;
		y1 = swap;

		inputPoint4f *swap4f = startT;
		startT = endT;
		endT = swap4f;
	}

	//dx = x1 - x0;
	dx = std::abs(x1 - x0);
	dy = y1 - y0;
	int yAdd = (dy < 0) ? -1 : 1;
	dy = std::abs(dy);

	int p = 2 * dy - dx;
	int c0 = 2 * dy;
	int c1 = c0 - 2 * dx;
	int tempY = y0;
	float lerpValue;


	if (drivingX)
	{
		if (x0 >= 0 && x0 < W && y0 >= 0 && y0 < H)
		{
			offset = y0*W * 3 + x0 * 3;

			*(colorBuffer + offset) = (*startT).r;
			*(colorBuffer + offset + 1) = (*startT).g;
			*(colorBuffer + offset + 2) = (*startT).b;
		}

	}
	else {
		if (y0 >= 0 && y0 < W && x0 >= 0 && x0 < H)
		{
			offset = x0*W * 3 + y0 * 3;

			*(colorBuffer + offset) = (*startT).r;
			*(colorBuffer + offset + 1) = (*startT).g;
			*(colorBuffer + offset + 2) = (*startT).b;
		}
	}


	for (int tempX = x0 + 1; tempX <= x1; tempX++)
	{
		lerpValue = (float)(tempX - x0) / (x1 - x0);
		if (p < 0)
		{
			p += c0;
			//lerpValue = 0;
		}else {
			p += c1;
			tempY += yAdd;
			//lerpValue = 1;
		}
		if(drivingX)
		{
			if (tempX >= 0 && tempX < W && tempY >= 0 && tempY < H)
			{
				offset = tempY*W * 3 + tempX * 3;

				*(colorBuffer + offset) = (1-lerpValue)*(*startT).r + (lerpValue)*(*endT).r;
				*(colorBuffer + offset + 1) = (1-lerpValue)*(*startT).g + (lerpValue)*(*endT).g;
				*(colorBuffer + offset + 2) = (1-lerpValue)*(*startT).b + (lerpValue)*(*endT).b;
			}

		}else{
			if (tempY >= 0 && tempY < W && tempX >= 0 && tempX < H)
			{
				offset = tempX*W * 3 + tempY * 3;

				*(colorBuffer + offset) = (1-lerpValue)*(*startT).r + (lerpValue)*(*endT).r;
				*(colorBuffer + offset + 1) = (1-lerpValue)*(*startT).g + (lerpValue)*(*endT).g;
				*(colorBuffer + offset + 2) = (1-lerpValue)*(*startT).b + (lerpValue)*(*endT).b;
			}
		}
	}

}

void drawMeALine(inputPoint4f& start, inputPoint4f& end)
{
	if(true)
	{
		drawMeALineBresenham(start, end);
	}else{
		drawMeALineNaive(start, end);
	}
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
	case sglEElementType::SGL_LINES:
		drawLines();
		break;
	case sglEElementType::SGL_LINE_STRIP:
		drawLineStrip();
		break;
	case sglEElementType::SGL_LINE_LOOP:
		drawLineLoop();
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
	(*point).w = w;

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
	(*point).w = 1;

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
	(*point).w = 1;

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

void sglPushMatrix(void) {}

void sglPopMatrix(void) {}

void sglLoadIdentity(void) {}

void sglLoadMatrix(const float *matrix) {}

void sglMultMatrix(const float *matrix) {}

void sglTranslate(float x, float y, float z) {}

void sglScale(float scalex, float scaley, float scalez) {}

void sglRotate2D(float angle, float centerx, float centery) {}

void sglRotateY(float angle) {}

void sglOrtho(float left, float right, float bottom, float top, float near, float far) {}

void sglFrustum(float left, float right, float bottom, float top, float near, float far) {}

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
