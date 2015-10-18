//---------------------------------------------------------------------------
// sgl.cpp
// Empty implementation of the SGL (Simple Graphics Library)
// Date:  2011/11/1
// Author: Jaroslav Krivanek, Jiri Bittner CTU Prague
//---------------------------------------------------------------------------

#include "sgl.h"
#include "sglcontext.h"

void ContextWrapper::add(SglContext* c) {
	int i = findFirstEmpty();
	if (i == -1) {
		setErrCode(sglEErrorCode::SGL_OUT_OF_RESOURCES);
	}
	else {
		contexts[i] = c;
		++count;
	}
}

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

void sglInit(void) {
	hasBegun = false;
	offsetX = offsetY = windowWidth = windowHeight = 0;
}

void sglFinish(void) {
	contexts.clear();
	contexts.activeContext = -1;
	hasBegun = false;
	offsetX = offsetY = windowWidth = windowHeight = 0;
}

int sglCreateContext(int width, int height) {
	SglContext* c = new SglContext(width, height);
	if (!c) {
		setErrCode(sglEErrorCode::SGL_OUT_OF_MEMORY);
		return -2;
	}
	
	contexts.add(c);
	return contexts.size() - 1;
}

void sglDestroyContext(int id) {
	contexts.clear(id);
}

void sglSetContext(int id) {
	if (id < contexts.size())
		contexts.activeContext = id;
	else
		setErrCode(sglEErrorCode::SGL_INVALID_VALUE);
}

int sglGetContext(void) {
	return contexts.activeContext;
}

float *sglGetColorBufferPointer(void) {
	if (contexts.activeContext == -1)
		return 0;

	return contexts[contexts.activeContext]->getColorBuffer();
}

//---------------------------------------------------------------------------
// Drawing functions
//---------------------------------------------------------------------------

void sglClearColor (float r, float g, float b, float alpha) {
	if (contexts.empty() || hasBegun || contexts.activeContext == -1)
		setErrCode(sglEErrorCode::SGL_INVALID_OPERATION);
	else
		contexts[contexts.activeContext]->setClearColor(r, g, b, alpha);
}

void sglClear(unsigned what) {
	if (contexts.empty() || hasBegun || contexts.activeContext == -1) {
		setErrCode(sglEErrorCode::SGL_INVALID_OPERATION);
		return;
	}

	if ((what & SGL_COLOR_BUFFER_BIT) == SGL_COLOR_BUFFER_BIT) {
			contexts[contexts.activeContext]->clearColor();
	}
	else if ((what & SGL_DEPTH_BUFFER_BIT) == SGL_DEPTH_BUFFER_BIT) {
			contexts[contexts.activeContext]->clearDepth();
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

void sglEnd(void) 
{
	if (!hasBegun) { setErrCode(sglEErrorCode::SGL_INVALID_OPERATION); return; }
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
	inputPoint3f *point = new inputPoint3f;
	(*point).x = x;
	(*point).y = y;
	(*point).z = z;

	(*point).r = colorVertexR;
	(*point).g = colorVertexG;
	(*point).b = colorVertexB;

	queue3f.push(*point);
}

void sglVertex2f(float x, float y) 
{
	inputPoint2f *point = new inputPoint2f;
	(*point).x = x;
	(*point).y = y;

	(*point).r = colorVertexR;
	(*point).g = colorVertexG;
	(*point).b = colorVertexB;

	queue2f.push(*point);
}

void sglCircle(float x, float y, float z, float radius) {}

void sglEllipse(float x, float y, float z, float a, float b) {}

void sglArc(float x, float y, float z, float radius, float from, float to) {}

//---------------------------------------------------------------------------
// Transform functions
//---------------------------------------------------------------------------

void sglMatrixMode( sglEMatrixMode mode ) {}

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
	else if (hasBegun || contexts.activeContext == -1 || contexts.empty()) {
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

void sglColor3f(float r, float g, float b) {}

void sglAreaMode(sglEAreaMode mode) {}

void sglPointSize(float size) 
{
	pointSize = size;
}

void sglEnable(sglEEnableFlags cap) 
{
	if (hasBegun) { setErrCode(sglEErrorCode::SGL_INVALID_OPERATION); return; }
	switch (cap) {
	case sglEEnableFlags::SGL_DEPTH_TEST:
		depthEnabled = true;
		break;
	default:
		setErrCode(sglEErrorCode::SGL_INVALID_ENUM);
		break;
	}
}

void sglDisable(sglEEnableFlags cap) {
	switch (cap)
	{
	case SGL_DEPTH_TEST:
		testDepth = false;
		break;
	default:
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
