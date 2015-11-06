//---------------------------------------------------------------------------
// sgl.cpp
// Empty implementation of the SGL (Simple Graphics Library)
// Date:  2011/11/1
// Author: Jaroslav Krivanek, Jiri Bittner CTU Prague
//---------------------------------------------------------------------------

#include "sgl.h"
#include "sglcontext.h"



//#define LINE_NAIVE
// decides which ellipse algoritm should be used
#define ELLIPSE

using namespace std;

void setPixel(float x0, float y0, float r, float g, float b);
void drawMeALine(inputPoint4f* start, inputPoint4f* end);
void drawPoints();
void drawLineStrip();
void drawLineLoop();

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
		setErrCode(SGL_OUT_OF_RESOURCES);
		return -1;
	}
	else {
		contexts[i] = c;
		++count;
		return i;
	}
}

void sglInit(void) {
	// init values to default
	hasBegun = false;
	viewportOffsetX = viewportOffsetY = viewportWidth = viewportHeight = 0;

	float *mv;
	float *proj;
	try
	{
		mv = new float[16];
		proj = new float[16];
	}
	catch (std::bad_alloc& ba)
	{
		setErrCode(SGL_OUT_OF_MEMORY);
		return;
	}
	// init stacks with identityMatrix matrices
	copyMatrix(mv, identityMatrix);
	copyMatrix(proj, identityMatrix);
	modelViewStack.push(mv);
	projectionStack.push(proj);
}

void sglFinish(void) {
	// reset values and clear
	contextWrapper.clear();
	contextWrapper.activeContext = -1;
	hasBegun = false;
	viewportOffsetX = viewportOffsetY = viewportWidth = viewportHeight = 0;
	// do some memory management on stacks
	for (;!modelViewStack.empty();) {
		delete[] modelViewStack.top();
		modelViewStack.pop();
}
	for (; !projectionStack.empty();) {
		delete[] projectionStack.top();
		projectionStack.pop();
	}
	for (;!queue4f.empty();) {
		delete &queue4f.front();
		queue4f.pop();
	}
}

int sglCreateContext(int width, int height) {
	SglContext* c = new SglContext(width, height);
	if (!c) {
		setErrCode(SGL_OUT_OF_MEMORY);
		return -2;
	}
	
	
	return contextWrapper.add(c);
}

void sglDestroyContext(int id) {
	if (!contextWrapper.clear(id))
	{
		setErrCode(SGL_INVALID_VALUE);
	}
}

void sglSetContext(int id) {
	if (id < contextWrapper.size() && id >= 0)
		contextWrapper.activeContext = id;
	else
		setErrCode(SGL_INVALID_VALUE);
}

int sglGetContext(void) {
	if (contextWrapper.empty())
	{
		setErrCode(SGL_INVALID_OPERATION);
		return -1;
	}
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
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	else {
		colorClearR = r;
		colorClearG = g;
		colorClearB = b;
	}
}

void sglClear(unsigned what) {
	if (contextWrapper.empty() || hasBegun) {
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}

	if ((what & SGL_COLOR_BUFFER_BIT) == SGL_COLOR_BUFFER_BIT) {
			contextWrapper[contextWrapper.activeContext]->clearColor(colorClearR, colorClearG, colorClearB );
	}
	else if ((what & SGL_DEPTH_BUFFER_BIT) == SGL_DEPTH_BUFFER_BIT) {
			contextWrapper[contextWrapper.activeContext]->clearDepth();
	}
	else {
			setErrCode(SGL_INVALID_VALUE);
	}
}

void sglBegin(sglEElementType mode) 
{
	if (hasBegun) 
	{ 
		setErrCode(SGL_INVALID_OPERATION);
		return; 
	}
	if (mode <= 0 || mode >= SGL_LAST_ELEMENT_TYPE) 
	{ 
		setErrCode(SGL_INVALID_ENUM); 
		return; 
	}
	hasBegun = true;
	drawingMethod = mode;
}

/*
Transforms input point - model, view, projection and viewport.
Input point is preserved and transformation is returned in its copy.
@param point input point to transform
@param output where transformed point is stored
*/
void transformThePoint(inputPoint4f* point, inputPoint4f& output)
{
	// test if between sglBegin and sglEnd
	if (hasBegun) {
		// for convenience, matrices are multiplied only once in sglEnd
		multiplyMatrixVector(multipliedMatrix, point, output);
	}
	else {
		// there is no sglBegin nor sglEnd, we can't be sure, where new 
		// transformation will appear so it needs to be done every time
		inputPoint4f tmp(*point);
		multiplyMatrixVector(modelViewStack.top(), point, output);
		multiplyMatrixVector(projectionStack.top(), &output, tmp);
		multiplyMatrixVector(viewportMatrix, &tmp, output);
	}
	// there should be perspective divide
}

// HERE UNUSED CODE STARTS
void transformScaleAndRotation(inputPoint4f* point, inputPoint4f& output) {
	float *mat = modelViewStack.top();
	output.x = mat[0] * point->x + mat[4] * point->y + mat[8] * point->z;
	output.y = mat[1] * point->x + mat[5] * point->y + mat[9] * point->z;
	output.z = mat[2] * point->x + mat[6] * point->y + mat[10] * point->z;
}

void invertTransformPoint(inputPoint4f* point, inputPoint4f& output) {
	if (!invertedForObject) {
		invertMatrix(viewportMatrix, inversedViewportMatrix);
		invertMatrix(projectionStack.top(), inversedProjectionMatrix);
		invertedForObject = true;
	}
	inputPoint4f tmp(*point);
	inputPoint4f tmp2(*point);
	multiplyMatrixVector(inversedViewportMatrix, point, tmp);
	multiplyMatrixVector(inversedProjectionMatrix, &tmp, tmp2);
	transformThePoint(&tmp2, output);
}

void transformPointModelView(inputPoint4f* point, inputPoint4f& output)
{
	multiplyMatrixVector(modelViewStack.top(), point, output);
	//inputPoint4f tmp (point);
	//multiplyMatrixVector(modelViewStack.top(), point, tmp);
	//multiplyMatrixVector(projectionStack.top(), tmp, output);
}

void drawPointNoTransform (inputPoint4f& point) {
	int W, H, x, y;

	SglContext *cont = contextWrapper.contexts[contextWrapper.activeContext];
	W = cont->getWidth();
	H = cont->getHeight();
	x = (int)round(point.x);
	y = (int)round(point.y);

	float *colorBuffer = cont->getColorBuffer();

	int offset;

	int size = (int)((pointSize - 1) / 2);
	int sizeCorrection = 1 - (int)pointSize % 2;

	for (int i = x - size; i <= x + (size + sizeCorrection); i++)
	{
		for (int j = y - size; j <= y + (size + sizeCorrection); j++)
		{
			if (i >= 0 && i < W && j >= 0 && j < H)
			{
				offset = j*W * 3 + i * 3;
				*(colorBuffer + offset) = point.r;
				*(colorBuffer + offset + 1) = point.g;
				*(colorBuffer + offset + 2) = point.b;
}
		}
	}
}

void drawPointRotatedAndScaled(inputPoint4f* point)
{
	inputPoint4f output;
	output.r = point->r;
	output.g = point->g;
	output.b = point->b;
	output.a = point->a;
	transformScaleAndRotation(point, output);

	int W, H, x, y;

	SglContext *cont = contextWrapper.contexts[contextWrapper.activeContext];
	W = cont->getWidth();
	H = cont->getHeight();
	x = (int)round(output.x);
	y = (int)round(output.y);

	float *colorBuffer = cont->getColorBuffer();

	int offset;

	int size = (int)((pointSize - 1) / 2);
	int sizeCorrection = 1 - (int)pointSize % 2;

	for (int i = x - size; i <= x + (size + sizeCorrection); i++)
	{
		for (int j = y - size; j <= y + (size + sizeCorrection); j++)
		{
			if (i >= 0 && i < W && j >= 0 && j < H)
			{
				offset = j*W * 3 + i * 3;
				*(colorBuffer + offset) = output.r;
				*(colorBuffer + offset + 1) = output.g;
				*(colorBuffer + offset + 2) = output.b;
			}
		}
	}

}
// HERE UNUSED CODE ENDS

/**
Method to draw point of specific size. Point is centered to middle (for odd size) with added right and bottom line (for even size).
*/
void drawMeAPoint(inputPoint4f* point) 
{
	//transform point
	inputPoint4f output;
	output.r = point->r;
	output.g = point->g;
	output.b = point->b;
	output.a = point->a;
	transformThePoint(point, output);

	//get color buffer and its properties
	int W, H, x, y;

	SglContext *cont = contextWrapper.contexts[contextWrapper.activeContext];
	W = cont->getWidth();
	H = cont->getHeight();
	x = (int)round(output.x);
	y = (int)round(output.y);

	float *colorBuffer = cont->getColorBuffer();

	int offset;

	//count size of pixel and right-bottom correction for even size.
	int size = (int)((pointSize-1) / 2);
	int sizeCorrection = 1 - (int) pointSize % 2;


	//set pixels in square
	for (int i = x - size; i <= x + (size + sizeCorrection); i++)
	{
		for (int j = y - size; j <= y + (size + sizeCorrection); j++)
		{
			if (i >= 0 && i < W && j >= 0 && j < H)
			{
				offset = j*W * 3 + i * 3;
				*(colorBuffer + offset) = output.r;
				*(colorBuffer + offset + 1) = output.g;
				*(colorBuffer + offset + 2) = output.b;
			}
		}
	}

}
/**
Naive method that was used to verify BresenhamLine work. It just samples line on X or Y axis and computes other coordinate.
*/
void drawMeALineNaive(inputPoint4f* start, inputPoint4f* end)
{
	int W, H, tempX, tempY, offset;

	SglContext *cont = contextWrapper.contexts[contextWrapper.activeContext];
	W = cont->getWidth();
	H = cont->getHeight();
	float *colorBuffer = cont->getColorBuffer();

	inputPoint4f startT;
	startT.r = start->r;
	startT.g = start->g;
	startT.b = start->b;
	transformThePoint(start, startT);
	inputPoint4f endT;
	endT.r = end->r;
	endT.g = end->g;
	endT.b = end->b;
	transformThePoint(end, endT);

	int x0, x1, y0, y1;
	float k;

	x0 = startT.x;
	y0 = startT.y;
	x1 = endT.x;
	y1 = endT.y;


	int endValue;
	float lerpValue;
	if (abs((float)(y1 - y0) / (x1 - x0)) < 1){
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
/**
Method drawing line using bresenham algoritm.
*/
void drawMeALineBresenham(inputPoint4f* start, inputPoint4f* end)

{
	int W, H, offset;
	int x0, x1, y0, y1, swap;
	//get collor buffer and its properties
	SglContext *cont = contextWrapper.contexts[contextWrapper.activeContext];
	W = cont->getWidth();
	H = cont->getHeight();
	float *colorBuffer = cont->getColorBuffer();

	//papply transformations to both ends
	inputPoint4f startT;
	startT.r = start->r;
	startT.g = start->g;
	startT.b = start->b;
	transformThePoint(start, startT);
	inputPoint4f endT;
	endT.r = end->r;
	endT.g = end->g;
	endT.b = end->b;
	transformThePoint(end, endT);

	//setup variables to work with
	x0 = (startT).x;
	y0 = (startT).y;
	x1 = (endT).x;
	y1 = (endT).y;

	int dx, dy;

	//decide which axis should be used for sampling (and make x be the one that grows faster)
	bool drivingX = abs(x1 - x0) > abs(y1 - y0);
	if (!drivingX)
	{
		swap = x0;
		x0 = y0;
		y0 = swap;

		swap = x1;
		x1 = y1;
		y1 = swap;
	}

	//make point 0 be the most left one
	if (x0 > x1)
	{
		swap = x0;
		x0 = x1;
		x1 = swap;

		swap = y0;
		y0 = y1;
		y1 = swap;

		inputPoint4f swap4f = startT;
		startT = endT;
		endT = swap4f;
	}

	//compute deltas
	dx = abs(x1 - x0);
	dy = y1 - y0;
	//decide it it will be quadrant 1,3 (value 1) or 2,4 (value -1)
	int yAdd = (dy < 0) ? -1 : 1;
	dy = abs(dy);

	//compute initial value of parameter p and constants c0, c1 for further modification of p
	int p = 2 * dy - dx;
	int c0 = 2 * dy;
	int c1 = c0 - 2 * dx;
	int tempY = y0;
	//variable for linear interpolation of color of line
	float lerpValue;

	//initial point (no p modifications)
	if (drivingX)
	{
		if (x0 >= 0 && x0 < W && y0 >= 0 && y0 < H)
		{
			offset = y0*W * 3 + x0 * 3;

			*(colorBuffer + offset) = startT.r;
			*(colorBuffer + offset + 1) = startT.g;
			*(colorBuffer + offset + 2) = startT.b;
		}

		}
	else {
		if (y0 >= 0 && y0 < W && x0 >= 0 && x0 < H)
		{
			offset = x0*W * 3 + y0 * 3;

			*(colorBuffer + offset) = startT.r;
			*(colorBuffer + offset + 1) = startT.g;
			*(colorBuffer + offset + 2) = startT.b;
		}
	}


	//rest of points (p is always modified)
	for (int tempX = x0 + 1; tempX <= x1; tempX++)
	{
		lerpValue = (float)(tempX - x0) / (x1 - x0);
		if (p < 0)
		{
			p += c0;
		}else {
			p += c1;
			tempY += yAdd;
		}
		if(drivingX)
		{
			if (tempX >= 0 && tempX < W && tempY >= 0 && tempY < H)
			{
				offset = tempY*W * 3 + tempX * 3;

				*(colorBuffer + offset) = (1-lerpValue)*startT.r + (lerpValue)*endT.r;
				*(colorBuffer + offset + 1) = (1-lerpValue)*startT.g + (lerpValue)*endT.g;
				*(colorBuffer + offset + 2) = (1-lerpValue)*startT.b + (lerpValue)*endT.b;
			}

		}else{
			if (tempY >= 0 && tempY < W && tempX >= 0 && tempX < H)
			{
				offset = tempX*W * 3 + tempY * 3;

				*(colorBuffer + offset) = (1 - lerpValue)*startT.r + (lerpValue)*endT.r;
				*(colorBuffer + offset + 1) = (1 - lerpValue)*startT.g + (lerpValue)*endT.g;
				*(colorBuffer + offset + 2) = (1 - lerpValue)*startT.b + (lerpValue)*endT.b;
			}
		}
	}

}

void fillLine(int xStart, int xEnd, int row, float colRStart, float colGStart, float colBStart, float colREnd, float colGEnd, float colBEnd)
{

	int W, H;
	SglContext *cont = contextWrapper.contexts[contextWrapper.activeContext];
	float *colorBuffer = cont->getColorBuffer();
	W = cont->getWidth();
	H = cont->getHeight();

	int y = 3*row*W;

	float lerpValue = 0;
	float lerpAdd = 1.0f / (xEnd - xStart);
	int offset = y + 3*xStart;

	for (int x = xStart; x <= xEnd; x++)
	{
		if (x >= 0 && x < W && row >= 0 && row < H)
		{
			*(colorBuffer + offset) = (1-lerpValue)*colRStart + (lerpValue)*colREnd;
			*(colorBuffer + offset + 1) = (1 - lerpValue)*colGStart + (lerpValue)*colGEnd;
			*(colorBuffer + offset + 2) = (1 - lerpValue)*colBStart + (lerpValue)*colBEnd;
		}
		offset += 3;
		lerpValue += lerpAdd;
	}
}

struct polyEdge
{
	int Y_upper, Y_lower, X_cross;
	float X_upper, X_step;
	polyEdge *next;
};

void drawMeAPolygon()
{
	switch (areaMode)
	{
	case SGL_POINT:
		drawPoints();
			break;
	case SGL_LINE:
		drawLineLoop();
		break;
	case SGL_FILL:
		polyEdge *root = new polyEdge;
		polyEdge *end = root;
		polyEdge *boundUpper;
		polyEdge *boundLower;

		inputPoint4f origin;
		inputPoint4f *tempPoint1;
		inputPoint4f *tempPoint2;

		inputPoint4f pointStart;
		inputPoint4f pointEnd;

		//store first point
		tempPoint2 = queue4f.front();
		origin = *tempPoint2;
		queue4f.pop();
		int counter = 1;
		float step;

		//build polyEdges
		while (!queue4f.empty())
		{
			tempPoint1 = tempPoint2;
			tempPoint2 = queue4f.front();
			queue4f.pop();
			counter++;

			if (tempPoint1->y > tempPoint2->y)
			{
				end->next = new polyEdge;
				end = end->next;

				end->Y_upper = tempPoint1->y;
				end->Y_lower = tempPoint2->y - 1;

				step = (tempPoint2->x - tempPoint1->x) / (tempPoint2->y - tempPoint1->y);
				end->X_step = step;
				end->X_upper = tempPoint1->x - ((tempPoint1->x) + (int)(tempPoint1->x))*step;
				end->X_cross = end->X_cross;

			}else if (tempPoint1->y > tempPoint2->y) {
				end->next = new polyEdge;
				end = end->next;

				end->Y_upper = tempPoint2->y;
				end->Y_lower = tempPoint1->y - 1;

				step = (tempPoint1->x - tempPoint2->x) / (tempPoint1->y - tempPoint2->y);
				end->X_step = step;
				end->X_upper = tempPoint2->x - ((tempPoint2->x) + (int)(tempPoint2->x))*step;
				end->X_cross = end->X_cross;

			}else {
				continue;
			}

			delete tempPoint1;
		}

		//draw line from last point to first point
		if (counter > 1)
		{
			drawMeALine(tempPoint2, &origin);
		}
		delete tempPoint2;


		break;
	}
	printf("drawMeAPolygon dont draw now \n");
}

void drawMeATriangleLineLoop(inputPoint4f* v1, inputPoint4f* v2, inputPoint4f* v3)
{
	printf("drawTriangleLineLoop not implemented yet! \n Transformation of points is still missing\n");

	drawMeALine(v1, v2);
	drawMeALine(v2, v3);
	drawMeALine(v3, v1);
}

void drawMeATriangle(inputPoint4f* v1, inputPoint4f* v2, inputPoint4f* v3)
{
	printf("trianglePrinting \n");
	inputPoint4f *p1, *p2, *p3;

	transformThePoint(v1, *v1);
	transformThePoint(v2, *v2);
	transformThePoint(v3, *v3);
	
	
	//order points
	if (v1->y >= v2->y)
	{
		// v1 >= v2
		if (v1->y >= v3->y)
		{
			// v1 >= v2
			// v1 >= v3
			if(v2->y >= v3->y)
			{
				// v1 >= v2
				// v1 >= v3
				// v2 >= v3
				p1 = v1;
				p2 = v2;
				p3 = v3;
			}else {
				// v1 >= v2
				// v1 >= v3
				// v3 > v2
				p1 = v1;
				p2 = v3;
				p3 = v2;
			}
		}else {
			// v1 >= v2
			// v3 > v1
			p1 = v3;
			p2 = v1;
			p3 = v2;
		}
	}else {
		// v2 > v1
		if (v2->y >= v3->y)
		{
			// v2 > v1
			// v2 >= v3
			if (v1->y >= v3->y)
			{
				// v2 > v1
				// v2 >= v3
				// v1 >= v3
				p1 = v2;
				p2 = v1;
				p3 = v3;
			}else {
				// c2 > v1
				// v2 >= v3
				// v3 > v1
				p1 = v2;
				p2 = v3;
				p3 = v1;
			}
		}else {
			// v2 > v1
			// v3 > v2
			p1 = v3;
			p2 = v2;
			p3 = v1;
		}
	}
	//points ordered

	inputPoint4f splittingPoint;
	float l;
	l = (p2->y - p1->y) / (p3->y - p1->y);

	splittingPoint.x = (1 - l)*p1->x + l*p3->x;
	splittingPoint.y = (1 - l)*p1->y + l*p3->y;
	splittingPoint.z = (1 - l)*p1->z + l*p3->z;
	splittingPoint.w = (1 - l)*p1->w + l*p3->w;

	splittingPoint.r = (1 - l)*p1->r + l*p3->r;
	splittingPoint.g = (1 - l)*p1->g + l*p3->g;
	splittingPoint.b = (1 - l)*p1->b + l*p3->b;

	//sglPointSize(5);
	//drawPointNoTransform(splittingPoint);

	int xInitLeft, xInitRight, yStart, yEnd, yCurrent, stepCount;
	float stepXLeft, stepXRight;

	float rL, rR, gL, gR, bL, bR;
	float rLd, rRd, gLd, gRd, bLd, bRd;

	if (p2->x <= splittingPoint.x)
	{
		//p2 is lefter
		{
			//upper triangle
			yStart = p1->y;
			yCurrent = yStart;
			yEnd = p2->y;
			stepCount = yStart - yEnd;
			stepXLeft = (p2->x - p1->x) / stepCount;
			stepXRight = (splittingPoint.x - p1->x) / stepCount;
			xInitLeft = xInitRight = p1->x;


			//prepare startingColor left and right
			rL = rR = p1->r;
			gL = gR = p1->g;
			bL = bR = p1->b;
			//prepare lerp increment
			//lerp left
			rLd = (p2->r - rL) / stepCount;
			gLd = (p2->g - gL) / stepCount;
			bLd = (p2->b - bL) / stepCount;
			//lerp right
			rRd = (splittingPoint.r - rR) / stepCount;
			gRd = (splittingPoint.g - gR) / stepCount;
			bRd = (splittingPoint.b - bR) / stepCount;

			
			for (int i = 0; i <= stepCount; i++)
			{
				//setPixel((int)(xInitLeft + i * stepXLeft), yCurrent, rL, gL, bL);
				//setPixel((int)(xInitRight + i * stepXRight), yCurrent, rR, gR, bR);
				fillLine((int)(xInitLeft + i * stepXLeft), (int)(xInitRight + i * stepXRight), yCurrent, rL, gL, bL, rR, gR, bR);
				//update left color
				rL += rLd;
				gL += gLd;
				bL += bLd;
				//update right color
				rR += rRd;
				gR += gRd;
				bR += bRd;
				//update row
				yCurrent--;
			}
		}

		{
			//lower triangle
			yStart = p3->y;
			yCurrent = yStart;
			yEnd = p2->y;
			stepCount = yEnd - yStart;
			stepXLeft = (p2->x - p3->x) / stepCount;
			stepXRight = (splittingPoint.x - p3->x) / stepCount;
			xInitLeft = xInitRight = p3->x;


			//prepare startingColor left and right
			rL = rR = p3->r;
			gL = gR = p3->g;
			bL = bR = p3->b;
			//prepare lerp increment
			//lerp left
			rLd = (p2->r - rL) / stepCount;
			gLd = (p2->g - gL) / stepCount;
			bLd = (p2->b - bL) / stepCount;
			//lerp right
			rRd = (splittingPoint.r - rR) / stepCount;
			gRd = (splittingPoint.g - gR) / stepCount;
			bRd = (splittingPoint.b - bR) / stepCount;

			for (int i = 0; i <= stepCount; i++)
			{
				fillLine((int)(xInitLeft + i * stepXLeft), (int)(xInitRight + i * stepXRight), yCurrent, rL, gL, bL, rR, gR, bR);
				//update left color
				rL += rLd;
				gL += gLd;
				bL += bLd;
				//update right color
				rR += rRd;
				gR += gRd;
				bR += bRd;
				//update row
				yCurrent++;
			}
		}
	}else {
		//splitting point is lefter
		{
			//upper triangle
			yStart = p1->y;
			yCurrent = yStart;
			yEnd = splittingPoint.y;
			stepCount = yStart - yEnd;
			stepXLeft = (splittingPoint.x - p1->x) / stepCount;
			stepXRight = (p2->x - p1->x) / stepCount;
			xInitLeft = xInitRight = p1->x;


			//prepare startingColor left and right
			rL = rR = p1->r;
			gL = gR = p1->g;
			bL = bR = p1->b;
			//prepare lerp increment
			//lerp left
			rLd = (splittingPoint.r - rL) / stepCount;
			gLd = (splittingPoint.g - gL) / stepCount;
			bLd = (splittingPoint.b - bL) / stepCount;
			//lerp right
			rRd = (p2->r - rR) / stepCount;
			gRd = (p2->g - gR) / stepCount;
			bRd = (p2->b - bR) / stepCount;


			for (int i = 0; i <= stepCount; i++)
			{
				setPixel((int)(xInitLeft + i * stepXLeft), yCurrent, rL, gL, bL);
				setPixel((int)(xInitRight + i * stepXRight), yCurrent, rR, gR, bR);
				fillLine((int)(xInitLeft + i * stepXLeft), (int)(xInitRight + i * stepXRight), yCurrent, rL, gL, bL, rR, gR, bR);
				//update left color
				rL += rLd;
				gL += gLd;
				bL += bLd;
				//update right color
				rR += rRd;
				gR += gRd;
				bR += bRd;
				//update row
				yCurrent--;
			}
		}

		{
			//lower triangle
			yStart = p3->y;
			yCurrent = yStart;
			yEnd = splittingPoint.y;
			stepCount = yEnd - yStart;
			stepXLeft = (splittingPoint.x - p3->x) / stepCount;
			stepXRight = (p2->x - p3->x) / stepCount;
			xInitLeft = xInitRight = p3->x;


			//prepare startingColor left and right
			rL = rR = p3->r;
			gL = gR = p3->g;
			bL = bR = p3->b;
			//prepare lerp increment
			//lerp left
			rLd = (splittingPoint.r - rL) / stepCount;
			gLd = (splittingPoint.g - gL) / stepCount;
			bLd = (splittingPoint.b - bL) / stepCount;
			//lerp right
			rRd = (p2->r - rR) / stepCount;
			gRd = (p2->g - gR) / stepCount;
			bRd = (p2->b - bR) / stepCount;

			for (int i = 0; i <= stepCount; i++)
			{
				fillLine((int)(xInitLeft + i * stepXLeft), (int)(xInitRight + i * stepXRight), yCurrent, rL, gL, bL, rR, gR, bR);
				//update left color
				rL += rLd;
				gL += gLd;
				bL += bLd;
				//update right color
				rR += rRd;
				gR += gRd;
				bR += bRd;
				//update row
				yCurrent++;
			}
		}

	}

	printf("drawMeATriangle not implemented yet! \n Transformation of points is still missing\n");
}

void drawTriangles()
{

	inputPoint4f *tempPoint1;
	inputPoint4f *tempPoint2;
	inputPoint4f *tempPoint3;
	switch (areaMode)
	{
	case SGL_POINT:
		drawPoints();
		break;
	case SGL_LINE:

		while (!queue4f.empty())
		{
			tempPoint1 = queue4f.front();
			queue4f.pop();

			if (queue4f.empty())
			{
				delete tempPoint1;
				break;
			}

			tempPoint2 = queue4f.front();
			queue4f.pop();

			if (queue4f.empty())
			{
				delete tempPoint1;
				delete tempPoint2;
				break;
			}

			tempPoint3 = queue4f.front();
			queue4f.pop();

			drawMeATriangleLineLoop(tempPoint1, tempPoint2, tempPoint3);

			delete tempPoint1;
			delete tempPoint2;
			delete tempPoint3;
		}
		break;
	case SGL_FILL:

		while (!queue4f.empty())
		{
			tempPoint1 = queue4f.front();
			queue4f.pop();

			if (queue4f.empty())
			{
				delete tempPoint1;
				break;
			}

			tempPoint2 = queue4f.front();
			queue4f.pop();

			if (queue4f.empty())
			{
				delete tempPoint1;
				delete tempPoint2;
				break;
			}

			tempPoint3 = queue4f.front();
			queue4f.pop();

			//printf("sending %f %f %f \n", tempPoint1->x, tempPoint1->y, tempPoint1->z);
			drawMeATriangle(tempPoint1, tempPoint2, tempPoint3);

			delete tempPoint1;
			delete tempPoint2;
			delete tempPoint3;
		}
		break;
	}
}


/**
Method to choose which algoritm will be used to draw line
*/
void drawMeALine(inputPoint4f* start, inputPoint4f* end)
{
	//choose which algorithm to use
	#ifdef LINE_NAIVE
		drawMeALineNaive(start, end);
	#else
		drawMeALineBresenham(start, end);
	#endif
}

/**
Method draw all points stored in vertex buffer.
*/
void drawPoints() 
{
	inputPoint4f *tempPoint;

	while (!queue4f.empty())
	{
		tempPoint = queue4f.front();
		drawMeAPoint(tempPoint);
		delete tempPoint;
		queue4f.pop();
	}
}

/**
draw lines
*/
void drawLines()
{
	inputPoint4f *tempPoint1;
	inputPoint4f *tempPoint2;

	//pick 2 following points and connect them with line
	while (!queue4f.empty())
	{
		tempPoint1 = queue4f.front();
		queue4f.pop();

		if (queue4f.empty())
		{
			delete tempPoint1;
			break;
		}

		tempPoint2 = queue4f.front();
		queue4f.pop();

		drawMeALine(tempPoint1, tempPoint2);
		delete tempPoint1;
		delete tempPoint2;
	}
}
/**
draw a strip of lines
*/
void drawLineStrip()
{
	//it is completely same as line loop but without line from end to start.
	inputPoint4f *tempPoint1;
	inputPoint4f *tempPoint2;

	if (queue4f.empty()) { return; }
	tempPoint2 = queue4f.front();
	queue4f.pop();

	while (!queue4f.empty())
	{
		tempPoint1 = tempPoint2;
		tempPoint2 = queue4f.front();
		queue4f.pop();

		drawMeALine(tempPoint1, tempPoint2);
		delete tempPoint1;
	}
	delete tempPoint2;
}
/**
draw a loop of lines
*/
void drawLineLoop()
{
	inputPoint4f origin;
	inputPoint4f *tempPoint1;
	inputPoint4f *tempPoint2;

	//store first point
	if (queue4f.empty()) { return; }
	tempPoint2 = queue4f.front();
	origin = *tempPoint2;
	queue4f.pop();
	int counter = 1;

	//draw lines form first point to last point
	while (!queue4f.empty())
	{
		tempPoint1 = tempPoint2;
		tempPoint2 = queue4f.front();
		queue4f.pop();
		counter++;

		drawMeALine(tempPoint1, tempPoint2);
		delete tempPoint1;
	}

	//draw line from last point to first point
	if (counter > 1)
	{
		drawMeALine(tempPoint2, &origin);
	}
	delete tempPoint2;
}

void sglEnd(void) 
{
	if (!hasBegun) { setErrCode(SGL_INVALID_OPERATION); return; }

	// so there is no need to compute multiplication every time, it is computed once here
	copyMatrix(multipliedMatrix, identityMatrix);
	// creates viewport (and there should be perspective divide)
	multipliedMatrix[0] = (viewportWidth - viewportOffsetX) / 2.0f;
	multipliedMatrix[5] = (viewportHeight - viewportOffsetY) / 2.0f;
	multipliedMatrix[10] = 0.5f;
	multipliedMatrix[12] = (viewportWidth / 2.0f) + viewportOffsetX;
	multipliedMatrix[13] = (viewportHeight / 2.0f) + viewportOffsetY;
	multipliedMatrix[14] = 0.5f;
	copyMatrix(viewportMatrix, multipliedMatrix);
	multiplyMatrix(multipliedMatrix, projectionStack.top());
	multiplyMatrix(multipliedMatrix, modelViewStack.top());

	invertMatrix(viewportMatrix, inversedViewportMatrix);

	switch (drawingMethod)
	{
	case SGL_POINTS:
		drawPoints();
		break;
	case SGL_LINES:
		drawLines();
		break;
	case SGL_LINE_STRIP:
		drawLineStrip();
		break;
	case SGL_LINE_LOOP:
		drawLineLoop();
		break;
	case SGL_TRIANGLES:
		drawTriangles();
		break;
	case SGL_POLYGON:
		drawMeAPolygon();
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

	queue4f.push(point);
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

	queue4f.push(point);
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

	queue4f.push(point);
}

/**
Used in Bresenham's algorithm for drawing circle. One point is copied on 8 
different places.
@param x untranslated position of point
@param y untranslated position of point
@param xs center of circle
@param ys center of circle
@param point input point with stored colors
*/
void setSymPoints(int x, int y, int xs, int ys, inputPoint4f& point) {
	point.x = x + xs;
	point.y = y + ys;
	//drawPointNoTransform(point);
	setPixel(point.x,point.y, point.r, point.g, point.b);

	point.x = xs - x;
	//drawPointNoTransform(point);
	setPixel(point.x, point.y, point.r, point.g, point.b);

	point.y = ys - y;
	//drawPointNoTransform(point);
	setPixel(point.x, point.y, point.r, point.g, point.b);

	point.x = x + xs;
	//drawPointNoTransform(point);
	setPixel(point.x, point.y, point.r, point.g, point.b);
	
	point.x = y + xs;
	point.y = x + ys;
	//drawPointNoTransform(point);
	setPixel(point.x, point.y, point.r, point.g, point.b);

	point.x = xs - y;
	//drawPointNoTransform(point);
	setPixel(point.x, point.y, point.r, point.g, point.b);

	point.y = ys - x;
	//drawPointNoTransform(point);
	setPixel(point.x, point.y, point.r, point.g, point.b);

	point.x = y + xs;
	//drawPointNoTransform(point);
	setPixel(point.x, point.y, point.r, point.g, point.b);
}

// HERE UNUSED CODE STARTS
void setSymPointsModified(int x, int y, int xs, int ys, inputPoint4f& point) {
	point.x = x + xs;
	point.y = y + ys;
	drawPointNoTransform(point);

	point.x = xs - x;
	drawPointNoTransform(point);

	point.y = ys - y;
	drawPointNoTransform(point);

	point.x = x + xs;
	drawPointNoTransform(point);
	/*
	point.x = y + xs;
	point.y = x + ys;
	drawPointNoTransform(point);

	point.x = xs - y;
	drawPointNoTransform(point);

	point.y = ys - x;
	drawPointNoTransform(point);

	point.x = y + xs;
	drawPointNoTransform(point);
	*/
}
// HERE UNUSED CODE ENDS

/**
Method to set single pixel in color buffer to specific color. (always size 1*1)
*/
void setPixel(float x0, float y0, float r, float g, float b)
{
	int W, H, x, y;

	SglContext *cont = contextWrapper.contexts[contextWrapper.activeContext];
	W = cont->getWidth();
	H = cont->getHeight();
	x = round(x0);
	y = round(y0);

	float *colorBuffer = cont->getColorBuffer();
	//printf("drawing: %d %d\n", x, y);

	if (x >= 0 && x < W && y >= 0 && y < H)
	{
		int offset = (y*W + x) * 3;
		*(colorBuffer + offset) = r;
		*(colorBuffer + offset + 1) = g;
		*(colorBuffer + offset + 2) = b;
	}
}

void sglCircle(float x, float y, float z, float radius) {
	if (hasBegun || contextWrapper.empty()) {
		setErrCode(SGL_INVALID_OPERATION); 
		return; 
	}

	if (radius < 0) {
		setErrCode(SGL_INVALID_VALUE);
		return;
	}

	float scaleFactor = sqrt(multipliedMatrix[0] * multipliedMatrix[5] - multipliedMatrix[1] * multipliedMatrix[4]);
	radius *= scaleFactor;

	inputPoint4f point;
	point.x = x;
	point.y = y;
	point.z = 0;
	point.w = 1;
	point.r = colorVertexR;
	point.g = colorVertexG;
	point.b = colorVertexB;
	point.a = 0;

	inputPoint4f output;
	transformThePoint(&point, output);
	x = output.x;
	y = output.y;

	switch (areaMode)
	{
	case SGL_POINT:
		sglBegin(SGL_POINTS);
		sglVertex3f(x, y, z);
		sglEnd();
		break;
	case SGL_LINE:
		// Bresenham's algorithm for drawing circle
		int xp, yp, p;
		xp = 0;
		yp = radius;
		p = 3 - 2 * radius;
		while (xp < yp) {
			setSymPoints(xp, yp, x, y, point);
			if (p < 0) {
				p = p + 4 * xp + 6;
			}
			else {
				p = p + 4 * (xp - yp) + 10;
				--yp;
			}
			++xp;
		}
		if (xp == yp)
			setSymPoints(xp, yp, x, y, point);

		break;
	case SGL_FILL:
		printf("No Circle filling algorithm implemented right now.\n");
		break;
	}
}

/**
Second algoritm to draw ellipses. 
*/
void sglEllipseSecond(float x, float y, float z, float a, float b) {
	if (contextWrapper.empty() || hasBegun) {
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	int psize = pointSize;
	pointSize = 1;

	if (a < 0 || b < 0) {
		setErrCode(SGL_INVALID_VALUE);
		return;
	}

	float scaleFactor = sqrt(multipliedMatrix[0] * multipliedMatrix[5] - multipliedMatrix[1] * multipliedMatrix[4]);
	a *= scaleFactor;
	b *= scaleFactor;

	inputPoint4f point;
	point.x = x;
	point.y = y;
	point.z = 0;
	point.w = 1;
	point.r = colorVertexR;
	point.g = colorVertexG;
	point.b = colorVertexB;
	point.a = 0;

	inputPoint4f output;
	transformThePoint(&point, output);
	x = output.x;
	y = output.y;
	//setPixel(x, y, 1.0f, 0.0f, 0.0f);

	float a2 = 2 * a * a;
	float b2 = 2 * b * b;
	float error = a*a*b;

	float tempX = 0;
	float tempY = b;

	float stopY = 0;
	float stopX = a2 * b;

	while (stopY <= stopX)
	{
		setPixel(x + tempX, y + tempY, colorVertexR, colorVertexG, colorVertexB);
		setPixel(x - tempX, y + tempY, colorVertexR, colorVertexG, colorVertexB);
		setPixel(x + tempX, y - tempY, colorVertexR, colorVertexG, colorVertexB);
		setPixel(x - tempX, y - tempY, colorVertexR, colorVertexG, colorVertexB);
		tempX++;
		error -= b2 * (tempX - 1);
		stopY += b2;
		if (error <= 0)
		{
			error += a2 * (tempY - 1);
			tempY--;
			stopX -= a2;
		}
	}

	error = b*b*a;
	tempX = a;
	tempY = 0;
	stopY = b2 * a;
	stopX = 0;

	while (stopY >= stopX)
	{
		setPixel(x + tempX, y + tempY, colorVertexR, colorVertexG, colorVertexB);
		setPixel(x - tempX, y + tempY, colorVertexR, colorVertexG, colorVertexB);
		setPixel(x + tempX, y - tempY, colorVertexR, colorVertexG, colorVertexB);
		setPixel(x - tempX, y - tempY, colorVertexR, colorVertexG, colorVertexB);
		tempY++;
		error -= a2 * (tempY - 1);
		stopX += a2;
		if (error < 0)
		{
			error += b2 * (tempX - 1);
			tempX--;
			stopY -= b2;
		}
	}
}

void sglEllipseSegmented(float x, float y, float z, float a, float b)
{
	if (contextWrapper.empty() || hasBegun) {
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	int psize = pointSize;
	pointSize = 1;

	if (a < 0 || b < 0) {
		setErrCode(SGL_INVALID_VALUE);
		return;
	}
	//decide how to draw it.

	switch (areaMode)
	{
	case SGL_POINT:
		sglBegin(SGL_POINTS); 
		sglVertex3f(x, y, z);
		sglEnd();
		break;
	case SGL_LINE:
		sglBegin(SGL_LINE_LOOP);
		break;
	case SGL_FILL:
		printf("check sgl.cpp sglEllipseSegmented fill branch.");
		sglBegin(SGL_POLYGON);
		break;
	}

	int segments = 40;
	float angle = 0.0f;
	float delta = 2.0f * 3.14159 / segments;
	for (int i = 0; i <= segments; i++)
	{
		sglVertex3f(x + a*cos(angle), y + b*sin(angle), z);
		angle += delta;
	}
	sglEnd();

}

/**
Bresenham's algorithm for drawing ellipses.
@param x untranslated position of point
@param y untranslated position of point
@param xs center of circle
@param ys center of circle
@param a horizontal squish
@param b vertical squish
*/
void sglEllipseFirst(float x, float y, float z, float a, float b) {
	invertedForObject = false;
	if (contextWrapper.empty() || hasBegun) {
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	int psize = pointSize;
	pointSize = 1;

	if (a < 0 || b < 0) {
		setErrCode(SGL_INVALID_VALUE);
		return;
	}

	float scaleFactor = sqrt(multipliedMatrix[0] * multipliedMatrix[5] - multipliedMatrix[1] * multipliedMatrix[4]);
	a *= scaleFactor;
	b *= scaleFactor;

	inputPoint4f point;
	point.x = x;
	point.y = y;
	point.z = 0;
	point.w = 1;
	point.r = colorVertexR;
	point.g = colorVertexG;
	point.b = colorVertexB;
	point.a = 0;

	inputPoint4f output;
	transformThePoint(&point, output);
	x = output.x;
	y = output.y;

	int a2 = a * a;
	int b2 = b * b;
	int fa2 = 4 * a2, fb2 = 4 * b2;
	float xp, yp, sigma;

	/* first half */
	for (xp = 0, yp = b, sigma = 2 * b2 + a2*(1 - 2 * b); b2*xp <= a2*yp; xp++)
	{
		point.x = xp + x;
		point.y = yp + y;
		drawPointRotatedAndScaled(&point);
		//drawPointNoTransform(point);

		point.y = y - yp;
		drawPointRotatedAndScaled(&point);
		//drawPointNoTransform(point);

		point.x = x - xp;
		drawPointRotatedAndScaled(&point);
		//drawPointNoTransform(point);

		point.y = y + yp;
		drawPointRotatedAndScaled(&point);
		//drawPointNoTransform(point);
		if (sigma >= 0)
		{
			sigma += fa2 * (1 - yp);
			yp--;
		}
		sigma += b2 * ((4 * xp) + 6);
	}

	/* second half */
	for (xp = a, yp = 0, sigma = 2 * a2 + b2*(1 - 2 * a); a2*yp <= b2*xp; yp++)
	{
		point.x = xp + x;
		point.y = yp + y;
		drawPointRotatedAndScaled(&point);
		//drawPointNoTransform(point);

		point.y = y - yp;
		drawPointRotatedAndScaled(&point);
		//drawPointNoTransform(point);

		point.x = x - xp;
		drawPointRotatedAndScaled(&point);
		//drawPointNoTransform(point);

		point.y = y + yp;
		drawPointRotatedAndScaled(&point);
		//drawPointNoTransform(point);
		if (sigma >= 0)
		{
			sigma += fb2 * (1 - xp);
			xp--;
		}
		sigma += a2 * ((4 * yp) + 6);
	}

	pointSize = psize;
}

void sglEllipse(float x, float y, float z, float a, float b) {
	#ifdef ELLIPSE
	sglEllipseSegmented(x, y, z, a, b);
	#elif ELLIPSE_SECOND
		sglEllipseSecond(x, y, z, a, b);
	#else
		sglEllipseFirst(x, y, z, a, b);
	#endif
	
}

/**
Bresenham's algorithm for drawing circle, though limited by minimal and maximal angle.
@param x untranslated position of point
@param y untranslated position of point
@param xs center of arc
@param ys center of arc
@param from Minimum required angle for arc.
@param to Maximum angle for arc.
*/
void setSymPointsLimit(int x, int y, int xs, int ys, inputPoint4f *point, float radius, float from, float to) {
	point->x = x + xs;
	point->y = y + ys;
	float angle = acos(x / radius);
	if (y < 0)
		angle = -angle + 2 * 3.14159274;
	if (angle >= from && angle <= to)
		setPixel(point->x, point->y, point->r, point->g, point->b);

	point->x = xs - x;
	angle = acos(-x / radius);
	if (y < 0)
		angle = -angle + 2 * 3.14159274;
	if (angle >= from && angle <= to)
		setPixel(point->x, point->y, point->r, point->g, point->b);

	point->y = ys - y;
	angle = acos(-x / radius);
	if (-y < 0)
		angle = -angle + 2 * 3.14159274;
	if (angle >= from && angle <= to)
		setPixel(point->x, point->y, point->r, point->g, point->b);

	point->x = x + xs;
	angle = acos(x / radius);
	if (-y < 0)
		angle = -angle + 2 * 3.14159274;
	if (angle >= from && angle <= to)
		setPixel(point->x, point->y, point->r, point->g, point->b);

	point->x = y + xs;
	point->y = x + ys;
	angle = acos(y / radius);
	if (x < 0)
		angle = -angle + 2 * 3.14159274;
	if (angle >= from && angle <= to)
		setPixel(point->x, point->y, point->r, point->g, point->b);

	point->x = xs - y;
	angle = acos(-y / radius);
	if (x < 0)
		angle = -angle + 2 * 3.14159274;
	if (angle >= from && angle <= to)
		setPixel(point->x, point->y, point->r, point->g, point->b);

	point->y = ys - x;
	angle = acos(-y / radius);
	if (-x < 0)
		angle = -angle + 2 * 3.14159274;
	if (angle >= from && angle <= to)
		setPixel(point->x, point->y, point->r, point->g, point->b);

	point->x = y + xs;
	angle = acos(y / radius);
	if (-x < 0)
		angle = -angle + 2 * 3.14159274;
	if (angle >= from && angle <= to)
		setPixel(point->x, point->y, point->r, point->g, point->b);
}

void sglArc(float x, float y, float z, float radius, float from, float to) {
	invertedForObject = false;
	if (hasBegun) {
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}

	if (radius < 0) {
		setErrCode(SGL_INVALID_VALUE);
		return;
	}


	int segments;
	float angle;
	float delta;

	switch (areaMode)
	{
	case SGL_POINT:
		sglBegin(SGL_POINTS);
		sglVertex3f(x, y, z);
		sglEnd();
		break;
	case SGL_LINE:
		sglBegin(SGL_LINE_STRIP);
		segments = 40;
		angle = from;
		delta = (to - from) / segments;
		for (int i = 0; i <= segments; i++)
		{
			sglVertex3f(x + radius*cos(angle), y + radius*sin(angle), z);
			angle += delta;
		}
		sglEnd();
		break;
	case SGL_FILL:
		printf("No Arc filling algorithm implemented right now. \n");
		break;
	}

	// per pixel arc, but didn't figure how to rotate
	/*
	while (from >= 2 * 3.14159274)
		from -= 2 * 3.14159274;

	while (to >= 2 * 3.14159274)
		to -= 2 * 3.14159274;

	if (from < 0.00001f)
		from = 0;

	if (to < 0.00001f)
		to = 0;

	float scaleFactor = sqrt(multipliedMatrix[0] * multipliedMatrix[5] - multipliedMatrix[1] * multipliedMatrix[4]);
	radius *= scaleFactor;
	
	inputPoint4f point;
	point.x = x;
	point.y = y;
	point.z = 0;
	point.w = 1;
	point.r = colorVertexR;
	point.g = colorVertexG;
	point.b = colorVertexB;
	point.a = 0;

	inputPoint4f output;
	transformThePoint(&point, output);
	x = output.x;
	y = output.y;

	float xp, yp, p;
	xp = 0;
	yp = radius;
	p = 3 - 2 * radius;
	//printf("drawing from %f to %f at [%f, %f] \n",from, to, x, y);
	while (xp < yp) {
		setSymPointsLimit(xp, yp, x, y, &point, radius, from, to);
		if (p < 0) {
			p = p + 4 * xp + 6;
		}
		else {
			p = p + 4 * (xp - yp) + 10;
			--yp;
		}
		++xp;
	}
	if (xp == yp)
		setSymPointsLimit(xp, yp, x, y, &point, radius, from, to);
		*/
}

//---------------------------------------------------------------------------
// Transform functions
//---------------------------------------------------------------------------

void sglMatrixMode( sglEMatrixMode mode ) {
	if (hasBegun || contextWrapper.empty()) {
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	if (mode != SGL_MODELVIEW && mode != SGL_PROJECTION) {
		setErrCode(SGL_INVALID_ENUM);
		return;
	}
	matrixMode = mode;
}

/**
Copies matrix and returns its deep copy.
@param matrix matrix 4x4 to be copied
@return copy of 4x4 matrix
*/
float* duplicateMatrix(const float* matrix) {
	float* ret = new float[16];
	for (int i = 0; i < 16; ++i) {
		ret[i] = matrix[i];
	}
	return ret;
}

void sglPushMatrix(void) {
	if (hasBegun || contextWrapper.empty()) {
		setErrCode(SGL_INVALID_OPERATION);
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
		setErrCode(SGL_INVALID_OPERATION);
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

void sglLoadIdentity(void) {
	if (hasBegun || contextWrapper.empty()) {
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}

	switch (matrixMode) {
	case SGL_MODELVIEW:
		copyMatrix(modelViewStack.top(), identityMatrix);
		break;
	case SGL_PROJECTION:
		copyMatrix(projectionStack.top(), identityMatrix);
		break;
	default:
		break;
	}
}

void sglLoadMatrix(const float *matrix) {
	if (hasBegun || contextWrapper.empty()) {
		setErrCode(SGL_INVALID_OPERATION);
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

void sglMultMatrix(const float *matrix) {
	if (hasBegun || contextWrapper.empty()) {
		setErrCode(SGL_INVALID_OPERATION);
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
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	float translate[16];
	copyMatrix(translate, identityMatrix);
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
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	float scale[16];
	copyMatrix(scale, identityMatrix);
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
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}

	float tmp[16];
	float rotate[16];
	copyMatrix(tmp, identityMatrix);
	copyMatrix(rotate, identityMatrix);
	tmp[12] = -centerx;
	tmp[13] = -centery;
	rotate[0] = cos(angle); rotate[4] = -sin(angle);
	rotate[1] = sin(angle); rotate[5] = cos(angle);

	multiplyMatrix(rotate, tmp);

	tmp[12] = centerx;
	tmp[13] = centery;

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
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}

	float rotate[16];
	copyMatrix(rotate, identityMatrix);
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
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}

	float ortho[16];
	copyMatrix(ortho, identityMatrix);
	ortho[0] = 2.0f / (right - left);
	ortho[5] = 2.0f / (top - bottom);
	ortho[10] = - 2.0f / (far - near);
	ortho[12] = - (right + left) / (right - left);
	ortho[13] = - (top + bottom) / (top - bottom);
	ortho[14] = -(far + near) / (far - near);
	/*for (int i = 0; i < 4; i++) {
		printf("\n");
		for (int j = i; j < 16; j += 4) {
			printf("%f ", ortho[j]);
		}
	}*/
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
		setErrCode(SGL_INVALID_VALUE);
		return;
	}
	else if (hasBegun || contextWrapper.empty()) {
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}

	float A = (right + left) / (right - left);
	float B = (top + bottom) / (top - bottom);
	float C = (far + near) / (far - near);
	float D = - 2.0f * far * near / (far - near);

	float persp[16];
	copyMatrix(persp, identityMatrix);
	persp[0] = 2.0f * near/ (right - left);
	persp[5] = 2.0f * near / (top - bottom);
	persp[8] = A;
	persp[9] = B;
	persp[10] = C;
	persp[11] = -1.0f;
	persp[14] = D;
	persp[15] = 0.0f;

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
		setErrCode(SGL_INVALID_VALUE);
	}
	else if (hasBegun || contextWrapper.empty()) {
		setErrCode(SGL_INVALID_OPERATION);
	}
	else {
		viewportOffsetX = x;
		viewportOffsetY = y;
		viewportWidth = width;
		viewportHeight = height;
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

void sglAreaMode(sglEAreaMode mode) 
{
	if (hasBegun || contextWrapper.empty())
	{
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	switch (mode)
	{
	case SGL_POINT:
		areaMode = SGL_POINT;
		break;
	case SGL_LINE:
		areaMode = SGL_LINE;
		break;
	case SGL_FILL:
		areaMode = SGL_FILL;
		break;
	default:
		setErrCode(SGL_INVALID_ENUM);
		break;
	}

}

void sglPointSize(float size) 
{
	if (size <= 0) 
	{
		setErrCode(SGL_INVALID_VALUE);
		return;
	}

	if (hasBegun || contextWrapper.empty() ) 
	{
		setErrCode(SGL_INVALID_OPERATION); 
		return;
	}
	pointSize = size;
}

void sglEnable(sglEEnableFlags cap) {
	if (hasBegun || contextWrapper.empty()) {
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	switch (cap)
	{
	case SGL_DEPTH_TEST:
		testDepth = true;
		break;
	default:
		setErrCode(SGL_INVALID_ENUM);
		break;
	}
}

void sglDisable(sglEEnableFlags cap) {
	if (hasBegun || contextWrapper.empty()) {
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	switch (cap)
	{
	case SGL_DEPTH_TEST:
		testDepth = false;
		break;
	default:
		setErrCode(SGL_INVALID_ENUM);
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
