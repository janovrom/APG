//---------------------------------------------------------------------------
// sgl.cpp
// Empty implementation of the SGL (Simple Graphics Library)
// Date:  2011/11/1
// Author: Jaroslav Krivanek, Jiri Bittner CTU Prague
//---------------------------------------------------------------------------

#include "sgl.h"
#include "sglcontext.h"
#include <limits>



//#define LINE_NAIVE
// decides which ellipse algoritm should be used
#define ELLIPSE

using namespace std;

inline void drawPixel(int offsetC, int offsetD, float z, float r, float g, float b, float *colorBuffer, float *depthBuffer);
inline void setPixel(float x0, float y0, float r, float g, float b, float z);
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

/**
Method for drawing into collor and depth buffers
*/
inline void drawPixel(int offsetC, int offsetD, float z, float r, float g, float b, float *colorBuffer, float *depthBuffer)
{
	if (testDepth == true)
	{
		if (*(depthBuffer + offsetD) > z)
		{
			*(colorBuffer + offsetC) = r;
			*(colorBuffer + offsetC + 1) = g;
			*(colorBuffer + offsetC + 2) = b;
			*(depthBuffer + offsetD) = z;
		}
	}
	else
	{
		*(colorBuffer + offsetC) = r;
		*(colorBuffer + offsetC + 1) = g;
		*(colorBuffer + offsetC + 2) = b;
		//printf("depth %d %d %f \n", offsetD, offsetC, z);
		*(depthBuffer + offsetD) = z;
	}
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

// DEBUG CODE STARTS (You shall not delete it again!!!) - I SHALL!
void drawPointNoTransform(inputPoint4f& point) {
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

// THIS DEBUG METHOD DOES NOT DO A DEPTH TEST
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
// DEBUG CODE ENDS

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
	bool noOperation = true;

	if ((what & SGL_COLOR_BUFFER_BIT) == SGL_COLOR_BUFFER_BIT) {
			contextWrapper[contextWrapper.activeContext]->clearColor(colorClearR, colorClearG, colorClearB );
			noOperation = false;
	}
	if ((what & SGL_DEPTH_BUFFER_BIT) == SGL_DEPTH_BUFFER_BIT) {
			contextWrapper[contextWrapper.activeContext]->clearDepth();
			noOperation = false;
	}
	if(noOperation) {
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
Input point is preserved and transformation is returned in output.
@param point input point to transform
@param output where transformed point is stored
*/
void transformThePointAndCopyColor(inputPoint4f* point, inputPoint4f& output)
{
	// test if between sglBegin and sglEnd
	if (hasBegun) {
		//printf("TRANSFORMING ----- has begun\n");
		// for convenience, matrices are multiplied only once in sglEnd
		multiplyMatrixVector(matrixMVP, point, output);
		// perspective divide
		output.x = output.x / output.w;
		output.y = output.y / output.w;
		output.z = output.z / output.w;
		output.w = 1.0f;
		// viewport multiplication
		output.x = output.x * viewportMatrix[0] + viewportMatrix[12];
		output.y = output.y * viewportMatrix[5] + viewportMatrix[13];
		output.z = output.z * viewportMatrix[10] + viewportMatrix[14];
	}
	else {
		//printf("TRANSFORMING ----- outside begin-end\n");
		// there is no sglBegin nor sglEnd, we can't be sure, where new 
		// transformation will appear so it needs to be done every time
		inputPoint4f tmp(*point);
		multiplyMatrixVector(modelViewStack.top(), point, output);
		multiplyMatrixVector(projectionStack.top(), &output, tmp);
		// perspective divide
		output.x = tmp.x / tmp.w;
		output.y = tmp.y / tmp.w;
		output.z = tmp.z / tmp.w;
		output.w = 1.0f;
		// viewport multiplication
		output.x = output.x * viewportMatrix[0] + viewportMatrix[12];
		output.y = output.y * viewportMatrix[5] + viewportMatrix[13];
		output.z = output.z * viewportMatrix[10] + viewportMatrix[14];
		//multiplyMatrixVector(viewportMatrix, &tmp, output);
	}
	output.r = point->r;
	output.g = point->g;
	output.b = point->b;
	output.a = point->a;
}

/**
Method to draw point of specific size. Point is centered to middle (for odd size) with added right and bottom line (for even size).
*/
void drawMeAPoint(inputPoint4f* point) 
{
	//transform point
	inputPoint4f output;
	/*output.r = point->r;
	output.g = point->g;
	output.b = point->b;
	output.a = point->a;*/
	transformThePointAndCopyColor(point, output);

	//get color buffer and its properties
	int W, H, x, y;

	SglContext *cont = contextWrapper.contexts[contextWrapper.activeContext];
	W = cont->getWidth();
	H = cont->getHeight();
	x = (int)round(output.x);
	y = (int)round(output.y);

	float *colorBuffer = cont->getColorBuffer();

	int offset, offsetD;

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
				offsetD = j*W + i;
				//printf("drawMeAPoint\n");
				drawPixel(offset, offsetD,
						  output.z,
						  output.r,
						  output.g,
						  output.b,
						  colorBuffer, cont->getDepthBuffer());
				//*(colorBuffer + offset) = output.r;
				//*(colorBuffer + offset + 1) = output.g;
				//*(colorBuffer + offset + 2) = output.b;
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
	transformThePointAndCopyColor(start, startT);
	inputPoint4f endT;
	endT.r = end->r;
	endT.g = end->g;
	endT.b = end->b;
	transformThePointAndCopyColor(end, endT);

	int x0, x1, y0, y1;
	float k;

	float z0, z1, dz;

	x0 = startT.x;
	y0 = startT.y;
	x1 = endT.x;
	y1 = endT.y;

	z0 = startT.z;
	z1 = endT.z;


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

			swap = z0;
			z0 = z1;
			z1 = swap;
		}
		dz = (z1 - z0) / (x1 - x0);
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
				int offsetD = tempY*W + tempX;
				//printf("line naive X\n");
				drawPixel(offset, offsetD,
						  z0,
						  (lerpValue)*startT.r + (1 - lerpValue)*endT.r,
						  (lerpValue)*startT.g + (1 - lerpValue)*endT.g,
						  (lerpValue)*startT.b + (1 - lerpValue)*endT.b,
						  colorBuffer, cont->getDepthBuffer());
				//*(colorBuffer + offset) = (lerpValue)*startT.r + (1-lerpValue)*endT.r;
				//*(colorBuffer + offset + 1) = (lerpValue)*startT.g + (1 - lerpValue)*endT.g;
				//*(colorBuffer + offset + 2) = (lerpValue)*startT.b + (1 - lerpValue)*endT.b;
			}
			z0 += dz;
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

			swap = z0;
			z0 = z1;
			z1 = swap;
		}
		dz = (z1 - z0) / (y1 - y0);
		k = (float)(x1 - x0) / (y1 - y0);
		endValue = y1 - y0;
		for (int i = 0; i <= endValue; i++) {
			//setPixel(x0 + i*k, y0 + i);
			lerpValue = (float)i / endValue;
			tempX = x0 + i*k;
			tempY = y0 + i;

			if (tempX >= 0 && tempX < W && tempY >= 0 && tempY < H)
			{
				//printf("line naive Y\n");
				offset = tempY*W * 3 + tempX * 3;
				int offsetD = tempY*W + tempX;
				drawPixel(offset, offsetD,
						  z0,
						  (lerpValue)*startT.r + (1 - lerpValue)*endT.r,
						  (lerpValue)*startT.g + (1 - lerpValue)*endT.g,
						  (lerpValue)*startT.b + (1 - lerpValue)*endT.b,
						  colorBuffer, cont->getDepthBuffer());
				//*(colorBuffer + offset) = (lerpValue)*startT.r + (1 - lerpValue)*endT.r;
				//*(colorBuffer + offset + 1) = (lerpValue)*startT.g + (1 - lerpValue)*endT.g;
				//*(colorBuffer + offset + 2) = (lerpValue)*startT.b + (1 - lerpValue)*endT.b;
			}
			z0 += dz;
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
	transformThePointAndCopyColor(start, startT);
	inputPoint4f endT;
	endT.r = end->r;
	endT.g = end->g;
	endT.b = end->b;
	transformThePointAndCopyColor(end, endT);

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
			int offsetD = y0*W + x0;
			//printf("line bress X init\n");
			drawPixel(offset, offsetD,
					  startT.z,
					  startT.r,
					  startT.g,
					  startT.b,
					  colorBuffer, cont->getDepthBuffer());
			//*(colorBuffer + offset) = startT.r;
			//*(colorBuffer + offset + 1) = startT.g;
			//*(colorBuffer + offset + 2) = startT.b;
		}

		}
	else {
		if (y0 >= 0 && y0 < W && x0 >= 0 && x0 < H)
		{
			offset = x0*W * 3 + y0 * 3;
			int offsetD = x0*W + y0;
			//printf("line bress Y init\n");
			drawPixel(offset, offsetD,
				      startT.z,
					  startT.r,
					  startT.g,
					  startT.b,
					  colorBuffer, cont->getDepthBuffer());
			//*(colorBuffer + offset) = startT.r;
			//*(colorBuffer + offset + 1) = startT.g;
			//*(colorBuffer + offset + 2) = startT.b;
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
				int offsetD = tempY*W + tempX;
				//printf("line bress X\n");
				drawPixel(offset, offsetD,
					      (1 - lerpValue)*startT.z + (lerpValue)*endT.z,
						  (1 - lerpValue)*startT.r + (lerpValue)*endT.r,
						  (1 - lerpValue)*startT.g + (lerpValue)*endT.g,
						  (1 - lerpValue)*startT.b + (lerpValue)*endT.b,
						  colorBuffer, cont->getDepthBuffer());
				//*(colorBuffer + offset) = (1-lerpValue)*startT.r + (lerpValue)*endT.r;
				//*(colorBuffer + offset + 1) = (1-lerpValue)*startT.g + (lerpValue)*endT.g;
				//*(colorBuffer + offset + 2) = (1-lerpValue)*startT.b + (lerpValue)*endT.b;
			}

		}else{
			if (tempY >= 0 && tempY < W && tempX >= 0 && tempX < H)
			{
				offset = tempX*W * 3 + tempY * 3;
				int offsetD = tempX*W + tempY;
				//printf("line bress Y\n");
				drawPixel(offset, offsetD,
					      (1 - lerpValue)*startT.z + (lerpValue)*endT.z,
						  (1 - lerpValue)*startT.r + (lerpValue)*endT.r,
						  (1 - lerpValue)*startT.g + (lerpValue)*endT.g,
						  (1 - lerpValue)*startT.b + (lerpValue)*endT.b,
						  colorBuffer, cont->getDepthBuffer());

				//*(colorBuffer + offset) = (1 - lerpValue)*startT.r + (lerpValue)*endT.r;
				//*(colorBuffer + offset + 1) = (1 - lerpValue)*startT.g + (lerpValue)*endT.g;
				//*(colorBuffer + offset + 2) = (1 - lerpValue)*startT.b + (lerpValue)*endT.b;
			}
		}
	}

}

/**
Method that draws pixels from start point to end point. It also interpolate colors.
*/
void fillLine(int xStart, int xEnd, float zStart, float zEnd, int row, float colRStart, float colGStart, float colBStart, float colREnd, float colGEnd, float colBEnd)
{

	int W, H;
	SglContext *cont = contextWrapper.contexts[contextWrapper.activeContext];
	float *colorBuffer = cont->getColorBuffer();
	W = cont->getWidth();
	H = cont->getHeight();

	int y = 3*row*W;

	float lerpValue = 0;
	float lerpAdd = 1.0f / (xEnd - xStart);
	int offsetC = y + 3*xStart;
	int offsetD = (row*W + xStart);

	for (int x = xStart; x <= xEnd; x++)
	{
		if (x >= 0 && x < W && row >= 0 && row < H)
		{
			//printf("fill\n");
			drawPixel(offsetC, offsetD, 
					  (1 - lerpValue)*zStart + (lerpValue)*zEnd,
					  (1 - lerpValue)*colRStart + (lerpValue)*colREnd, 
					  (1 - lerpValue)*colGStart + (lerpValue)*colGEnd, 
					  (1 - lerpValue)*colBStart + (lerpValue)*colBEnd, 
					  colorBuffer, cont->getDepthBuffer());
			//*(colorBuffer + offset) = (1-lerpValue)*colRStart + (lerpValue)*colREnd;
			//*(colorBuffer + offset + 1) = (1 - lerpValue)*colGStart + (lerpValue)*colGEnd;
			//*(colorBuffer + offset + 2) = (1 - lerpValue)*colBStart + (lerpValue)*colBEnd;
		}
		offsetC += 3;
		offsetD += 1;
		lerpValue += lerpAdd;
	}
}
/**
Structure representing one line for polygon filling algorithm.
*/
struct polyEdge
{
	int Y_upper, Y_lower, X_cross;
	float X_upper, X_step;
	float Z_upper, Z_step;
	polyEdge *next;

	float R, G, B;
	float RD, GD, BD;
};

/**
Method to set up edge for polygon filling.
*/
void setPolyEdge(polyEdge *end, inputPoint4f *high, inputPoint4f *low)
{
	int stepCount;
	float step;

	end->Y_upper = high->y;
	end->Y_lower = low->y + 1;
	stepCount = end->Y_upper - end->Y_lower + 1;

	step = (low->x - high->x) / (high->y - low->y);
	end->X_step = step;
	//printf("correction %f \n", ((float)(high->y) - (int)(high->y))*step);
	end->X_upper = high->x;// +((float)(high->y) - (int)(high->y))*step;
	end->X_cross = end->X_upper;

	end->R = high->r;
	end->G = high->g;
	end->B = high->b;

	end->RD = (low->r - high->r) / stepCount;
	end->GD = (low->g - high->g) / stepCount;
	end->BD = (low->b - high->b) / stepCount;

	end->Z_upper = high->z;
	end->Z_step= (low->z - high->z) / stepCount;
	/*
	printf("\npolyEdge preparing\n");
	printf("Yup %d Ylow %d \n", end->Y_upper, end->Y_lower);
	printf("Xcro %d \n", end->X_cross);
	printf("Xtrue %f Xstep %f \n", end->X_upper, end->X_step);
	printf("polyEdge prepared\n\n");
	*/
}
/**
Order list of edges by Y_upper
*/
void listOrderByY_Upper(polyEdge *root, polyEdge *end)
{
	polyEdge *currentPred;
	polyEdge *tempEdge1;
	polyEdge *tempEdge2;
	int changes = 1;
	while (changes != 0)
	{
		changes = 0;
		currentPred = root;
		if (currentPred->next == end) { return; }
		while (true)
		{
			if (currentPred->next->next == end){break;}
			if (currentPred->next->Y_upper < currentPred->next->next->Y_upper)
			{
				tempEdge1 = currentPred->next;
				tempEdge2 = tempEdge1->next;
				currentPred->next = tempEdge2;
				tempEdge1->next = tempEdge2->next;
				tempEdge2->next = tempEdge1;
				changes++;
			}
			currentPred = currentPred->next;
		}
	}
}
/**
Print all exdes from list
*/
/*
void printList(polyEdge *root, polyEdge *end)
{
	polyEdge *current;
	current = root->next;
	while (current != end)
	{
		printf("-----Yu %d Yl %d XC %d Xt %f S %f\n",current->Y_upper, current->Y_lower, current->X_cross, current->X_upper, current->X_step);
		current = current->next;
	}
}
*/
/**
Order list of edges by X_Cross 
*/
void listOrderByX_Cross(polyEdge *root, polyEdge *end)
{
	polyEdge *currentPred;
	polyEdge *tempEdge1;
	polyEdge *tempEdge2;
	int changes = 1;
	while (changes != 0)
	{
		changes = 0;
		currentPred = root;
		if (currentPred->next == end) { return; }
		while (true)
		{
			if (currentPred->next->next == end) { break; }
			if (currentPred->next->X_cross > currentPred->next->next->X_cross)
			{
				tempEdge1 = currentPred->next;
				tempEdge2 = tempEdge1->next;
				currentPred->next = tempEdge2;
				tempEdge1->next = tempEdge2->next;
				tempEdge2->next = tempEdge1;
				changes++;
			}
			currentPred = currentPred->next;
		}
	}
}
/**
Move all edges from old list to new list, if their Y_upper is greater than threshold
*/
void listChangeList(polyEdge *rootNew, polyEdge *rootOld, polyEdge *endOld, int threshold)
{
	//polyEdge *currentPred;
	polyEdge *current;
	while (true)
	{
		current = rootOld->next;
		if (current == endOld) { break; }
		if (current->Y_upper >= threshold)
		{
			rootOld->next = current->next;
			current->next = rootNew->next;
			rootNew->next = current;
		}else {
			break;
		}
	}
}
/**
Fill segments defined by lines in list
*/
void drawActiveList(polyEdge *root, polyEdge *end)
{
	if (root->next == end) { return; }
	polyEdge *current = root;
	polyEdge *first;
	polyEdge *second;
	while (current != end)
	{
		//
		current = current->next;
		first = current;
		if (current == end) { break; }
		current = current->next;
		if (current == end) { break; }
		second = current;
		fillLine(first->X_cross, second->X_cross, first->Z_upper, second->Z_upper, first->Y_upper, first->R, first->G, first->B, second->R, second->G, second->B);
	}
}
/**
Decrement Y_upper and update colors, depth and x intersection of lines. Lines where Y_lower is greater than Y_upper are deleted
*/
void listDecrementActiveAndRemove(polyEdge *root, polyEdge *end)
{
	polyEdge *currentPred = root;
	polyEdge *current;
	while (currentPred->next != end)
	{
		current = currentPred->next;
		current->Y_upper--;
		if (current->Y_upper >= current->Y_lower)
		{
			current->X_upper += current->X_step;
			current->X_cross = current->X_upper;

			current->R += current->RD;
			current->G += current->GD;
			current->B += current->BD;

			current->Z_upper += current->Z_step;

			currentPred = currentPred->next;
		}
		else {
			currentPred->next = current->next;
			delete current;
		}
	}
}
/**
Draw polygon
*/
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
		int debugCount = 0;
		polyEdge *rootActive = new polyEdge;
		rootActive->next = new polyEdge;
		polyEdge *endActive = rootActive->next;

		polyEdge *rootPrepared = new polyEdge;
		rootPrepared->next = new polyEdge;
		polyEdge *endPrepared = rootPrepared->next;

		inputPoint4f origin;
		inputPoint4f *tempPoint1;
		inputPoint4f *tempPoint2;

		inputPoint4f *tempTempPoint;

		inputPoint4f pointStart;
		inputPoint4f pointEnd;
		
		int top = -std::numeric_limits<int>::max();
		int bottom = std::numeric_limits<int>::max();
		//printf("limits - top %d bottom %d\n", top, bottom);

		//store first point
		tempTempPoint = queue4f.front();
		tempPoint2 = new inputPoint4f();
		transformThePointAndCopyColor(tempTempPoint, *tempPoint2);
		delete tempTempPoint;
		
		/*tempPoint2 = queue4f.front();
		transformThePointAndCopyColor(tempPoint2, *tempPoint2);*/

		//DEBUG
		/*sglPointSize(5);
		drawPointNoTransform(*tempPoint2);*/
		//DEBUG ENDS

		origin = *tempPoint2;
		queue4f.pop();
		int counter = 1;

		//printf("ORIGIN %f %f\n",origin.x, origin.y);

		//build polyEdges
		while (!queue4f.empty())
		{
			tempPoint1 = tempPoint2;
			tempPoint2 = new inputPoint4f();
			tempTempPoint = queue4f.front();
			transformThePointAndCopyColor(tempTempPoint, *tempPoint2);
			delete tempTempPoint;

			//DEBUG
			/*sglPointSize(5);
			drawPointNoTransform(*tempPoint2);*/
			//DEBUG ENDS

			queue4f.pop();
			counter++;

			int tempInt1 = tempPoint1->y;
			int tempInt2 = tempPoint2->y;

			if (tempInt1 > tempInt2)
			{

				//printf("POINT %f %f\n", tempPoint2->x, tempPoint2->y);
				setPolyEdge(endPrepared, tempPoint1, tempPoint2);

				endPrepared->next = new polyEdge;
				endPrepared = endPrepared->next;

				if (tempPoint1->y > top) { top = tempPoint1->y; }
				if (tempPoint2->y < bottom) { bottom = tempPoint2->y; }
				delete tempPoint1;

				debugCount++;
			}else if (tempInt2 > tempInt1) {

				//printf("POINT %f %f\n", tempPoint2->x, tempPoint2->y);
				setPolyEdge(endPrepared, tempPoint2, tempPoint1);

				endPrepared->next = new polyEdge;
				endPrepared = endPrepared->next;

				if (tempPoint2->y > top) { top = tempPoint2->y; }
				if (tempPoint1->y < bottom) { bottom = tempPoint1->y; }
				delete tempPoint1;

				debugCount++;
			}else {
				continue;
			}
		}

		//final segment


		if (tempPoint2->y > origin.y)
		{
			setPolyEdge(endPrepared, tempPoint2, &origin);
			delete tempPoint2;

			endPrepared->next = new polyEdge;
			endPrepared = endPrepared->next;
			debugCount++;
		}
		else if (origin.y > tempPoint2->y) {
			setPolyEdge(endPrepared, &origin, tempPoint2);
			delete tempPoint2;

			endPrepared->next = new polyEdge;
			endPrepared = endPrepared->next;
			debugCount++;
		}
		else {
		}
		//printf("\nDEBUG segments %d \n", debugCount);

		//order tail
		//printf("orderind\n");
		listOrderByY_Upper(rootPrepared, endPrepared);
		//printf("ordered\n");
		//draw
		//printf("drawing from %d to %d\n", top, bottom);
		/*printf("\n");
		printf("initialized lists start\n");
		printList(rootActive, endActive);
		printf("\n");
		printList(rootPrepared, endPrepared);
		printf("initialized lists end\n");*/
		while (top >= bottom)
		{

			listDecrementActiveAndRemove(rootActive, endActive);
			listChangeList(rootActive, rootPrepared, endPrepared, top);
			listOrderByX_Cross(rootActive, endActive);
			/*printf("top %d\n",top);
			printf("lists start\n");
			printList(rootActive, endActive);
			printf("\n");
			printList(rootPrepared, endPrepared);
			printf("lists end\n");*/
			
			drawActiveList(rootActive, endActive);
			top--;
		}
		delete endPrepared;
		delete rootActive;
		delete endActive;
		delete rootPrepared;
		//printf("drawn\n");

		break;
	}
	//printf("drawMeAPolygon dont draw now \n");
}

/**
Draw triangle outlines
*/
void drawMeATriangleLineLoop(inputPoint4f* v1, inputPoint4f* v2, inputPoint4f* v3)
{
	//printf("drawTriangleLineLoop not implemented yet! \n Transformation of points is still missing\n");

	drawMeALine(v1, v2);
	drawMeALine(v2, v3);
	drawMeALine(v3, v1);
}

/**
Fill triangle defined by points t1, t2 and t3
*/
void drawMeATriangle(inputPoint4f* t1, inputPoint4f* t2, inputPoint4f* t3)
{
	//printf("trianglePrinting \n");
	inputPoint4f *p1, *p2, *p3;
	inputPoint4f *v1 = new inputPoint4f();
	inputPoint4f *v2 = new inputPoint4f();
	inputPoint4f *v3 = new inputPoint4f();

	transformThePointAndCopyColor(t1, *v1);
	transformThePointAndCopyColor(t2, *v2);
	transformThePointAndCopyColor(t3, *v3);
	
	
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
	float zInitLeft, zInitRight;
	float stepXLeft, stepXRight, stepZLeft, stepZRight;

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

			stepZLeft = (p2->z - p1->z) / stepCount;
			stepZRight = (splittingPoint.z - p1->z) / stepCount;
			zInitLeft = zInitRight = p1->z;


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
				fillLine((int)(xInitLeft + i * stepXLeft), (int)(xInitRight + i * stepXRight), zInitLeft, zInitRight, yCurrent, rL, gL, bL, rR, gR, bR);
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
				//update depth
				zInitLeft += stepZLeft;
				zInitRight += stepZRight;
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

			stepZLeft = (p2->z - p3->z) / stepCount;
			stepZRight = (splittingPoint.z - p3->z) / stepCount;
			zInitLeft = zInitRight = p3->z;


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
				fillLine((int)(xInitLeft + i * stepXLeft), (int)(xInitRight + i * stepXRight), zInitLeft, zInitRight, yCurrent, rL, gL, bL, rR, gR, bR);
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
				//update depth
				zInitLeft += stepZLeft;
				zInitRight += stepZRight;
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

			stepZLeft = (splittingPoint.z - p1->z) / stepCount;
			stepZRight = (p2->z - p1->z) / stepCount;
			zInitLeft = zInitRight = p1->z;


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
				fillLine((int)(xInitLeft + i * stepXLeft), (int)(xInitRight + i * stepXRight), zInitLeft, zInitRight, yCurrent, rL, gL, bL, rR, gR, bR);
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
				zInitLeft += stepZLeft;
				zInitRight += stepZRight;
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

			stepZLeft = (splittingPoint.z - p3->z) / stepCount;
			stepZRight = (p2->z - p3->z) / stepCount;
			zInitLeft = zInitRight = p3->z;


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
				fillLine((int)(xInitLeft + i * stepXLeft), (int)(xInitRight + i * stepXRight), zInitLeft, zInitRight, yCurrent, rL, gL, bL, rR, gR, bR);
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
				zInitLeft += stepZLeft;
				zInitRight += stepZRight;
			}
		}

	}
	delete v1, v2, v3;

	//printf("drawMeATriangle not implemented yet! \n Transformation of points is still missing\n");
}
/**
draw triangles
*/
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
	copyMatrix(matrixMVP, identityMatrix);
	// creates viewport (and there should be perspective divide)
	/*multipliedMatrix[0] = (viewportWidth - viewportOffsetX) / 2.0f;
	multipliedMatrix[5] = (viewportHeight - viewportOffsetY) / 2.0f;
	multipliedMatrix[10] = 0.5f;
	multipliedMatrix[12] = (viewportWidth / 2.0f) + viewportOffsetX;
	multipliedMatrix[13] = (viewportHeight / 2.0f) + viewportOffsetY;
	multipliedMatrix[14] = 0.5f;*/
	copyMatrix(viewportMatrix, identityMatrix);
	viewportMatrix[0] = (viewportWidth - viewportOffsetX) / 2.0f;
	viewportMatrix[5] = (viewportHeight - viewportOffsetY) / 2.0f;
	viewportMatrix[10] = 0.5f;
	viewportMatrix[12] = (viewportWidth / 2.0f) + viewportOffsetX;
	viewportMatrix[13] = (viewportHeight / 2.0f) + viewportOffsetY;
	viewportMatrix[14] = 0.5f;
	multiplyMatrix(matrixMVP, projectionStack.top());
	multiplyMatrix(matrixMVP, modelViewStack.top());

	copyMatrix(multipliedMatrix, viewportMatrix);
	multiplyMatrix(multipliedMatrix, matrixMVP);

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
	setPixel(point.x,point.y, point.r, point.g, point.b, point.z);

	point.x = xs - x;
	//drawPointNoTransform(point);
	setPixel(point.x, point.y, point.r, point.g, point.b, point.z);

	point.y = ys - y;
	//drawPointNoTransform(point);
	setPixel(point.x, point.y, point.r, point.g, point.b, point.z);

	point.x = x + xs;
	//drawPointNoTransform(point);
	setPixel(point.x, point.y, point.r, point.g, point.b, point.z);
	
	point.x = y + xs;
	point.y = x + ys;
	//drawPointNoTransform(point);
	setPixel(point.x, point.y, point.r, point.g, point.b, point.z);

	point.x = xs - y;
	//drawPointNoTransform(point);
	setPixel(point.x, point.y, point.r, point.g, point.b, point.z);

	point.y = ys - x;
	//drawPointNoTransform(point);
	setPixel(point.x, point.y, point.r, point.g, point.b, point.z);

	point.x = y + xs;
	//drawPointNoTransform(point);
	setPixel(point.x, point.y, point.r, point.g, point.b, point.z);
}

/**
Used in Bresenham's algorithm for drawing circle. One point is copied on 8
different places. Modified for filling.
@param x untranslated position of point
@param y untranslated position of point
@param xs center of circle
@param ys center of circle
@param point input point with stored colors
*/
void setSymPointsFillLine(int x, int y, int xs, int ys, float z, inputPoint4f& point) {
	float r = point.r;
	float g = point.g;
	float b = point.b;
	//printf("x %d y %d ||| XS %d XE %d Y %d \n",xs, ys, xs -x, xs + x, ys + y);
	fillLine(xs - x, xs + x, z, z, ys + y, r, g, b, r, g, b);
	fillLine(xs - x, xs + x, z, z, ys - y, r, g, b, r, g, b);
	fillLine(xs - y, xs + y, z, z, ys - x, r, g, b, r, g, b);
	fillLine(xs - y, xs + y, z, z, ys + x, r, g, b, r, g, b);
}

/**
Method to set single pixel in color buffer to specific color. (always size 1*1)
*/
inline void setPixel(float x0, float y0, float r, float g, float b, float z)
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
		int offsetC = (y*W + x) * 3;
		int offsetD = (y*W + x);
		//printf("setPixel\n");
		drawPixel(offsetC, offsetD, z, r, g, b, colorBuffer, cont->getDepthBuffer());
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
	//printf("radius %f scalefactor %f \n", radius, scaleFactor);
	//scaleFactor = 100;
	radius *= scaleFactor;

	inputPoint4f point;
	point.x = x;
	point.y = y;
	point.z = z;
	point.w = 1;
	point.r = colorVertexR;
	point.g = colorVertexG;
	point.b = colorVertexB;
	point.a = 0;

	inputPoint4f output;
	transformThePointAndCopyColor(&point, output);
	x = output.x;
	y = output.y;
	z = output.z;

	int xp, yp, p;
	switch (areaMode)
	{
	case SGL_POINT:
		sglBegin(SGL_POINTS);
		sglVertex3f(x, y, z);
		sglEnd();
		break;
	case SGL_LINE:
		// Bresenham's algorithm for drawing circle
		//printf("circle line\n");
		xp = 0;
		yp = radius;
		p = 3 - 2 * radius;
			//printf("%f %f \n", xp, yp);
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
		// Bresenham's algorithm for drawing circle
		xp = 0;
		yp = radius;
		p = 3 - 2 * radius;
		while (xp < yp) {
			setSymPointsFillLine(xp, yp, x, y, z, point);
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
			setSymPointsFillLine(xp, yp, x, y, z, point);

		//printf("No Circle filling algorithm implemented right now.\n");
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

	float scaleFactor = sqrt(matrixMVP[0] * matrixMVP[5] - matrixMVP[1] * matrixMVP[4]);
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
	transformThePointAndCopyColor(&point, output);
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
		setPixel(x + tempX, y + tempY, colorVertexR, colorVertexG, colorVertexB, output.z);
		setPixel(x - tempX, y + tempY, colorVertexR, colorVertexG, colorVertexB, output.z);
		setPixel(x + tempX, y - tempY, colorVertexR, colorVertexG, colorVertexB, output.z);
		setPixel(x - tempX, y - tempY, colorVertexR, colorVertexG, colorVertexB, output.z);
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
		setPixel(x + tempX, y + tempY, colorVertexR, colorVertexG, colorVertexB, output.z);
		setPixel(x - tempX, y + tempY, colorVertexR, colorVertexG, colorVertexB, output.z);
		setPixel(x + tempX, y - tempY, colorVertexR, colorVertexG, colorVertexB, output.z);
		setPixel(x - tempX, y - tempY, colorVertexR, colorVertexG, colorVertexB, output.z);
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
		//printf("check sgl.cpp sglEllipseSegmented fill branch.\n");
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

void sglEllipse(float x, float y, float z, float a, float b) {
	#ifdef ELLIPSE
	sglEllipseSegmented(x, y, z, a, b);
	#elif ELLIPSE_SECOND
		sglEllipseSecond(x, y, z, a, b);
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
		setPixel(point->x, point->y, point->r, point->g, point->b, point->z);

	point->x = xs - x;
	angle = acos(-x / radius);
	if (y < 0)
		angle = -angle + 2 * 3.14159274;
	if (angle >= from && angle <= to)
		setPixel(point->x, point->y, point->r, point->g, point->b, point->z);

	point->y = ys - y;
	angle = acos(-x / radius);
	if (-y < 0)
		angle = -angle + 2 * 3.14159274;
	if (angle >= from && angle <= to)
		setPixel(point->x, point->y, point->r, point->g, point->b, point->z);

	point->x = x + xs;
	angle = acos(x / radius);
	if (-y < 0)
		angle = -angle + 2 * 3.14159274;
	if (angle >= from && angle <= to)
		setPixel(point->x, point->y, point->r, point->g, point->b, point->z);

	point->x = y + xs;
	point->y = x + ys;
	angle = acos(y / radius);
	if (x < 0)
		angle = -angle + 2 * 3.14159274;
	if (angle >= from && angle <= to)
		setPixel(point->x, point->y, point->r, point->g, point->b, point->z);

	point->x = xs - y;
	angle = acos(-y / radius);
	if (x < 0)
		angle = -angle + 2 * 3.14159274;
	if (angle >= from && angle <= to)
		setPixel(point->x, point->y, point->r, point->g, point->b, point->z);

	point->y = ys - x;
	angle = acos(-y / radius);
	if (-x < 0)
		angle = -angle + 2 * 3.14159274;
	if (angle >= from && angle <= to)
		setPixel(point->x, point->y, point->r, point->g, point->b, point->z);

	point->x = y + xs;
	angle = acos(y / radius);
	if (-x < 0)
		angle = -angle + 2 * 3.14159274;
	if (angle >= from && angle <= to)
		setPixel(point->x, point->y, point->r, point->g, point->b, point->z);
}

void sglArc(float x, float y, float z, float radius, float from, float to) {
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
		sglBegin(SGL_POLYGON);
		segments = 40;
		angle = from;
		delta = (to - from) / segments;
		sglVertex3f(x, y, z);
		for (int i = 0; i <= segments; i++)
		{
			sglVertex3f(x + radius*cos(angle), y + radius*sin(angle), z);
			angle += delta;
		}
		sglEnd();
		//printf("No Arc filling algorithm implemented right now. \n");
		break;
	}
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
	float C = -(far + near) / (far - near);
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
