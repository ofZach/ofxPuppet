#ifndef __RMS_DEBUG_H__
#define __RMS_DEBUG_H__

//#define USE_WINDOWS_DEBUGPRINT

#ifdef USE_WINDOWS_DEBUG
#include "windows.h"
#pragma warning(disable:4505)
#endif

#include <string>
using namespace std;

#include <stdio.h>
#include <stdarg.h>

#ifndef ASSERT
#ifdef _DEBUG
#define ASSERT(x)  { if(! (x) ) abort(); }
#else
#define ASSERT(x) 
#endif
#endif

static void _RMSInfo(char * str, ...)
{	
	static char buf[1024];
	va_list args;

	va_start(args, str);
	sprintf(buf, str, args);
#ifdef USE_WINDOWS_DEBUG
	OutputDebugString(buf);
#else
	fprintf(stderr, "%s", buf);
#endif
	va_end(args);
}

static string _RMSInfoString(char * str, ...)
{	
	static char buf[1024];
	va_list args;

	va_start(args, str);
	sprintf(buf, str, args);
	va_end(args);
	return string(buf);
}


#ifndef USE_WINDOWS_DEBUG
static void Debugbreak() {
	abort();
}
#endif


/*
 *
 * OpenGL debug utility functions
 *
 */


#ifdef __RMS_DEBUG_WANT_GL


#include <GL/gl.h>


static char _RMSDEBUG_GLErrorStrings[][32] = {  "GL_NO_ERROR",
												"GL_INVALID_ENUM", 
												"GL_INVALID_VALUE", 
												"GL_INVALID_OPERATION", 
												"GL_STACK_OVERFLOW",
												"GL_STACK_UNDERFLOW",
												"GL_OUT_OF_MEMORY",
												"Unknown GL Error"};

static char * glErrorString(GLenum error) 
{
	switch(error) {
		case GL_NO_ERROR: return _RMSDEBUG_GLErrorStrings[0];
		case GL_INVALID_ENUM: return _RMSDEBUG_GLErrorStrings[1];
		case GL_INVALID_VALUE: return _RMSDEBUG_GLErrorStrings[2];
		case GL_INVALID_OPERATION: return _RMSDEBUG_GLErrorStrings[3];
		case GL_STACK_OVERFLOW: return _RMSDEBUG_GLErrorStrings[4];
		case GL_STACK_UNDERFLOW: return _RMSDEBUG_GLErrorStrings[5];
		case GL_OUT_OF_MEMORY: return _RMSDEBUG_GLErrorStrings[6];
		default: return _RMSDEBUG_GLErrorStrings[7];
	}
}

static GLenum _RMSGLError = GL_NO_ERROR;

#define _RMSHaveGLError() (_RMSGLError != GL_NO_ERROR)

static bool glError()
{
	_RMSGLError = glGetError();
	int nextErr = _RMSGLError;
	while (nextErr != GL_NO_ERROR)
		nextErr = glGetError();
	return _RMSGLError != GL_NO_ERROR;
}



#endif // __RMS_DEBUG_WANT_GL

#endif   // __RMS_DEBUG_H__