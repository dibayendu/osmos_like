#ifndef UTILITY_H
#define UTILITY_H

#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <sys/time.h>

#ifdef _WIN32
#include <windows.h>
#endif
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#define ARENA_RADIUS 9.0

typedef double Real;

/* Arena. */
typedef struct
{
  	Real min[2], max[2];
  	Real momentum[2];
} Arena;

typedef struct {
	float x, y, z;
} vec3f;

typedef struct {
	float x, y;
} vec2f;

typedef struct {
	float x, y, z, w;
} vec4f;

void panic(const char *m); 

#endif
