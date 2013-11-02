#include "arena.h"
#include "texture.h"
#include <GL/gl.h>
#include <GL/glut.h>

GLuint circular_arena_image;
GLUquadric *quadric;

/* drawing a quad which defines the square arena where the motes should be bound within */
void draw_square_arena(Arena *arena)
{
  	glBegin(GL_LINE_LOOP);
  	glVertex3f(arena->min[0], arena->min[1], 0.0);
  	glVertex3f(arena->max[0], arena->min[1], 0.0);
  	glVertex3f(arena->max[0], arena->max[1], 0.0);
  	glVertex3f(arena->min[0], arena->max[1], 0.0);
  	glEnd();
}

/* drawing a circle which defines the circular arena where the motes should be bound within */
void draw_circular_arena()
{
	// only draw perimeter lines of circle
	gluQuadricDrawStyle(quadric, GLU_SILHOUETTE);
  	gluDisk(quadric, 0.0, ARENA_RADIUS, 64, 64);
}

/* initialises both the square and the circular arena */
void initialiseArena(Arena *arena)
{
  	arena->min[0] = -HALF_LENGTH_ARENA;
  	arena->min[1] = -HALF_LENGTH_ARENA;
  	arena->max[0] = HALF_LENGTH_ARENA;
  	arena->max[1] = HALF_LENGTH_ARENA;

  	arena->momentum[0] = 0.0;
  	arena->momentum[1] = 0.0;
	
	circular_arena_image = texture_load("images/background.jpg");
 	quadric = gluNewQuadric();
}

/* draws the ciruclar or square arena depending on the parameter passed */
void displayArena(Arena *arena, bool square_arena)
{
	glPushAttrib(GL_ENABLE_BIT);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glPushMatrix();	
	if(square_arena)
		draw_square_arena(arena);
	else
		draw_circular_arena();
	glPopMatrix();
	glPopAttrib();
}
