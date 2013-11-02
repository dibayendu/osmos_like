/* sdl-base contains opengl/SDL init code and provides
 * a main loop. this file implements expected callback
 * functions. see sdl-base.h */
#include "sdlbase.h"
#include "utility.h"
#include "arena.h"
#include "particle.h"
#include "collision.h"

Real total_mass;
bool game_over;
bool player_dead;
bool player_won;
int number_of_particles;
bool square_arena;
Particle *particle;
Particle *player;
Arena arena;

/* Rendering info. */
enum renderMode 
{ 
	wire, 
	solid 
};
static int renMode = wire; 
static Real elapsedTime = 0.0, startTime = 0.0;
static Real time = 0.0, h;
static const int milli = 1000;
static bool go = false;
float orthoValue;

/* Collision detection method. */
enum CollisionDetectionMethod 
{
  	bruteForce,
  	uniformGrid
};

int CDmethod = bruteForce;

/* Basic camera struct */
typedef struct 
{
	int zooming;
	float rotX, rotY;
	float zoom;
	float sensitivity;
} Camera;

/* Scene globals */
Camera camera;
float currentFramerate;
float currentFrametime;
float aspect;
int windowWidth;
int windowHeight;
int lastMouseX = 0;
int lastMouseY = 0;

static const float pi = 3.14159265f;

void setRenderMode(int rm)
{
  	/* Example of GNU C/C++ brace indentation style.  */
  	if (rm == wire) 
    {
      	glDisable(GL_LIGHTING);
      	glDisable(GL_DEPTH_TEST);
      	glDisable(GL_NORMALIZE);
      	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    }
  	else if (rm == solid)
    {
      	glEnable(GL_NORMALIZE);
      	glShadeModel(GL_SMOOTH);
      	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }
}

/* move the particles in the direction of its velocity */
void integrateMotionParticles(Real h)
{
  	for (int i = 0; i < number_of_particles; i++)
    	eulerStepSingleParticle(&particle[i], h);
}

void updateParticles(void)
{
  	/* Calculate time increment. */
  	h = elapsedTime - time;
  	time = elapsedTime;
    
	//printf("updateParticles: time %f %f\n", time, h);

  	/* Compute new positions of particles. */
  	integrateMotionParticles(h);

  	/* Collisions against walls. */
  	collideParticlesWall(particle, number_of_particles, &arena, square_arena);

	int number_of_dead_particles;

  	/* Collisions amongst particles. */
  	if (CDmethod == bruteForce)
    	number_of_dead_particles = collideParticlesBruteForce(particle, number_of_particles, &arena, h, square_arena);
  	else
	{
    	number_of_dead_particles = collideParticlesUniformGrid(particle, number_of_particles, &arena, h, square_arena);
		/* number_of_dead_particles can return -1
		 * this is incase gridCellSize is smaller than the particle diameter
		 * in such a situation brute force will be used
		 */
		number_of_dead_particles = number_of_dead_particles < 0 ? collideParticlesBruteForce(particle, number_of_particles, &arena, h, square_arena) : number_of_dead_particles;
	}

	if(!player->alive)
	{
		game_over = true;
		player_dead = true;
		player->position[0] = player->position[1] = 0.0;
		player->velocity[0] = player->velocity[1] = 0.0;
	}
	if(number_of_dead_particles > 0)
	{
		// all the dead particles need to be removed from the particle array
		Particle *new_particle = reinitialise_particles(particle, &number_of_particles,  number_of_dead_particles);
		free(particle);
		number_of_particles -= number_of_dead_particles;
		particle = new_particle;
		player = &particle[0];		// the first particle in the array is always the player
	}
	
	if(player->mass >= (total_mass / 2.0))
	{
		game_over = true;
		player_won = true;
	}

	/*	used for debugging
	 	
    	Real p[2];
    	sumMomentum(p, particle, number_of_particles, &arena);
    	printf("p = %f %f\n", p[0], p[1]);
    	
		//printf("K = %f\n", sumKineticEnergy(particle, number_of_particles));
  	*/
}

void initialise_all_game_variables_and_methods()
{
	glClearColor(0, 0, 0, 0);
	
	memset(&camera, 0, sizeof(Camera));
	camera.sensitivity = 0.3f;
	camera.zoom = 2.0f;
	
	total_mass = 0.0;
	time = h = 0.0;
	game_over = false;
	player_dead = false;
	player_won = false;
	number_of_particles = 100;
	renMode = solid; 
	elapsedTime = startTime = 0.0;
	go = false;
	orthoValue = 0.0;

  	setRenderMode(renMode);
  	initialiseArena(&arena);
   	particle = initialiseParticlesRandomly(number_of_particles, &arena, &total_mass, square_arena);
	player = &particle[0];
}

/* Called once at program start */
void init(bool square)
{
	int argc = 0;  /* fake glutInit args */
	char *argv = "";
	glutInit(&argc, &argv);
	
	square_arena = square;
	initialise_all_game_variables_and_methods();
}

/* Called once at start and again on window resize */
void reshape(int width, int height)
{
	windowWidth = width;
	windowHeight = height;
	
	/* Portion of viewport to render to */
	glViewport(0, 0, width, height);
	
	/* Calc aspect ratio */
	aspect = width / (float)height;
	
	/* Begin editing projection matrix */
	glMatrixMode(GL_PROJECTION);
	
	/* Clear previous projection */
	glLoadIdentity();
	
	orthoValue = 5.0 * camera.zoom;
	float left = -orthoValue * aspect;
	float right = orthoValue * aspect;
	float bottom = -orthoValue;
	float top = orthoValue;
  	glOrtho(left, right, bottom, top, -1.0, 1.0);
	
	/* used for debugging
		printf("Ortho value : %f\n", right);
		printf("width : %d\n", width);
		printf("height : %d\n", height);
	*/

	/* Restore modelview as current matrix */
	glMatrixMode(GL_MODELVIEW);
}

/* used to display the characters on the screen */
void draw_game_details(float *colour, int pos_x, int pos_y, char *bufp, char *buffer, void *font)
{
	glPushMatrix();
	glLoadIdentity();
	glColor3fv(colour);
	glRasterPos2i(pos_x, pos_y);
	for (bufp = buffer; *bufp; bufp++)
		glutBitmapCharacter(font, *bufp);
	glPopMatrix();
}

/* not only draws the OSD but also other game details which the user needs to know */
void drawOSD()
{
	char *bufp;
	char buffer[32];
	float colour[3];
	
	/* Backup previous "enable" state */
	glPushAttrib(GL_ENABLE_BIT);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);

	/* Create a temporary orthographic projection, matching
	 * window dimensions, and push it onto the stack */
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0, windowWidth, 0, windowHeight, -1, 1);
	glMatrixMode(GL_MODELVIEW);
	
	/* draw framerate */
	snprintf(buffer, sizeof(buffer), "Frame_Rate: %d", (int)currentFramerate);
	colour[0] = 1.0f; 	colour[1] = 1.0f;	colour[2] = 0.0f;	// yellow
	draw_game_details(colour, 10, 10, bufp, buffer, GLUT_BITMAP_9_BY_15);
	
	
	/* draw frametime */
	snprintf(buffer, sizeof(buffer), "Frame_Time: %d", (int)currentFrametime);
	colour[0] = 0.0f; 	colour[1] = 1.0f;	colour[2] = 1.0f;	// turquoise
	draw_game_details(colour, 200, 10, bufp, buffer, GLUT_BITMAP_9_BY_15);
	
	/* draw game over */
	if(game_over)
	{
		colour[0] = 1.0f; 	colour[1] = 1.0f;	colour[2] = 1.0f;	// white
		snprintf(buffer, sizeof(buffer), " GAME OVER!");
		draw_game_details(colour, windowWidth/2, windowHeight/2, bufp, buffer, GLUT_BITMAP_HELVETICA_18);
		
		if(player_dead)
		{
			snprintf(buffer, sizeof(buffer), " You Lost! press 'u' to restart the game.");
			draw_game_details(colour, 400, 10, bufp, buffer, GLUT_BITMAP_HELVETICA_18);
		}
		else if(player_won)
		{
			snprintf(buffer, sizeof(buffer), " You WON! press 'u' to restart the game.");
			draw_game_details(colour, 400, 10, bufp, buffer, GLUT_BITMAP_HELVETICA_18);
		}
	}
	
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();  /* pop projection */
	glMatrixMode(GL_MODELVIEW);

	/* Restore "enable" state */
	glPopAttrib();
}

/* called every from to display the environment and all its variables */
void display()
{
  	static int frameNo = 0;
  	GLenum err;
	
	/* Clear the colour and depth buffer */
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  	glColor3f (0.8, 0.8, 0.8);
	
	/* Camera transformations */
	glLoadIdentity();
	
	/* if the player is dead, the camera comes to the center of the arena */
	if(player_dead)
  		glTranslatef(0.0f, 0.0f, 0.0f);
	else
  		glTranslatef(-player->position[0], -player->position[1], 0.0f);

  	glPushMatrix();
  		glPushMatrix();
  			displayArena(&arena, square_arena);
  			displayParticles(particle, number_of_particles, player_dead);
		glPopMatrix();
	
  		glPushMatrix();
			drawOSD();
  		glPopMatrix();
  	glPopMatrix();

  	/* Check for errors. */
  	while ((err = glGetError()) != GL_NO_ERROR)
    	printf("%s\n",gluErrorString(err));
}

/* Called continuously. dt is time between frames in seconds */
void update(float dt)
{
	static float fpsTime = 0.0f;
	static int fpsFrames = 0;
	fpsTime += dt;
	fpsFrames += 1;
	if (fpsTime > 1.0f)
	{
		currentFramerate = fpsFrames / fpsTime;
		currentFrametime = 1000.0f / currentFramerate;
		fpsTime = 0.0f;
		fpsFrames = 0;
	}
	if (!go) 
    	return;

  	elapsedTime = glutGet(GLUT_ELAPSED_TIME) / (Real)milli - startTime;
  	updateParticles();
}

/* motes the player particle when mouse clicked */
void move_player(int x, int y)
{
	/* no mouse interaction can be received when the game is over */
	if(game_over)
		return;
	
	/* relative position of the click on the 2D screen is calculated */
	float rel_x = ( ( ((float)x / (float)windowWidth) * 2 * (float)orthoValue ) - (float)orthoValue ) * (float) aspect;
	float rel_y = -( ( ((float)y / (float)windowHeight) * 2 * (float)orthoValue ) - (float)orthoValue);
	
	/* calculating the relative angle of the click and figuring out the angle relative to the player */
	float d_x = 0.0 - rel_x;	// since the player particle is always in the center of the display
	float d_y = 0.0 - rel_y;
	float value = d_x / d_y;
	value = value < 0 ? -value : value;
	float radian = atan(value);

	/* velocity of the player particle in the opposite direction of the click */
	double vel_x = rel_x < 0.0 ? (double)sinf(radian) : -(double)sinf(radian);
	double vel_y = rel_y < 0.0 ? (double)cosf(radian) : -(double)cosf(radian);
	
	/* calculating the mass and the relative position of the releasing particle */
	double new_particle_mass = player->mass / 100.0;
	double np_pos_x = player->position[0] + (player->radius* -vel_x) + (sqrt(new_particle_mass)* -vel_x);
	double np_pos_y = player->position[1] + (player->radius * -vel_y) + (sqrt(new_particle_mass)* -vel_y);

	/* decreasing the velocity of the released particle by a certain constant */
	vel_x /= 10.0;
	vel_y /= 10.0;

	/* calculating the velocity of the released particle in the direction of its parent */
	Real n_v_x = ( (player->mass - new_particle_mass) / new_particle_mass ) * vel_x;
	Real n_v_y = ( (player->mass - new_particle_mass) / new_particle_mass ) * vel_y;
	
	/* adding new particle to the array by passing the mass, position and velocity of the new particle in the opposite direction of the player */
	Particle *new_particle = add_particle(particle, number_of_particles, new_particle_mass, np_pos_x, np_pos_y, (player->velocity[0]-n_v_x), (player->velocity[1]-n_v_y) );
	free(particle);
	number_of_particles++;
	particle = new_particle;
	player = &particle[0];
  	
	/* since the player released some mass, it now has a new mass and velocity so that the momentum of the environment is conserved */
	player->mass -= new_particle_mass;
	player->radius = (player->mass <= 0) ? 0 : sqrt(player->mass);
	player->alive = !(player->mass <= 0);
	player->velocity[0] += vel_x;
  	player->velocity[1] += vel_y;

	/*	the code below is commented as this is used for debugging

		printf("click_x : %f\n", rel_x);
		printf("click_y : %f\n", rel_y);
		printf("player->position[0] : %f\n", player->position[0]);
		printf("player->position[1] : %f\n", player->position[1]);
		float degree = (radian * 180.0) / 22 * 7;
		degree = degree < 0 ? -degree : degree;
		printf("degree: %f\n", degree);
		printf("velocity[0]: %f\n", player->velocity[0]);
		printf("velocity[1]: %f\n", player->velocity[1]);
		printf("d_x : %f\n", d_x);
		printf("d_y : %f\n", d_y);
		printf("\n\n");
	*/
}

void mouseDown(int button, int state, int x, int y)
{
	if (button == SDL_BUTTON_LEFT)
	{
		if(state == 1)
			move_player(x, y);
	}
	if (button == SDL_BUTTON_RIGHT)
		camera.zooming = (state == 1);
}

void mouseMove(int x, int y)
{
	int dx = x - lastMouseX;
	int dy = y - lastMouseY;
	if (camera.zooming)
	{
		camera.zoom -= dy * camera.zoom * camera.sensitivity * 0.03f;
		/* in ortho project if zoomed in and zoomed out, calling all the parameters to resize the projection */
		reshape(windowWidth, windowHeight);
	}
	lastMouseX = x;
	lastMouseY = y;
}

void keyDown(int key)
{
	switch (key) 
  	{
		/* press q or esc to quit the application */
  		case SDLK_q:
		case SDLK_ESCAPE:
			free(particle);
    		exit(EXIT_SUCCESS);
    		break;
  		case SDLK_w:
			renMode = (renMode == wire) ? solid : wire;
			setRenderMode(renMode);
    		break;
  		case SDLK_d:
			if(CDmethod == uniformGrid)
			{
				CDmethod = bruteForce;
				printf("Collision detection: brute force\n");
			}
			else
			{
				CDmethod = uniformGrid;
				printf("Collision detection: uniform grid\n");
			}
    		break;
  		case SDLK_s:
    		if (!go) 
			{
      			startTime = glutGet(GLUT_ELAPSED_TIME) / (Real)milli;
      			go = true;
    		}
			else
			{
				go = false;
				elapsedTime = 0.0, startTime = 0.0;
				time = 0.0, h;
			}
    		break;
		case SDLK_u:
			free(particle);
			initialise_all_game_variables_and_methods();
			reshape(windowWidth, windowHeight);
			break;	
	}
}

void keyUp(int key)
{
}

void cleanup()
{
}

