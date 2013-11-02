#include "particle.h"
#include "texture.h"

#define SLICES_COUNT 32
#define LOOPS_COUNT 32

GLuint particle_image;
GLuint player_image;

Real random_uniform() 
{
  	return rand()/(float)RAND_MAX;
}

/* Small number to handle numerical imprecision. */
const Real epsilon = 1.0e-6;

/* copy particle from source to destination */
void copy_particle(Particle source, Particle *destination)
{
	destination->position[0] = source.position[0];
	destination->position[1] = source.position[1];
	destination->velocity[0] = source.velocity[0];
	destination->velocity[1] = source.velocity[1];
    destination->mass = source.mass;
	destination->radius = source.radius;
	destination->elasticity = source.elasticity ;
    destination->quadric = source.quadric ;
    destination->slices = source.slices ;
    destination->loops = source.loops ;
	destination->alive = source.alive ;
	destination->colour[0] = source.colour[0];
	destination->colour[1] = source.colour[1];
	destination->colour[2] = source.colour[2];
}

/* setting up the player particle in the center of the screen */
void setup_player_particle(Particle *particle, GLUquadric *quadric)
{
  	const Real maxVelocity = 1.0;
	
	particle->position[0] = particle->position[1] = 0.0;
	particle->velocity[0] = (random_uniform() - 0.5) * maxVelocity;
	particle->velocity[1] = (random_uniform() - 0.5) * maxVelocity;
    particle->mass = 0.05;
	particle->radius = sqrt(particle->mass);
	particle->elasticity = 1.0;
    particle->quadric = quadric;
    particle->slices = SLICES_COUNT;
    particle->loops = LOOPS_COUNT;
	particle->alive = true;
	particle->colour[0] = particle->colour[1] = particle->colour[2] = 1.0;
}

/* when a mote wants to gain velocity at a particular direction, it releases a new mote which gets added in the environment */
Particle* add_particle(Particle *particle, int number_of_particles, double mass, double pos_x, double pos_y, double vel_x, double vel_y)
{
	// allocating memory for all the particles
	Particle *new_particle = (Particle*)malloc(sizeof(Particle) * (number_of_particles+1));
	
	// copying all the particles in a new array so that the new particle can be added at the end.
	int i, j = 0;
  	for (i = 0; i < number_of_particles; i++)
		copy_particle(particle[i], &new_particle[j++]);

	// adding the new particle at the end of the array
	new_particle[number_of_particles].position[0] = pos_x;
	new_particle[number_of_particles].position[1] = pos_y;
	new_particle[number_of_particles].velocity[0] = vel_x;
	new_particle[number_of_particles].velocity[1] = vel_y;
    new_particle[number_of_particles].mass = mass;
	new_particle[number_of_particles].radius = sqrt(mass);
	new_particle[number_of_particles].elasticity = 1.0;
    new_particle[number_of_particles].quadric = gluNewQuadric();
    new_particle[number_of_particles].slices = SLICES_COUNT;
    new_particle[number_of_particles].loops = LOOPS_COUNT;
	new_particle[number_of_particles].alive = true;
	new_particle[number_of_particles].colour[0] = particle[0].colour[0];
	new_particle[number_of_particles].colour[1] = particle[0].colour[1];
	new_particle[number_of_particles].colour[2] = particle[0].colour[2];
	
	return new_particle;
}

/* removes all the dead motes in the environment and reinitialises the motes array */
Particle* reinitialise_particles(Particle *particle, int *number_of_particles, int number_of_dead_particles)
{
	// allocating memory for the particles since some of the particles died
	Particle *new_particle = (Particle*)malloc(sizeof(Particle) * (*number_of_particles - number_of_dead_particles));
	
	// copying all the paticles in the new array as some of them died
	int i, j = 0;
  	for (i = 0; i < *number_of_particles; i++)
	{
		// do not include the particle if it is dead
		if(!particle[i].alive)
			continue;

		copy_particle(particle[i], &new_particle[j++]);
	}

	return new_particle;
}

/* returns true if the particle collided with the walls of the circular arena */
bool collided_with_circular_arena(Particle *p)
{
	Real dist_from_center = sqrt( (p->position[0]*p->position[0]) + (p->position[1]*p->position[1]) );
	dist_from_center += p->radius;	// the diameter of the particle should be inside the arena

	return (dist_from_center < 9.0) ? false : true;
}

/* particles are initialised in the arena and given random positions in the environment */
Particle* initialiseParticlesRandomly(int number_of_particles, Arena *arena, Real *total_mass, bool square_arena)
{
  	GLUquadric *quadric = gluNewQuadric();
  	const Real maxVelocity = 1.0;
  	Real n[2], n_mag_sq, sum_radii, sum_radii_sq;
  	bool collision, done;
  	int i, j;

	// giving a texture to the particles
	particle_image = texture_load("images/sparkle.jpg");
	player_image = texture_load("images/player.png");

	// allocating memory for the particles
	Particle *particle = (Particle*)malloc(sizeof(Particle) * number_of_particles);

  	for (i = 0; i < number_of_particles; i++)
	{
		if(i == 0)
		{
			setup_player_particle(&particle[0], quadric);
			continue;
		}

   		particle[i].velocity[0] = (random_uniform() - 0.5) * maxVelocity;
    	particle[i].velocity[1] = (random_uniform() - 0.5) * maxVelocity;
    	particle[i].mass = random_uniform() * 0.05;
    	particle[i].radius = sqrt(particle[i].mass);
    	particle[i].elasticity = 1.0;
    	particle[i].quadric = quadric;
    	particle[i].slices = SLICES_COUNT;
    	particle[i].loops = LOOPS_COUNT;
		particle[i].alive = true;
		// giving the particles some random colour
		float c = (float)random_uniform();
		particle[i].colour[0] = c < 0.5 ? (1.0 - c) : c;
		c = (float)random_uniform();
		particle[i].colour[1] = c < 0.5 ? (1.0 - c) : c;
		c = (float)random_uniform();
		particle[i].colour[2] = c < 0.5 ? (1.0 - c) : c;

		*total_mass += particle[i].mass;

		/* making sure all the particles are spaced out and no particle is overlapping */
    	done = false;
		while (!done) 
		{
      		particle[i].position[0] = random_uniform() *
        		(arena->max[0] - arena->min[0] - 2.0 * particle[i].radius) + 
        			arena->min[0] + particle[i].radius + epsilon;
      		particle[i].position[1] = random_uniform() *
        		(arena->max[1] - arena->min[1] - 2.0 * particle[i].radius) + 
        			arena->min[1] + particle[i].radius + epsilon;
			
			/* if the arena is circular checking if the particle is created inside the circular arena */
			if(!square_arena && collided_with_circular_arena(&particle[i]))
				continue;		// if the particles is not inside the ciruclar arena, try initialising its position again

      		/* Check for collision with existing particles. */
      		collision = false;
      		j = 0;
			while (!collision && j < i) 
			{
        		sum_radii = particle[i].radius + particle[j].radius;
        		sum_radii_sq = sum_radii * sum_radii;
        		n[0] = particle[j].position[0] - particle[i].position[0];
        		n[1] = particle[j].position[1] - particle[i].position[1];
        		n_mag_sq = n[0] * n[0] + n[1] * n[1];
        		if (n_mag_sq < sum_radii_sq)
          			collision = true;
        		else
          			j++;
      		}
      		
			if (!collision) 
        		done = true;
    	}
      	// printf ("initialiseParticles: x %f y %f\n", particle[i].position[0], particle[i].position[1]);
  	}

	return particle;
}

/* using inbild functions in the GLU library to draw the particles in the environment */
void displayParticle(Particle *p, float sx, float sy, float sz)
{
  	glPushMatrix();
  	glScalef(sx, sy, sz);
	gluQuadricDrawStyle(p->quadric, GLU_FILL);
	gluQuadricTexture(p->quadric, GL_TRUE);
  	gluDisk(p->quadric, 0.0, p->radius, p->slices, p->loops);
  	glPopMatrix();
}

/* drawing all the live particles in the environment giving them some texture and some random colour */
void displayParticles(Particle *particle, int number_of_particles, bool player_dead)
{
  	int i;

  	/* Display particles. */
  	for (i = 0; i < number_of_particles; i++) 
	{
		if(!particle[i].alive)
			continue;

      	// printf ("displayParticles: x %f y %f\n", particle[i].position[0], particle[i].position[1]);
		
		glPushAttrib(GL_ENABLE_BIT);
    	glPushMatrix();
    		glEnable(GL_BLEND);
			glColor4f(particle[i].colour[0], particle[i].colour[1], particle[i].colour[2], 1.0f);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE);
   			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    		glEnable(GL_DEPTH_TEST);
    		glEnable(GL_TEXTURE_2D);
    		glTranslatef(particle[i].position[0], particle[i].position[1], 0.0); 
			
			if(i == 0 && !player_dead) 		// player gets a different texture
				glBindTexture(GL_TEXTURE_2D, player_image);
    		else
				glBindTexture(GL_TEXTURE_2D, particle_image);
				
			displayParticle(&particle[i], 1.0, 1.0, 1.0);
			glDisable(GL_TEXTURE_2D);
    		glDisable(GL_DEPTH_TEST);
			glDisable(GL_BLEND);
    	glPopMatrix();
		glPopAttrib();
  	}
}

/* kinetic energy is not conserved in the applicaition but the method is there to debug */
float sumKineticEnergy(Particle *particle, int number_of_particles) 
{
  	Real v_sq, K;
  	K = 0;
  	
	for (int i = 0; i < number_of_particles; i++) 
	{
    	v_sq = particle[i].velocity[0] * particle[i].velocity[0] +
      			particle[i].velocity[1] * particle[i].velocity[1];
    	K += 0.5 * particle[i].mass * v_sq;
  	}

  	return K;
}

/* calcuate the total momentum of the environment */
void sumMomentum(Real *p, Particle *particle, int number_of_particles, Arena *arena) 
{
  	p[0] = p[1] = 0;
  	for (int i = 0; i < number_of_particles; i++) 
	{
   		p[0] += particle[i].mass * particle[i].velocity[0];
    	p[1] += particle[i].mass * particle[i].velocity[1];
  	}
  	p[0] += arena->momentum[0];
  	p[1] += arena->momentum[1];
}
 
