#ifndef PARTICLE_H
#define PARTICLE_H

#include "utility.h"

/* Particles (particles). */
typedef struct
{
  	Real position[2];
  	Real velocity[2];
  	Real radius;
  	Real mass;
  	Real elasticity;
  	GLUquadric *quadric;  /* For rendering. */
  	int slices, loops;    /* For rendering. */
	bool alive;
	float colour[3];
} Particle;

Particle* add_particle(Particle *particle, int number_of_particles, double mass, double pos_x, double pos_y, double vel_x, double vel_y);
Particle* reinitialise_particles(Particle *particle, int *number_of_particles, int number_of_dead_particles);
Particle* initialiseParticlesRandomly(int number_of_particles, Arena *arena, Real *total_mass, bool square_arena);
void displayParticles(Particle *particle, int number_of_particles, bool player_dead);
float sumKineticEnergy(Particle *particle, int number_of_particles);
void sumMomentum(Real *p, Particle *particle, int number_of_particles, Arena *arena);
bool collided_with_circular_arena(Particle *p);

#endif
