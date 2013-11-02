#ifndef COLLISION_H
#define COLLISION_H

#include "utility.h"
#include "particle.h"

/* Control collision reaction calculation */
enum ReactionCalculation 
{
  	basisChange,
  	projNormal
};

void collideParticlesWall(Particle *particle, int number_of_particles, Arena *arena, bool square_arena);
void eulerStepSingleParticle(Particle *p, Real h);
int collideParticlesBruteForce(Particle *particle, int number_of_particles, Arena *arena, Real h, bool square_arena);
int collideParticlesUniformGrid(Particle *particle, int number_of_particles, Arena *arena, Real h, bool square_arena);

#endif
