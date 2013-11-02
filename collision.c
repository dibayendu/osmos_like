#include "collision.h"

void eulerStepSingleParticle(Particle *p, Real h)
{
  	p->position[0] += h * p->velocity[0];
  	p->position[1] += h * p->velocity[1];
}

/* calculate the particles position, velocity and the momentum of the arena is conserved */
void collide_square_wall(Particle *p, Arena *a)
{
  	float dp[2];
  	dp[0] = dp[1] = 0.0;
  
  	if ((p->position[0] - p->radius) < a->min[0]) 
	{
    	p->position[0] += 
      		2.0 * (a->min[0] - (p->position[0] - p->radius));
    	p->velocity[0] *= -1.0;
    	dp[0] += p->mass * -2.0 * p->velocity[0];
  	}
  	if ((p->position[1] - p->radius) < a->min[1]) 
	{
    	p->position[1] += 
      		2.0 * (a->min[1] - (p->position[1] - p->radius));
    	p->velocity[1] *= -1.0;
    	dp[1] += p->mass * -2.0 * p->velocity[1];
  	}
  	if ((p->position[0] + p->radius) > a->max[0]) 
	{
    	p->position[0] -= 
      		2.0 * (p->position[0] + p->radius - a->max[0]);
    	p->velocity[0] *= -1.0;
    	dp[0] += p->mass * -2.0 * p->velocity[0];
  	}
  	if ((p->position[1] + p->radius) > a->max[1]) 
	{
   		p->position[1] -= 
      		2.0 * (p->position[1] + p->radius - a->max[1]);
    	p->velocity[1] *= -1.0;
    	dp[1] += p->mass * -2.0 * p->velocity[1];
  	}
  	a->momentum[0] += dp[0];
  	a->momentum[1] += dp[1];
}

Real get_distance(Real x1, Real y1, Real x2, Real y2)
{
	Real x2_x1 = x2 - x1;
	Real y2_y1 = y2 - y1;
	Real distance = sqrt( (x2_x1 * x2_x1) + (y2_y1 * y2_y1) );
	return distance;
}

/* using 2D basis change, calculating the position of the particle if they hit the wall */
void collide_circular_wall(Particle *p, Arena *a)
{
  	Real initial_vel[2];
	Real n[2], t[2], n_mag;
  	Real v_nt[2];
  	Real m, vi, vf;
	
	/* storing the initial velocities of the particle first before we change it */
	initial_vel[0] = p->velocity[0];
	initial_vel[1] = p->velocity[1];

  	/* Normal vector n between centres. */
  	n[0] = p->position[0];
  	n[1] = p->position[1];
  
  	/* Normalise n. */
  	n_mag = sqrt(n[0] * n[0] + n[1] * n[1]);
  	n[0] /= n_mag;
  	n[1] /= n_mag;

  	/* Tangent vector t. */
  	t[0] = -n[1];
  	t[1] = n[0];

  	/* Change basis for velocities from standard basis to nt basis. */
  	v_nt[0] = n[0] * p->velocity[0] + n[1] * p->velocity[1];
  	v_nt[1] = t[0] * p->velocity[0] + t[1] * p->velocity[1];

  	/* Use 1D equations to calculate final velocities in n direction. */
  	m = p->mass;
  	vi = v_nt[0];
  	vf = -vi;

  	/* Update the 2D velocity. Force is in n direction, so in nt basis,
  	 * velocity change is in n direction only, no change in t direction.
   	 */
  	v_nt[0] = vf;

	Real length_outside_arena = get_distance(p->position[0], p->position[1], 0.0, 0.0) + p->radius - ARENA_RADIUS;
	/* moving backward the length which is outside the arena */
	p->position[0] += p->position[0] < 0 ? length_outside_arena : -length_outside_arena;
	p->position[1] += p->position[1] < 0 ? length_outside_arena : -length_outside_arena;
  	
	/*	the code below is used for debugging
	 	
		length_outside_arena = get_distance(p->position[0], p->position[1], 0.0, 0.0) + p->radius - ARENA_RADIUS;
		if(length_outside_arena > 0.0)
		{
			printf("length_outside_arena = %f\n", length_outside_arena);
		}
	*/

	/* Change back to standard basis. */
  	p->velocity[0] = n[0] * v_nt[0] + t[0] * v_nt[1];
  	p->velocity[1] = n[1] * v_nt[0] + t[1] * v_nt[1];
	
	/* taking a step ahead */
	p->position[0] += length_outside_arena *  p->velocity[0];
	p->position[1] += length_outside_arena *  p->velocity[1];

	/* conserving momentum of the environment */
	a->momentum[0] += p->mass * (initial_vel[0] - p->velocity[0]);
	a->momentum[1] += p->mass * (initial_vel[1] - p->velocity[1]);

}

/* used for determining a square arena or circular arena */
void collide_particle_wall(Particle *p, Arena *a, bool square_arena)
{
	if(square_arena)
		collide_square_wall(p, a);
	else
		if(collided_with_circular_arena(p))
			collide_circular_wall(p, a);
}

void collideParticlesWall(Particle *particle, int number_of_particles, Arena *arena, bool square_arena)
{
  	for (int i = 0; i < number_of_particles; i++) 
	{
      	//printf("%d %f %f\n", i, particle[i].position[0], particle[i].position[1]);
		collide_particle_wall(&particle[i], arena, square_arena);
  	}
}

bool collisionDetectionParticles(Particle *p1, Particle *p2)
{
	Real sum_radii, sum_radii_sq, n[2], n_mag_sq;

  	sum_radii = p1->radius + p2->radius;
  	sum_radii_sq = sum_radii * sum_radii;
  	n[0] = p2->position[0] - p1->position[0];
  	n[1] = p2->position[1] - p1->position[1];
  	n_mag_sq = n[0] * n[0] + n[1] * n[1];
  	if (n_mag_sq <= sum_radii_sq) 
    	return true;
  	else 
    	return false;
}

/* returns 1 if a dead particle is found so that the array can be reinitialised */
int mote_absorbing_mote(Particle *big, Particle *small, Real new_mass)
{
	new_mass = (small->mass - new_mass) < 0 ? small->mass : new_mass;
	small->alive = (small->mass - new_mass) > 0.0 ? true : false;
	
	if(small->alive)
	{
		big->velocity[0] = ( (big->mass * big->velocity[0]) + (new_mass * small->velocity[0]) ) / (big->mass + new_mass);
		big->velocity[1] = ( (big->mass * big->velocity[1]) + (new_mass * small->velocity[1]) ) / (big->mass + new_mass);
	}
	else
	{
		big->velocity[0] = ( (big->mass * big->velocity[0]) + (small->mass * small->velocity[0]) ) / (big->mass + new_mass);
		big->velocity[1] = ( (big->mass * big->velocity[1]) + (small->mass * small->velocity[1]) ) / (big->mass + new_mass);
		small->velocity[0] = small->velocity[1] = 0.0;
	}

	/* bigger particle absorbs the smaller one, so the mass and the radius of both the particles change */
	big->mass += new_mass;
    big->radius = sqrt(big->mass);
	small->mass -= new_mass;
    small->radius = sqrt(small->mass);
	return small->alive ? 0 : 1;
}

/* since there is a collision some mass from the small mote will be absorbed by the bigger mote */
int absorb_mote(Particle *big, Particle *small)
{
	Real sum_radii = (big->radius + small->radius);
	Real dist_between_two_particles = get_distance(big->position[0], big->position[1], small->position[0], small->position[1]);
	Real overlapping_length = sum_radii -  dist_between_two_particles;
	Real new_mass;
	
	/* calculating amount of mass that can be sucked in */
	if(dist_between_two_particles <= big->radius)
		new_mass = small->mass;
	else
	{
		Real final_radius_of_small_particle = small->radius - overlapping_length;
		Real final_mass_of_small_particle = final_radius_of_small_particle * final_radius_of_small_particle;
		new_mass = small->mass - final_mass_of_small_particle;
	}
	return mote_absorbing_mote(big, small, new_mass);
}

/* 
 * brute force algorithm to check if any particle collided with any other
 * return total number of dead particles 
 */
int collideParticlesBruteForce(Particle *particle, int number_of_particles, Arena *arena, Real h, bool square_arena)
{
  int i, j;
  int number_of_dead_particles = 0;
  
  for (i = 0; i < number_of_particles - 1; i++) 
  {
    	for (j = i + 1; j < number_of_particles; j++) 
		{

			if(!particle[i].alive || !particle[j].alive)
				continue;

      		if (collisionDetectionParticles(&particle[i], &particle[j])) 
			{
				if(particle[i].mass > particle[j].mass)
					number_of_dead_particles += absorb_mote(&particle[i], &particle[j]);
				else
					number_of_dead_particles += absorb_mote(&particle[j], &particle[i]);
				
        		/* Check walls. */
				collide_particle_wall(&particle[i], arena, square_arena);
				collide_particle_wall(&particle[j], arena, square_arena);
      		}
    	}
  	}
	return number_of_dead_particles;
}

void calcGridIndex(Particle *p, Arena *a, 
                   Real *gridCellSize, int *gridNumCells,
                   int *index)
{
  	index[0] = (int)((p->position[0] - a->min[0]) / gridCellSize[0]);
  	index[1] = (int)((p->position[1] - a->min[1]) / gridCellSize[1]);
	
	/*	code below used for debugging 
	 
    	if (index[0] < 0 || index[0] > gridNumCells[0] - 1)
      		panic("gridIndex: index out of range\n");
    	if (index[1] < 0 || index[1] > gridNumCells[1] - 1)
      		panic("gridIndex: index out of range\n");
  	*/
}

int collideParticlesUniformGrid(Particle *particle, int number_of_particles, Arena *arena, Real h, bool square_arena)
{
  	Real gridCellSize[2];
  	int **gridCellParticleCount, **gridCellParticleListEnd, *gridCellParticleList;
  	int gridNumCells[2], gridSize, gridIndex[2], gridCellParticleListStart;
  	int gridIndexMin[2], gridIndexMax[2];
  	int i, j, k, s, t, p1, p2, total;
  	int number_of_dead_particles = 0;
  
  	/* Work out grid dimensions and allocate. */
  	gridNumCells[0] = (int)(sqrt(number_of_particles) + 1);
  	gridNumCells[1] = (int)(sqrt(number_of_particles) + 1);
  	gridCellSize[0] = (arena->max[0] - arena->min[0]) / gridNumCells[0];
  	gridCellSize[1] = (arena->max[1] - arena->min[1]) / gridNumCells[1];
  	gridSize = gridNumCells[0] * gridNumCells[1];

  	/* Assumption. */
  	for (i = 0; i < number_of_particles; i++)
    	if (particle[i].radius * 2.0 > gridCellSize[0] ||
        		particle[i].radius * 2.0 > gridCellSize[1])
		{
			printf("collideParticlesUniformGrid: particle diameter > cellSize\t Hence using Brute Force\n");
      		return -1;
		}
 	
  	/* Allocate arrays. */
  	gridCellParticleCount = (int **)malloc(gridNumCells[0] * sizeof(int *));
  	if (gridCellParticleCount == 0)
    	panic("collideParticlesUniformGrid: malloc failed\n");
  	gridCellParticleListEnd = (int **)malloc(gridNumCells[0] * sizeof(int *));
  	if (gridCellParticleListEnd == 0)
    	panic("collideParticlesUniformGrid: malloc failed\n");
  	for (i = 0; i < gridNumCells[0]; i++) 
	{
    	gridCellParticleCount[i] = (int *)malloc(gridNumCells[1] * sizeof(int));
    	if (gridCellParticleCount[i] == 0)
      		panic("collideParticlesUniformGrid: malloc failed\n");
    	gridCellParticleListEnd[i] = (int *)malloc(gridNumCells[1] * sizeof(int));
    	if (gridCellParticleListEnd[i] == 0)
      		panic("collideParticlesUniformGrid: malloc failed\n");
  	}
  	gridCellParticleList = (int *)malloc(number_of_particles * sizeof(int));

  	/* Initialise grid particle count. */
  	for (i = 0; i < gridNumCells[0]; i++)
 		for (j = 0; j < gridNumCells[1]; j++)
    		gridCellParticleCount[i][j] = 0;

  	/* Cell counts. */
  	for (i = 0; i < number_of_particles; i++) 
	{
     	calcGridIndex(&particle[i], arena, gridCellSize, gridNumCells, gridIndex);
     	gridCellParticleCount[gridIndex[0]][gridIndex[1]] += 1;
  	}

	/*
    	printf("collideParticlesUniformGrid: gridCellParticleCount\n");
    	for (i = 0; i < gridNumCells[0]; i++)
      		for (j = 0; j < gridNumCells[1]; j++)
        			printf("%d %d %d\n", i, j, gridCellParticleCount[i][j]);
  	*/

  	/* Work out end of cell lists by accumulating cell counts. */
  	for (i = 0; i < gridNumCells[0]; i++)
    	for (j = 0; j < gridNumCells[1]; j++)
      		gridCellParticleListEnd[i][j] = 0;
  
  	total = 0;
  	for (i = 0; i < gridNumCells[0]; i++)
	{
    	for (j = 0; j < gridNumCells[1]; j++) 
		{
      		total = total + gridCellParticleCount[i][j];
      		gridCellParticleListEnd[i][j] = total - 1;
    	}
	}

	/*
    	printf("collideParticlesUniformGrid: gridCellParticleListEnd\n");
    	for (i = 0; i < gridNumCells[0]; i++)
      		for (j = 0; j < gridNumCells[1]; j++)
        		printf("%d %d %d\n", i, j, gridCellParticleListEnd[i][j]);
  	*/

  	/* Build particle lists. */
  	for (i = 0; i < gridNumCells[0]; i++)
    	for (j = 0; j < gridNumCells[1]; j++)
      		gridCellParticleCount[i][j] = 0;

  	for (i = 0; i < number_of_particles; i++) 
	{
    	calcGridIndex(&particle[i], arena, gridCellSize, gridNumCells, gridIndex);
    	gridCellParticleList[gridCellParticleListEnd[gridIndex[0]][gridIndex[1]] - 
      	gridCellParticleCount[gridIndex[0]][gridIndex[1]]] = i;
    	gridCellParticleCount[gridIndex[0]][gridIndex[1]] += 1;
  	}

	/*	code below used for debugging

    	printf("collideParticlesUniformGrid: gridCellParticleList\n");
    	for (i = 0; i < gridNumCells[0]; i++) 
		{
      		for (j = 0; j < gridNumCells[1]; j++) 
			{
        		gridCellParticleListStart = 
          			gridCellParticleListEnd[i][j] - gridCellParticleCount[i][j] + 1;
        		printf("particle list %d %d\n", i, j);
        		for (k = gridCellParticleListStart; k < gridCellParticleListEnd[i][j]; k++)
          			printf("%d\n", gridCellParticleList[k]);
        		printf("\n");
      		}
    	}
  	*/

  	/* Collision detection. */
  	for (i = 0; i < number_of_particles; i++) 
	{
    	calcGridIndex(&particle[i], arena, gridCellSize, gridNumCells, gridIndex);

    	/* Grid index bounds for this particle. */
    	gridIndexMin[0] = gridIndex[0] - 1;
    	if (gridIndexMin[0] < 0) 
      		gridIndexMin[0] = 0;
   		gridIndexMin[1] = gridIndex[1] - 1;
    	if (gridIndexMin[1] < 0) 
      		gridIndexMin[1] = 0;
    	gridIndexMax[0] = gridIndex[0] + 1;
    	if (gridIndexMax[0] > gridNumCells[0] - 1) 
      		gridIndexMax[0] = gridNumCells[0] - 1;
    	gridIndexMax[1] = gridIndex[1] + 1;
    	if (gridIndexMax[1] > gridNumCells[1] - 1) 
      		gridIndexMax[1] = gridNumCells[1] - 1;

    	p1 = i;

	    for (s = gridIndexMin[0]; s <= gridIndexMax[0]; s++) 
		{
    		for (t = gridIndexMin[1]; t <= gridIndexMax[1]; t++) 
			{
        		gridCellParticleListStart = 
          			gridCellParticleListEnd[s][t] - gridCellParticleCount[s][t] + 1;
        		for (j = gridCellParticleListStart; j <= gridCellParticleListEnd[s][t]; j++) 
				{
          			p2 = gridCellParticleList[j];

	  				/* Don't test particle against itself. */
          			if (p2 == p1)
            			continue;

          			/* Only test pairs once. */
	  				if (p2 < p1)
	    				continue;

	    			//printf("collideParticlesUniformGrid: testing %d %d\n", p1, p2);

	  				if (collisionDetectionParticles(&particle[p1], &particle[p2])) 
					{
              			//printf("collision: %d %d\n", p1, p2);

						if(particle[p1].mass > particle[p2].mass)
							number_of_dead_particles += absorb_mote(&particle[p1], &particle[p2]);
						else
							number_of_dead_particles += absorb_mote(&particle[p2], &particle[p1]);

            			/* Check walls. */
						collide_particle_wall(&particle[p1], arena, square_arena);
						collide_particle_wall(&particle[p2], arena, square_arena);
	  				}
				}
      		}
    	}
	}
 
  	/* Free arrays. */
  	for (i = 0; i < gridNumCells[0]; i++) 
	{
    	free(gridCellParticleCount[i]);
    	free(gridCellParticleListEnd[i]);
	}
  	free(gridCellParticleCount);
  	free(gridCellParticleListEnd);
  	free(gridCellParticleList);

	return number_of_dead_particles;
}

