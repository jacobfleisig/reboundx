/**
 * @file    modify_orbits_direct.c
 * @brief   Update orbital element every 100 orbital periods as described by interactions w/ Neptune
 * @author  Dan Tamayo <tamayo.daniel@gmail.com>
 * 
 * @section     LICENSE
 * Copyright (c) 2015 Dan Tamayo, Hanno Rein
 *
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

static struct reb_particle rebx_calculate_modify_orbits_direct(struct reb_simulation* const sim, struct rebx_effect* const effect, struct reb_particle* p, struct reb_particle* primary, const double dt){
	int err=0;
	struct reb_orbit o = reb_tools_particle_to_orbit_err(sim->G, *p, *primary, &err);
	double time = sim->t;

	if(err){        // mass of primary was 0 or p = primary.  Return same particle without doing anything.
		return *p;
	}
	
	// change in perihelion
	double q0 = o.a*(1-o.e);
	double d_a = 0.02*(1.0/o.a);

	// every 100 orbital periods, give a particle a kick in semi-major axis	
	if (time > pow(o.a,3.0/2.0)*2.0*M_PI){
		if (reb_output_check(sim,100.0*2.0*M_PI*pow(o.a,3.0/2.0))){
			if ( q0 <= 0.35){    		
				o.a += d_a;
			}
		}
	}
	return reb_tools_orbit_to_particle(sim->G, *primary, p->m, o.a, o.e, o.inc, o.Omega, o.omega, o.f);
}

void rebx_modify_orbits_direct(struct reb_simulation* const sim, struct rebx_effect* const effect, const double dt, enum rebx_timing timing){
	const int* const ptr = rebx_get_param_check(effect, "coordinates", REBX_TYPE_INT);
	enum REBX_COORDINATES coordinates = REBX_COORDINATES_JACOBI;
	if (ptr != NULL){
		coordinates = *ptr;
	}
	const int back_reactions_inclusive = 1;
	const char* reference_name = "primary";
	rebxtools_com_ptm(sim, effect, coordinates, back_reactions_inclusive, reference_name, rebx_calculate_modify_orbits_direct, dt);
}
