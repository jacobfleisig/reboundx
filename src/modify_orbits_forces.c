/**
 * @file    modify_orbits_forces.c
 * @brief   Update orbital elements with prescribed timescales using forces.
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
#include "modify_orbits_forces.h"
#include "rebound.h"
#include "reboundx.h"
#include "rebxtools.h"

#define TWOPI 6.2831853071795862

struct rebx_params_modify_orbits_forces* rebx_add_modify_orbits_forces(struct rebx_extras* rebx){
	struct rebx_params_modify_orbits_forces* params = malloc(sizeof(*params));
    params->coordinates = JACOBI;
    int force_is_velocity_dependent = 1;
    rebx_add_force(rebx, params, "modify_orbits_forces", rebx_modify_orbits_forces, force_is_velocity_dependent);
    return params;
}

reb_vec3d rebx_calculate_modify_orbits_forces(void* params, struct reb_particle* p, struct reb_particle* source){
    double tau_a = rebx_get_param_double(p, "tau_a");
    double tau_e = rebx_get_param_double(p, "tau_e");
    double tau_inc = rebx_get_param_double(p, "tau_inc");

    if(isnan(tau_a)){
        rebx_set_param_double(p, "tau_a", INFINITY);
        tau_a = INFINITY;
    }
    if(isnan(tau_e)){
        rebx_set_param_double(p, "tau_e", INFINITY);
        tau_e = INFINITY;
    }
    if(isnan(tau_inc)){
        rebx_set_param_double(p, "tau_inc", INFINITY);
        tau_inc = INFINITY;
    }

    const double dvx = p->vx - source.vx;
    const double dvy = p->vy - source.vy;
    const double dvz = p->vz - source.vz;

    reb_vec3d a = {0};

    a.x =  dvx/(2.*tau_a);
    a.y =  dvy/(2.*tau_a);
    a.z =  dvz/(2.*tau_a);

    if (tau_e < INFINITY || tau_inc < INFINITY){   // need h and e vectors for both types
        const double dx = p->x-source.x;
        const double dy = p->y-source.y;
        const double dz = p->z-source.z;
        const double r = sqrt ( dx*dx + dy*dy + dz*dz );
        const double vr = (dx*dvx + dy*dvy + dz*dvz)/r;
        const double prefac = 2*vr/r;
        a.x += prefac*dx/tau_e;
        a.y += prefac*dy/tau_e;
        a.z += prefac*dz/tau_e + 2.*dvz/tau_inc;
    }
    return a;
}

void rebx_modify_orbits_forces(struct reb_simulation* const sim, struct rebx_effect* const effect){
    const struct rebx_params_modify_orbits_forces* const params = effect->paramsPtr; // if coordinates hardcoded can remove this and set NULL below
    const int first_index = 1;
    const int last_index = -1;
    const int back_reactions_inclusive = 1;

    rebx_ghost_effect(sim, params, first_index, last_index, REBX_JACOBI, back_reactions_inclusive, rebx_calculate_modify_orbits_forces);
}

void rebx_modify_orbits_around_particle_forces(struct reb_simulation* const sim, struct rebx_effect* const effect){
    const int N_real = sim->N - sim_>N_var;
    struct reb_vec3d a;
    const int back_reactions_inclusive = 1;

    for(int i=0; i < N_real; i++){
        struct reb_particle* central_body = &sim->particles[i];
        int modify_orbits = rebx_get_param_int(central_body, "modify_orbits");
        if(modify_orbits){
            int first_index = rebx_get_param_int(central_body, "modify_orbits_first_index");
            int last_index= rebx_get_param_int(central_body, "modify_orbits_last_index");
            rebx_particle_effect(sim, NULL, first_index, last_index, central_body, i, back_reactions_inclusive, rebx_calculate_modify_orbits_forces);
        }
    }
}
