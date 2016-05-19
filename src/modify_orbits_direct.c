/**
 * @file    modify_orbits_direct.c
 * @brief   Update orbital with prescribed timescales by directly changing orbital elements after each timestep.
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
#include "modify_orbits_direct.h"
#include "rebound.h"
#include "reboundx.h"

#define TWOPI 6.2831853071795862

struct rebx_params_modify_orbits_direct* rebx_add_modify_orbits_direct(struct rebx_extras* rebx){
	struct rebx_params_modify_orbits_direct* params = malloc(sizeof(*params));
	params->p = 0;
    params->coordinates = JACOBI;
    rebx_add_post_timestep_modification(rebx, params, "modify_orbits_direct", rebx_modify_orbits_direct);
    return params;
}

void rebx_modify_orbits_direct(struct reb_simulation* const sim, struct rebx_effect* const effect){
    const struct rebx_params_modify_orbits_direct* const params = effect->paramsPtr;
    const int N_real = sim->N - sim->N_var;

    struct reb_particle com = reb_get_com(sim);
    double Mtot = com.m;
    for(int i=N_real-1;i>0;--i){
        struct reb_particle p = sim->particles[i];
        struct reb_particle* pptr = &sim->particles[i];
        rebxtools_update_com_without_particle(&com, pptr);
        /*if(params->coordinates == JACOBI){
            rebxtools_update_com_without_particle(&com, p);
        }*/
        struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, p, com);
        const double dt = sim->dt_last_done;
        //TODO step through linked list to populate timescales in one pass, and don't calculate orbit above if none are set 
        const double tau_a = rebx_get_param_double(pptr,"tau_a");
        const double tau_e = rebx_get_param_double(pptr,"tau_e");
        const double tau_inc = rebx_get_param_double(pptr,"tau_inc");
        const double tau_omega = rebx_get_param_double(pptr,"tau_omega");
        const double tau_Omega = rebx_get_param_double(pptr,"tau_Omega");
        
        const double a0 = o.a;
        const double e0 = o.e;
        const double inc0 = o.inc;

        if(!isnan(tau_a)){
            o.a += a0*(dt/tau_a);
        }
        if(!isnan(tau_e)){
            o.e += e0*(dt/tau_e);
        }
        if(!isnan(tau_inc)){
            o.inc += inc0*(dt/tau_inc);
        }
        if(!isnan(tau_omega)){
            o.omega += TWOPI*dt/tau_omega;
        }
        if(!isnan(tau_Omega)){
            o.Omega += TWOPI*dt/tau_Omega;
        }
        if(!isnan(tau_e)){
            o.a += 2.*a0*e0*e0*params->p*dt/tau_e; // Coupling term between e and a
        }

        rebxtools_orbit2p(sim->G, pptr, &com, &o);

        // TODO write function that takes particle dp, sim, and index to do the below.  Write subtract_particles function in rebound
        if(p.m > 0.){   // change in p will have shifted COM. Add shift to all particles with index j <= i to keep COM fixed. 
            double massratio = p.m/(com.m+p.m);
            struct reb_particle dp = {0};
            dp.x = massratio*(p.x - pptr->x);
            dp.y = massratio*(p.y - pptr->y);
            dp.z = massratio*(p.z - pptr->z);
            dp.vx = massratio*(p.vx - pptr->vx);
            dp.vy = massratio*(p.vy - pptr->vy);
            dp.vz = massratio*(p.vz - pptr->vz);
            for (int j=0;j<=i;j++){
                sim->particles[j].x += dp.x;
                sim->particles[j].y += dp.y;
                sim->particles[j].z += dp.z;
                sim->particles[j].vx += dp.vx;
                sim->particles[j].vy += dp.vy;
                sim->particles[j].vz += dp.vz;
            }
        }
    }
}


