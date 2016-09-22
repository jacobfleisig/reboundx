/**
 * @file    input.c
 * @brief   Input functions.
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
#include <string.h>
#include "rebound.h"
#include "reboundx.h"
#include "core.h"

void rebx_load_param(void* const object, struct rebx_binary_field* field, struct rebx_extras* rebx, FILE* inf, enum reb_input_binary_messages* warnings){
    struct rebx_binary_param_metadata metadata;
    printf("before metadata %lu\n", ftell(inf));
    fread(&metadata, sizeof(metadata), 1, inf);
    printf("%d\t%d\t%lu\n", metadata.type, metadata.ndim, metadata.namelength);
    
    char* name = malloc(metadata.namelength);
    printf("before name %lu\n", ftell(inf));
    fread(name, sizeof(*name), metadata.namelength, inf);
    printf("%s\t%lu\t%lu\n", name, strlen(name), metadata.namelength);
    int* shape = malloc(metadata.ndim*sizeof(*shape));
    printf("ndim %d\n", metadata.ndim);
    printf("before shape %lu\n", ftell(inf));
    fread(shape, sizeof(*shape), metadata.ndim, inf);
    printf("after shape %lu\n", ftell(inf));
    struct rebx_param* param = rebx_add_param_node(object, name, metadata.type, metadata.ndim, shape);
    free(name);
    free(shape);
    fread(param->contents, rebx_sizeof(param->param_type), param->size, inf);
    printf("after contents %lu\n", ftell(inf));
}
    
static size_t rebx_load_effect(struct rebx_binary_field* field, struct rebx_extras* rebx, FILE* inf, enum reb_input_binary_messages* warnings){
    size_t namelength;
    //printf("aL %lu\n", ftell(inf));
    fread(&namelength, sizeof(namelength), 1, inf);
    char* name = malloc(namelength);
    fread(name, sizeof(*name), namelength, inf);
    printf("%s\t%lu\t%lu\n", name, strlen(name), namelength);
    struct rebx_effect* effect = rebx_add(rebx, name);
    free(name);
    /*fseek(inf, field.size, SEEK_CUR);
     size_t elements_read = fread(&field, sizeof(field), 1, inf);
     printf("%lu\n", field.size);
     */
    /*double a,b;
    if(!fread(&a, sizeof(a), 1, inf)){
        *warnings |= REB_INPUT_BINARY_WARNING_INCOMPATIBLE_FORMAT;
        return 0;
    }
    if(!fread(&b, sizeof(b), 1, inf)){
        *warnings |= REB_INPUT_BINARY_WARNING_INCOMPATIBLE_FORMAT;
        return 0;
    }
    
    printf("%f\t%f\n", a, b);*/
    //printf("bef next read %lu\n", ftell(inf));
    size_t field_successfully_read = fread(field, sizeof(*field), 1, inf);
    while(field_successfully_read && field->type == REBX_BINARY_FIELD_TYPE_PARAM){
        rebx_load_param(effect, field, rebx, inf, warnings);
        field_successfully_read = fread(field, sizeof(*field), 1, inf);
    }
    //printf("%lu", sizeof(field));
    //printf("after read %lu\n", ftell(inf));
    return field_successfully_read;
}

static void rebx_create_extras_from_binary_with_messages(struct rebx_extras* rebx, const char* const filename, enum reb_input_binary_messages* warnings){
    FILE* inf = fopen(filename,"rb");
    if (!inf){
        *warnings |= REB_INPUT_BINARY_ERROR_NOFILE;
        return;
    }
    
    long objects = 0;
    // Input header.
    const char str[] = "REBOUNDx Binary File. Version: ";
    char readbuf[65], curvbuf[65];
    sprintf(curvbuf,"%s%s",str,rebx_version_str);
    for(size_t j=strlen(curvbuf);j<63;j++){
        curvbuf[j] = ' ';
    }
    curvbuf[63] = '\0';
    objects += fread(readbuf,sizeof(*str),64,inf);
    if(strcmp(readbuf,curvbuf)!=0){
        *warnings |= REB_INPUT_BINARY_WARNING_VERSION;
    }
    printf("%s\n", readbuf);
    struct rebx_binary_field field;
    size_t field_successfully_read = fread(&field, sizeof(field), 1, inf); // fread will return 1 if it successfully read 1 element
    printf("%lu\n", field.size);
    while(field_successfully_read) {
        switch (field.type){
            case REBX_BINARY_FIELD_TYPE_EFFECT:
            {
                field_successfully_read = rebx_load_effect(&field, rebx, inf, warnings);
                break;
            }
            case REBX_BINARY_FIELD_TYPE_PARTICLE:
            {
                break;
            }
            case REBX_BINARY_FIELD_TYPE_PARAM:
            {
                break;
            }
            default:
                *warnings |= REB_INPUT_BINARY_WARNING_FIELD_UNKOWN;
                fseek(inf,field.size,SEEK_CUR); // type unrecognized (diff version?) try skipping
                break;
        }
    }
    
    fclose(inf);
    return;
}

struct rebx_extras* rebx_create_extras_from_binary(struct reb_simulation* sim, const char* const filename){
    enum reb_input_binary_messages warnings = REB_INPUT_BINARY_WARNING_NONE;
    struct rebx_extras* rebx = rebx_init(sim);
    rebx_create_extras_from_binary_with_messages(rebx, filename, &warnings);
    
    if (warnings & REB_INPUT_BINARY_WARNING_VERSION){
        reb_warning(sim,"REBOUNDx: Binary file was saved with a different version of REBOUNDx. Binary format might have changed and corrupted the loading. Check that effects and parameters are loaded as expected.");
    }
    if (warnings & REB_INPUT_BINARY_WARNING_FIELD_UNKOWN){
        reb_warning(sim,"REBOUNDx: Unknown field found in binary file. Any unknown fields not loaded.  This can happen if the binary was created with a later version of REBOUNDx than the one used to read it.");
    }
    if (warnings & REB_INPUT_BINARY_ERROR_NOFILE){
        reb_error(sim,"REBOUNDx: Cannot read binary file. Check filename and file contents.");
        rebx_free(rebx);
        rebx = NULL;
    }
    return rebx;
}