/*
 *  Copyright (C) 2010 by Tim Massingham
 *  tim.massingham@ebi.ac.uk
 *
 *  This file is part of the AYB base-calling software.
 *
 *  AYB is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AYB is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with AYB.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <inttypes.h>
#include "utility.h"
#include "coordinates.h"

typedef struct {
    uint16_t lane,tile,x,y;
    bool okay;
} LINE;

COORD new_COORD(const uint32_t ncluster){
    COORD coord = calloc(1,sizeof(*coord));
    if(NULL==coord){ return NULL;}
    coord->ncluster = ncluster;
    coord->x = calloc(ncluster,sizeof(uint16_t));
    coord->y = calloc(ncluster,sizeof(uint16_t));
    if(NULL==coord->x || NULL==coord->y){
        goto cleanup;
    }
    
    return coord;
    
cleanup:
    xfree(coord->x);
    xfree(coord->y);
    xfree(coord);
    
    return NULL;
}


void * free_COORD(COORD coord){
    if(NULL==coord) return NULL;
    xfree(coord->x);
    xfree(coord->y);
    xfree(coord);
    return NULL;
}


COORD copy_COORD( const COORD coord){
    if(NULL==coord){ return NULL;}

    COORD coord_copy = new_COORD(coord->ncluster);
    if(NULL==coord_copy){ return NULL;}

    // Copy values, ncluster already copied via new_COORD
    coord_copy->lane = coord->lane;
    coord_copy->tile = coord->tile;
    for ( uint32_t cl=0 ; cl<coord->ncluster ; cl++){
        coord_copy->x[cl] = coord->x[cl];
        coord_copy->y[cl] = coord->y[cl];
    }

    return coord_copy;
}

void show_COORD( FILE * fp, const COORD coord){
    if(NULL==fp || NULL==coord){ return;}
    fprintf(fp,"%" SCNu32 " clusters\n",coord->ncluster);
    fprintf(fp,"Lane %" SCNu16 "   Tile %" SCNu16 "\n", coord->lane, coord->tile);
    for ( uint32_t i=0 ; i<5 ; i++){
        fprintf(fp,"%4" SCNu32 ":\t%" SCNu16 "\t%" SCNu16 "\n", i+1, coord->x[i], coord->y[i]);
    }
}


LINE read_coordinate_line( FILE * fp){
    LINE l;
    int ret = fscanf(fp,"%" SCNu16 "%" SCNu16 "%" SCNu16 "%" SCNu16, &l.lane,&l.tile,&l.x,&l.y);
    l.okay = (ret==4) ? true : false;

    return l;
}


COORD read_coordinates(const CSTRING str, const uint32_t ncluster){
    if(NULL==str){ return NULL;}

    // Open file and return if fail
    FILE * fp = fopen(str,"r");
    if(NULL==fp){ return NULL;}

    // Memory for coordinates (also LINE store first line)
    COORD coord = new_COORD(ncluster);
    LINE l1;
    for ( uint32_t cl=0 ; cl<ncluster ; cl++){
        LINE l = read_coordinate_line(fp);
        if(!l.okay){ errx(EXIT_FAILURE,"Failed to parse coordinate file, line %u",cl);}
        // Set initial to check other lines' lane and tile against
        if(0==cl){ l1 = l;}
        // Basic checks of lane and tile.
        if(l1.lane!=l.lane){
            errx(EXIT_FAILURE,"Failed to parse coordinate file.\n"
                              "Lane number of line %u does not match, got %u expecteding %u",cl,l.lane,l1.lane);
        }
        if(l1.tile!=l.tile){
            errx(EXIT_FAILURE,"Failed to parse coordinate file.\n"
                              "Tile number of line %u does not match, got %u expecteding %u",cl,l.tile,l1.tile);
        }
        // Set coordinates
        coord->x[cl] = l.x;
        coord->y[cl] = l.y;
    }
    coord->lane = l1.lane;
    coord->tile = l1.tile;
    
    fclose(fp);
    return coord;
}

