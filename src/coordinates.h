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

#ifndef _COORDINATE_H
#define _COORDINATE_H

#include "utility.h"

typedef struct {
    uint16_t lane,tile;
    uint32_t ncluster;
    uint16_t *x, *y;
} * COORD;

COORD new_COORD(const uint32_t ncluster);
void * free_COORD(COORD coord);
COORD copy_COORD( const COORD coord);
void show_COORD( FILE * fp, const COORD coord);

COORD read_coordinates(const CSTRING str, const uint32_t ncluster);

#endif /* _COORDINATE_H */
