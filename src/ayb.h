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
 
#ifndef _AYB_H
#define _AYB_H

#include "cif.h"
#include "matrix.h"
#include "nuc.h"
#include "coordinates.h"

#define X(A) A ## int16_t
#include "array.def"
#undef X

#define CALIBRATION_FILE "tables/76cyL2.h"


typedef struct {
    uint32_t ncycle,ncluster;
    ARRAY(int16_t) intensities;
    ARRAY(NUC) bases;
    ARRAY(PHREDCHAR) quals;
    MAT M,P,N,lambda,lss;
    MAT At;
    MAT we, cycle_var;
    // Information about filtering
    bool filtered;
    bool * passed_filter;
    // coordinates
    COORD coordinates;
    uint16_t index, readnum;
} * AYB;

/* Basic functions */
AYB new_AYB(const uint32_t ncycle, const uint32_t ncluster);
void free_AYB(AYB ayb);
AYB copy_AYB(const AYB ayb);

void show_AYB(XFILE * fp, const AYB ayb);

/* More complex initialisation */
AYB initialise_AYB(const CIFDATA cif);

real_t estimate_MPC( AYB ayb );
real_t estimate_Bases(AYB ayb, const bool lastIt);

#endif /* _AYB_H */

