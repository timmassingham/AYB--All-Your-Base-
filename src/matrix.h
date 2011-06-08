/*
 *  Copyright (C) 2010 by Tim Massingham, European Bioinformatics Institute
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


#ifndef _MATRIX_H
#define _MATRIX_H


#include <stdint.h>
#include <stdbool.h>
#include "xio.h"
#include "utility.h"

/*  Specification did not specify what type of input should be taken. Use a
 * typedef so it can be easily changed, although care must be taken with the
 * input and output routines (printf etc) to make sure types match.
 */ 

struct _matrix_str {
    int    nrow, ncol;
    real_t    * x;
    };

// Make future abstraction easier
typedef struct _matrix_str * MAT;

// create/free
MAT new_MAT( const int nrow, const int ncol );
void free_MAT( MAT mat );
void shallow_free_MAT( MAT mat);
MAT copy_MAT( const MAT mat);
MAT copyinto_MAT( MAT matout, const MAT matin);
MAT new_MAT_from_array( const uint32_t nrow, const uint32_t ncol, const real_t * x);
MAT coerce_MAT_from_array(const uint32_t nrow, const uint32_t ncol, real_t * x);

// Input, output
void show_MAT( XFILE * fp, const MAT mat, const uint32_t mrow, const uint32_t mcol);

// Identities
bool is_square(const MAT mat);

// Special matrices
MAT identity_MAT( const int nrow);
MAT set_MAT( MAT mat, const real_t x);

// Operations
MAT reshape_MAT( MAT mat, const int nrow);
MAT trim_MAT( MAT mat, const int mrow, const int mcol, const bool forwards);
MAT * block_diagonal_MAT( const MAT mat, const int n);
MAT scale_MAT(MAT mat, const real_t f);
MAT transpose_inplace( MAT mat);
MAT transpose( const MAT mat);
MAT invert(const MAT mat);
MAT invertSym(const MAT mat);

real_t xMy( const real_t * x, const MAT M, const real_t * y);
real_t normalise_MAT(MAT mat, const real_t delta_diag);
MAT cholesky( MAT mat);

struct structLU { MAT mat; int * piv; };
struct structLU LUdecomposition( const MAT mat);
#endif

