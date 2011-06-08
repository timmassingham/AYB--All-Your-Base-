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

#ifndef _ESTIMATE_LAMBDA_H
#define _ESTIMATE_LAMBDA_H

#include "nuc.h"
#include "matrix.h"
#include "utility.h"

real_t estimate_lambdaWLS( const MAT p, const NUC * base, const real_t oldlambda, const real_t * v);
real_t estimate_lambdaGWLS( const MAT p, const NUC * base, const real_t oldlambda, const real_t * v, const MAT * V);
real_t estimate_lambdaOLS( const MAT p, const NUC * base);
real_t estimate_lambda_A ( const int16_t * intensity, const MAT N, const MAT lamN, const MAT At, const NUC * base, const int ncycle);


#endif /* _ESTIMATE_LAMBDA_H */
