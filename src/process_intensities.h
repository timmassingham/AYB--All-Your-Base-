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

#ifndef _PROCESS_INTENSITIES_H
#define _PROCESS_INTENSITIES_H

#include <stdint.h>
#include "matrix.h"

MAT process_intensities(const int16_t * intensities, const MAT Minv_t, const MAT Pinv_t, const MAT N, MAT p);
MAT expected_intensities( const real_t lambda, const NUC * bases, const MAT M, const MAT P, const MAT N, MAT e);

#endif /* _PROCESS_INTENSITIES_H */

