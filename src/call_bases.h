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
 
#ifndef _CALL_BASES_H
#define _CALL_BASES_H

struct basequal { NUC base; PHREDCHAR qual;};

NUC call_base_simple( const real_t * restrict p);
struct basequal call_base( const real_t * restrict p, const real_t lambda, const MAT omega);

MAT * calculate_covariance(const AYB ayb);

#endif /* _CALL_BASES_H */

