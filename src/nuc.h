/*
 *  Copyright (C) 2010 by Tim Massingham, European Bioinformatics Institute
 *  tim.massingham@ebi.ac.uk
 *
 *  This file is part of the AYB software for simulating likelihoods
 *  for next-generation sequencing machines.
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

#ifndef _NUC_H
#define _NUC_H

#include <stdio.h>
#include <stdint.h>

typedef char NUC;
#define NBASE       4

#define NUC_AMBIG   4
#define NUC_A       0
#define NUC_C       1
#define NUC_G       2
#define NUC_T       3

typedef char PHREDCHAR;

void show_NUC(FILE * fp, const NUC nuc);
void show_PHREDCHAR(FILE * fp, const PHREDCHAR nuc);

NUC read_NUC(FILE * fp);
PHREDCHAR read_PHREDCHAR(FILE * fp);
#define MIN_PHRED   33
#define MAX_PHRED   126
#define ERR_PHRED   0

#define X(A) A ## NUC
#include "array.def"
#undef X
#define X(A) A ## PHREDCHAR
#include "array.def"
#undef X


    
NUC nuc_from_char( const char c);
char char_from_nuc(const NUC nuc);
ARRAY(NUC) nucs_from_string( const char * nucstr );
NUC complement(const NUC nuc);
ARRAY(NUC) reverse_complement(const ARRAY(NUC) nucs);
real_t qual_from_prob( const real_t prob);
PHREDCHAR phredchar_from_quality( const real_t qual);
PHREDCHAR phredchar_from_char( const char c);
PHREDCHAR phredchar_from_prob( const real_t p);
real_t qualadd(const real_t q1, const real_t q2);

#endif

