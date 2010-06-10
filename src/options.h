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

#ifndef _AYB_OPTIONS_H
#define _AYB_OPTIONS_H

#include <stdint.h>
#include <stdio.h>
#include "utility.h"

extern unsigned int NOUTPUTFORMAT;
extern CSTRING output_format_str[];

enum output_format_enum { OUTPUT_FASTA, OUTPUT_FASTQ, OUTPUT_QSEQ, OUTPUT_LIKE, OUTPUT_INVALID};

typedef struct {
    CSTRING machine_name;
    uint32_t run_number,lane,tile;
    enum output_format_enum output_format;
    uint32_t niter;
    real_t mu;
    FILE * output_fp;
    CSTRING coordinate_file, dump_file_prefix;
    real_t spike_threshold;
    bool remove_negative;
    real_t min_lambda;
} AYBOPT;


#endif /* _AYB_OPTIONS_H */

