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

#include <stdio.h>
#include <tgmath.h>
#include "utility.h"
#include "nuc.h"

void show_NUC(FILE * fp, const NUC nuc){
    validate(NULL!=fp,);
    fputc(char_from_nuc(nuc),fp);
}

void show_PHREDCHAR(FILE *fp, const PHREDCHAR pc){
	validate(NULL!=fp,);
	fputc(pc,fp);
}

NUC read_NUC(FILE *fp){
	return nuc_from_char(fgetc(fp));
}

PHREDCHAR read_PHREDCHAR(FILE * fp){
	return phredchar_from_char(fgetc(fp));
}

NUC nuc_from_char( const char c){
    switch(c){
        case 'A':
        case 'a': return NUC_A;
        case 'C':
        case 'c': return NUC_C;
        case 'G':
        case 'g': return NUC_G;
        case 'T':
        case 't': return NUC_T;
        case 'N':
        case 'n': return NUC_AMBIG;
        default:
            fprintf(stderr,"Unrecognised nucleotide \"%c\". Returning NUC_AMBIG\n",c);
    }
    return NUC_AMBIG;
}

char char_from_nuc(const NUC nuc){
    switch(nuc){
    case NUC_A: return 'A';
    case NUC_C: return 'C';
    case NUC_G: return 'G';
    case NUC_T: return 'T';
    }
    return 'N';
}

ARRAY(NUC) nucs_from_string( const char * nucstr ){
    validate(NULL!=nucstr,null_ARRAY(NUC));
    const uint32_t len = strlen(nucstr);
    ARRAY(NUC) nucs = new_ARRAY(NUC)(len);
    validate(0!=nucs.nelt,nucs);
    for( uint32_t i=0 ; i<len ; i++){
        nucs.elt[i] = nuc_from_char(nucstr[i]);
    }
    return nucs;
}
    

NUC complement(const NUC nuc){
    switch(nuc){
    case NUC_A: return NUC_T;
    case NUC_C: return NUC_G;
    case NUC_G: return NUC_C;
    case NUC_T: return NUC_A;
    }
    return NUC_AMBIG;
}

ARRAY(NUC) reverse_complement(const ARRAY(NUC) nucs){
    validate(NULL!=nucs.elt,null_ARRAY(NUC));
    ARRAY(NUC) new_nuc = new_ARRAY(NUC)(nucs.nelt);
    validate(NULL!=new_nuc.elt,new_nuc);
    for ( uint32_t i=0 ; i<nucs.nelt ; i++){
        new_nuc.elt[i] = complement(nucs.elt[nucs.nelt-i-1]);
    }
    return new_nuc;
}

real_t qual_from_prob( real_t prob){
    return -10.*log1p(-prob)/log(10.);
}

PHREDCHAR phredchar_from_quality( real_t qual){
   real_t c = 33+qual;
   if(!isfinite(qual)){ c = MIN_PHRED; }
   if(c<MIN_PHRED){c=MIN_PHRED;}
   if(c>MAX_PHRED){c=MAX_PHRED;}
   return (PHREDCHAR)(c);
}

PHREDCHAR phredchar_from_char( const char c){
    validate(c>=MIN_PHRED,MIN_PHRED);
    validate(c<=MAX_PHRED,MAX_PHRED);
    return c;
}

PHREDCHAR phredchar_from_prob( real_t p){
    //validate(isprob(p),ERR_PHRED);
    real_t c = 33.0-10.0*log1p(-p)/log(10.0);
    if(c<MIN_PHRED){c=MIN_PHRED;}
    if(c>MAX_PHRED){c=MAX_PHRED;}
    return (PHREDCHAR)(c+0.5);
}

real_t qualadd (const real_t q1, const real_t q2){
	return (q1>q2)?
	       q1 - 10.0 * log1p( pow(10.0,(q1-q2)/10.0))/M_LN10 :
	       q2 - 10.0 * log1p( pow(10.0,(q2-q1)/10.0))/M_LN10;
}


#ifdef TEST
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main ( int argc, char * argv[]){
    if(argc!=2){
        fputs("Usage: test sequence\n",stderr);
        return EXIT_FAILURE;
    }
    uint32_t nnuc = strlen(argv[1]);
    ARRAY(NUC) nucs = new_ARRAY(NUC)(nnuc);
    for ( uint32_t i=0 ; i<nucs.nelt ; i++){
        nucs.elt[i] = nuc_from_char(argv[1][i]);
    }
    
    fputs("Read sequence:     ",stdout);
    for ( uint32_t i=0 ; i<nucs.nelt ; i++){
        fputc(char_from_nuc(nucs.elt[i]),stdout);
    }
    fputc('\n',stdout);
    
    ARRAY(NUC) rc = reverse_complement(nucs);
    fputs("Reversed sequence: ",stdout);
    for ( uint32_t i=0 ; i<rc.nelt ; i++){
        fputc(char_from_nuc(rc.elt[i]),stdout);
    }
    fputc('\n',stdout);
    
}
#endif

