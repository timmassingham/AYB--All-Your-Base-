/*
 *  Copyright (C) 2010 by Tim Massingham
 *  tim.massingham@ebi.ac.uk
 *
 *  This file is part of AYB.
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
#include <getopt.h>
#include <math.h>
#include "cif.h"
#include "ayb.h"
#include "utility.h"
#include "options.h"

#define PROGNAME "AYB"

void fprint_usage( FILE * fp){
    validate(NULL!=fp,);
    fputs(
"\t\"" PROGNAME "\"\n"
"Call bases for Illumina data\n"
"\n"
"Usage:\n"
"\t" PROGNAME " data.cif\n"
"\t" PROGNAME " --help\n"
"\t" PROGNAME " --licence\n"
"\n"
PROGNAME " from CIF file and prints FASTQ formated sequence to stdout.\n"
"\n"
"Example:\n"
"\t" PROGNAME " s_4_0033.cif > s_4_0033.fastq\n"
,fp);
}

void fprint_licence (FILE * fp){
    validate(NULL!=fp,);
    fputs(
"  " PROGNAME ": Call bases for Illumina data\n"
#include "copyright.inc"
    ,fp);
}

void fprint_help( FILE * fp){
    validate(NULL!=fp,);
    fputs(
/*
12345678901234567890123456789012345678901234567890123456789012345678901234567890
*/
"\n"
"-h, --help\n"
"\tDisplay information about usage and exit.\n"
"\n"
"--licence\n"
"\tDisplay licence terms and exit.\n"
, fp);
}

static struct option longopts[] = {
    { "mu",         required_argument, NULL, 'm'},
    { "niter",      required_argument, NULL, 'n'},
    { "help",       no_argument,       NULL, 'h' },
    { "licence",    no_argument,       NULL, 0 },
};


AYBOPT aybopt = {
    .niter = 5,
    .mu = 1e-5
};


unsigned int parse_uint( const CSTRING str){
    validate(NULL!=str,0);
    unsigned int n=0;
    sscanf(str,"%u",&n);
    return n;
}


real_t parse_real( const CSTRING str){
    validate(NULL!=str,NAN);
    real_t x = NAN;
    sscanf(str,real_format_str,&x);
    return x;
}


void parse_arguments( const int argc, char * const argv[] ){
        int ch;
        while ((ch = getopt_long(argc, argv, "m:n:h", longopts, NULL)) != -1){
        switch(ch){
            case 'm':  aybopt.mu = parse_real(optarg);
                       if(aybopt.mu<0){ errx(EXIT_FAILURE,"Mu mustbe greater than zero.");}
                       break;
            case 'n':  aybopt.niter = parse_uint(optarg);
                       break;
            case 'h':
                fprint_usage(stderr);
                fprint_help(stderr);
                exit(EXIT_SUCCESS);
            case 0:
                fprint_licence(stderr);
                exit(EXIT_SUCCESS);
            default:
                fprint_usage(stderr);
                exit(EXIT_FAILURE);
            }
        }
}

void dump_fastq(FILE * fp, const AYB ayb){
    validate(NULL!=fp,);
    validate(NULL!=ayb,);
    const uint32_t ncycle = ayb->ncycle;
    const uint32_t ncluster = ayb->ncluster;
    
    for ( uint32_t cl=0 ; cl<ncluster ; cl++){
        fprintf(fp,"@cluster_%u\n",cl+1);
        for ( uint32_t cy=0 ; cy<ncycle ; cy++){
            show_NUC(fp,ayb->bases.elt[cl*ncycle+cy]);
        }
        fputs("\n+\n",fp);
        for ( uint32_t cy=0 ; cy<ncycle ; cy++){
            show_PHREDCHAR(fp,ayb->quals.elt[cl*ncycle+cy]);
        }
        fputc('\n',fp);
    }
}

void analyse_tile( XFILE * fp){
    CIFDATA cif = NULL;
    AYB ayb = NULL;
    
    validate(NULL!=fp && xnotnull_file(fp),);
    cif = readCIFfromStream(fp);
    if(NULL==cif){
        xfputs("Failed to read tile.\n",xstderr);
        goto cleanup;
    }
    xfprintf(xstderr,"Read tile: %u cycles from %u clusters.\n",cif_get_ncycle(cif),cif_get_ncluster(cif));

    ayb = initialise_AYB(cif);
    if(NULL==ayb){
        xfputs("Failed to initialise.\n",xstderr);
        goto cleanup;
    }
    free_cif(cif);
    timestamp("Initialised\n",stderr);

    for ( int i=0 ; i<aybopt.niter ; i++){
        xfprintf(xstderr,"\tIteration %d\n",i+1);
        timestamp(" * parameter estimation\n",stderr);
        estimate_MPC(ayb);
        timestamp(" * base calling\n",stderr);
        estimate_Bases(ayb);
    }
    
    timestamp("Writing fastq\n",stderr);
    
    dump_fastq(stdout,ayb);
    free_AYB(ayb);
    return;
    
cleanup:
    free_AYB(ayb);
    free_cif(cif);
}

int main(int argc, char * argv[]){
    XFILE * fp=NULL;
    parse_arguments(argc,argv);
    argc -= optind;
    argv += optind;
    
    xfputs("Starting.\n",xstderr);
    if( 0==argc ){
        xfputs("No CIF file given on command line. Reading from stdin\n",xstderr);
        fp = xstdin;
        analyse_tile(fp);
    } else {
        for ( int i=0 ; i<argc ; i++){
            fp = xfopen(argv[i],XFILE_UNKNOWN,"r");
            if(NULL==fp){
                xfprintf(xstderr,"Failed to open file \"%s\" for input. Skipping.\n",argv[i]);
                continue;
            }
            analyse_tile(fp);
            xfclose(fp); fp=NULL;
        }
    }
    xfputs("No more tiles left. Finished.\n",xstderr);
    return EXIT_SUCCESS;
}
