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
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include "cif.h"
#include "ayb.h"
#include "utility.h"
#include "options.h"
#include "call_bases.h"
#include "process_intensities.h"

#define PROGNAME "AYB"


void fprint_usage( FILE * fp){
    validate(NULL!=fp,);
    fputs(
"\t\"" PROGNAME "\"\n"
"Call bases for Illumina data\n"
"\n"
"Usage:\n"
"\t" PROGNAME " [-c coordinate_file] [-d prefix] [-f format] [-i iterations] [-l lane] [-L minlam]\n"
"\t    [-m mu] [-n machine_name] [-o outfile] [-r run_number] [-s threshold] [-t tile] -z\n"
"\t    data.cif\n"
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
"-c, --coordinates file [default: none]\n"
"\tFile to read coordinates of clusters. Required for qseq output\n"
"-d, --dump prefix [default: don't dump]\n"
"\tDump information about run parameters.\n"
"-f, --format format [default: fastq]\n"
"\tFormat to output results in. Choices are fasta, fastq and qseq.\n"
"-i, --iterations niterations [default: 5]\n"
"\tNumber of iterations of refinement to use.\n"
"-l, --lane lane [default: 0]\n"
"\tLane number for output\n"
"-L, --lambda lambda [default: 0]\n"
"\tMinimum allowed lambda\n"
"-m, --mu value [default: 1e-5]\n"
"\tMu for quality adjustment\n"
"-n, --name machine_name [default: unknown]\n"
"\tName of machine for qseq output\n"
"-o, --outfile [default: stdout]\n"
"\tFile to write output to.\n"
"-r, --run run_number [default: 0]\n"
"\tRun number for qseq output\n"
"-s , --spike threshold [default: don't remove]\n"
"\tMove spikes by zeroing intensities. If any base in a cycle for\n"
"a given cluster is greater than threshold then set all intensities\n"
"to zero.\n"
"-t, --tile tile [default: 0]\n"
"\tTile number for output.\n"
"-z, --zero [default: don't zero]\n"
"\tZero observations with a negative intensity.\n"
"-h, --help\n"
"\tDisplay information about usage and exit.\n"
"\n"
"--licence\n"
"\tDisplay licence terms and exit.\n"
, fp);
}

static struct option longopts[] = {
    { "coordinates",required_argument, NULL, 'c'},
    { "dump",       required_argument, NULL, 'd'},
    { "format",     required_argument, NULL, 'f'},
    { "iterations", required_argument, NULL, 'i'},
    { "lane",       required_argument, NULL, 'l'},
    { "lambda",	    required_argument, NULL, 'L'},
    { "mu",         required_argument, NULL, 'm'},
    { "name",       required_argument, NULL, 'n'},
    { "outfile",    required_argument, NULL, 'o'},
    { "run",        required_argument, NULL, 'r'},
    { "spike",      required_argument, NULL, 's'},
    { "tile",       required_argument, NULL, 't'},
    { "zero",       no_argument,       NULL, 'z'},
    { "help",       no_argument,       NULL, 'h' },
    { "licence",    no_argument,       NULL, 0 },
};


AYBOPT aybopt = {
    .machine_name = "unknown",
    .run_number = 0,
    .lane = 0,
    .tile = 0,
    .output_format = OUTPUT_FASTQ,
    .niter = 5,
    .mu = 1e-5,
    .coordinate_file = NULL,
    .dump_file_prefix = NULL,
    .spike_threshold = 0.,
    .remove_negative = false,
    .min_lambda = 0.
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
        aybopt.output_fp = stdout;
        
        int ch;
        while ((ch = getopt_long(argc, argv, "c:d:f:i:l:L:m:n:o:r:s:t:zh", longopts, NULL)) != -1){
        switch(ch){
            case 'c':  aybopt.coordinate_file = copy_CSTRING(optarg);
                       break;
	    case 'd':  aybopt.dump_file_prefix = copy_CSTRING(optarg);
			if(NULL==aybopt.dump_file_prefix){
			   errx(EXIT_FAILURE,"Failed to parse prefix for dumping");
			}
		       break;
            case 'f':  aybopt.output_format = OUTPUT_INVALID;
                       for ( int i=0 ; i<NOUTPUTFORMAT ; i++){
                           if(!strcasecmp(optarg,output_format_str[i])){
                               aybopt.output_format = i;
                               break;
                           }
                       }
                       if(aybopt.output_format==OUTPUT_INVALID){
                           fprintf(stderr,"Unrecognised output format %s.\n"
                                          "Choices are:\n", optarg);
                           for ( int i=0 ; i<NOUTPUTFORMAT ; i++){
                               fprintf(stderr,"  %s\n", output_format_str[i]);
                           }
                           exit(EXIT_FAILURE);
                       }
                       break;
            case 'i':  aybopt.niter = parse_uint(optarg);
                       break;
            case 'l':  aybopt.lane = parse_uint(optarg);
                       break;
	    case 'L':  aybopt.min_lambda = parse_real(optarg);
		       if(aybopt.min_lambda<0){ errx(EXIT_FAILURE,"Minimum lambda must be greater than zero.");}
		       break;
            case 'm':  aybopt.mu = parse_real(optarg);
                       if(aybopt.mu<0){ errx(EXIT_FAILURE,"Mu must be greater than zero.");}
                       break;
            case 'n':  aybopt.machine_name = copy_CSTRING(optarg);
                       break;
            case 'o':  aybopt.output_fp = fopen(optarg,"w");
                       // Revert to stdout if failed
                       if(NULL==aybopt.output_fp){
                           fprintf(stderr,"Failed to open %s for output, reverting to stdout\n",optarg);
                           aybopt.output_fp = stdout;
                       }
            case 'r':  aybopt.run_number = parse_uint(optarg);
                       break;
	    case 's':  aybopt.spike_threshold = parse_real(optarg);
		       if(aybopt.spike_threshold<0){
			    errx(EXIT_FAILURE,"Threshold for zeroing intensities must be greater than zero");
		       }
		       break;
            case 't':  aybopt.tile = parse_uint(optarg);
                       break;
	    case 'z':  aybopt.remove_negative = true;
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

void dump_fasta(FILE * fp, const AYB ayb){
    validate(NULL!=fp,);
    validate(NULL!=ayb,);
    const uint32_t ncycle = ayb->ncycle;
    const uint32_t ncluster = ayb->ncluster;
    
    for ( uint32_t cl=0 ; cl<ncluster ; cl++){
        fprintf(fp,">cluster_%u\n",cl+1);
        for ( uint32_t cy=0 ; cy<ncycle ; cy++){
            show_NUC(fp,ayb->bases.elt[cl*ncycle+cy]);
        }
        fputc('\n',fp);
    }
}

void dump_qseq( FILE * fp, const AYB ayb){
    validate(NULL!=fp,);
    validate(NULL!=ayb,);
    const uint32_t ncycle = ayb->ncycle;
    const uint32_t ncluster = ayb->ncluster;

    // Check that coordinates exist, should not be able to happen due to previous checks
    if(!ayb->coordinates){
        errx(EXIT_FAILURE,"Attempting to write qseq output but coordinates not available");
    }
    
    const uint32_t lane = (aybopt.lane==0)?ayb->coordinates->lane:aybopt.lane;
    const uint32_t tile = (aybopt.tile==0)?ayb->coordinates->tile:aybopt.tile;
    for ( uint32_t cl=0 ; cl<ncluster ; cl++){
        // General information about machine
        fprintf(fp,"%s\t%" SCNu32 "\t%" SCNu32 "\t%" SCNu32 "\t",aybopt.machine_name,aybopt.run_number,lane,tile);
        // Coordinates of cluster
        fprintf(fp,"%" SCNu16 "\t%" SCNu16 "\t", ayb->coordinates->x[cl], ayb->coordinates->y[cl]); 
        // Write sequence
        for ( uint32_t cy=0 ; cy<ncycle ; cy++){
            show_NUC(fp,ayb->bases.elt[cl*ncycle+cy]);
        }
        fputc('\t',fp);
        // Write qualities
        for ( uint32_t cy=0 ; cy<ncycle ; cy++){
            show_PHREDCHAR(fp,ayb->quals.elt[cl*ncycle+cy]);
        }
        fputc('\t',fp);
        // Index and read number
        fprintf( fp, "%" SCNu16 "\t%" SCNu16, ayb->index, ayb->readnum);
        // Was filtered
        if(ayb->filtered){
            fputc('\t',fp);
            fputc( ayb->passed_filter[cl] ?'1':'0',fp);
        }

        fputc('\n',fp);
    }
}

/* Refactor: this function repeats much of estimate_Bases in ayb.c */
void dump_likelihood( FILE * fp, const AYB ayb){
    validate(NULL!=fp,);
    validate(NULL!=ayb,);
    const uint32_t ncycle = ayb->ncycle;
    const uint32_t ncluster = ayb->ncluster;

    if(!ayb->coordinates){
	fprintf(stderr,"Attempting to write likelihood output but coordinates not available\nCoordinates will be zero'd.\n");
    }

    // Calculate variance matrices and invert
    MAT * V = calculate_covariance(ayb);
    XFILE * vout = xfopen("information.txt",XFILE_RAW,"w");
    for ( uint32_t cy=0 ; cy<ncycle ; cy++){
       MAT a = V[cy];
       show_MAT(vout,V[cy],0,0);
       V[cy] = invert(a);
       free_MAT(a);
    }
    xfclose(vout);
    
    MAT pcl_int = NULL;
    MAT Minv_t = transpose_inplace(invert(ayb->M));
    MAT Pinv_t = transpose_inplace(invert(ayb->P));
    real_t like[4] = {NAN,NAN,NAN,NAN};
    for ( uint32_t cl=0 ; cl<ncluster ; cl++){
       // Process intensities
       const real_t lambda = ayb->lambda->x[cl];
       pcl_int = process_intensities(ayb->intensities.elt+cl*ncycle*NBASE,Minv_t,Pinv_t,ayb->N,pcl_int);
       // Coordinates
       uint16_t x=0,y=0;
       if(ayb->coordinates){
	  x = ayb->coordinates->x[cl];
	  y = ayb->coordinates->y[cl];
       }
       fprintf(fp,"%" SCNu32 "\t%" SCNu32 "\t%" SCNu16 "\t%" SCNu16, aybopt.lane,aybopt.tile,x,y);
       for ( uint32_t cy=0 ; cy<ncycle ; cy++){
	  call_likelihoods(pcl_int->x+cy*NBASE,lambda,V[cy],like);
          fprintf(fp,"\t%3.2f %3.2f %3.2f %3.2f",like[0],like[1],like[2],like[3]);
       }
       fputc('\n',fp);
    }
    free_MAT(Pinv_t);
    free_MAT(Minv_t);
    free_MAT(pcl_int);
    for(uint32_t cy=0 ; cy<ncycle ; cy++){
       free_MAT(V[cy]);
    }
    free(V);
}



void dump_results(FILE * fp, const AYB ayb){
    switch(aybopt.output_format){
    case OUTPUT_FASTA: dump_fasta(fp,ayb); break;
    case OUTPUT_FASTQ: dump_fastq(fp,ayb); break;
    case OUTPUT_QSEQ:  dump_qseq(fp,ayb);  break;
    case OUTPUT_LIKE:  dump_likelihood(fp,ayb); break;
    default:
        errx(EXIT_FAILURE,"Unrecognised output format %s in %s (%s:%d)\n",
            output_format_str[aybopt.output_format],__func__,__FILE__,__LINE__);
    }
}

size_t xstrlen( const CSTRING str){
   return (NULL!=str)?strlen(str):0;
}

CSTRING xstrcpy( CSTRING dst, const CSTRING src){
   return (NULL!=src) ? strcpy(dst,src) : NULL;
}

CSTRING concat_CSTRING( CSTRING str1, const CSTRING str2){
   if( NULL==str1 && NULL==str2 ){ return NULL; }
   const size_t len1 = xstrlen(str1);
   const size_t tlen = len1 + xstrlen(str2) + 1;
   str1 = reallocf(str1,tlen*sizeof(char));
   if(NULL==str1){ return str1;}
   xstrcpy(str1+len1,str2);
   return str1;
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

    // Remove spikes of high intensity
    if(aybopt.spike_threshold!=0){
       uint32_t nfailed = 0;
       for ( uint32_t cl=0 ; cl<ayb->ncluster ; cl++){
	  for ( uint32_t cy=0 ; cy<ayb->ncycle ; cy++){
	     const uint32_t offset = (cl*ayb->ncycle+cy)*NBASE;
	     bool fail_threshold = false;
	     for ( uint32_t b=0 ; b<NBASE ; b++){
		if(ayb->intensities.elt[offset+b]>aybopt.spike_threshold){
		   fail_threshold = true;
		   break;
		}
	     }
	     if(fail_threshold){
		for ( uint32_t b=0 ; b<NBASE ; b++){
		   ayb->intensities.elt[offset+b] = 0.;
		}
		nfailed++;
	     }
	  }
       }
       fprintf(stderr,"Number failing filter = %" SCNu32 "\n",nfailed);
    }
    // Remove negative observations
    if(aybopt.remove_negative){
       const uint32_t n = ayb->ncluster * ayb->ncycle * NBASE;
       for ( uint32_t i=0 ; i<n ; i++){
	  if(ayb->intensities.elt[i]<0.){
	     ayb->intensities.elt[i] = 0.;
	  }
       }
    }
 
    // Read coordinates if necessary
    if(aybopt.coordinate_file!=NULL){
        ayb->coordinates = read_coordinates(aybopt.coordinate_file,ayb->ncluster);
        if(NULL==ayb->coordinates){
            errx(EXIT_FAILURE,"Failed to read coordinates from %s",aybopt.coordinate_file);
        }
    	if(ayb->coordinates->ncluster != ayb->ncluster){
        	errx(EXIT_FAILURE,"Number of coordinates (%" SCNu32 ") disagrees with number of clusters (%" SCNu32, ayb->coordinates->ncluster, ayb->ncluster);
	}
    }
    
    
    timestamp("Initialised\n",stderr);

    for ( int i=0 ; i<aybopt.niter ; i++){
        xfprintf(xstderr,"\tIteration %d\n",i+1);
        timestamp(" * parameter estimation\n",stderr);
        estimate_MPC(ayb);
        timestamp(" * base calling\n",stderr);
        estimate_Bases(ayb);
    }

    
    timestamp("Writing results\n",stderr);
    dump_results(aybopt.output_fp,ayb);

    // Want to dump parameters
    if(aybopt.dump_file_prefix){
        // Cross-talk
        CSTRING filename = NULL;
	filename = concat_CSTRING(filename,aybopt.dump_file_prefix);
	filename = concat_CSTRING(filename,".M.txt");
        XFILE * fp = xfopen(filename,XFILE_RAW,"w");
	if(fp){
            show_MAT(fp,ayb->M,0,0);
	    xfclose(fp);
	} else {
	    fprintf(stderr,"Failed to open %s for output\n",filename);
	}
        // Phasing
        xstrcpy(filename,aybopt.dump_file_prefix);
	filename = concat_CSTRING(filename,".P.txt");
	fp = xfopen(filename,XFILE_RAW,"w");
	if(fp){
	    show_MAT(fp,ayb->P,0,0);
	    xfclose(fp);
	} else {
	    fprintf(stderr,"Failed to open %s for output\n",filename);
        }
	// Noise
        xstrcpy(filename,aybopt.dump_file_prefix);
	filename = concat_CSTRING(filename,".N.txt");
        fp = xfopen(filename,XFILE_RAW,"w");
        if(fp){
            show_MAT(fp,ayb->N,0,0);
            xfclose(fp);
        } else {
	   fprintf(stderr,"Failed to open %s for output\n",filename);
        }
	// lambda
        xstrcpy(filename,aybopt.dump_file_prefix);
	filename = concat_CSTRING(filename,".lambda.txt");
	fp = xfopen(filename,XFILE_RAW,"w");
        if(fp){
            show_MAT(fp,ayb->lambda,0,0);
            xfclose(fp);
        } else {
	    fprintf(stderr,"Failed to open %s for output\n",filename);
        }
	// Weights
	xstrcpy(filename,aybopt.dump_file_prefix);
	filename = concat_CSTRING(filename,".we.txt");
        fp = xfopen(filename,XFILE_RAW,"w");
        if(fp){
            show_MAT(fp,ayb->we,0,0);
            xfclose(fp);
        } else {
	    fprintf(stderr,"Failed to open %s for output\n",filename);
        }
	xfree(filename);
    }
    free_AYB(ayb);
    return;
    
cleanup:
    free_AYB(ayb);
    free_cif(cif);
}

int main(int argc, char * argv[]){
    srandom(23413);
    XFILE * fp=NULL;
    parse_arguments(argc,argv);
    argc -= optind;
    argv += optind;
    
    // Validate arguments
    if( aybopt.output_format==OUTPUT_QSEQ && aybopt.coordinate_file==NULL){
        errx(EXIT_FAILURE,"Format qseq specified but no coordinate file given");
    }
        
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
