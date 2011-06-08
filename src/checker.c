#include <stdio.h>
#include <stdint.h>

int main ( int argc, char * argv[]){
   FILE * fp = fopen(argv[1],"r");
   uint32_t count=0,lcount=1,lccount=1;
   int c = 0;

   fprintf(stderr,"Checks if file contains null characters\n");
   while( (c=fgetc(fp)) != EOF){
      if(c==0){ fprintf(stderr,"Error at char %u, line %u\n",lccount,lcount);}
      lccount++;
      if(c=='\n'){
	 lcount++;
	 lccount=1;
      }
   }
}
