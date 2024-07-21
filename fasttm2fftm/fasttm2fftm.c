/*
;
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 		HRC tm2ftm.c
		S. Murray
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Read in HRC FAST TM file (xxx.rm) and screen at the
      highest level for proper sync pattern.
      If necessary slip the data stream until a complete good
      major frame is found. Write that out, and continue checking

      No attempt is made to recover partial minor frames of data
      this is a future enhancement that can be added to this process

      This version is for FAST (TM) Data

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
$Header: /usr/home/ssm/hrc/src/ssm/work/RCS/fasttm2fftm.c,v 3.1 1999/01/28 17:20:30 ssm Exp $

$Log: fasttm2fftm.c,v $
Revision 3.1  1999/01/28 17:20:30  ssm
Major 1999 Cleanup

Revision 1.1  1997/06/01 01:01:33  ssm
Initial revision

 * Revision 1.1  1996/12/30  23:34:07  ssm
 * Initial revision
 *
*/

/*** include files ***/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "nhrc.h"

#define RCS "$Revision: 3.1 $"

extern char *optarg;
extern int optind;

int debug = 1;

/*============================================================*/

int main(argc,argv)
int argc;
char *argv[];

{
  int c, evt, i, j, k;
  int ftel_frame, stat;
  char raw_bytes[8192];

  struct ftm_data *ftel;

  ftel = CALLOC(1,struct ftm_data);

  while ((c = getopt(argc, argv, "D:h?")) != EOF){
    switch (c) {
    case 'D':
      debug = atoi(optarg);
      break;
    case 'h':
    case '?': /* print the usage */
      fprintf(stderr, RCS);
      fprintf(stderr,"\nusage:tm2ftm -[D] <rd >frd");
      fprintf(stderr,"\n\tD[0]:\tDebug level\n");
      exit(0);
    }
  }
  
  ftel_frame = 0;
  ftel->frame_count = 0;
  evt = 0;
  stat = BAD;
  /* read in some fast tm  data */
  while( (c=fread( raw_bytes, 8192, 1, stdin)) == 1 ) {
    for ( i=0; i<8192; i=i+2 ) {
      if( raw_bytes[i] ==1 ) {
	stat = GOOD;
	/* event sync */
	k = 16*evt;
	ftel->tm[k++] = raw_bytes[i+1];
	for ( j=0; j<15; j++) {
	  i = i + 2;
	  if ( raw_bytes[i] == 0 ) {
	    ftel->tm[k++] = raw_bytes[i+1];
	  } else {
	    stat = BAD;
	    break;
	  }
	}
	if( stat == GOOD ) {
	  evt++;
	  ftel_frame++;
	  if( evt == 512 ) {
	    fwrite( ftel->tm, FAST_HRC_FRAME_SIZE, 1, stdout);
	    evt = 0;
	  }
	}
      }
    }
  }
  fwrite( ftel->tm, FAST_HRC_FRAME_SIZE, 1, stdout );

  return EXIT_SUCCESS;
}

