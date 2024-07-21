/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 fftm2l0_fits
                 M. Juda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

created 1997 August 11

Converts HRC IPI team "fixed fast telem" format files (fftm) into 
ASC level 0 FITS files.

$header$

*/

/*** include files ***/
#include <stdio.h>
#include <string.h>

#include "fitsio.h"

#include "nhrc.h"

#define CLK_PERIOD 15.625e-6
#define FRM_PERIOD 0.205

void printerror( int status);

/*============================================================*/
int main(argc,argv)
int argc;
char *argv[];

{
  FILE *fftmfile;
  fitsfile *fitsfile;
  char *progname;

  /* FITSIO variables */
  int bitpix = 16, hdutype, status = 0;
  long naxis = 0;
  long naxes[2] = {0, 0};
  long firstrow, firstelem;
  int colnum;
  int tfields = 20;
  long nrows = 0, numrows;
  char *inst[] = { "HRC-I", "HRC-S" };
  char extname[] = "EVENTS";
  char *ttype[] = {"TIME", "MJF", "MNF", "ENDMNF", "CRSV", "CRSU", "AMP_SF", 
		   "AV1", "AV2", "AV3", "AU1", "AU2", "AU3", "PHA", "E_TRIG", 
		   "VETOSTT", "DET_ID", "SUB_MJF", "CLKTICKS", "QUALITY"};
  char *tform[] = { "1D", "1J", "1J", "1J", "1B", "1B", "1B", "1I", "1I", 
		    "1I", "1I", "1I", "1I", "1B", "2X", "8X", "1X", "1B", 
		    "1J", "20X" };
  char *tunit[] = { "s", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", 
		    " ", " ", "chan", " ", " ", " ", " ", " ", " " };

  int tfields2 = 2;
  char *ttype2[] = {"START", "STOP"};
  char *tform2[] = {"1D", "1D"};
  char *tunit2[] = {"s", "s"};

  int c, ctr, jj, ncycle, init_row;
  int ftel_frame = 0;
  unsigned short int lsbyte, msbyte;
  int p, byte2, byte1, byte0;
  int majfc = 0, minfc = 0, evt = 0;

  struct ftm_data *ftel;

  /* these are the output columns of the EVENTS extension of the FITS file */
  double *time;
  long int *mjf, *startmnf, *stopmnf, *clkticks;
  unsigned char  *crsv, *crsu, *amp_sf, *pha, *e_trig, *vetostt, *det_id, 
    *sub_mjf, *quality;
  short int *av1, *av2, *av3, *au1, *au2, *au3;

  /* these are the output columns of the TIMES extension of the FITS file */
  double tstart[1], tstop[1];

  short int last_sub = 0, sub_roll = 0;

  ftel = CALLOC(1, struct ftm_data);

  progname = strrchr(argv[0], '/');
  if(progname)
    progname++;
  else
    progname = argv[0];

  if( argc != 2 )
    {
      fprintf(stderr, "Usage:\n\t%s <FITS_file> < fftm\n", progname);
      exit(1);
    }

  /*  if( (fftmfile = fopen(argv[1],"rb")) == NULL)
    {
      fprintf(stderr, "File named: %s can not be opened\n", argv[1]);
      exit(1);
    } */

  if (fits_create_file(&fitsfile, argv[1], &status)) /* create new FITS file */
    printerror( status );           /* call printerror if error occurs */

  /* write the required keywords for the primary array image */
  if ( fits_create_img(fitsfile, bitpix, naxis, naxes, &status) )
    printerror( status );          


  /* write a MISSION keyword to the header */
  if( fits_write_key_str( fitsfile, "MISSION", "AXAF", 
			  "Advanced X-Ray Astrophysics Facility", &status))
    printerror( status );
  /* write a TELESCOP keyword to the header */
  if( fits_write_key_str( fitsfile, "TELESCOP", "NONE",
			  "No telescope used - Laboratory data", &status))
    printerror( status );

  /* write an INSTRUME keyword to the header */
  if( fits_write_key_str( fitsfile, "INSTRUME", "HRC", "Instrument", 
			  &status))
    printerror( status );

  if( fits_write_key_str( fitsfile, "ORIGIN", "ASC", "", &status))
    printerror( status );
  if( fits_write_key_str( fitsfile, "CREATOR", progname, 
			  "Program that created this file", &status))
    printerror( status );

  /* append a new empty binary table onto the FITS file */
  if ( fits_create_tbl( fitsfile, BINARY_TBL, nrows+1, tfields, ttype, tform,
			tunit, extname, &status) )
    printerror( status );

  /* write keywords to extension header */
  if( fits_write_key_str( fitsfile, "ORIGIN", "ASC", "", &status))
    printerror( status );
  if( fits_write_key_str( fitsfile, "CREATOR", progname, 
			  "Program that created this file", &status))
    printerror( status );
  if( fits_write_key_str( fitsfile, "MISSION", "AXAF", 
			  "Advanced X-Ray Astrophysics Facility", &status))
    printerror( status );
  if( fits_write_key_str( fitsfile, "TELESCOP", "NONE",
			  "No telescope used - Laboratory data", &status))
    printerror( status );
  if( fits_write_key_str( fitsfile, "INSTRUME", "HRC", "Instrument", 
			  &status))
    printerror( status );

  fits_get_rowsize(fitsfile, &numrows, &status);

  time = CALLOC(numrows, double);
  mjf = CALLOC(numrows, long);
  startmnf = CALLOC(numrows, long);
  stopmnf = CALLOC(numrows, long);
  clkticks = CALLOC(numrows, long);
  crsv = CALLOC(numrows, unsigned char);
  crsu = CALLOC(numrows, unsigned char);
  amp_sf = CALLOC(numrows, unsigned char);
  pha = CALLOC(numrows, unsigned char);
  e_trig = CALLOC(numrows, unsigned char);
  vetostt = CALLOC(numrows, unsigned char);
  det_id = CALLOC(numrows, unsigned char);
  sub_mjf = CALLOC(numrows, unsigned char);
  quality = CALLOC(3*numrows, unsigned char);
  av1 = CALLOC(numrows, short);
  av2 = CALLOC(numrows, short);
  av3 = CALLOC(numrows, short);
  au1 = CALLOC(numrows, short);
  au2 = CALLOC(numrows, short);
  au3 = CALLOC(numrows, short);

  firstrow  = 1;  /* first row in table to write   */
  firstelem = 1;  /* first element in row  (ignored in ASCII tables) */

  jj = 0;
  ncycle = 0;

  /* read in a fixed raw data frame */
  while( (c=fread( ftel->tm, FAST_HRC_FRAME_SIZE, 1, stdin)) == 1 )
    {

    /* increment the total input frame counter and other counters */
    ftel_frame++;
    for (ctr=0; ctr<512; ctr++) 
      {
	evt++;
	if (evt == 64 ) 
	  {
	    evt = 0;
	    minfc++;
	    if (minfc == 64 ) 
	      {
		minfc = 0;
		majfc++;
	      }
	  }
	/* initialize byte pointer */
	p = 16 * ctr;

	/* in the "fast telem" format there are no real telemetry major or
	   minor frames - we will use these arbitrary counters as the 
	   HRC IPI team does in their processing */
	mjf[jj] = majfc;
	startmnf[jj] = minfc;
	stopmnf[jj] = minfc;

	/* look for event trigger bits for good event; they're in Byte 12 */
	e_trig[jj] = (ftel->tm[p+12] & 0xc0) >> 6;
	if(e_trig[jj] != 0) {
	  /* have an event */
	  nrows++;

	  crsv[jj] = ftel->tm[p];
	  crsu[jj] = ((ftel->tm[p+1] >> 2) & 0x003f);

	  /* get the amplier (u and v ) scale factor Byte 1 */
	  p++;
	  amp_sf[jj] = (ftel->tm[p] & 0x0003);

	  /* get amplifier av1 Bytes 2 & 3 */
	  p++;
	  msbyte = (USI)((ftel->tm[p] << 4) & 0x0ff0);
	  lsbyte = (USI)((ftel->tm[p+1] >>4) & 0x000f);
	  av1[jj] = msbyte ^ lsbyte;

	  /* get amplifier av2 Bytes 3 & 4 */
	  p++;
	  msbyte = (USI)((ftel->tm[p] << 8) & 0x0f00);
	  lsbyte = (USI)(ftel->tm[p+1] & 0x00ff);
	  av2[jj] =  msbyte ^ lsbyte;

	  /* get amplifier av3 Bytes 5 & 6 */
	  p++;
	  p++;
	  msbyte = (USI)((ftel->tm[p] << 4) & 0x0ff0);
	  lsbyte = (USI)((ftel->tm[p+1] >>4) & 0x000f);
	  av3[jj] = msbyte ^ lsbyte;

	  /* get amplifier au1 Bytes 6 & 7 */
	  p++;
	  msbyte = (USI)((ftel->tm[p] << 8) & 0x0f00);
	  lsbyte = (USI)(ftel->tm[p+1] & 0x00ff);
	  au1[jj] = msbyte ^ lsbyte;

	  /* get amplifier au2 Bytes 8 & 9 */
	  p++;
	  p++;
	  msbyte = (USI)((ftel->tm[p] << 4) & 0x0ff0);
	  lsbyte = (USI)((ftel->tm[p+1] >>4) & 0x000f);
	  au2[jj] = msbyte ^ lsbyte;

	  /* get amplifier au3 Bytes 9 & 10 */
	  p++;
	  msbyte = (USI)((ftel->tm[p] << 8) & 0x0f00);
	  lsbyte = (USI)(ftel->tm[p+1] & 0x00ff);
	  au3[jj] = msbyte ^ lsbyte;

	  /* get the total pha Byte 11 */
	  p++;
	  p++;
	  pha[jj] = (ftel->tm[p]);

	  /* get the veto status byte Byte 12 */
	  p++;
	  vetostt[jj] = (((ftel->tm[p] << 2) & 0xfc) ^ 
			    ((ftel->tm[p+1] >> 6) & 0x03) );

	  /* get the detector ID bit from Byte 13 */
	  p++;
	  det_id[jj] = ((ftel->tm[p] >> 5) & 0x0001);

	  /* the subframe counter is in Byte 13 also */
	  sub_mjf[jj] = ((ftel->tm[p] >> 2) & 0x0007);
	  /* keep track of sub_mjf roll over for time calculation */
	  if(sub_mjf[jj] < last_sub) sub_roll++;
	  last_sub = sub_mjf[jj];

	  /* get the number of clock ticks from Bytes 13, 14 & 15 */
	  byte2 = (int)((ftel->tm[p] <<  16) & 0x00030000);
	  byte1 = (int)((ftel->tm[p+1] << 8) & 0x0000ff00);
	  byte0 = (int)(ftel->tm[p+2]);
	  clkticks[jj] = byte2 ^ byte1 ^ byte0;

	  /* we can calculate the event time - period of the internal clock 
	     is CLK_PERIOD - the clkticks counter is reset every FRM_PERIOD */
	  time[jj] = CLK_PERIOD * clkticks[jj] 
	    + (sub_mjf[jj] + 8 * sub_roll) * FRM_PERIOD;

	  /* time of first event is the start time */
	  if(nrows == 1) tstart[0] = time[jj];

	  /* time of last event read is the stop time */
	  tstop[0] = time[jj];

	  jj++;

	  if(jj == numrows)
	    {
	      /* At end of our pre-allocated arrays so write to FITS file */
	      init_row = 1 + ncycle * numrows;

              fits_get_colnum(fitsfile, CASEINSEN, "TIME", &colnum, &status);
              fits_write_col(fitsfile, TDOUBLE, colnum, init_row, 1, numrows, 
			     time, &status);
              fits_get_colnum(fitsfile, CASEINSEN, "MJF", &colnum, &status);
              fits_write_col(fitsfile, TLONG, colnum, init_row, 1, numrows, 
			     mjf, &status);
              fits_get_colnum(fitsfile, CASEINSEN, "MNF", &colnum, 
			      &status);
              fits_write_col(fitsfile, TLONG, colnum, init_row, 1, numrows, 
			     startmnf, &status);
              fits_get_colnum(fitsfile, CASEINSEN, "ENDMNF", &colnum, 
			      &status);
              fits_write_col(fitsfile, TLONG, colnum, init_row, 1, numrows, 
			     stopmnf, &status);
              fits_get_colnum(fitsfile, CASEINSEN, "CRSV", &colnum, &status);
              fits_write_col(fitsfile, TBYTE, colnum, init_row, 1, numrows, 
			     crsv, &status);
              fits_get_colnum(fitsfile, CASEINSEN, "CRSU", &colnum, &status);
              fits_write_col(fitsfile, TBYTE, colnum, init_row, 1, numrows, 
			     crsu, &status);
              fits_get_colnum(fitsfile, CASEINSEN, "AMP_SF", &colnum, &status);
              fits_write_col(fitsfile, TBYTE, colnum, init_row, 1, numrows, 
			     amp_sf, &status);
              fits_get_colnum(fitsfile, CASEINSEN, "AV1", &colnum, &status);
              fits_write_col(fitsfile, TSHORT, colnum, init_row, 1, numrows, 
			     av1, &status);
              fits_get_colnum(fitsfile, CASEINSEN, "AV2", &colnum, &status);
              fits_write_col(fitsfile, TSHORT, colnum, init_row, 1, numrows, 
			     av2, &status);
              fits_get_colnum(fitsfile, CASEINSEN, "AV3", &colnum, &status);
              fits_write_col(fitsfile, TSHORT, colnum, init_row, 1, numrows, 
			     av3, &status);
              fits_get_colnum(fitsfile, CASEINSEN, "AU1", &colnum, &status);
              fits_write_col(fitsfile, TSHORT, colnum, init_row, 1, numrows, 
			     au1, &status);
              fits_get_colnum(fitsfile, CASEINSEN, "AU2", &colnum, &status);
              fits_write_col(fitsfile, TSHORT, colnum, init_row, 1, numrows, 
			     au2, &status);
              fits_get_colnum(fitsfile, CASEINSEN, "AU3", &colnum, &status);
              fits_write_col(fitsfile, TSHORT, colnum, init_row, 1, numrows, 
			     au3, &status);
              fits_get_colnum(fitsfile, CASEINSEN, "PHA", &colnum, &status);
              fits_write_col(fitsfile, TBYTE, colnum, init_row, 1, numrows, 
			     pha, &status);
              fits_get_colnum(fitsfile, CASEINSEN, "E_TRIG", &colnum, &status);
              fits_write_col(fitsfile, TBYTE, colnum, init_row, 1, numrows, 
			     e_trig, &status);
              fits_get_colnum(fitsfile, CASEINSEN, "VETOSTT", &colnum, 
			      &status);
              fits_write_col(fitsfile, TBYTE, colnum, init_row, 1, numrows, 
			     vetostt, &status);
              fits_get_colnum(fitsfile, CASEINSEN, "DET_ID", &colnum, &status);
              fits_write_col(fitsfile, TBYTE, colnum, init_row, 1, numrows, 
			     det_id, &status);
              fits_get_colnum(fitsfile, CASEINSEN, "SUB_MJF", &colnum, 
			      &status);
              fits_write_col(fitsfile, TBYTE, colnum, init_row, 1, numrows, 
			     sub_mjf, &status);
              fits_get_colnum(fitsfile, CASEINSEN, "CLKTICKS", &colnum, 
			      &status);
              fits_write_col(fitsfile, TLONG, colnum, init_row, 1, numrows, 
			     clkticks, &status);
              fits_get_colnum(fitsfile, CASEINSEN, "QUALITY", &colnum, 
			      &status);
              fits_write_col(fitsfile, TBYTE, colnum, init_row, 1, numrows, 
			     quality, &status);

	      ncycle++;
	      jj = 0;
	    }
	}
      }
    }

  /* If we hadn't written the last block of events do it now */
  if(jj != numrows)
    {
      init_row = 1 + ncycle * numrows;

      fits_get_colnum(fitsfile, CASEINSEN, "TIME", &colnum, &status);
      fits_write_col(fitsfile, TDOUBLE, colnum, init_row, 1, jj, 
		     time, &status);
      fits_get_colnum(fitsfile, CASEINSEN, "MJF", &colnum, &status);
      fits_write_col(fitsfile, TLONG, colnum, init_row, 1, jj, 
		     mjf, &status);
      fits_get_colnum(fitsfile, CASEINSEN, "MNF", &colnum, 
		      &status);
      fits_write_col(fitsfile, TLONG, colnum, init_row, 1, jj, 
		     startmnf, &status);
      fits_get_colnum(fitsfile, CASEINSEN, "ENDMNF", &colnum, 
		      &status);
      fits_write_col(fitsfile, TLONG, colnum, init_row, 1, jj, 
		     stopmnf, &status);
      fits_get_colnum(fitsfile, CASEINSEN, "CRSV", &colnum, &status);
      fits_write_col(fitsfile, TBYTE, colnum, init_row, 1, jj, 
		     crsv, &status);
      fits_get_colnum(fitsfile, CASEINSEN, "CRSU", &colnum, &status);
      fits_write_col(fitsfile, TBYTE, colnum, init_row, 1, jj, 
		     crsu, &status);
      fits_get_colnum(fitsfile, CASEINSEN, "AMP_SF", &colnum, &status);
      fits_write_col(fitsfile, TBYTE, colnum, init_row, 1, jj, 
		     amp_sf, &status);
      fits_get_colnum(fitsfile, CASEINSEN, "AV1", &colnum, &status);
      fits_write_col(fitsfile, TSHORT, colnum, init_row, 1, jj, 
		     av1, &status);
      fits_get_colnum(fitsfile, CASEINSEN, "AV2", &colnum, &status);
      fits_write_col(fitsfile, TSHORT, colnum, init_row, 1, jj, 
		     av2, &status);
      fits_get_colnum(fitsfile, CASEINSEN, "AV3", &colnum, &status);
      fits_write_col(fitsfile, TSHORT, colnum, init_row, 1, jj, 
		     av3, &status);
      fits_get_colnum(fitsfile, CASEINSEN, "AU1", &colnum, &status);
      fits_write_col(fitsfile, TSHORT, colnum, init_row, 1, jj, 
		     au1, &status);
      fits_get_colnum(fitsfile, CASEINSEN, "AU2", &colnum, &status);
      fits_write_col(fitsfile, TSHORT, colnum, init_row, 1, jj, 
		     au2, &status);
      fits_get_colnum(fitsfile, CASEINSEN, "AU3", &colnum, &status);
      fits_write_col(fitsfile, TSHORT, colnum, init_row, 1, jj, 
		     au3, &status);
      fits_get_colnum(fitsfile, CASEINSEN, "PHA", &colnum, &status);
      fits_write_col(fitsfile, TBYTE, colnum, init_row, 1, jj, 
		     pha, &status);
      fits_get_colnum(fitsfile, CASEINSEN, "E_TRIG", &colnum, &status);
      fits_write_col(fitsfile, TBYTE, colnum, init_row, 1, jj, 
		     e_trig, &status);
      fits_get_colnum(fitsfile, CASEINSEN, "VETOSTT", &colnum, 
		      &status);
      fits_write_col(fitsfile, TBYTE, colnum, init_row, 1, jj, 
		     vetostt, &status);
      fits_get_colnum(fitsfile, CASEINSEN, "DET_ID", &colnum, &status);
      fits_write_col(fitsfile, TBYTE, colnum, init_row, 1, jj, 
		     det_id, &status);
      fits_get_colnum(fitsfile, CASEINSEN, "SUB_MJF", &colnum, 
		      &status);
      fits_write_col(fitsfile, TBYTE, colnum, init_row, 1, jj, 
		     sub_mjf, &status);
      fits_get_colnum(fitsfile, CASEINSEN, "CLKTICKS", &colnum, 
		      &status);
      fits_write_col(fitsfile, TLONG, colnum, init_row, 1, jj, 
		     clkticks, &status);
      fits_get_colnum(fitsfile, CASEINSEN, "QUALITY", &colnum, 
		      &status);
      fits_write_col(fitsfile, TBYTE, colnum, init_row, 1, jj, 
		     quality, &status);
    }

  /* since we've finished writing to the EVENTS extension let's update the 
     number of rows in the table */

  if(fits_modify_key_lng(fitsfile, "NAXIS2", nrows, "Nunber of rows in table",
			 &status))
    printerror( status );

  /* now to put on the TIMES extension */

  /* append a new empty binary table onto the FITS file */

  if ( fits_create_tbl( fitsfile, BINARY_TBL, 1, tfields2, ttype2, tform2,
                        tunit2, "TIMES", &status) )
    printerror( status );

  fits_write_col(fitsfile, TDOUBLE, 1, 1, 1, 1, tstart, &status);
  fits_write_col(fitsfile, TDOUBLE, 2, 1, 1, 1, tstop, &status);

  if (fits_close_file(fitsfile, &status))
    printerror( status );

  exit(0);
}
/*--------------------------------------------------------------------------*/
void printerror( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/

    char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];
  
    if (status)
      fprintf(stderr, "\n*** Error occurred during program execution ***\n");

    fits_get_errstatus(status, status_str);   /* get the error description */
    fprintf(stderr, "\nstatus = %d: %s\n", status, status_str);

    /* get first message; null if stack is empty */
    if ( fits_read_errmsg(errmsg) ) 
    {
         fprintf(stderr, "\nError message stack:\n");
         fprintf(stderr, " %s\n", errmsg);

         while ( fits_read_errmsg(errmsg) )  /* get remaining messages */
             fprintf(stderr, " %s\n", errmsg);
    }

    exit( status );       /* terminate the program, returning error status */
}
