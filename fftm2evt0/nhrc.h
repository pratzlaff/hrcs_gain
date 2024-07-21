/* HRC Software System Header File

$Header: /juda1.real/juda1/juda/asc/hrc/code/evt_tools/RCS/nhrc.h,v 1.1 2001/06/25 14:10:28 juda Exp $

$Log: nhrc.h,v $
Revision 1.1  2001/06/25 14:10:28  juda
Initial revision

Revision 3.1  1999/01/28 17:20:30  ssm
Major 1999 Cleanup

Revision 2.1  1998/06/02 15:52:36  ssm
Major clean up to include FITS file formats

Revision 1.1  1998/06/02 01:32:12  ssm
Initial revision


$Revision: 1.1 $

*/

/*** General Defines ***/
#define CALLOC(n,x)  ((x *) calloc(n,sizeof(x)))
#define DEBUG(x)	if( debug >= x )

#define YES 0
#define NO -1
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#define ON 0
#define OFF -1
#define OK (0)
#define DONE -1
#define BAD -1
#define GOOD 0
#define NOP 0
#define OKAY 0

#define EOS '\0'
#define SPACE ' '
#define NL '\n'
#define BUFSIZE 80
#define EOD -1
#define LARGE 2147483647

/*** HRC DATA defines ***/

#define UAMPS 0
#define VAMPS 1

#define DATA 1		
#define NOTYET 0
#define HOUSE 2
#define BADPOS -1
#define INVALID -1

#define DET_CMD 0
#define DET_DATA 1

#define PGOOD 0x00
#define ZERO_SUM 0x01
#define GRID_RATIO 0x02
#define ZERO_PSUM 0x04
#define PHA_RATIO 0x08
#define U_CENTER 0x10
#define V_CENTER 0x20
#define FINE_POS 0x40
#define HOT_SPOT 0x80

#define EVENTS 0
#define SEC 1
#define FRAMES 2

#define CTUE_HRC_FRAME_SIZE 1029 /* (1019 bytes + 6 CCSDS + 4 SYNC) */
#define LAB_HRC_FRAME_SIZE 8192 /* Somewhat arbitrary */
#define FAST_HRC_FRAME_SIZE 8192 /* A singe event * 512*/

#define USI unsigned short int
#define UC unsigned char

/* raw data structure LAB */
struct  raw_data {
  int read_count;	/* number of bytes read on read_data call */
  unsigned char tm[LAB_HRC_FRAME_SIZE];
  int ID;
  int frame_count;
  long buf_time;
};
typedef struct raw_data RawRec, *Raw;

/* tm data structure CTUE */
struct  tm_data {
  int read_count;	/* number of bytes read on read_data call */
  unsigned char tm[CTUE_HRC_FRAME_SIZE];
  int ID;
  int frame_count;
  long buf_time;
};
typedef struct tm_data TmRec, *Tm;

/* fast tm data structure FAST */
struct  ftm_data {
  int read_count;	/* number of bytes read on read_data call */
  unsigned char tm[FAST_HRC_FRAME_SIZE];
  int ID;
  int frame_count;
  long buf_time;
};
typedef struct ftm_data FTmRec, *FTm;

/* define the processed raw data  data structure */
/* this data contains primary science data only  */
struct  processed_raw_data {
  int major_frame;
  int time;
  int day_sec;
  short axaf_day;
  short day_msec;
  short amps[2][3];
  unsigned char minor_frame;
  unsigned char event;
  unsigned char cp[2];
  unsigned char pha;
  unsigned char scale;
  unsigned short uraw;
  unsigned short vraw;
  unsigned char veto_status;
  unsigned char event_status;
};
typedef struct processed_raw_data ProcessedRec, *Processed;

/* define secondary science data structure */
struct secondary_science_data {
  int major_frame;
  int day_sec;
  short axaf_day;
  short day_msec;
  unsigned char minor_frame;
  unsigned short tevt_rate_1;
  unsigned short tevt_rate_2;
  unsigned short vevt_rate_1;
  unsigned short vevt_rate_2;
  unsigned short shield_rate_1;
  unsigned short shield_rate_2;
};
typedef struct secondary_science_data SecondaryRec, *Secondary;

#define PRD_XDIM  64
#define PRD_YDIM  64

/* define expanded event file data structure */
struct epr_event {
  /* these are identical to values in the prd record */
  int major_frame;
  int time;
  int day_sec;
  short axaf_day;
  short day_msec;
  short amps[2][3];
  unsigned short uraw;
  unsigned short vraw;
  unsigned short upos;
  unsigned short vpos;
  unsigned short fupos;
  unsigned short fvpos;
  unsigned char minor_frame;
  unsigned char event;
  unsigned char cp[2];
  unsigned char pha;
  unsigned char scale;
  unsigned char veto_status;
  unsigned char event_status;
  /* these are new values */
  unsigned char smu;
  unsigned char smv;
  unsigned char fmu;
  unsigned char fmv;
  unsigned char sum_amps;
  unsigned char process_status;
};
typedef struct epr_event EprRec, *Epr;

struct fb_event {
  double time;
  double detx;
  double dety;
  double x;
  double y;
  int rawx;
  int rawy;
  int tdetx;
  int tdety;
  short crsv;
  short crsu;
  short amp_sf;
  short av1;
  short av2;
  short av3;
  short au1;
  short au2;
  short au3;
  short chipx;
  short chipy;
  short pha;
  short sumamps;
  short chip_id;
  unsigned char status[4];
};

typedef struct fb_event FbRec, *Fb;

#define EPR_XDIM  16384
#define EPR_YDIM  16384

struct statistics {
	long total_events_in;
	long total_events_out;
	long bad_grid_ratio;
	long bad_pha_ratio;
	long bad_u_dist;	
	long bad_v_dist;	
	long bad_bot;
	long fixed_mfinpos;
	long fixed_pfinpos;
	long vetoed;
	long saturation;
	long hot_spot_events;
};
typedef struct statistics StatRec, *Stat;

struct secondaryscience {
        int dummy;
};
typedef struct secondaryscience SSRec, *SS;

struct rd_sel_input_parameters {
  int debug;
  int d;
  int start;
  int stop;
};

struct frd2prd_input_parameters {
  int debug;
  int d;
  char hot_file[64];
};

struct ftm2prd_input_parameters {
  int debug;
  int d;
  char hot_file[64];
};


struct prd_sel_input_parameters {
  int debug;
  int cul;
  int cuu;
  int cvl;
  int cvu;
  long estart;
  long estop;
};

struct epr_sel_input_parameters {
  int debug;
  int uposl;
  int uposu;
  int vposl;
  int vposu;
  int lsmu;
  int usmu;
  int lsmv;
  int usmv;
  int phal;
  int phau;
  unsigned char veto;
  unsigned char event;
  unsigned char process;
  short saturated;
  long estart;
  long estop;
};
typedef struct epr_sel_input_parameters EprSelIPRec, *EprSelIP;


struct prd2epr_input_parameters {
  int debug;
  int d;
  int degap_sel;
  int spill;
  int wire_charge;
  double cf[2][2];
  double pha_ratio;
  double grid_ratio;
  double amp_gain;
  char degap_file[64];
  char hot_file[64];
};
typedef struct prd2epr_input_parameters Prd2EprIPRec, *Prd2EprIP;

struct epr2epr_input_parameters {
  int debug;
  int d;
  int degap_sel;
  int spill;
  int wire_charge;
  double cf[2][2];
  double pha_ratio;
  double grid_ratio;
  double amp_gain;
  char file_name[64];
  char degap_file[64];
  char hot_file[64];
};
typedef struct epr2epr_input_parameters Epr2EprIPRec, *Epr2EprIP;

struct det_cell {
  int u;
  int v;
  int cts;
  double cu;
  double cv;
};
typedef struct det_cell DetCell, *Cell;

struct fdet_cell {
  int u;
  int v;
  double cts;
  double cu;
  double cv;
};
typedef struct fdet_cell FDetCell, *FCell;

struct fine_pos {
  double fu;
  double usig;
  double fv;
  double vsig;
  int sum;
  double pha;
  double psig;
  double period;
  double persig;
};
typedef struct fine_pos FinePos, *Fine;

struct degap {
  double la;
  double lb;
  double ra;
  double rb;

};
typedef struct degap DeGap, *Degap;

struct hotspot_map {
  short spot;
  unsigned short ul;
  unsigned short vl;
  unsigned short uu;
  unsigned short vu;
};
typedef struct hotspot_map HotSpot, *Hotspot;

void *calloc();


