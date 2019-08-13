/************************************************************************

  Demo for  TimeHarp 200   Software version 6.0, Format version 6.0  
  
  T3R File Access

  Michael Wahl, PicoQuant GmbH, April 2006

  NEW: - Special records for external markers
       - Time Tag overflow is encoded in a special record

  Tested with the following compilers:

  - MinGW 2.0.0-3 (free compiler for Win 32 bit)
  - MS Visual C++ 4.0/5.0/6.0 (Win 32 bit)
  - Borland C++ 5.3 / C++Builder 3.0 (Win 32 bit)

  it should work with most others.

  This is demo code. Use at your own risk. No warranties.

************************************************************************/

#include <stdio.h>
#include <string.h>

/*
The following structures are used to hold the TimeHarp file data
They reflect the file structure.
The data types used here to match the file structure are correct
for MSVC++. They may have to be changed for other compilers.
*/

typedef struct{ float Start;
                float Step;
				float End;  } tParamStruct;

typedef struct{ long int MapTo;
				long int Show; } tCurveMapping;

/* The following represents the readable ASCII file header portion} */

struct {
	char Ident[16];
	char FormatVersion[6];	
	char CreatorName[18];		
	char CreatorVersion[12];	
	char FileTime[18];
	char CRLF[2];
	char CommentField[256]; } TxtHdr;

/* The following is binary header information */

struct {
	long int Channels;
	long int Curves;
	long int BitsPerChannel;
	long int RoutingChannels;
	long int NumberOfBoards;
	long int ActiveCurve;
	long int MeasMode;
	long int SubMode;
	long int RangeNo;
	long int Offset;			/* in ns */
	long int Tacq;				/* in ms */
	long int StopAt;
	long int StopOnOvfl;
	long int Restart;
	long int DispLinLog;
	long int DispTimeFrom;
	long int DispTimeTo;
	long int DispCountsFrom;
	long int DispCountsTo;
	tCurveMapping DispCurves[8];
	tParamStruct Params[3];
	long int RepeatMode;
	long int RepeatsPerCurve;
	long int RepeatTime;
	long int RepeatWaitTime;
	char ScriptName[20];} BinHdr;

struct {
	char HardwareIdent[16];  
	char HardwareVersion[8];  
	long int BoardSerial;
	long int CFDZeroCross;
	long int CFDDiscrMin;
	long int SyncLevel;
	long int CurveOffset;
	float Resolution;} BoardHdr;

struct {
	long int Globclock;
	long int ExtDevices; 
	long int Reserved1;
	long int Reserved2;
	long int Reserved3;
	long int Reserved4;
	long int Reserved5;
	long int SyncRate;
	long int TTTRCFDRate;
	long int TTTRStopAfter;
	long int TTTRStopReason;
	long int NoOfRecords;
	long int SpecialHeaderSize;} TTTRHdr;

/*The following data records will appear in the file NoOfRecords times*/

struct {
	unsigned TimeTag	:10;
	unsigned Channel	:15;
	unsigned Route		:6;
	unsigned Valid		:1; }  TTTRrecord;


extern "C" __declspec(dllexport) int _stdcall ht3(char* fpinpath, int *tac, int *channels, double *t, int ovfl, int pos, int nrec)
{

 FILE *fpin;
 long int result,ii;
 unsigned long counts=0, overflows=0;
 double ofltime=0, truetime=0;
 int mode=1;





fpin=fopen(fpinpath,"rb");

/*
 result = fread( &TxtHdr, 1, sizeof(TxtHdr) ,fpin);
 if (result!= sizeof(TxtHdr))
 {
   goto close;
 }


 result = fread( &BinHdr, 1, sizeof(BinHdr) ,fpin);


 if (result!= sizeof(BinHdr))
 {
   goto close;
 }


 result = fread( &BoardHdr, 1, sizeof(BoardHdr) ,fpin);

 result = fread( &TTTRHdr, 1, sizeof(TTTRHdr) ,fpin);
*/

 /* skip the special header (you may need to read it if you
    want to interpret an imaging file */
// fseek(fpin,TTTRHdr.SpecialHeaderSize*4,SEEK_CUR);
fseek(fpin,pos,0);
 /* Now read and interpret the TTTRrecords */


 for(ii=0; ii<nrec; ii++)
 {
	result = fread( &TTTRrecord, 1, sizeof(TTTRrecord) ,fpin);
	if (result!= sizeof(TTTRrecord))
	{
		goto close;
	}


//	truetime = (ofltime + TTTRrecord.TimeTag) * TTTRHdr.Globclock * 1e-9; /* convert to seconds */

	if(TTTRrecord.Valid)
	{


		t[counts]=ofltime + TTTRrecord.TimeTag;
		tac[counts]=TTTRrecord.Channel;
		channels[counts]=TTTRrecord.Route;
		
		counts++;
	}
	else /* this means we have a 'special record' */
	{
		if(TTTRrecord.Channel & 0x800)  /* evaluate this AFTER evaluating the valid record. */
		{								          
			ofltime += 1024; /* unwrap the time tag overflow */
			overflows++;
		}
		
	}
	
 }



close:
 fclose(fpin);



 ovfl=overflows;
 return(counts);

}
