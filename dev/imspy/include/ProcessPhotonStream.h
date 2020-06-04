//Author: Nicolaas van der Voort
//AG Seidel, HHU Dusseldorf
//June 3, 2020
#ifndef PROCESSPHOTONSTREAM
#define PROCESSPHOTONSTREAM

#include<vector>

typedef struct{
	//e.g. line_ids = [1,2]. 1 codes for FRET, 2 for PIE (by convention)
	std::vector<int> line_ids; 
	float dwelltime;
	float counttime;//laser rep rate
	long long NumRecords;
	float linestep; //e.g. 10nm px, 2 lines makes linstep 5e-9m
	float pxsize;
} imOpts;

struct ph{
	int tac;
	long long t;
	unsigned char can;
	long long n; 
	float x;//in m
	float y;//in m
	int frame;
} ;

typedef struct {
	std::vector<unsigned char> can;
	int tacmin;//microtime range
	int tacmax;
	long long tmin;// macrotime range
	long long tmax; 
	int line_id;//e.g. line_id = 1
	std::vector<ph> phstream;
} imChannel;

void ProcessPhotonStream(
		int * tac,
		long long * t,
		unsigned char * can,
		imOpts ImOpts,
		std::vector<imChannel> Channels
		);

#endif
