//author: Nicolaas van der Voort
//AG Seidel, HHU Dusseldorf
//May 2019
#include "ProcessPhotonStream.h"
#include <stdio.h>
//#include <Eigen/Eigen>
#include <vector>
#include <algorithm>


void ProcessPhotonStream(
		int * tac,
		long long * t,
		unsigned char * can,
		imOpts ImOpts,
		std::vector<imChannel> Channels
		)
{
	int mode_iter;
	long long startline, i;
	int currentline = 0;
	int currentframe = 0;
	int currentline_id = 0;
	bool incan, intac, in_t, inmode, inscan;
	//scanspeed = pxxsize / dwelltime.
	// time_since_line_start = macrot_since_line_start * counttime
	// pos = scanspeed * time_since_line_start
	float macrot2pos = ImOpts.pxsize * ImOpts.counttime / ImOpts.dwelltime;
	for (i = 0; i < ImOpts.NumRecords; i++){
		if (can[i] == 65) {
			inscan = true;
			startline = t[i];
			continue;
		}
		else if (can[i] == 66) {
			inscan = false;
			currentline++;//line numbers start with 0
			mode_iter = currentline % ImOpts.line_ids.size();
			currentline_id = ImOpts.line_ids[mode_iter];
			continue;
		}
		//The first frame marker comes before the first frame
		//and the last frame marker comes after the last frame
		//for n frames, there are n+1 frame markers
		//the first frame has label '1'
		else if (can[i] == 68){
			currentframe++;
			currentline = 0;
			continue;
		}
	}
	
	for (imChannel Ch : Channels){
		intac = tac[i] > Ch.tacmin and tac[i] < Ch.tacmax;
		in_t = t[i] > Ch.tmin and t[i] < Ch.tmax; 
		inmode = Ch.line_id == currentline_id;
		incan = std::find(Ch.can.begin(),Ch.can.end(),can[i]) 
				!= Ch.can.end();
		if (intac and in_t and incan and inscan and inmode){
			ph Ph;
			Ph.tac = tac[i];
			Ph.t = t[i];
			Ph.can = can[i];
			Ph.n = 	i;
			Ph.x = (t[i]-startline) * macrot2pos;
			Ph.y = currentline * ImOpts.linestep;
			Ph.frame = currentframe;	
			
			Ch.phstream.push_back(Ph);
			//current behaviour is that a photon can only be in a single
			//channel
			break;
		}
	}
}



extern "C"
{
	__declspec(dllexport) int SplitOnTacs(long long * eventN,
		int * tac,
		long long * t,
		unsigned char * can,
		int dimX,
		int dimY,
		float dwelltime,
		float counttime,
		long long NumRecords,
		int gate,
		int nlines,
		char * uselines,
		int * imA,
		int * imB) {
		int i, startline, pixel;
		int j = 0, line = 0;
		bool inscan(false);
		float macrotime2pixel = counttime / dwelltime;
		int timedif;
		for (i = 0; i < NumRecords; i++) {
			if (can[i] == 65) {
				inscan = true;
				startline = t[i];
			}
			else if (can[i] == 66) {
				inscan = false;
				j++;
				if (j >= nlines) {
					line++;
					j = 0;
				}
			}
			else if ((can[i] == 1 or can[i] == 3) and inscan and uselines[j] and tac[i] > gate) {
				//image is stored as 1D array
				//alculated once for computational efficiency
				timedif = t[i] - startline; //cast long long to int
				pixel = timedif * macrotime2pixel;
				if (pixel < dimX) {
					if (tac[i] % 2 == 0) imA[line * dimX + pixel]++;
					else if (tac[i] % 2 == 1) imB[line * dimX + pixel]++;
				}

			}
		}
		return 0;
	}
	
//this function will be deprecated in favour of ProcessPhotonStrean function
	__declspec(dllexport) int genGRYlifetime(
		long long * eventN,
		int * tac,
		long long * t,
		unsigned char * can,
		int dimX,
		int dimY,
		int ntacs,
		int tac_range,
		float dwelltime,
		float counttime,
		long long NumRecords,
		int nlines,
		int framestop,
		char * uselines,
		char * Gchan,
		char * Rchan,
		char * Ychan,
		int * imG,
		int * imR,
		int * imY) {
		long long index;
		int i, startline, pixel;
		int j = 0, line = 0, frame = 0;
		bool inscan(false);
		float macrotime2pixel = counttime / dwelltime;
		int tacconvert = tac_range / ntacs;
		int timedif;
		for (i = 0; i < NumRecords; i++) {
			if (can[i] == 65) {
				inscan = true;
				startline = t[i];
			}
			else if (can[i] == 66) {
				inscan = false;
				j++;
				if (j >= nlines) {
					line++;
					j = 0;
				}
			}
			//The first frame marker comes before the first frame
			//and the last frame marker comes after the last frame
			//for n frames, there are n+1 frame markers
			//the first frame has label '1'
			//currently only all or untill slicestop frames are taken
			//slicestop = -1 means all frames are taken (default)
			//multiple frame selection may be implemented later.
			else if (can[i] == 68) {
				line = 0;
				frame++;
				if (frame == framestop) {
					break;
				}
			}
			//check if pixels goes into G, R or Y channel.
			//uselines[j] = 0 means skip line
			//uselines[j] = 1 means FRET sensitised (G or R)
			//uselines[j] = 2 means Acceptor sensitiyed (Y)
			//Gchan, Rchan and Ychan are used in addition to select correct image
			else if (inscan) {
				//image is stored as 1D array
				//macrotime2pixel is calculated once for computational efficiency
				timedif = t[i] - startline; //cast long long to int
				pixel = timedif * macrotime2pixel;
				index = line * ntacs * dimX + pixel * ntacs + 
					tac[i] / tacconvert;
				if (pixel < dimX and 
					(can[i] == Gchan[0] or can[i] == Gchan[1]) and uselines[j] == 1) {
					imG[index]++;
				}
				else if (pixel < dimX and
					(can[i] == Rchan[0] or can[i] == Rchan[1]) and uselines[j] == 1) {
					imR[index]++;
				}
				else if (pixel < dimX and
					(can[i] == Ychan[0] or can[i] == Ychan[1]) and uselines[j] == 2) {
					imY[index]++;
				}

			}
		}
		return 0;
	}
}

