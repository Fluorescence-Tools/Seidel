//Eigen wrappers of existing c code for fit2DGaussian function
//Eigen code will facilitate wrapping by pybind to numpy arrays
//Author: Nicolaas van der Voort, AG Seidel, HHU Düsseldorf
//Created December 2, 2019

#include "ProcessPhotonStream.h"
#include <Eigen/Core>
#include <Eigen/LU>
//#include <pybind11/pybind11.h>
//#include <pybind11/iostream.h>// for py::print debugging function
//namespace py = pybind11;
//#include <iostream>


//using namespace Eigen;

// takes numpy array as input and returns another numpy array
Eigen::MatrixXd inv(Eigen::MatrixXd xs) {
	return xs.inverse();
}

int subtract(int i, int j)
{
	return (i - j);
}
//write to separate file later
using namespace Eigen;
int ProcessPhotonStream_pywrap(
		ArrayX<int> tac,
		ArrayX<long long> t,
		ArrayX<unsigned char> can,
		ImOpts imOpts,
		std::vector<ImChannel> Channels){
			int * tac_p = tac.data();
			long long * t_p = t.data();
			unsigned char * can_p = can.data();
			return ProcessPhotonStream(tac_p, t_p, can_p, imOpts, Channels);
		}
int ProcessPhotonStream(
		//tac can be caught in 16 bit short, but my PQ_ptu_sf function expects an int.
		//potentially this already occurs in Surens code, but maybe also my wrapper is at fault.
		//I will except this inefficiency for now.
		int * tac,
		long long * t,
		unsigned char * can,
		ImOpts imOpts,
		std::vector<ImChannel> Channels
		)
{
	int outflag = 0;
	int mode_iter, currentline_id;
	long long startline, i, timedif;
	int lineIndex, colIndex, tacIndex;
	long long index;
	//int currentframe = 0; saving per frame is too computationally intensive.
	int currentframe; //frames not used but kept as legacy
	bool incan, intac, in_t, inmode, inscan;
	//scanspeed = pxxsize / dwelltime.
	// time_since_line_start = macrot_since_line_start * counttime
	// pos = scanspeed * time_since_line_start
	//float macrot2pos = ImOpts.pxsize * ImOpts.counttime / ImOpts.dwelltime;
	float macrotime2pixel = imOpts.counttime / imOpts.dwelltime;
	int tacconvert = imOpts.TAC_range / imOpts.ntacs;
	int chsize = Channels.size();
	
	
	for (i = 0; i < imOpts.NumRecords; i++){
		if (can[i] == 65) {
			inscan = true;
			startline = t[i];
			continue;
		}
		else if (can[i] == 66) {
			inscan = false;
			lineIndex++;//line numbers start with 0
			mode_iter = lineIndex % imOpts.line_ids.size();
			currentline_id = imOpts.line_ids[mode_iter];
			continue;
		}
		//not using frame counting
		//The first frame marker comes before the first frame
		//and the last frame marker comes after the last frame
		//for n frames, there are n+1 frame markers
		//the first frame has label '1'
		else if (can[i] == 68){
			currentframe++;
			lineIndex = 0;
			continue;
		}
		for(int j = 0; j < chsize ; j++) {
		//one could use the line below instead to loop, but it does not work for some reason
		//for (ImChannel Ch : Channels){ 
			intac = tac[i] > Channels[j].tacmin and tac[i] < Channels[j].tacmax;
			in_t = t[i] > Channels[j].tmin and t[i] < Channels[j].tmax; 
			inmode = Channels[j].line_id == currentline_id;
			
			incan = std::find(Channels[j].can.begin(),Channels[j].can.end(),can[i]) 
					!= Channels[j].can.end();
			
			//it seems very inefficient to me to make an object for every photon, as we have so many.
			//I guess the channel sorting is OK, but I should just cast it in a 3D array as before.
			//can use the bookkeeping from previous function
			//TODO: how do I know the dimensions of the array?
			//answer: dimX and dimY are read by the header,
			//if tac binning is used, that must be specified by the user
			//TODO: create a container for the imChannel data, probably an Eigen array that translates no numpy and that I have to give the proper dimensions.
			
			if (intac and in_t and incan and inscan and inmode){
				timedif = t[i] - startline; //no cast
				colIndex = timedif * macrotime2pixel;
				tacIndex = tac[i] / tacconvert;
				
				index = lineIndex * imOpts.ntacs * imOpts.dimX 
					+ colIndex * imOpts.ntacs 
					+ tacIndex;
				//check that there is no index overflow
				if ( index > Channels[j].ltImage.size() ) {
					//outflag += 1; //bit code to signify index overflow
					break;
				} 
				outflag += 1;
				//issue: somehow this does not really change the object. Pointer problem?
				//yes -> objects are being copied, hence any modification here is not refleted in python.
				//this is also not possible, making this whole pybind11 module pretty sucky.
				Channels[j].ltImage[index] = 10;
				//current behaviour is that a photon can only be in a single
				break;

			}
			
		}
	Channels[0].ltImage[0] = 200;
	Channels[0].name = "lala";
	Channels[0].ltImage.push_back(20);
	tac[0] = 100;
		
	}
	return outflag;
}
//old stuff
//old photon struct, discarded
				/*
				ph Ph;
				Ph.tac = tac[i];
				Ph.t = t[i];
				Ph.can = can[i];
				Ph.n = 	i;
				Ph.x = (t[i]-startline) * macrot2pos;
				Ph.y = lineIndex * ImOpts.linestep;
				Ph.frame = currentframe;	
				
				Ch.phstream.push_back(Ph);
				
				//channel
				break;
				*/
/*
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
				//bug: Chan must have length 2, but this is not enforced.
				//use std::<vector> type instead.
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
		*/