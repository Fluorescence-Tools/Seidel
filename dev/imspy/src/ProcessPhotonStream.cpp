//author: Nicolaas van der Voort
//AG Seidel, HHU Dusseldorf
//May 2019

#include "ProcessPhotonStream.h"
#include <stdio.h>
#include <Eigen/Core>
#include <vector>
#include <algorithm>

void imChannel::gentacdecay(int ntacs) {
	//unsigned char bitshift;
	int divider = 32768 / ntacs;
	int i, tac;
	//bitshift = log2(32768 / ntacs);
	tacdecay = Eigen::ArrayXi::Zero(ntacs);
	
	for (i = 0; i < phstream.size(); i++) {
		//need to declar all vars unsigned for this to work
		//tac = phstream[i] >>= bitshift;
		tac = phstream[i].tac / divider;
		tacdecay[tac]++;
	}
	
}

int imChannel::log2(int n) {
	int count = 0;
	while (n) {
		count += n > 1;
		n >>= 1;
	}
	return count;
}

//idea: instead of sorting the photons directly, keep index lists
//pass the original tac, t, can by reference or store in object.
//the sorting function appends the photon number to the index file.
//functions exist to return a tac, t or can based on an index list
//other functions exist to generate these indexes. E.g. to select
//only gated data, intensity threshold, macrotime selection.
//also add function that calculated posx and posy and puts in array
void imspy::ProcessPhotonStream()
{
	ph Ph;
	int mode_iter, Ch_iter;
	long long startline, i;
	//long long startline = 0, i;
	int currentline = 0;
	int currentframe = 0;
	int currentline_id = 0;
	bool incan, intac, in_t, inmode, inscan = false;
	//scanspeed = pxxsize / dwelltime.
	// time_since_line_start = macrot_since_line_start * counttime
	// pos = scanspeed * time_since_line_start
	float macrot2pos = ImOpts.pxsize * ImOpts.counttime / ImOpts.dwelltime;


	for (i = 0; i < ImOpts.NumRecords - 1; i++) {

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
		else if (can[i] == 68) {
			currentframe++;
			currentline = 0;
			continue;
		}
		//I experimented with imChannel Ch : Channels, but then I could no longer
		//acces the class variable for writing.
		for (Ch_iter = 0; Ch_iter < Channels.size(); Ch_iter++) {
			imChannel Ch = Channels[Ch_iter];
					intac = tac[i] > Ch.tacmin && tac[i] < Ch.tacmax;
					in_t = t[i] > Ch.tmin && t[i] < Ch.tmax;
					inmode = Ch.line_id == currentline_id;
					incan = std::find(Ch.can.begin(), Ch.can.end(), can[i])
						!= Ch.can.end();
			if (intac && in_t && incan && inscan && inmode) {
				Ph.tac = tac[i];
				Ph.t = t[i];
				Ph.can = can[i];
				Ph.n = i;
				Ph.x = (t[i] - startline) * macrot2pos;
				Ph.y = currentline * ImOpts.linestep;
				Ph.frame = currentframe;

				Channels[Ch_iter].phstream.push_back(Ph);
				//current behaviour is that a photon can only be in a single
				//channel
				break;
			}

		}
	}
}
