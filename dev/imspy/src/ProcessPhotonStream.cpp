//author: Nicolaas van der Voort
//AG Seidel, HHU Dusseldorf
//May 2019

#include "ProcessPhotonStream.h"
#include <stdio.h>
//#include <Eigen/Eigen>
#include <vector>
#include <algorithm>

//class constructors
//imOpts::imOpts() {}
//ph::ph() {}
//imChannel::imChannel() :tacmin(0), tacmax(32768), tmin(0), tmax(9223372036854775807) {}

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
	for (i = 0; i < ImOpts.NumRecords; i++) {
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
	}

	for (imChannel Ch : Channels) {
		intac = tac[i] > Ch.tacmin && tac[i] < Ch.tacmax;
		in_t = t[i] > Ch.tmin && t[i] < Ch.tmax;
		inmode = Ch.line_id == currentline_id;
		incan = std::find(Ch.can.begin(), Ch.can.end(), can[i])
			!= Ch.can.end();
		if (intac && in_t && incan && inscan && inmode) {
			ph Ph;
			Ph.tac = tac[i];
			Ph.t = t[i];
			Ph.can = can[i];
			Ph.n = i;
			Ph.x = (t[i] - startline) * macrot2pos;
			Ph.y = currentline * ImOpts.linestep;
			Ph.frame = currentframe;

			Ch.phstream.push_back(Ph);
			//current behaviour is that a photon can only be in a single
			//channel
			break;
		}
	}
}
