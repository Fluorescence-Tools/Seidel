//author: Nicolaas van der Voort
//AG Seidel, HHU Dusseldorf
//May 2019

#include "ProcessPhotonStream.h"
#include <stdio.h>
#include <Eigen/Core>
#include <vector>
#include <algorithm>

//typedef Eigen::Array<unsigned char, Eigen::Dynamic, 1> ArrayXu;

ph Eigen_array(
	Eigen::Array<long long, Eigen::Dynamic, 1> t
	, Eigen::ArrayXi tac
	, Eigen::Array<unsigned char, Eigen::Dynamic, 1>  can
) {
	ph Ph;
	tac[0] = 2;
	if (t[0] == 1) {
		tac[0] = 3;
	}
	Ph.tac = tac[0];
	Ph.t = t[0];
	Ph.can = can[0];
	printf("hellow world! \n");
	return Ph;
}

std::vector<imChannel> ProcessPhotonStream(
	Eigen::ArrayXi tac
	, Eigen::Array<long long, Eigen::Dynamic, 1> t
	, Eigen::Array<unsigned char, Eigen::Dynamic, 1>  can
	, imOpts ImOpts
	, std::vector<imChannel> Channels
)
{
	int mode_iter;
	long long startline, i;
	int currentline = 0;
	int currentframe = 0;
	int currentline_id = 0;
	bool incan, intac, in_t, inmode, inscan = false;
	//bool incan = true, intac = true, in_t = true, inscan = true, inmode = true;
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
	return Channels;
}
