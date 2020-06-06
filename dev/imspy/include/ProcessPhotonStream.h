//Author: Nicolaas van der Voort
//AG Seidel, HHU Dusseldorf
//June 3, 2020
//TODO: build constructors to take the needed arguments
#ifndef PPS_H
#define PPS_H

#include<vector>
#include <Eigen/Core>

class imOpts{
public:
	imOpts() {};
	//e.g. line_ids = [1,2]. 1 codes for FRET, 2 for PIE (by convention)
	std::vector<int> line_ids; 
	float dwelltime;
	float counttime;//laser rep rate
	long long NumRecords;
	float linestep; //e.g. 10nm px, 2 lines makes linstep 5e-9m
	float pxsize;
} ;

class ph{
public:
	ph() {};
	int tac;
	long long t;
	unsigned char can;
	long long n; 
	float x;//in m
	float y;//in m
	int frame;
};


class imChannel {
public:
	imChannel() :tacmin(0), tacmax(32768), tmin(0), tmax(9223372036854775807) {};
	std::vector<unsigned char> can;
	int tacmin;//microtime range
	int tacmax;
	long long tmin;// macrotime range
	long long tmax; 
	int line_id;//e.g. line_id = 1
	//for some reason this variable gets really slow when accessed from Python
	std::vector<ph> phstream;
};

class imspy {
public:
	imspy() {};
	void ProcessPhotonStream();

	Eigen::ArrayXi tac;
	Eigen::Array<long long, Eigen::Dynamic, 1> t;
	Eigen::Array<unsigned char, Eigen::Dynamic, 1>  can;
	imOpts ImOpts;
	std::vector<imChannel> Channels;
};

#endif
