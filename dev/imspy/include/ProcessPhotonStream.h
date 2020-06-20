//Author: Nicolaas van der Voort
//AG Seidel, HHU Dusseldorf
//June 3, 2020
//TODO: build constructors to take the needed arguments
#ifndef PPS_H
#define PPS_H

#include<vector>
#include <Eigen/Core>

typedef std::vector<long long int> int64vec;

class imOpts{
public:
	imOpts() {};
	//e.g. line_ids = [1,2]. 1 codes for FRET, 2 for PIE (by convention)
	std::vector<int> line_ids; 
	float counttime;//laser rep rate
	long long NumRecords;
	float linestep; //e.g. 10nm px, 2 lines makes linstep 5e-9m
	float dwelltime; //not an image constant, move?
	float pxsize; //not an image const, move?
} ;

/*
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
*/


class imChannel {
public:
	void gentacdecay(int ntacs);
	imChannel() :tacmin(0), tacmax(32768), tmin(0), tmax(9223372036854775807) {};
	std::vector<unsigned char> can;
	int tacmin;//microtime range
	int tacmax;
	long long tmin;// macrotime range
	long long tmax; 
	int line_id;//e.g. line_id = 1
	//for some reason this variable gets really slow when accessed from Python
	//Idea 1: change to Eigen vector like so:
	//Eigen::Vector<ph, Eigen::Dynamic, 1> phstream;
	//Idea 2: get rid of the photon class and just replace it with a series
	//of arrays for each tac, t, can, x, y etc.
	//Then, a protected mechanism is needed to delete from all arrays 
	//simultaneously
private:
	int log2(int n);

};

class imspy {
private:
	//class variables
	//variables are passed in initialization list
	const std::vector<int> tac;
	const int64vec t;
	const std::vector<unsigned char> can;
	const imOpts ImOpts;
	//x, y, frame are calculated in constructor and cannot be const
	std::vector<float> x;
	std::vector<float> y;
	std::vector<int> slice;//t-slice or z-slice

	//member functions
	void calcxyslice();

public:
	//constructor
	imspy(
		std::vector<int> tac,
		int64vec t,
		std::vector<unsigned char> can,
		imOpts Imopts
		) 
		: tac(tac), t(t), can(can), ImOpts(Imopts) {
		x.reserve(ImOpts.NumRecords);
		y.reserve(ImOpts.NumRecords);
		slice.reserve(ImOpts.NumRecords);
		void calcxyslice();
	};


	 //build indexes for each imChannel
	int64vec indexchannels(std::vector<imChannel> Channels);
	void ProcessPhotonStream(); //obsolete


	//get functions
	//Eigen::ArrayXi gettac(int64vec index);
	Eigen::ArrayXi gettac();
	Eigen::Array< long long, Eigen::Dynamic, 1> gett(int64vec index);
	Eigen::Array< unsigned char, Eigen::Dynamic, 1> getcan(int64vec index);
	Eigen::ArrayXf getx(int64vec index);
	Eigen::ArrayXf gety(int64vec index);
	Eigen::ArrayXi getslice(int64vec index);



};

#endif
