#ifndef PYWRAP
#define PYWRAP

//get import needed in function declarations
#include <Eigen/Core>
//#include <stdio.h>
#include<vector>
using namespace Eigen;


/*! inverse a numpy matrix
*/
MatrixXd inv(MatrixXd xs);

/*! Substract one integer from another
	\param i an integer
	\param j an integer to subtract from \p i
*/
int subtract(int i, int j);

struct Pet {
	Pet(const std::string &name) : name(name) { }
	void setName(const std::string &name_) { name = name_;}
	const std::string &getName() const { return name; }
	
	std::string name;
};


class imOpts {
	public:
		//e.g. line_ids = [1,2]. 1 codes for FRET, 2 for PIE (by convention)
		std::vector<int> line_ids; 
		float dwelltime;
		float counttime;//laser rep rate
		long long NumRecords;
		float linestep; //e.g. 10nm px, 2 lines makes linstep 5e-9m
		float pxsize;
	imOpts(
		std::vector<int> line_ids,
		float dwelltime,
		float counttime,
		long long NumRecords,
		float linestep,
		float pxsize) 
		: line_ids(line_ids),
		dwelltime(dwelltime),
		counttime(counttime),
		NumRecords(NumRecords),
		linestep(linestep),
		pxsize(pxsize)
		{ };
};

class ph {
	public:
		int tac;
		long long t;
		unsigned char can;
		long long n; 
		float x;//in m
		float y;//in m
		int frame;
	ph (
		int tac,
		long long t,
		unsigned char can,
		long long n,
		float x,
		float y,
		int frame)
		: tac(tac),
		t(t),
		can(can),
		n(n),
		x(x),
		y(y),
		frame(frame)
		{ };
};

class imChannel {
	public:
		std::vector<unsigned char> can;
		int tacmin;//microtime range
		int tacmax;
		long long tmin;// macrotime range
		long long tmax; 
		int line_id;//e.g. line_id = 1
		std::vector<ph> phstream;
	imChannel (
		std::vector<unsigned char> can,
		int tacmin,//microtime range
		int tacmax,
		long long tmin,// macrotime range
		long long tmax, 
		int line_id,//e.g. line_id = 1
		std::vector<ph> phstream)
		: can(can),
		tacmin(tacmin),
		tacmax(tacmax),
		tmin(tmin),
		tmax(tmax),
		line_id(line_id),
		phstream(phstream)
		{ };
};

/*
typedef struct{
	//e.g. line_ids = [1,2]. 1 codes for FRET, 2 for PIE (by convention)
	std::vector<int> line_ids; 
	float dwelltime;
	float counttime;//laser rep rate
	long long NumRecords;
	float linestep; //e.g. 10nm px, 2 lines makes linstep 5e-9m
	float pxsize;
} imOpts;
*/
/*
typedef struct{
	int tac;
	long long t;
	unsigned char can;
	long long n; 
	float x;//in m
	float y;//in m
	int frame;
} ph;

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
*/

#endif