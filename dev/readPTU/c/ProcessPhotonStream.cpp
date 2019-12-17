#include <stdio.h>

extern "C"
{
	__declspec(dllexport) int SplitOnTacs(long long * eventN,
		int * tac,
		long long * t,
		unsigned __int8 * can,
		int dimX,
		int dimY,
		float dwelltime,
		float counttime,
		int NumRecords,
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
				//macrotime2pixel is calculated once for computational efficiency
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

	__declspec(dllexport) int genGRYlifetime(
		long long * eventN,
		int * tac,
		long long * t,
		unsigned __int8 * can,
		int dimX,
		int dimY,
		int ntacs,
		int tac_range,
		float dwelltime,
		float counttime,
		int NumRecords,
		int nlines,
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
			//currently no frame separation is implemented
			else if (can[i] == 68) {
				line = 0;
				frame++;
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

