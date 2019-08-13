
#include <stdio.h>

extern "C" int SplitOnTacs(long long * eventN,
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
                          int * imB){
    int i, startline, pixel;
    int j = 0, line = 0;
    int inscanpxl = 0, unaccounted = 0, linemarkers = 0, outscan = 0, GChannel = 0;
    int badcalculatedpxl=0, gatedpxl =0;
    bool inscan(false);
    float macrotime2pixel = counttime / dwelltime;
    int timedif;
    //printf("macrotime2pixel is %d \n", (int)(macrotime2pixel*timedif));
    
    for(i = 0; i < NumRecords; i++){
        if(can[i] == 65){
            inscan = true;
            startline = t[i];
            linemarkers++;
        } else if( can[i] == 66){
            inscan = false;
            j++;
            linemarkers ++;
            if (j >= nlines) {
                line ++;
                j = 0;
            }
        } else if ((can[i] == 1 or can[i] == 3) and inscan and uselines[j] and tac[i] > gate) {
            //image is stored as 1D array
            //macrotime2pixel is calculated once for computational efficiency
            
            timedif = t[i] - startline; //cast long long to int
            pixel = timedif * macrotime2pixel; 
            if (pixel < dimX){
                if (tac[i] % 2 == 0) imA[line * dimX + pixel]++;
                else if (tac[i] % 2 == 1) imB[line * dimX + pixel]++;
                inscanpxl ++;
            } else badcalculatedpxl++;
        } else if ((can[i] == 0 or can[i] == 2)) GChannel++;
        else if (not inscan) outscan++;
        else if (tac[i] < gate) gatedpxl ++;
        else unaccounted ++;
    }
    printf("inscan counts is %d \n", inscanpxl);
    printf("linemarker count is %d \n", linemarkers);
    printf("pixels recorded while outside of scan is %d \n", outscan);
    printf("pixels recorded in the Green channel is %d \n", GChannel);
    printf("bad calculated pixels is %d \n", badcalculatedpxl);
    printf("number of gated pixels is %d \n", gatedpxl);
    printf("unaccounted pixels is %d \n", unaccounted);
    printf("photon bookkeeping: %d photons remaining\n", NumRecords - inscanpxl - linemarkers - unaccounted - 
           outscan - GChannel - badcalculatedpxl - gatedpxl);
    return 0;
}