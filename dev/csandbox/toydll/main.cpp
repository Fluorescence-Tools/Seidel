#include <stdio.h>
#include <string.h>

// fitting options
typedef struct Options{
	char mode[50];
	int threshold; 			// threshold for peak search
} Option;

int printfunc(
	char * filename,
	long long length,
	long long * eventN, 
	int * tac, 
	unsigned __int8 * can);

int printStruct ( Option * option);

int main(int argc, char **argv){
	Option option1;

	strcpy(option1.mode, "normal");
	option1.threshold = 5;
	int NumRecords = 5921583;
	printf("hello world\n");
	printStruct ( &option1 );
	char filename[100] = "mijnnaam.ptu";
	char * Pfilename = filename;
	long long length = 100;
//	long long* pteventN = new long long int[NumRecords];
//	printfunc(Pfilename, length, pteventN);
//	printf("eventN entry is %lld \n", pteventN[NumRecords-1]);
	return 0;
}

int printStruct ( Option * option){
	printf("mode is %s and threshold is %d \n",option->mode, option->threshold);
	return 0;
}

int printfunc(char * filename, long long length, long long * eventN, int * tac, unsigned __int8 * can)
{
	long long int i = 0;
	printf("filename is: %s\n", filename);
	printf("length is: %lld \n", length);
	for(i = 0; i<length; i++)
	{
		eventN[i] = i+10;
		tac[i] = i;
		can[i] = i;
	}
	return i;
}