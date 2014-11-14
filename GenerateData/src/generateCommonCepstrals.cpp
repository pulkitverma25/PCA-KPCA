/**********************************************************
 * @author  Pulkit Verma
 * @email   technopreneur[dot]pulkit[at]gmail[dot]com
 **********************************************************/

#include <stdio.h>
#include <string.h>
#include <cstdlib>
#include <fstream>
#include <dirent.h>
#include <sys/stat.h>

using namespace std;

#include "startEndDetection.h"
#include "cepstral.h"

#define FRAME_SIZE 320
#define FILENAME_SIZE 50
#define CONFIG_FILE "../IOFiles/createCodebook.config"
#define PATH "../IOFiles/Samples"

void readConfigFile(char* allCepsFile,char* inputFile,char* dcFile,char* riFile,char* lpcFile,char* cepsFile,char* startEndLogFile,int* frameSize, int* noOfP, double* piValue);
char* mystrsep(char** stringp, const char* delim);
void initializeFile(char* vqcFile, char* distOpFile,char* logFile);
void renameFile(char* tempFile,char* fileName,char* token,char* extension);

int main()
{
	char inputFile[FILENAME_SIZE],dcFile[FILENAME_SIZE], riFile[FILENAME_SIZE], lpcFile[FILENAME_SIZE], cepsFile[FILENAME_SIZE],startEndLogFile[FILENAME_SIZE],allCepsFile[FILENAME_SIZE];
	char tempInputFile[FILENAME_SIZE],tempDcFile[FILENAME_SIZE], tempRiFile[FILENAME_SIZE], tempLpcFile[FILENAME_SIZE], tempCepsFile[FILENAME_SIZE],tempStartEndLogFile[FILENAME_SIZE];
	int frameSize, noOfP;
	double piValue;


	char vqcFile[FILENAME_SIZE],tvFile[FILENAME_SIZE],distOpFile[FILENAME_SIZE],logFile[FILENAME_SIZE];
	int vqcSize,z,*distance;
	float sParameter,distortionPercent;
	distance=(int*)malloc(sizeof(int)*12);

	readConfigFile(allCepsFile,inputFile,dcFile, riFile,lpcFile, cepsFile,startEndLogFile, &frameSize,  &noOfP,  &piValue);

	long double hammingWindow[FRAME_SIZE];
	int zcrArray[100000],starting=0,ending=0,i=0;
	unsigned long energyArray[100000]; //unsigned int used to accomodate large values of energy
	for(int i=0;i<FRAME_SIZE;i++)
		hammingWindow[i]=0;

	FILE *fp;
	if( (fp = fopen(allCepsFile,"w")) == NULL ) //open the file
	{
		perror("Error while opening the input file.\n");
		exit(-1);
	}
	fclose(fp);
	
	ifstream fin;
	string dir, filepath;
	int num;
	DIR *dp;
	struct dirent *dirp;
	struct stat filestat;

	dir=PATH;

	dp = opendir( dir.c_str() );
	if (dp == NULL)
	{
		printf("Error opening %s",dir.c_str()); 
		return -1;
    }

	while ((dirp = readdir( dp )))
	{
		filepath = dir + "/" + dirp->d_name;

		// If the file is a directory (or is in some way invalid) we'll skip it 
		if (stat( filepath.c_str(), &filestat )) continue;
		if (S_ISDIR( filestat.st_mode ))         continue;

		char* token;
		char* string;
		char* tofree;

		string = strdup(dirp->d_name);

		if (string != NULL) 
		{
			while ((token = mystrsep(&string, ".")) != NULL)
			{
				renameFile(tempInputFile,inputFile,token,".txt");
				renameFile(tempDcFile,dcFile,token,"DC.txt");
				renameFile(tempRiFile,riFile,token,"RI.csv");
				renameFile(tempLpcFile,lpcFile,token,"LPC.csv");
				renameFile(tempCepsFile,cepsFile,token,"CEPS.csv");
				renameFile(tempStartEndLogFile,startEndLogFile,token,"SE.txt");
				
				starting=ending=0;
				adjustDcShift(tempInputFile,tempDcFile);
				calculateHammingWindow(hammingWindow,piValue,frameSize);
				findStartEnd(tempDcFile,zcrArray,energyArray,&starting,&ending,tempStartEndLogFile);
				performAutocorrelation(hammingWindow,starting,ending,tempDcFile,tempRiFile);
				applyDurbins(tempRiFile,tempLpcFile);
				calculateCepstrals(tempLpcFile,tempCepsFile,tempRiFile,allCepsFile);
				
				break;
			}
		}
    }

	closedir( dp );
	free(distance);

	return 0;
}


void readConfigFile(char* allCepsFile,char* inputFile,char* dcFile,char* riFile,char* lpcFile,char* cepsFile,char* startEndLogFile,int* frameSize, int* noOfP, double* piValue)
{
	FILE *fp;
	if( (fp = fopen(CONFIG_FILE,"r")) == NULL ) //open the file
	{
		perror("Error while opening the input file.\n");
		exit(-1);
	}
	char varName[FILENAME_SIZE];

	while(fscanf(fp,"%s ",varName)>0)
	{
		if(!strcmp(varName,"INPUT_FILE"))
			fscanf(fp,"%s\n",inputFile);

		else if(!strcmp(varName,"ALL_CEPS"))
			fscanf(fp,"%s\n",allCepsFile);

		else if(!strcmp(varName,"DC_FILE"))
			fscanf(fp,"%s\n",dcFile);

		else if(!strcmp(varName,"RI_FILE"))
			fscanf(fp,"%s\n",riFile);

		else if(!strcmp(varName,"LPC_FILE"))
			fscanf(fp,"%s\n",lpcFile);

		else if(!strcmp(varName,"CEPS_FILE"))
			fscanf(fp,"%s\n",cepsFile);

		else if(!strcmp(varName,"SELOG_FILE"))
			fscanf(fp,"%s\n",startEndLogFile);

		else if(!strcmp(varName,"FRAME_SIZE"))
			fscanf(fp,"%d\n",frameSize);

		else if(!strcmp(varName,"NO_OF_P"))
			fscanf(fp,"%d\n",noOfP);

		else if(!strcmp(varName,"PI"))
			fscanf(fp,"%d\n",piValue);
	}

	fclose(fp);
}

char* mystrsep(char** stringp, const char* delim)
{
	char* start = *stringp;
	char* p;
 	p = (start != NULL) ? strpbrk(start, delim) : NULL;

	if (p == NULL)
    	*stringp = NULL;
	else
  	{
 		*p = '\0';
		*stringp = p + 1;
	}

	return start;
}

//Ensures that files are empty when initialized
void initializeFile(char* vqcFile, char* distOpFile,char* logFile)
{
	FILE *fp,*fp2,*fp3;
	if( (fp = fopen(vqcFile,"w")) == NULL ) //open the file
	{
		perror("Error while opening the input file.\n");
		exit(-1);
	}
	fclose(fp);

	if( (fp2 = fopen(distOpFile,"w")) == NULL ) //open the file
	{
		perror("Error while opening the input file.\n");
		exit(-1);
	}
	fclose(fp2);

	if( (fp3 = fopen(logFile,"w")) == NULL ) //open the file
	{
		perror("Error while opening the input file.\n");
		exit(-1);
	}
	fclose(fp3);
}

void renameFile(char* tempFile,char* fileName,char* token,char* extension)
{
	memset(tempFile, 0, FILENAME_SIZE);
	strcat(tempFile,fileName);
	strcat(tempFile,token);
	strcat(tempFile,extension);
	return;
}
