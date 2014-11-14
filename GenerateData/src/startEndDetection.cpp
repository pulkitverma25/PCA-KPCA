/**********************************************************
 * @author  Pulkit Verma
 * @email   technopreneur[dot]pulkit[at]gmail[dot]com
 **********************************************************/

#include<stdio.h>
#include<stdlib.h>

#define CROSSCHECK_RANGE 10
#define FRAME_SIZE 100

void adjustDcShift(char *inputFile,char *outputFile)
{
	FILE *fp,*fp2,*fp3;
	double sampleAvg=0;
	int temp;
	int count=0;
	char s[10]; //9 digit amplitude can be accomodated in this string
	if( (fp = fopen(inputFile,"r")) == NULL ) //open the file
	{
		perror("Error while opening the input file.\n");
		exit(-1);
	}
	while( fgets(s,9,fp)  != NULL ) //read the file line by line
	{
		sampleAvg+=atoi(s); //convert the string into integer
		count++;
	}
	sampleAvg/=count;

	fclose(fp);
	if( (fp3 = fopen(inputFile,"r")) == NULL ) //open the file
	{
		perror("Error while opening the input file.\n");
		exit(-1);
	}
	if( (fp2 = fopen(outputFile,"w")) == NULL ) //open the file
	{
		perror("Error while opening the input file.\n");
		exit(-1);
	}
		while( fgets(s,9,fp3)  != NULL ) //read the file line by line
	{
		temp=atoi(s);
		fprintf(fp2,"%d\n",(int)(temp-sampleAvg));
	}
	fclose(fp3);
	fclose(fp2);
}

void findStartEnd (char *inputFile,int *zcrArray, unsigned long* energyArray,int* starting,int* ending, char* finalOutputFile)
{
	FILE *fp;
	char s[10]; //9 digit amplitude can be accomodated in this string
	long sCurr,sPrev;
	unsigned long energy=0;
	int count=0,i=0;
	

	if( (fp = fopen(inputFile,"r")) == NULL ) //open the file
	{
		perror("Error while opening the input file.\n");
		exit(-1);
	}

	int z=0; //z keeps track of number of samples
	while( fgets(s,9,fp)  != NULL ) //read the file line by line
	{
		if(i>0)
			sPrev=sCurr; //store the previous value to calculate the ZCR

		sCurr=atoi(s); //convert the string into integer
		energy+=(atoi(s)*atoi(s)); //not the energy, will be used to calculate enrgy later
		if(i>0)
		{
			if(sCurr*sPrev<0) //if change of sign then increase the count for ZCR
				count++;
		}
		
		if(i==FRAME_SIZE-1) //if end of frame reached
		{
			energy/=FRAME_SIZE; //energy of the frame
			energyArray[z]=energy; //store energy into the corresponding index on array
			zcrArray[z]=count; //store ZCR into the corresponding index on array
			z++;
			count=0;energy=0;i=0; //initialize all variables back to zero for the next frame
		}
		else{
			i++; //if end of frame not reached then read next sample
		}
	}
	fclose(fp); //close the file

	int p; //p and n will be used in for loops
	unsigned long totalEnergy=0; 

	for(i=CROSSCHECK_RANGE+1;i<z;i++) //to calculate the starting point
	{
		if(i==CROSSCHECK_RANGE+1) //to calculate average energy till the starting point of this loop
			for(int b=0;b<CROSSCHECK_RANGE-1;b++)
				totalEnergy+=energyArray[b];

		totalEnergy+=energyArray[i-1]; 

		//if for any sample there is rise in energy by 20% from previous sample and average energy also
		if((energyArray[i]>1.2*energyArray[i-1]) && i<z-2 && energyArray[i]>(1.2*totalEnergy/i))
		{
			unsigned long tempMax=0; 
			for(int b=i-1; b>i-CROSSCHECK_RANGE;b--) //calculates the maximum energy in previous n samples, here n=CROSSCHECK_RANGE
				if(energyArray[b]>tempMax)
					tempMax=energyArray[b];

			for(p=i;p<i+CROSSCHECK_RANGE;p++) //detects if for any of next n samples energy increase is not 20% from maximum of previous n samples
				if(energyArray[p+1]<1.2*tempMax || energyArray[p+1]<(1.2*totalEnergy/(i)))
					break;

			if(p==i+CROSSCHECK_RANGE) //if energy of each of next n samples > energy of previous n samples by 20%
			{
				*starting=i; //starting point is detected
				break;
			}
		}

	}

	totalEnergy=0; //reset to be used again for calculating ending point
	//to calculate the ending point
	for(i=z-2-CROSSCHECK_RANGE;i>0;i--)
	{
		if(i==z-2-CROSSCHECK_RANGE) //to calculate average energy till the starting point of this loop
			for(int b=z-2-CROSSCHECK_RANGE;b<z-1;b++)
				totalEnergy+=energyArray[b];

		totalEnergy+=energyArray[i+1];

		//if for any sample there is rise in energy by 20% from next sample and average energy also
		if((energyArray[i]>1.2*energyArray[i+1]) && i>0 && energyArray[i]>1.2*totalEnergy/(z-i-1))
		{
			unsigned long tempMax=0; //calculates the maximum energy in next n samples, here n=CROSSCHECK_RANGE
			for(int b=i+1; b<i+CROSSCHECK_RANGE;b++)
				if(energyArray[b]>tempMax)
					tempMax=energyArray[b];

			for(p=i;p>i-CROSSCHECK_RANGE;p--) //detects if for any of previous n samples energy increase is not 20% from maximum of next n samples 
				if(energyArray[p-1]<1.2*tempMax && energyArray[p-1]<(1.2*totalEnergy/(z-i-1)))
					break;

			if(p==i-CROSSCHECK_RANGE) //if energy of each of next n samples > energy of previous n samples by 20%
			{
				*ending=i;
				break;
			}
		}	
	}

	FILE *fp3;
	if( (fp3 = fopen(finalOutputFile,"w")) == NULL ) //open the file
	{
		perror("Error while opening the input file.\n");
		exit(-1);
	}
	fprintf(fp3,"Starting=%d \nEnding=%d\n",*starting,*ending);
	//printf("Starting=%d \nEnding=%d\n",*starting,*ending);
	fclose(fp3);
	return;
}