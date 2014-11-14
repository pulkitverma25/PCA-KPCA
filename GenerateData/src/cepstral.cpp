/**********************************************************
 * @author  Pulkit Verma
 * @email   technopreneur[dot]pulkit[at]gmail[dot]com
 **********************************************************/

#include<stdio.h>
#include<stdlib.h>
#include <math.h>

#define NO_OF_P 12
#define FRAME_SIZE 320

void calculateHammingWindow(long double* hammingWindow,double piValue,int frameSize)
{
	for(int i=0;i<FRAME_SIZE;i++)
		hammingWindow[i]=0.54-(0.46*cos(2*piValue*i/(frameSize-1))); //Calculates the hamming distance

	return;
}

void performAutocorrelation(long double* hammingWindow,int starting,int ending,char* dcFile, char* riFile)
{
	FILE *fp,*fp2;
	char s[10]; //9 digit amplitude can be accomodated in this string
	int count=0,i=0,linecount=0;
	long double rArray[NO_OF_P+1],sampleValues[FRAME_SIZE],windowedValues[FRAME_SIZE];
	
	if( (fp = fopen(dcFile,"r")) == NULL ) //open the file
	{
		perror("Error while opening the input file.\n");
		exit(-1);
	}	

	if( (fp2 = fopen(riFile,"w")) == NULL ) //open the file
	{
		perror("Error while opening the input file.\n");
		exit(-1);
	}
	
	while( fgets(s,9,fp)  != NULL ) //read the file line by line
	{
		linecount++; //calculates the number of lines read
		if(linecount<=starting*100 || linecount>ending*100)
			continue; //skips the lines other than between starting and ending

		if(count>0) //will shift the values by 80 for next frame
		{ 
			for(i=0;i<FRAME_SIZE*0.75;i++)
					sampleValues[i]=sampleValues[i+80];
			count=0;
		}	
		for(int i=0;i<320;i++)
		{
			hammingWindow[i]=0.54-(0.46*cos(2*3.14*i/(FRAME_SIZE-1)));
		}

		sampleValues[i]=(long double)(atoi(s)); //reads sample value from the file
		
		if(i==319) //if end of frame reached
		{
			for(int n=0;n<FRAME_SIZE;n++) //applies the window on the sample
				windowedValues[n]=sampleValues[n]*hammingWindow[n];

			for(int j=0;j<=NO_OF_P;j++) 
			{
				rArray[j]=0; //initializes the value to zero
				for(int k=0;k+j<FRAME_SIZE;k++) //calculate the R[i]
					rArray[j]+=windowedValues[k]*windowedValues[k+j];
				fprintf(fp2,"%Lf,",rArray[j]); //writes R[i] to the file
			}
			fprintf(fp2,"\n");
			i=0;count=1; //initialize all variables back to zero for the next frame
		}
		i++;
	}

	fclose(fp); //close the file
	fclose(fp2); //close the file
}

void applyDurbins(char *riFile, char *lpcFile) //apply durbins to get alpha values
{
	FILE *fp,*fp2;
	char s[10]; //9 digit amplitude can be accomodated in this string
	int count=0,i=0,sampleValues[320];
	long double rArray[NO_OF_P+1],eArray[NO_OF_P+1],p,kArray[NO_OF_P+1],rootArray[NO_OF_P+2][NO_OF_P+2];      

	for (int k=0;k<=NO_OF_P;k++) //Initialize the array to zero
	{
		rArray[k]=kArray[k]=eArray[k]=0;
		for (int j=0;j<=NO_OF_P;j++)
			rootArray[k][j]=0;
	}

	if( (fp = fopen(riFile,"r")) == NULL ) //open the file
	{
		perror("Error while opening the input file.\n");
		exit(-1);
	}       
	if( (fp2 = fopen(lpcFile,"w")) == NULL ) //open the file
	{
		perror("Error while opening the input file.\n");
		exit(-1);
	}

	while(fscanf(fp,"%Lf,",&rArray[0])>0) //read the file line by line
	{
		for (int k=1;k<=NO_OF_P;k++) //read the R[i] values
			fscanf(fp,"%Lf,",&rArray[k]);
        
		eArray[0]=rArray[0]; //E0=R0
        for (int k=1;k<=NO_OF_P;k++) //Apply the Durbin's algo
        {
			kArray[k]=rArray[k]; 
			if(k>1)
			{
				 for(int j=1;j<=k-1;j++)
					kArray[k]-=(rootArray[j][k-1]*rArray[k-j]);
			}

			kArray[k]/=eArray[k-1];
            rootArray[k][k]=kArray[k];
			eArray[k]=(1-(kArray[k]*kArray[k]))*eArray[k-1];

			if(k>1)
			for (int j=1;j<=k-1;j++)
				rootArray[j][k]=rootArray[j][k-1]-(kArray[k]*rootArray[k-j][k-1]);
        }
        
        for (int k=1;k<=NO_OF_P;k++) //Write the values of alpha in file
                fprintf(fp2,"%Lf,",rootArray[k][NO_OF_P]);
        fprintf(fp2,"\n");
	}
        
	fclose(fp); //close the file
	fclose(fp2); //close the file

	return;
}

void calculateCepstrals(char* lpcFile,char* cepsFile,char* riFile,char* allCepsFile) //calculate the cepstral values
{
	FILE *fp,*fp2,*fp3,*fp4;
	int count=0,i=0,sampleValues[320];
	long double lpcArray[NO_OF_P+1];
	long double cepsArray[NO_OF_P+1],p;
	long double r0Value;
	char s[1000];
			
	for (int k=0;k<=NO_OF_P;k++) //Initialize the array
		lpcArray[k]=cepsArray[k]=0;

	if( (fp = fopen(lpcFile,"r")) == NULL ) //open the file
	{
		perror("Error while opening the input file.\n");
		exit(-1);
	}
	if( (fp2 = fopen(cepsFile,"w")) == NULL ) //open the file
	{
		perror("Error while opening the input file.\n");
		exit(-1);
	}
	if( (fp3 = fopen(riFile,"r")) == NULL ) //open the file
	{
		perror("Error while opening the input file.\n");
		exit(-1);
	}
	if( (fp4 = fopen(allCepsFile,"a")) == NULL ) //open the file
	{
		perror("Error while opening the input file.\n");
		exit(-1);
	}

	while(fscanf(fp,"%Lf,",&lpcArray[0])>0) //read the file line by line
	{
		for (int k=1;k<NO_OF_P;k++) //Read the alpha values
		{
			if(k<NO_OF_P-1)
				fscanf(fp,"%Lf,",&lpcArray[k]);
			else if(k==NO_OF_P-1)
				fscanf(fp,"%Lf,\n",&lpcArray[k]);
		}

		fscanf(fp3,"%Lf,",&r0Value); //read R[0] value
		fgets(s,999,fp3);

		cepsArray[0]=log(r0Value)/ log(2.0); //c0=log(R[0])
	
		for (int k=1;k<=NO_OF_P;k++) //calculate the cepstral values
		{
			cepsArray[k]=lpcArray[k-1];
			if(k>1)
			{
				for(int j=1;j<=k-1;j++)
					cepsArray[k]+=(cepsArray[j]*lpcArray[k-j-1]*j/k);
			}
		}
			
		for (int k=0;k<=NO_OF_P;k++) //write cepstral values to file
		{
			fprintf(fp2,"%Lf,",cepsArray[k]);
			fprintf(fp4,"%Lf,",cepsArray[k]);
		}
		fprintf(fp2,"\n");
		fprintf(fp4,"\n");
	}

	fclose(fp); //close the files
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);

	return;
}