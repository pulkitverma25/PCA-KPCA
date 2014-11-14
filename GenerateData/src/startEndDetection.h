/**********************************************************
 * @author  Pulkit Verma
 * @email   technopreneur[dot]pulkit[at]gmail[dot]com
 **********************************************************/

#ifndef startEndDetection_H_
#define startEndDetection_H_

void adjustDcShift(char *inputFile,char *outputFile);
void findStartEnd (char *inputFile,int *zcrArray, unsigned long* energyArray,int* starting,int* ending, char* finalOutputFile);

#endif