/**********************************************************
 * @author  Pulkit Verma
 * @email   technopreneur[dot]pulkit[at]gmail[dot]com
 **********************************************************/

#ifndef cepstral_H_
#define cepstral_H_

void calculateHammingWindow(long double* hammingWindow,double piValue,int frameSize);
void performAutocorrelation(long double* hammingWindow,int starting,int ending,char* dcFile, char* riFile);
void applyDurbins(char *riFile, char *lpcFile);
void calculateCepstrals(char* lpcFile,char* cepsFile,char* riFile,char* allCepsFile);

#endif