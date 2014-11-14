/**********************************************************
 * @author  Pulkit Verma
 * @email   technopreneur[dot]pulkit[at]gmail[dot]com

 * Dependency: Requires Eigen Library
			   Requires gnuplot for graph plotting
 **********************************************************/

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../lib/Eigen/SVD"
#include "../lib/Eigen/Eigenvalues"
#include "../lib/unsupported/Eigen/MatrixFunctions"

using namespace Eigen;
using namespace std;

void getMatrix(char* fileName,int dimension, int exampleCount, float* matrixX, int* exampleClass);
VectorXf getEigenValues(MatrixXf M, int dimension, int exampleCount) ;
MatrixXf getEigenVectors(MatrixXf M) ;
MatrixXf generateKernelMatrix (MatrixXf mT, int kernelNo, int dimension, int exampleCount, int power, int sigma);
void plotEigenValueVariation(VectorXf eigValue,int dimension, int kernelNo);
MatrixXf centerKernelMatrix(MatrixXf kernelMatrix, int exampleCount);
MatrixXf calculateAlpha(MatrixXf eigVector, VectorXf eigValue, int exampleCount,int dimension);
MatrixXf findTransformedData(MatrixXf kernelMatrix, MatrixXf alpha);

int main(int argc, char *argv[])
{
	if ( argc != 8 )
	{
		cout<<"usage: "<<argv[0]<<" filename dimesnion exampleCount alpha noOfDim power sigma"<<endl;
		exit(-1);
	}
			
	int exampleCount=atoi(argv[3]);
	int dimension=atoi(argv[2]);
	float alpha = atof(argv[4]);
	int errorDimension=atoi(argv[5]);
	int power=atoi(argv[6]);
	int sigma=atoi(argv[7]);
	float *matrixX= new float[exampleCount*dimension];
	int* exampleClass=new int[exampleCount];

	getMatrix(argv[1],dimension,exampleCount,matrixX,exampleClass);
	
	Map<MatrixXf> m(matrixX,dimension,exampleCount); 
	MatrixXf mT = m.transpose();
	
	MatrixXf kernelMatrix1 = generateKernelMatrix(mT,1,dimension,exampleCount,power,sigma);
	MatrixXf kernelMatrix2 = generateKernelMatrix(mT,2,dimension,exampleCount,power,sigma);

	cout<<"\n\nKernel 1: "<<endl<<kernelMatrix1<<endl;
	cout<<"\n\nKernel 2: "<<endl<<kernelMatrix2<<endl;

	MatrixXf centeredKernelMatrix1 = centerKernelMatrix(kernelMatrix1,exampleCount);
	MatrixXf centeredKernelMatrix2 = centerKernelMatrix(kernelMatrix2,exampleCount);

	cout<<"\n\nCentered Kernel Matrix 1: "<<endl<<centeredKernelMatrix1<<endl;
	cout<<"\n\nCentered Kernel Matrix 2: "<<endl<<centeredKernelMatrix2<<endl;

	VectorXf eigValues1 = getEigenValues(centeredKernelMatrix1, dimension, exampleCount);
	//eigValues1 /= exampleCount;
	VectorXf eigValues2 = getEigenValues(centeredKernelMatrix2, dimension, exampleCount);
	//eigValues2 /= exampleCount;

	cout<<"\n\nEigValues 1: "<<endl<<eigValues1<<endl;
	cout<<"\n\nEigValues 2: "<<endl<<eigValues2<<endl;

	plotEigenValueVariation(eigValues1, exampleCount,1);
	plotEigenValueVariation(eigValues2, exampleCount,2);

	MatrixXf eigVector1 = getEigenVectors(centeredKernelMatrix1);
	MatrixXf eigVector2 = getEigenVectors(centeredKernelMatrix2);

	cout<<"\n\nEigen Vector 1: "<<endl<<eigVector1<<endl;
	cout<<"\n\nEigen Vector 2: "<<endl<<eigVector2<<endl;

	MatrixXf alpha1 = calculateAlpha( eigVector1, eigValues1, exampleCount,dimension);
	MatrixXf alpha2 = calculateAlpha( eigVector2, eigValues2, exampleCount,dimension);

	MatrixXf principalComponentMatrix1 = findTransformedData( centeredKernelMatrix1,  alpha1);
	MatrixXf principalComponentMatrix2 = findTransformedData( centeredKernelMatrix2,  alpha2);

	cout<<"Principal Component 1: "<<endl<<principalComponentMatrix1<<endl;
	cout<<"Principal Component 2: "<<endl<<principalComponentMatrix2<<endl;

	return 0;
}

MatrixXf findTransformedData(MatrixXf kernelMatrix, MatrixXf alpha)
{
	MatrixXf principalComponentMatrix = (alpha.transpose() * kernelMatrix).transpose();

	return principalComponentMatrix;
}

MatrixXf calculateAlpha(MatrixXf eigVector, VectorXf eigValue, int exampleCount,int dimension)
{
	MatrixXf alpha(exampleCount,dimension);
	for(int i=0;i<dimension;i++)
		alpha.col(i)=(1/sqrt(eigValue(i)*exampleCount))*eigVector.col(i);

	return alpha;
}

MatrixXf centerKernelMatrix(MatrixXf kernelMatrix, int exampleCount)
{
		MatrixXf matrixD = kernelMatrix.colwise().mean();
		float matrixE = matrixD.mean();

		MatrixXf oneMatrix(exampleCount,1);
		oneMatrix.fill(1);

		MatrixXf oneMatrix2(exampleCount,exampleCount);
		oneMatrix2.fill(1);

		MatrixXf matrixJ = oneMatrix * matrixD;
		MatrixXf centeredKernelMatrix = kernelMatrix - matrixJ - matrixJ.transpose() + (matrixE*oneMatrix2);

		return centeredKernelMatrix;
}

void getMatrix(char* fileName,int dimension, int exampleCount, float* matrixX, int* exampleClass)
{
	FILE *filePointer;
    if ((filePointer = fopen(fileName,"r"))== NULL ) 
    {
		cout<<"Could not open file"<<endl;
		exit(-1);	
    }

	for(int i=0;i<exampleCount;i++)
		for(int j=0;j<dimension+1;j++)
			if(j<dimension)
				fscanf(filePointer,"%f,",&matrixX[(i*dimension)+j]);
			else
				fscanf(filePointer,"%d\n",&exampleClass[i]);

	fclose(filePointer);

	return;	
}

MatrixXf generateKernelMatrix (MatrixXf mT, int kernelNo, int dimension, int exampleCount, int power, int sigma)
{
	MatrixXf kernelMatrix(exampleCount,exampleCount);
	kernelMatrix.fill(0);
	//MatrixXf m = mT.transpose();
	for(int k=0;k<exampleCount;k++)
	{
		for(int i=0;i<exampleCount;i++)
		{
			
			if(kernelNo==1)
			{
				long double dp=(mT.row(k))*(mT.row(i).transpose());
				kernelMatrix(k,i)=pow(dp+1,power);
			}
			else if (kernelNo==2)
			{
				long double dp=(mT.row(k)-mT.row(i)).norm();
				kernelMatrix(k,i)=exp(-1*pow(dp,2)/(2*pow(sigma,2)));
			}


		}
	}

	return kernelMatrix;
}

MatrixXf getEigenVectors(MatrixXf M) 
{
	JacobiSVD<MatrixXf> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
	MatrixXf eigVector = svd.matrixV();

	cout<<"\n\nReturning"<<endl;
	return eigVector;
}


VectorXf getEigenValues(MatrixXf M, int dimension, int exampleCount) 
{
	JacobiSVD<MatrixXf> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
	VectorXf eigValue = svd.singularValues();
	
	for(int i=0;i<dimension;i++)
		eigValue(i)= eigValue(i)/(exampleCount-1);

	return eigValue;
}

void plotEigenValueVariation(VectorXf eigValue,int dimension, int kernelNo)
{
	VectorXf v1(dimension);
	for(int i=0;i<dimension;i++)
		v1(i)=i+1;

	for(int i=0;i<dimension;i++)
		eigValue(i) = log(eigValue(i));

	MatrixXf mat1(dimension,2);
	mat1<< v1,eigValue;
	

	ofstream myfile ("../IOFiles/eigenVal.dat");
	if (myfile.is_open())
	{
		myfile << mat1 << endl;
		myfile.close();
	}
	else 
		cout << "Unable to open file";

	FILE *pipe = popen("gnuplot -persist 2>/dev/null", "w");
	fprintf(pipe, "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5\n");
	if(kernelNo==1)
		fprintf(pipe,"set title \"Eigenvalue Variation : Polynomial Kernel\" \n");
	else
		fprintf(pipe,"set title \"Eigenvalue Variation : Radial Basis Kernel\" \n");
	fprintf(pipe,"set xlabel \"Eigenvalue Number\" \n");
	fprintf(pipe,"set ylabel \"log(eigenvalue) (base 10)\" \n");
	fprintf(pipe, "plot \"../IOFiles/eigenVal.dat\" using 1:2::xtic(1) with linespoints ls 1\n");
	pclose(pipe);

	return;
}
