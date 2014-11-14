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
#include "../lib/Eigen/SVD"
#include "../lib/Eigen/Eigenvalues"
#include "../lib/unsupported/Eigen/MatrixFunctions"


using namespace Eigen;
using namespace std;

void getMatrix(char* fileName,int dimension, int exampleCount, float* matrixX, int* exampleClass);
void plotBiplotMethod1(MatrixXf mT, float power, int dimension, int exampleCount);
int getMatrixRank (MatrixXf A);
MatrixXf getEigenVectorMatrix(MatrixXf mT);
void reconstructDataPoints(MatrixXf matrixV, MatrixXf mT, int exampleCount, int dimension, int errorDimension);
void plotEigenValueVariation(MatrixXf mT,int dimension,int exampleCount);
VectorXf getEigenValues(MatrixXf M);

int main(int argc, char *argv[])
{
	if ( argc != 6 )
	{
		cout<<"usage: "<<argv[0]<<" filename dimesnion exampleCount alpha noOfDim"<<endl;
		exit(-1);
	}
			
	int exampleCount=atoi(argv[3]);
	int dimension=atoi(argv[2]);
	float power = atof(argv[4]);
	int errorDimension=atoi(argv[5]);
	float *matrixX= new float[exampleCount*dimension];
	int* exampleClass=new int[exampleCount];

	getMatrix(argv[1],dimension,exampleCount,matrixX,exampleClass);
	
	
	Map<MatrixXf> m(matrixX,dimension,exampleCount); 
	MatrixXf mT = m.transpose();
	
	plotEigenValueVariation(mT,dimension,exampleCount);
	plotBiplotMethod1(mT,power, dimension,exampleCount);
	MatrixXf vMat = getEigenVectorMatrix(mT);
	reconstructDataPoints(vMat, mT, exampleCount, dimension, errorDimension);
	

	return 0;
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

MatrixXf getEigenVectorMatrix(MatrixXf mT)
{
	MatrixXf centered = mT.rowwise() - mT.colwise().mean();
	MatrixXf cov = (centered.adjoint() * centered);

	JacobiSVD<MatrixXf> svd(cov, Eigen::ComputeThinU | Eigen::ComputeThinV);
	return svd.matrixV();

}
void plotBiplotMethod1(MatrixXf mT, float power, int dimension, int exampleCount)
{

	MatrixXf centered = mT.rowwise() - mT.colwise().mean();

	JacobiSVD<MatrixXf> svd(centered, Eigen::ComputeThinU | Eigen::ComputeThinV);	
	VectorXf powAlphaSingVal = svd.singularValues().array().pow(power);
	VectorXf powAlphaSingVal2 = svd.singularValues().array().pow(1-power);

	MatrixXf l1 = powAlphaSingVal.asDiagonal();
	MatrixXf l2 = powAlphaSingVal2.asDiagonal();
	MatrixXf g = svd.matrixU()*(l1);
	MatrixXf h = (l2)*(svd.matrixV().transpose());

	cout << "\n\n[Method 1] G matrix :" << endl << g << endl;
	cout << "\n\n[Method 1] H matrix :" << endl << h << endl;

	VectorXf v1(dimension);
	v1.fill(0);
	VectorXf v2(dimension);
	v2.fill(0);

	MatrixXf mat1(exampleCount,2);
	MatrixXf mat2(dimension,4);
	mat1<< g.col(0),g.col(1);
	mat2<< v1,v2,h.col(0),h.col(1);	

	ofstream myfile ("../IOFiles/gMat.dat");
	if (myfile.is_open())
	{
		myfile << mat1 << endl;
		myfile.close();
	}
	else 
		cout << "Unable to open file";

	ofstream myfile2 ("../IOFiles/hMat.dat");
	if (myfile2.is_open())
	{
		myfile2 <<mat2<<endl;
		myfile2.close();
	}
	else 
		cout << "Unable to open file";	

	FILE *pipe = popen("gnuplot -persist 2>/dev/null", "w");
	fprintf(pipe, "set term x11 enhanced \n");
	fprintf(pipe, "set style arrow 1 filled head lw 1 linecolor rgb \"red\" \n");
	fprintf(pipe, "set style line 2 linecolor rgb \"blue\" \n");
	fprintf(pipe, "plot \"../IOFiles/hMat.dat\" using 1:2:3:4 with vectors arrowstyle 1, \"../IOFiles/gMat.dat\" using 1:2 ls 2\n");
	pclose(pipe);

	return;
}

int getMatrixRank (MatrixXf A)
{
	FullPivLU<MatrixXf> luA(A);
	return (luA.rank());
}

void reconstructDataPoints(MatrixXf matrixV, MatrixXf mT, int exampleCount, int dimension, int errorDimension)
{
	cout<<"\n\nError array:"<<endl;
	MatrixXf centered = mT.rowwise() - mT.colwise().mean();

	int rank = getMatrixRank(matrixV);
	long  double* errorArray=new long  double[exampleCount];
	int* tempError = new int [dimension];
	for(int j=0;j<exampleCount;j++)
	{
		VectorXf tempVect (dimension);
		tempVect.fill(0);
		for(int i=0;i<errorDimension;i++)
		{
			VectorXf normalizedVect = matrixV.col(i).normalized();
			double dp=(normalizedVect.transpose())*(centered.row(j).transpose());

			tempVect = tempVect + (dp*normalizedVect);
		}

		errorArray[j]= (tempVect - centered.row(j).transpose()).norm();
		cout<<"error["<<j<<"] = "<<errorArray[j]<<endl;
	}


	ofstream myfile ("../IOFiles/error.dat");
	if (myfile.is_open())
	{
		for(int i=0; i<exampleCount;i++)
			myfile << i<<" "<<i+1<<" "<< errorArray[i]<< endl;
		myfile.close();
	}
	else 
		cout << "Unable to open file";

	FILE *pipe = popen("gnuplot -persist 2>/dev/null", "w");
	fprintf(pipe, "set term x11 enhanced \n");
	fprintf(pipe, "set boxwidth 0.5 \n");
	fprintf(pipe, "set style fill solid \n");
	fprintf(pipe, "plot \"../IOFiles/error.dat\" using 1:3:xtic(2) with boxes \n");
	pclose(pipe);

	return;
}

VectorXf getEigenValues(MatrixXf M) {
    SelfAdjointEigenSolver<MatrixXf> es(M);
    return es.eigenvalues();
}

void plotEigenValueVariation(MatrixXf mT,int dimension,int exampleCount)
{

	VectorXf v1(dimension);
	for(int i=0;i<dimension;i++)
		v1(i)=i+1;


	MatrixXf centered = mT.rowwise() - mT.colwise().mean();
	MatrixXf cov = (centered.adjoint() * centered);// / double(mT.rows());

	JacobiSVD<MatrixXf> svd(cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

	VectorXf eigValue = svd.singularValues();
	
	
	for(int i=0;i<dimension;i++)
		eigValue(i) = eigValue(i)/(exampleCount-1);
	cout<<"\n\nEigenvalues:"<<endl<<eigValue<<endl;
	MatrixXf mat1(dimension,2);
	mat1<< v1,eigValue;
	

	ofstream myfile ("../IOFiles/eigenValPCA.dat");
	if (myfile.is_open())
	{
		myfile << mat1 << endl;
		myfile.close();
	}
	else 
		cout << "Unable to open file";

	FILE *pipe = popen("gnuplot -persist 2>/dev/null", "w");
	fprintf(pipe, "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5\n");
	fprintf(pipe, "plot \"../IOFiles/eigenValPCA.dat\" using 1:2::xtic(1) with linespoints ls 1\n");
	pclose(pipe);

	return;
}
