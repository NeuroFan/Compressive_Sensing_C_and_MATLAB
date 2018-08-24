// ConsoleApplication4.cpp : Defines the entry point for the console application.
//  if number_of_iteration is increased the quality of the recovery increased

#include "math.h"
#include <stdio.h>
#include <stdlib.h>

#define M 64
#define N 224
#define lambda 1.35 // lambda is determinded from exprimental data from RICE DSP github
#define number_of_iterations 120

int matMul(float *mat1, float *mat2, float *result, int m, int n, int q);
int matMul_M(float *mat1, float *mat2, float *result, int m, int n, int q);
int matMul_N(float *mat1, float *mat2, float *result, int m, int n, int q);
int transpose(float *mat1, float *mat2, int m, int n);
int matSUB(float *mat1, float *mat2, float *result, int m, int n);
int matSUB(float *mat1, float *mat2, float *result, int m, int n);
int matADD(float *mat1, float *mat2, float *result, int m, int n);
float norm2(float *vec);
int print(float *mat, int m, int n);
float sum_array(float *arr, int m);
float SNR(float *a, float *b, int Length);
float MSE(float *a, float *b, int Length);


float sign(float input)
{
	return (float)((0 < input) - (input < 0));
}

int main()
{
	/*
	Original Matlab Code
	r = y;
	s = zeros(n, 1);
	for iter = 1:iters
		pseudo_data = At*(r)+s;
	sigma_hat = sqrt(1 / m*sum(fabs(r). ^ 2));
	s = (fabs(pseudo_data)> lambda*sigma_hat).*(fabs(pseudo_data) - lambda*sigma_hat).*sign(pseudo_data);
	r = y - A*(s)+1 / m.*r.*length(find(fabs(s)>0));
	end
		x_hat = s;
	end*/

	float A[M][N];
	int i,j,k,K;
	int d;
	for (i = 0; i < M; i++)
		for (j = 0; j < N; j++){
		int temp=0;
		//_TCE_RAND(1,temp);
		A[i][j] = (float)(rand()%3)-1;
		}

	float x[N],  Y[M], At[N][M];
	float middle_data[N], s[N], r[M], temp[M];
	transpose((float *)A, (float *)At, M, N);
	float sigma;
	float Minve = 1 / (float)(M);
	float inv_sqrt_M = 1 / sqrtf(M);
	float inv_sqrt_N = 1 / sqrtf(N);
	float factor=lambda / sqrtf(M);

	//for (i = 1; i < M; i++)
	//	for (j = 1; j < N; j++)
	//		A[i][j] = (rand() % 2); // DO NOT USE THIS FUNCTION

	for (i = 0; i < M; i++)
		for (j = 0; j<N; j++)
			A[i][j] = A[i][j] * inv_sqrt_N;


	for (i = 0; i < N; i++) // only 6 out of N coeff of x are nonzero
		x[i] = 0;


	x[6] = (float)(102);
	x[2] = (float)(12);
	x[3] = (float)(32);



	matMul_N((float *)A,(float *) x,(float *) Y, M, N, 1);

	for (i = 0; i < N; i++)
		s[i] = 0;

	transpose((float *)A, (float *)At, M, N);
	for (i = 0; i < M; i++)
		r[i] = Y[i];

int ii,iii;
printf("Algorithm Starts here:\n");
for (iii=0;iii<10;iii++)
	for (ii = 0; ii < number_of_iterations; ii++){ //AMP
	{
		matMul_M((float *)At, (float *)r, (float *)middle_data, N, M, 1);
		matADD(s, middle_data, s, N, 1);
		sigma = norm2(r) ;
		float astane= sigma* factor;
		int b = 0; 	//b = sum(fabs(s)>0) / m;
		for (i = 0; i < N; i+=8)
		{
			float mediatevalue1=fabs(s[i]) - astane;
			float mediatevalue2=fabs(s[i+1]) - astane;
			float mediatevalue3=fabs(s[i+2]) - astane;
			float mediatevalue4=fabs(s[i+3]) - astane;
			float mediatevalue5=fabs(s[i+4]) - astane;
			float mediatevalue6=fabs(s[i+5]) - astane;
			float mediatevalue7=fabs(s[i+6]) - astane;
			float mediatevalue8=fabs(s[i+7]) - astane;

			if (mediatevalue1>0) {	s[i] =    (mediatevalue1)*sign(s[i]);//s = (abs(pseudo_data)> lambda*sigma_hat).*(abs(pseudo_data) - lambda*sigma_hat).*sign(pseudo_data);
				b = b +1; }else s[i] = 0;

			if (mediatevalue2>0) {	s[i+1] =  (mediatevalue2)*sign(s[i+1]);//s = (abs(pseudo_data)> lambda*sigma_hat).*(abs(pseudo_data) - lambda*sigma_hat).*sign(pseudo_data);
				b = b +1;}else	s[i+1] = 0;

			if (mediatevalue3>0){ 	s[i+2] =  (mediatevalue3)*sign(s[i+2]);//s = (abs(pseudo_data)> lambda*sigma_hat).*(abs(pseudo_data) - lambda*sigma_hat).*sign(pseudo_data);
				b = b +1;}else 	s[i+2] = 0;

			if (mediatevalue4>0) {	s[i+3] =  (mediatevalue4)*sign(s[i+3]);//s = (abs(pseudo_data)> lambda*sigma_hat).*(abs(pseudo_data) - lambda*sigma_hat).*sign(pseudo_data);
				b = b +1;}else 	s[i+3] = 0;

			if (mediatevalue5>0) {	s[i+4] =  (fabs(s[i+4]) - astane)*sign(s[i+4]);//s = (abs(pseudo_data)> lambda*sigma_hat).*(abs(pseudo_data) - lambda*sigma_hat).*sign(pseudo_data);
				b = b +1;}else 	s[i+4] = 0;

			if (mediatevalue6>0) {	s[i+5] =  (fabs(s[i+5]) - astane)*sign(s[i+5]);//s = (abs(pseudo_data)> lambda*sigma_hat).*(abs(pseudo_data) - lambda*sigma_hat).*sign(pseudo_data);
				b = b +1;}else 	s[i+5] = 0;

			if (mediatevalue7>0) {	s[i+6] =  (fabs(s[i+6]) - astane)*sign(s[i+6]);//s = (abs(pseudo_data)> lambda*sigma_hat).*(abs(pseudo_data) - lambda*sigma_hat).*sign(pseudo_data);
				b = b +1;}else 	s[i+6] = 0;

			if (mediatevalue8>0) {	s[i+7] =  (fabs(s[i+7]) - astane)*sign(s[i+7]);//s = (abs(pseudo_data)> lambda*sigma_hat).*(abs(pseudo_data) - lambda*sigma_hat).*sign(pseudo_data);
				b = b +1;}else 	s[i+7] = 0;

		}
		//	r = y - H*s + b.*r;
		b = b * Minve; ;
		matMul_N((float *)A, s, temp, M, N, 1);
		for (i = 0; i < M; i+=4){
			r[i] = Y[i] - temp[i] + b*r[i];
			r[i+1] = Y[i+1] - temp[i+1] + b*r[i+1];
			r[i+2] = Y[i+2] - temp[i+2] + b*r[i+2];
			r[i+3] = Y[i+3] - temp[i+3] + b*r[i+3];
			}

		//if (MSE(s,x,N)<0.1) /////// BREAKING CONDITION
		//	break;

	}
}
	printf("Iteration Number: %d",ii);
	SNR(x,s,N);
 // print(s,1,N);
 // SNR(x,s,N);
    return(-1- 0);
}





int matSUB(float *mat1, float *mat2, float *result, int m, int n) {
	//Calculates mat1_m*n - mat2_m*n and reutrns the resultant matrix to result
	//Matrices should be the same size
	int i,j,k,K;
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++)
			*(result + i*n + j) = *(mat1 + i*n + j) - *(mat2 + i*n + j);
	}
	return(0);
}


int matADD(float *mat1, float *mat2, float *result, int m, int n) {
	//Calculates mat1_m*n + mat2_m*n and reutrns the resultant matrix to result
	//Matrices should be the same size
	int i,j,k,K;
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++)
			*(result + i*n + j) = *(mat1 + i*n + j) + *(mat2 + i*n + j);
	}
	return(0);
}







int matMul(float *mat1, float *mat2, float *result, int m, int n, int q) {
	//Multiplication  for Mat1_(m*n) * Mat2_(n*k)
	// M,N and K are dimentions of matrices
	// Output is going to sit in the matrix that result represents here
	float sum;
	int i,j,k,K;
	for (i = 0; i < m; i++) {
		for (j = 0; j < q; j++) {
			sum = 0;
			for (k = 0; k < n; k++) {
				sum = sum + *(mat1 + i*n + k) * *(mat2 + k*q + j);
			}
			*(result + i*q + j) = sum;
		}
	}
	return(0);
}

int matMul_M( float *mat1, float *mat2, float *result, int m, int n, int q) {
	float sum1,sum2,sum3,sum4;
	int i,j,k,K;
	for (i = 0; i < m; i++) {
		for (j = 0; j < q; j++) {
			sum1 = 0;sum2 = 0;sum3 = 0;sum4 = 0;
			float* index_t1=mat1 + i*n;
			float* index_t2=mat2 + j;
			for (k = 0; k < M; k+=4) {
				sum1 = sum1 + *(index_t1 + k) * *(index_t2 + k*q);	//sum = sum + *(mat1 + i*n + k) * *(mat2 + k*q + j);
				sum2 = sum2 + *(index_t1 + k+1) * *(index_t2 + (k+1)*q);
				sum3 = sum3 + *(index_t1 + k+2) * *(index_t2 + (k+2)*q);
				sum4 = sum4 + *(index_t1 + k+3) * *(index_t2 + (k+3)*q);
			}
			*(result + i*q + j) = sum1 + sum2+sum3+sum4;
		}
	}
	return(0);
}

int matMul_N( float *mat1, float *mat2, float *result, int m, int n, int q) {
	float sum1,sum2,sum3,sum4;
	int i,j,k,K;
	for (i = 0; i < m; i++) {
		for (j = 0; j < q; j++) {
			sum1 = 0;sum2 = 0;sum3 = 0;sum4 = 0;
			float* index_t1=mat1 + i*n;
			float* index_t2=mat2 + j;
			for (k = 0; k < N; k+=4) {
				sum1 = sum1 + *(index_t1 + k) * *(index_t2 + k*q);	//sum = sum + *(mat1 + i*n + k) * *(mat2 + k*q + j);
				sum2 = sum2 + *(index_t1 + k+1) * *(index_t2 + (k+1)*q);
				sum3 = sum3 + *(index_t1 + k+2) * *(index_t2 + (k+2)*q);
				sum4 = sum4 + *(index_t1 + k+3) * *(index_t2 + (k+3)*q);
			}
			*(result + i*q + j) = sum1+sum2+sum3+sum4;
		}
	}
	return(0);
}

int transpose(float *mat1, float *mat2, int m, int n)
{
	//transposes mat1_m*n which is input to mat2_n*m which is out
	// M,N  dimentions of matrices
	int i,j,k,K;
	for (i = 0; i<m; i++)
		for (j = 0; j < n; j++)
		{
			*(mat2 + j*m + i) = *(mat1 + i*n + j);
		}
	return(0);
}



float norm2(float *vec) { //Calculate norm 2 of input vector, since all vectors are M length here, for loop unrolling loops are defined over constant M
	float sum=0;
	int i,j,k,K;
	for (i = 0; i < M; i++)
	{
		sum += *(vec + i) * *(vec + i);
	}
	return sqrtf(sum);
}


int print(float *mat, int m, int n) {
	//  Prints mat with dimentions of m and n
	int i,j,k,K;
	for (i = 0;i < m;i++) {
		for (j = 0;j < n;j++){
		//_TCE_STREAM_OUT((char) *(mat+i*n+j));
		printf("%f  ",*(mat+i*n+j));
				  }
	printf("\n");
	}
	return(1);
}



float sum_array(float *arr, int m) {
	float sum = 0;
	int i,j,k,K;
	for (i = 0; i < m; i++)
		sum += *(arr + i);
	return(sum);
}

float MSE(float *a, float *b, int Length)
{
	int l;
	float temp;
	float mse = 0;
	for (l = 0; l<Length; l++)
	{
		temp = a[l] - b[l];
		mse += temp*temp;
	}
	return mse;
}

float SNR(float *a, float *b, int Length)
{
	int l,i;
	int temp;
	float mse = MSE(a, b, Length);
	float signal_power = 0;
	for (i = 0; i<Length; i++)
		signal_power += a[i] * a[i];

	printf("mse : %d \n", (int)mse);
	printf("signal_power : %d \n",(int) signal_power);
	float snr = 10 * log10(signal_power / mse);
	printf("SNR : %d \n", (int)snr);
	return snr;
}
