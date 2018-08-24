// ConsoleApplication6.cpp : Defines the entry point for the console application.
//



#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#define N 256
#define M 64
#define number_of_iteration 100 //Number of Iteration is usually higher in IHT
#define sparsity 4

int swap(int *xp, int *yp);
int bubbleSort(int *arr, int n);
int matMul(float *mat1, float *mat2, float *result, int m, int n, int q);
int transpose(float *mat1, float *mat2, int m, int n);
int matSUB_M(float *mat1, float *mat2, float *result, int m, int n);
int matADD_N(float *mat1, float *mat2, float *result, int m, int n);
//int print(float *mat, int m, int n);
int sorted_order(float *arr, int *sorted_ind, int n, int k);
int matMul_N_in_Middle(float *mat1, float *mat2, float *result, int m, int n, int q);
int matMul_M_in_Middle(float *mat1, float *mat2, float *result, int m, int n, int q);
float SNR(float *a, float *b, int Length);
float MSE(float *a, float *b, int Length);

int main()
{

	float A[M][N], x[N], y[M], x_hat[N];
	float step = 0.5;
	int i,j;
	for (i = 0; i < M; i++)
		for (j = 0; j < N; j++)
			A[i][j] =2*(rand() % 2)-1;
	for (i = 0; i < N; i++) // only 6 out of N coeff of x are nonzero
		x[i] = 0;

	/* for i = 1:iteration
		g = A'*(r);
		w = (g'*g)/(g'*(A'*A)*g);
			x = xnew + w*g;
	[v, s] = sort(abs(x), 'descend');
	T = s(1:k);
	xnew(:) = 0;
	xnew(T) = x(T);
	r = y - A*xnew;
	end */

	x[1] = 5;
	x[6] = (float)(23);
	x[2] = (float)(12);





	//printf("x:\n");
	//print((float *)x, 1, N);

	matMul_N_in_Middle((float *)A, (float *)x, y, M, N, 1);

	//sorted_order((float *)x, (*)ind, 5, 5);
	//print_int((*)ind, 5, 1);

	//IHT((float *)A, (float *)y, (float *)x_hat, M, N, 5, number_of_iterations, 0.001);



	float xnew[N], xnew_Intermediate[N], r[M], r_Intermediate[M];
	for (j = 0; j < N; j++)//xnew = zeros(N, 1);
		xnew[j] = 0;
	float At[N][M];
	float A_in_xnew[N];
	float xnew_in_A[N];
	transpose((float *)A, (float *)At, M, N);

	matMul((float *)A, (float *)xnew, (float *)r_Intermediate, M, N, 1);//r = y - A*xnew;
	matSUB_M((float *)y, (float *)r_Intermediate, (float *)r, M, 1);//r = y - A*xnew;

	for (i = 0; i < N; i++)
	{
		xnew_in_A[i] = 0;
		A_in_xnew[i] = 0;
	}

	int ii,iii;
	printf("Begining of Algorithms:");
//for (iii=0;iii<10;iii++)
	for (ii = 0; ii < number_of_iteration; ii++)
	{
		matMul_M_in_Middle((float *)At, (float *)r, (float *)xnew_Intermediate, N, M, 1); //A*g

		/*	g = A'*(r);
			w = (g'*g)/(g'*(A'*A)*g);
				x = xnew + w*g;
		*/
		matMul_N_in_Middle((float *)xnew_Intermediate, (float *)At, xnew_in_A, 1, N, M);   //g'*At
		matMul_N_in_Middle((float *)A, (float *)xnew_Intermediate, (float *)A_in_xnew,M, N, 1); //A*g
		float step =0;
		float Num = 0; float Denum = 0;
		for (i = 0; i < N; i++)
		{
			Num += (xnew_Intermediate[i] * xnew_Intermediate[i]);
			Denum += (A_in_xnew[i] * xnew_in_A[i]);
		}
		step = Num / Denum;
		matMul((float *)xnew_Intermediate, &step, (float *)xnew_Intermediate, N, 1, 1);
		matADD_N((float *)xnew, (float *)xnew_Intermediate, (float *)x, N, 1);//x = xnew + w*A'*(r);

		int ind[N];
		sorted_order((float *)x, (int *)ind, N, sparsity);//[v, s] = sort(abs(x), 'descend');

		for (j = 0; j < N; j++)//xnew = zeros(N, 1);
		{
			xnew[j] = 0;
		}
		for (j = 0; j < sparsity; j++)//xnew(T) = x(T);
		{
			xnew[ind[j]] = x[ind[j]];
		}
		matMul_N_in_Middle((float *)A, (float *)xnew, (float *)r_Intermediate, M, N, 1);//r = y - A*xnew;
		matSUB_M((float *)y, (float *)r_Intermediate, (float *)r, M, 1);//r = y - A*xnew;


	}
	printf("End of Algorithms:");
	SNR(x,xnew,N);
	return 0;
}

int swap(int *xp, int *yp)
{
	int temp = *xp;
	*xp = *yp;
	*yp = temp;
	return(0);
}

int bubbleSort(int *arr, int n)
{
	int i, j;
	for (i = 0; i < n - 1; i++)

		// Last i elements are already in place
		for (j = 0; j < n - i - 1; j++)
			if (arr[j] > arr[j + 1])
				swap(&arr[j], &arr[j + 1]);
	return(0);
}


int matMul(float *mat1, float *mat2, float *result, int m, int n, int q) {
	//Multiplication  for Mat1_(m*n) * Mat2_(n*k)
	// M,N and K are dimentions of matrices
	// Output is going to sit in the matrix that result represents here
	float sum;	int i,j,k;
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



int matSUB_M(float *mat1, float *mat2, float *result, int m, int n) {
	//Calculates mat1_m*n - mat2_m*n and reutrns the resultant matrix to result
	//Matrices should be the same size
	int i,j;
	for (i = 0; i < M; i++) {
		for (j = 0; j < n; j++)
			*(result + i*n + j) = *(mat1 + i*n + j) - *(mat2 + i*n + j);
	}
	return(0);
}

int matADD_N(float *mat1, float *mat2, float *result, int m, int n) {
	//Calculates mat1_m*n + mat2_m*n and reutrns the resultant matrix to result
	//Matrices should be the same size
	int i,j;
	for (i = 0; i < N; i++) {
		for (j = 0; j < n; j++)
			*(result + i*n + j) = *(mat1 + i*n + j) + *(mat2 + i*n + j);
	}
	return(0);
}

int transpose(float *mat1, float *mat2, int m, int n)
{
	//transposes mat1_m*n which is input to mat2_n*m which is out
	// M,N  dimentions of matrices
	int i,j;
	for (i = 0; i<m; i++)
		for (j = 0; j < n; j++)
		{
			*(mat2 + j*m + i) = *(mat1 + i*n + j);
		}
	return(0);
}
/*int print(float *mat, int m, int n) {
	//  Prints mat with dimentions of m and n
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			//_TCE_STREAM_OUT((char) *(mat+i*n+j));
			printf("%f  ", *(mat + i*n + j));
		}
		printf("\n");
	}
	return(1);
}*/
int sorted_order(float *arr, int *sorted_ind, int n, int k)
// arr is input array, n is length of input array, k is number of returned indeces,sorted_int is array of sorted indices
{
	int  i, j;

	for (i = 0; i<n; i++)
	{
		*(sorted_ind + i) = i;
	}
	for (i = 0; i<k + 1; i++) //Only k largests are required
	{
		for (j = i + 1; j<n; j++)
		{
			if (fabs(*(arr + *(sorted_ind + i))) < fabs(*(arr + *(sorted_ind + j))))
			{
				swap(&sorted_ind[i], &sorted_ind[j]);
			}
		}
	}

	return (0);
}




int matMul_N(float *mat1, float *mat2, float *result, int m, int n, int q) {
	float sum1, sum2, sum3, sum4;	int i,j,k;
	for (i = 0; i < m; i++) {
		for (j = 0; j < q; j++) {
			sum1 = 0; sum2 = 0; sum3 = 0; sum4 = 0;
			float* index_t1 = mat1 + i*n;
			float* index_t2 = mat2 + j;
			for (k = 0; k < N; k += 4) {
				sum1 = sum1 + *(index_t1 + k) * *(index_t2 + k*q);	//sum = sum + *(mat1 + i*n + k) * *(mat2 + k*q + j);
				sum2 = sum2 + *(index_t1 + k + 1) * *(index_t2 + (k + 1)*q);
				sum3 = sum3 + *(index_t1 + k + 2) * *(index_t2 + (k + 2)*q);
				sum4 = sum4 + *(index_t1 + k + 3) * *(index_t2 + (k + 3)*q);
			}
			*(result + i*q + j) = sum1 + sum2 + sum3 + sum4;
		}
	}
	return(0);
}

int matMul_N_in_Middle(float *mat1, float *mat2, float *result, int m, int n, int q) {
	float sum1, sum2, sum3, sum4;	int i,j,k;
	for (i = 0; i < m; i++) {
		for (j = 0; j < q; j++) {
			sum1 = 0; sum2 = 0; sum3 = 0; sum4 = 0;
			float* index_t1 = mat1 + i*n;
			float* index_t2 = mat2 + j;
			for (k = 0; k < N; k += 4) {
				sum1 = sum1 + *(index_t1 + k) * *(index_t2 + k*q);	//sum = sum + *(mat1 + i*n + k) * *(mat2 + k*q + j);
				sum2 = sum2 + *(index_t1 + k + 1) * *(index_t2 + (k + 1)*q);
				sum3 = sum3 + *(index_t1 + k + 2) * *(index_t2 + (k + 2)*q);
				sum4 = sum4 + *(index_t1 + k + 3) * *(index_t2 + (k + 3)*q);
			}
			*(result + i*q + j) = sum1 + sum2 + sum3 + sum4;
		}
	}
	return(0);
}

int matMul_M_in_Middle(float *mat1, float *mat2, float *result, int m, int n, int q) {
	float sum1, sum2, sum3, sum4;int i,j,k;
	for (i = 0; i < m; i++) {
		for (j = 0; j < q; j++) {
			sum1 = 0; sum2 = 0; sum3 = 0; sum4 = 0;
			float* index_t1 = mat1 + i*n;
			float* index_t2 = mat2 + j;
			for (k = 0; k < M; k += 4) {
				sum1 = sum1 + *(index_t1 + k) * *(index_t2 + k*q);	//sum = sum + *(mat1 + i*n + k) * *(mat2 + k*q + j);
				sum2 = sum2 + *(index_t1 + k + 1) * *(index_t2 + (k + 1)*q);
				sum3 = sum3 + *(index_t1 + k + 2) * *(index_t2 + (k + 2)*q);
				sum4 = sum4 + *(index_t1 + k + 3) * *(index_t2 + (k + 3)*q);
			}
			*(result + i*q + j) = sum1 + sum2 + sum3 + sum4;
		}
	}
	return(0);
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
	int l;
	int i;
	int temp;
	float mse = MSE(a, b, Length);
	float signal_power = 0;
	for (i = 0; i<Length; i++)
		signal_power += (a[i] * a[i]);

	printf("mse : %d \n",(int) mse);
	printf("signal_power : %d \n",(int) signal_power);
	float snr = 10 * log10(signal_power / mse);
	printf("SNR : %d \n", (int)snr);
	return snr;
}
