/* Orthogonal Matching Persuits Implementation By Mehdi Safarpour,
The QR decompoision and Jacobi based solvered have been used in this version
This implementation intended to be used in embedded systems,
The current version is not suitible yet and can be improved to reducing redundancies in calculations
2017
In this version, the functions that require definition of arrays inside are inlinded in the main
*/


#include "stdio.h"
#include <stdlib.h>
#include "math.h"

#define M 64 // Number of measurements
#define N 256//Number of samples
#define number_of_finded_Elements 15
#define S 6 //Sparsity degree
#define number_of_iterations 10


int number_of_finded = 0;

int matMul(float *mat1, float *mat2, float *result, int m, int n, int q);
int matMul_M(float *mat1, float *mat2, float *result, int m, int n, int q);
int matMul_N(float *mat1, float *mat2, float *result, int m, int n, int q);
int transpose(float *mat1, float *mat2, int m, int n);
float innerMulColumn_N(float *mat1, float *vector, int numberof_of_rows_Mat1, int numberof_of_columns_Mat1, int C);
int matSUB(float *mat1, float *mat2, float *result, int m, int n);
int matADD(float *mat1, float *mat2, float *result, int m, int n);
int QR(float *A, float *Q, float *R, int n);
int swap(int *xp, int *yp);
int bubbleSort(int *arr, int n);
int print(float *mat, int m, int n);
float sum_array(float *arr, int m);
float innerMatColumnMAT(float *mat, int n, int C1, int C2);
int Union(int *vec, int newval);
int backsubstitotion(float *R, float *y_Qt, float *x_hat, int n);
float norm_Col(float *A, int m, int n, int C);
int max_index(float *vector, int size);
int find_max_index(float *A, float *r, int m, int n, int *finded_index, int *remaining_index, int number_of_remaining, int number_of_finded);
float SNR(float *a, float *b, int Length);
float MSE(float *a, float *b, int Length);


static struct timeval tm1;
static inline void start()
{
    gettimeofday(&tm1, NULL);
}

static inline void stop()
{
    struct timeval tm2;
    gettimeofday(&tm2, NULL);

    unsigned long long t = 1000 * (tm2.tv_sec - tm1.tv_sec) + (tm2.tv_usec - tm1.tv_usec) / 1000;
    printf("%llu ms\n", t);
}

int main() {

	float A[M][N];
	int i,j;
	for (i = 0; i < M; i++)
		for ( j = 0; j < N; j++){
		A[i][j] =(float)(rand()%10);
		}
	float x[N];
	float  y[M];
	float x_hat[N];

/*	for (int i = 0;i < M;i++)
		for (int j = 0;j < N/4;j=j+4)
{
				int data;
			_TCE_STREAM_IN("fifo_s16_stream_in_fifo_s16_stream_in_status",	data );
			A[i][j]=(float)data;
			_TCE_STREAM_IN("fifo_s16_stream_in_fifo_s16_stream_in_status",	data );
			A[i][j+1]=(float)data;
			_TCE_STREAM_IN("fifo_s16_stream_in_fifo_s16_stream_in_status",	data );
			A[i][j+2]=(float)data;
			_TCE_STREAM_IN("fifo_s16_stream_in_fifo_s16_stream_in_status",	data );
			A[i][j+3]=(float)data;
}*/

/*
for (int i = 0; i < M; i++)  // WHEN I USE this part all programs inlcuding this one behave strangely!!!!!!!!
		for (int j = 0; j < N; j++)
			A[i][j] = (2*(rand() % 2)-1)   ;
*/

	for (i = 0;i <N;i=i+1) // only 6 out of N coeff of x are nonzero
		x[i] = i/(i+10);

	x[6] = (float)(11);
	x[2] = (float)(12);
	x[3] = (float)(32);


	//for (int i = 0;i < S;i++) // only 4 out of N coeff of x are nonzero
		//x[(int)(rand() % N - 1)] = (float)(rand() % 25);


	matMul_N(( float *)A, (float *)x, (float *)y, M, N, 1);

//print(y,M,1);
	/* BEGINNIG OF OMP ALGORITHM */

		 {

		float r[M ], norm_of_colum[N ], correlation[N];
		int finded_index[N];

		number_of_finded = 0;

		for ( i = 0;i < N;i++)
		{
			*(x_hat + i)= 0;
			finded_index[i] = 999999; //Initialize with some value less than all
			norm_of_colum[i] = norm_Col((float *)A, M, N, i); //Calculate norm of column of C
		}


		for ( i = 0;i < M;i++) { // in matlab code x meant to be y
			r[i] = *(y + i); //r = x;
		}

		printf("Algorim Begins \n");
		int ii,iii;
		start();
		for (iii=0;iii<100;iii++)
		for ( ii = 0;ii < number_of_iterations;ii++) { //main loop,  q is number of expected coefficients or iterations

			for ( i = 0;i < N;i+=4) //// IMPORTANT LOOP
			{
				correlation[i] = fabs(((float)innerMulColumn_N((float *)A, (float *)r, M, N, i)) / ((float)norm_of_colum[i]));
				correlation[i+1] = fabs(((float)innerMulColumn_N((float *)A, (float *)r, M, N, i+1)) / ((float)norm_of_colum[i+1]));
				correlation[i+2] = fabs(((float)innerMulColumn_N((float *)A, (float *)r, M, N, i+2)) / ((float)norm_of_colum[i+2]));
				correlation[i+3] = fabs(((float)innerMulColumn_N((float *)A, (float *)r, M, N, i+3)) / ((float)norm_of_colum[i+3]));
			}
			int max_value_index = max_index((float *)correlation, N);// find(f == max(f));
			Union((int *)finded_index, max_value_index);

			{
				//Solve_LS_with_finding_columns((float *)A, y, x_hat, (int *)finded_index, number_of_finded, m, n);
				/* Start of Solve Least Sq. Algorithm */

				float Acopy[M*number_of_finded_Elements], At[number_of_finded_Elements*M], A_square[number_of_finded_Elements*number_of_finded_Elements],
					Q[number_of_finded_Elements*number_of_finded_Elements], R[number_of_finded_Elements*number_of_finded_Elements], Qt[number_of_finded_Elements*number_of_finded_Elements],
					yQt[number_of_finded_Elements];
				// Defining a two D array
				/*float *Acopy = (float *)calloc(M * number_of_finded, sizeof(float));
				float *At = (float *)calloc(M * number_of_finded, sizeof(float)); //Transpose of A copy
				float *A_square = (float *)calloc(number_of_finded *number_of_finded, sizeof(float)); //Transpose of A copy
				float *Q = (float *)calloc(number_of_finded * number_of_finded, sizeof(float)); //Transpose of A copy
				float *R = (float *)calloc(number_of_finded *number_of_finded, sizeof(float)); //Transpose of A copy
				float *Qt = (float *)calloc(number_of_finded*number_of_finded, sizeof(float)); //Transpose of Q copy
				float *yQt = (float *)calloc(number_of_finded, sizeof(float)); //Transpose of A copy
				*/
				float y_p[number_of_finded_Elements]; // will be used as At * y (before Jacobi)
				float x_hat_temp[M];

				for ( i = 0;i < M;i++) // Get a copy of A with only selected index
					for ( j = 0;j < number_of_finded;j++) {
						Acopy[i*number_of_finded+ j] = A[i][finded_index[j]];
					}



				transpose((float *)Acopy, (float *)At, M, number_of_finded);

				matMul_M((float *)At, y, (float *)y_p, number_of_finded, M, 1); //At*y
				matMul_M((float *)At, (float *)Acopy, (float *)A_square, number_of_finded, M, number_of_finded);

				QR((float *)A_square, (float *)Q, (float *)R, number_of_finded);

				transpose((float *)Q, (float *)Qt, number_of_finded, number_of_finded);
				/*for (int i = 0;i<number_of_finded;i++)
					for (int j = 0; j < number_of_finded;j++)
					{
						Qt[j][i] = Q[i][j];
					}*/
				matMul((float *)Qt, (float *)y_p, (float *)yQt, number_of_finded, number_of_finded, 1);

				backsubstitotion((float *)R, (float *)yQt, (float *)x_hat_temp, number_of_finded);

				//Jaocabi(A_square, y_p, x_hat_temp, number_of_finded, 100); // Solve the least squar for and get the results in x, of course length of x_hat_temp is more than necessary but we know how many of the are valid with number_of_finded

				for ( i = 0;i < number_of_finded;i++)
					x_hat[finded_index[i]] = x_hat_temp[i]; //put the right value in the right positions



			float A_x_h[M];
			matMul_N((float *)A, x_hat, A_x_h, M, N, 1);

			matSUB(y, A_x_h, r, M, 1); //r = y - A*s_hat;

			}
		}/*END OF OMP SECTION*/
	}
	stop();
	//print((float *)A, 1, N);

	//print((float *)x_hat, 1, N);
	printf("Finished");

	SNR(x,x_hat,N);
		return (0);

}


int matMul(float *mat1, float *mat2, float *result, int m, int n, int q) {
	//Multiplication  for Mat1_(m*n) * Mat2_(n*k)
	// M,N and K are dimentions of matrices
	// Output is going to sit in the matrix that result represents here
	float sum;
	int i,j,k;
	for ( i = 0; i < m; i++) {
		for ( j = 0; j < q; j++) {
			sum = 0;
			for ( k = 0; k < n; k++) {
				sum = sum + *(mat1 + i*n + k) * *(mat2 + k*q + j);
			}
			*(result + i*q + j) = sum;
		}
	}
	return(0);
}

int matMul_M( float *mat1, float *mat2, float *result, int m, int n, int q) {
	float sum1,sum2,sum3,sum4;
	int i,j,k;
	for ( i = 0; i < m; i++) {
		for ( j = 0; j < q; j++) {
			sum1 = 0;sum2 = 0;sum3 = 0;sum4 = 0;
			float* index_t1=mat1 + i*n;
			float* index_t2=mat2 + j;
			for ( k = 0; k < M; k+=4) {
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
	int i,j,k;
	for ( i = 0; i < m; i++) {
		for ( j = 0; j < q; j++) {
			sum1 = 0;sum2 = 0;sum3 = 0;sum4 = 0;
			float* index_t1=mat1 + i*n;
			float* index_t2=mat2 + j;
			for ( k = 0; k < N; k+=4) {
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
	int i,j;
	for ( i = 0;i<m;i++)
		for ( j = 0; j < n;j++)
		{
			*(mat2 + j*m + i) = *(mat1 + i*n + j);
		}
	return(0);
}

float innerMulColumn_N(float *mat1, float *vector, int numberof_of_rows_Mat1, int numberof_of_columns_Mat1, int C) {
	//Returns inner porduct of a column vector of a and b matrix, colum number 1,2,3,...
	//d C is the column of interest
	float sum1= 0,sum2 = 0,sum3 = 0,sum4 = 0;
	float*	commont_index=mat1+C;
	int i;
	for ( i = 0;i<M;i+=4){
		sum1 = sum1 + *(commont_index + i*numberof_of_columns_Mat1) * *(vector + i);
		sum2 = sum2 + *(commont_index + (i+1)*numberof_of_columns_Mat1 ) * *(vector + i+1);
		sum3 = sum3 + *(commont_index + (i+2)*numberof_of_columns_Mat1 ) * *(vector + i+2);
		sum4 = sum4 + *(commont_index + (i+3)*numberof_of_columns_Mat1 ) * *(vector + i+3);
		}
	return(sum1+sum2+sum3+sum4); // The scalar value of inner product
}

int max_index(float *vector, int size) {
	int index_of_largest_valyue = 0;
	float max = *(vector + index_of_largest_valyue);
	int i;
	for ( i = 0; i < size; i++)
	{
		if (*(vector + i) > max)
		{
			max = *(vector + i);
			index_of_largest_valyue = i;
		}
	}
	return(index_of_largest_valyue);
}
int QR(float *A, float *Q, float *R, int n) {

	/*Q=A;

	Q(:,1) = Q(:,1);
	R(1,1) = norm(Q(:,1));
	Q(:,1) = Q(:,1)/R(1,1);

	for k = 2:n,
	R(1:k-1,k) = Q(:,1:k-1)'*Q(:,k);
	Q(:,k) = Q(:,k)- Q(:,1:k-1)*R(1:k-1,k);
	R(k,k) = norm(Q(:,k));
	Q(:,k) = Q(:,k)/R(k,k);
	end
	end
	*/
	int i,j,K;
	for ( i = 0;i < n;i++) //Q = A;
		for ( j = 0;j < n;j++) {
			*(Q + i*n + j) = *(A + i*n + j);
			*(R + i*n + j) = 1e-20;
		}
	/*	R(1:k-1, k) = Q(1..n, 1 : k-1)'*Q(1..n,k);
	Q(1..n, k) = Q(1..n, k) - Q(1..n, 1 : k-1)*R(1:k-1,k);
	R(k, k) = norm(Q(1..n, k));
	Q(1..n, k) = Q(1..n, k) / R(k, k);
	*/
	for ( K = 0;K < n;K++) {
		for ( i = 0;i < K;i++)
			*(R + i*n + K) = innerMatColumnMAT(Q, n, i, K); //	R(1:k-1,k) = Q(:,1:k-1)'*Q(:,k);

		for ( i = 0;i < n;i++) {
			float QinR = 0;
			for ( j = 0;j < K;j++)
				QinR = QinR + *(Q + i*n + j) * *(R + j*n + K);
			*(Q + i*n + K) = *(Q + i*n + K) - QinR;
		}
		*(R + K*n + K) = norm_Col(Q, n, n, K);
		for ( i = 0;i < n;i++) {

			*(Q + i*n + K) = *(Q + i*n + K) / *(R + K*n + K);
		}
	}



	return(1);
}
int matSUB(float *mat1, float *mat2, float *result, int m, int n) {
	//Calculates mat1_m*n - mat2_m*n and reutrns the resultant matrix to result
	//Matrices should be the same size
	int i,j;
	for ( i = 0;i < m;i++) {
		for ( j = 0;j < n;j++)
			*(result + i*n + j) = *(mat1 + i*n + j) - *(mat2 + i*n + j);
	}
	return(0);
}


int matADD(float *mat1, float *mat2, float *result, int m, int n) {
	//Calculates mat1_m*n + mat2_m*n and reutrns the resultant matrix to result
	//Matrices should be the same size
	int i,j;
	for ( i = 0;i < m;i++) {
		for ( j = 0;j < n;j++)
			*(result + i*n + j) = *(mat1 + i*n + j) + *(mat2 + i*n + j);
	}
	return(0);
}


int print(float *mat, int m, int n) {
	//  Prints mat with dimentions of m and n
	int i,j;
	for ( i = 0; i < m; i++) {
		for ( j = 0; j < n; j++) {
			//_TCE_STREAM_OUT((char) *(mat+i*n+j));
				printf("%d ", (int)  *(mat + i*n + j));
		}
		printf("\n");
	}
	return(1);
}
/*
int print_int(int *mat, int m, int n) {
	//  Prints mat with dimentions of m and n
	for (int i = 0;i < m;i++) {
		for (int j = 0;j < n;j++)
			std::cout << *(mat + i*n + j) << " ";
		std::cout << std::endl;
	}
	return(1);
}*/
float norm_Col(float *A, int m, int n, int C)
{
	int i;
	float sum = 0;
	float *fixed_part_of_addrss=A + C ;
	for ( i = 0;i < m;i++)
		sum = sum + *(fixed_part_of_addrss + i*n ) * *(fixed_part_of_addrss + i*n);
	sum=sqrt(sum);
	return sum;
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
int backsubstitotion(float *R, float *y_Qt, float *x_hat, int n) {
	int i,j,d;
	for ( d = 0;d < n;d++)
		*(x_hat + d) = 0;
	for ( i = n - 1;i >= 0;i--)
	{
		*(x_hat + i) = *(y_Qt + i) / *(R + i*n + i);
		for ( j = 0;j < i;j++)
			*(y_Qt + j) = *(y_Qt + j) - *(R + j*n + i) * *(x_hat + i);
	}
	return(1);
}
float innerMatColumnMAT(float *mat, int n, int C1, int C2)
{// n is size, C1 and C2 are columns that are correlated
	float sum = 0;
	int i;
	for ( i = 0;i < n;i++)
		sum += *(mat + i*n + C1) * *(mat + i*n + C2);
	return(sum);
}
int Union(int *vec, int newval) {
	// Using this function values of index are sorted and there is no need to check for repitition
	int i;
	for ( i = 0;i<number_of_finded;i++)
		if (*(vec + i) == newval) // If there is repeatition, ignore
			return(0);

	//if (*(vec + number_of_finded - 1) != newval)
	{
		*(vec + number_of_finded) = newval;  //If it is bigger than all put it in the righmost place
		number_of_finded++;
		bubbleSort(vec, number_of_finded);
	}
	return(0);
}

float sum_array(float *arr, int m) {
	float sum = 0;
	int i;
	for ( i = 0;i < m;i++)
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
	int l;
	int i;
	int temp;
	float mse = MSE(a, b, Length);
	float signal_power = 0;
	for ( i = 0; i<Length; i++)
		signal_power += a[i] * a[i];

	printf("mse : %d \n",(int) mse);
	printf("signal_power : %d \n",(int) signal_power);
	float snr = 10 * log10(signal_power / mse);
	printf("SNR : %d \n",(int) snr);
	return snr;
}



























