//QR_Decompostion

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