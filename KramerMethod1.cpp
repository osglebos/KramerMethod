#include <stdio.h>
#include <windows.h>
#include <math.h>

#define ARR_SIZE 12

void GetMatr(double **mas, double **p, int a, int b, int m) {

	for (int k = 0; k < m - 1; k++)
		for (int d = 0; d < m - 1; d++)
			p[k][d] = 0;
	int di = 0, dj = 0;
	for (int i = 0; i < m - 1; i++)
	{
		if (i == a)di = 1;
		for (int j = 0; j < m - 1; j++)
		{
			if (j == b)dj = 1;
			p[i][j] = mas[i + di][j + dj];
		}
	}
}

double Determinant(double **mas, int  m)
{
	if (m < 1) printf("ERROR");
	if (m == 1) return mas[0][0];
	if (m == 2)	return mas[0][0] * mas[1][1] - mas[1][0] * mas[0][1];
	if (m > 2) {
		int i, j = 0, k = 1;
		double d = 0;

		double **p = new double*[m - 1];
		for (int i = 0; i < 4; i++) p[i] = new double[m - 1];
		for (i = 0; i < m; i++) {
			GetMatr(mas, p, i, 0, m);
			d = d + k * mas[i][0] * Determinant(p, m - 1);
			k = -k;
		}
		return(d);
	}
	return -1;
}

double CountRoot(double **matrix_x, double *matrix_y, int column)
{
	for (int i = 0; i < 4; i++)
		matrix_x[column][i] = matrix_y[i];

	return Determinant(matrix_x, 4);
}

void CopyMatrix(double **matrix, double **matrix_copy)
{
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			matrix[i][j] = matrix_copy[i][j];
}

int main(void)
{
	double Coords[2][ARR_SIZE];
	double **Powers_X = new double*[4];
	double **Powers_X_copy = new double*[4];
	double *Powers_Y = new double[4];

	for (int i = 0; i < 4; i++)
	{
		Powers_X[i] = new double[4];
		Powers_X_copy[i] = new double[4];
		Powers_Y[i] = 0;
		for (int j = 0; j < 4; j++) Powers_X[i][j] = Powers_X_copy[i][j] = 0;
	}

	Coords[0][0] = -3.10;
	Coords[0][1] = -2.30;
	Coords[0][2] = -1.50;
	Coords[0][3] = -0.70;
	Coords[0][4] = 0.10;
	Coords[0][5] = 0.90;
	Coords[0][6] = 1.70;
	Coords[0][7] = 2.50;
	Coords[0][8] = 3.30;
	Coords[0][9] = 4.10;
	Coords[0][10] = 4.90;
	Coords[0][11] = 5.70;

	Coords[1][0] = 6.07;
	Coords[1][1] = 5.88;
	Coords[1][2] = 5.61;
	Coords[1][3] = 5.34;
	Coords[1][4] = 5.17;
	Coords[1][5] = 4.89;
	Coords[1][6] = 4.65;
	Coords[1][7] = 4.41;
	Coords[1][8] = 4.12;
	Coords[1][9] = 2.91;
	Coords[1][10] = 3.65;
	Coords[1][11] = 3.47;
	
	for (int i = 0; i < ARR_SIZE; i++)
		printf("x: %5.2lf | y: %5.2lf\n", Coords[0][i], Coords[1][i]);
	printf("\n");

	for (int i = 0; i < ARR_SIZE; i++)
	{
		Powers_X[0][0] += 1;
		Powers_X[0][1] = Powers_X[1][0] += Coords[0][i];
		Powers_X[0][2] = Powers_X[1][1] = Powers_X[2][0] += Coords[0][i] * Coords[0][i];
		Powers_X[0][3] = Powers_X[1][2] = Powers_X[2][1] = Powers_X[3][0] += Coords[0][i] * Coords[0][i] * Coords[0][i];
		Powers_X[1][3] = Powers_X[2][2] = Powers_X[3][1] += Coords[0][i] * Coords[0][i] * Coords[0][i] * Coords[0][i];
		Powers_X[2][3] = Powers_X[3][2] += Coords[0][i] * Coords[0][i] * Coords[0][i] * Coords[0][i] * Coords[0][i];
		Powers_X[3][3] += Coords[0][i] * Coords[0][i] * Coords[0][i] * Coords[0][i] * Coords[0][i] * Coords[0][i];

		Powers_Y[0] += Coords[1][i];
		Powers_Y[1] += Coords[0][i] * Coords[1][i];
		Powers_Y[2] += Coords[0][i] * Coords[0][i] * Coords[1][i];
		Powers_Y[3] += Coords[0][i] * Coords[0][i] * Coords[0][i] * Coords[1][i];
	}

	CopyMatrix(Powers_X_copy, Powers_X);

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++)
			printf("%11.5lf ", Powers_X[j][i]);
		printf("\n");
	}
	printf("\n");


	for (int i = 0; i < 4; i++)
		printf("Powers_Y[%i]: %10.6lf\n", i, Powers_Y[i]);
	printf("\n");

	double det = Determinant(Powers_X, 4);

	if (det == 0)
	{
		printf("Determinant is zero, so there is no roots for this system.\n\n");
		system("pause");
		return 0;
	}

	double d = CountRoot(Powers_X, Powers_Y, 0) / det;

	CopyMatrix(Powers_X, Powers_X_copy);
	double c = CountRoot(Powers_X, Powers_Y, 1) / det;

	CopyMatrix(Powers_X, Powers_X_copy);
	double b = CountRoot(Powers_X, Powers_Y, 2) / det;

	CopyMatrix(Powers_X, Powers_X_copy);
	double a = CountRoot(Powers_X, Powers_Y, 3) / det;

	printf("A = %5.5lf\n", a);
	printf("B = %5.5lf\n", b);
	printf("C = %5.5lf\n", c);
	printf("D = %5.5lf\n\n", d);

	double summ = 0;
	for (int i = 0; i < ARR_SIZE; i++)
	{
		double temp = a * Coords[0][i] * Coords[0][i] * Coords[0][i] + b * Coords[0][i] * Coords[0][i] + c * Coords[0][i] + d;
		printf("x: %5.2lf | y: %5.2lf | P(x): %5.2lf | y-P(x): %5.2lf\n", Coords[0][i], Coords[1][i], temp, Coords[1][i] - temp);
		summ += (Coords[1][i] - temp)*(Coords[1][i] - temp);
	}
	printf("\nDelta: %5.4lf%%\n\n", sqrt(summ / (ARR_SIZE + 1)) * 100);

	system("pause");
	return 0;
}

