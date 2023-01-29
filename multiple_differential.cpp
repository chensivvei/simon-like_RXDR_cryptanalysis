#include <stdio.h>
#include <stdint.h>
#include <vector>
#include <iostream>
#include <math.h>
#include <string.h>
#include <time.h>
using namespace std;
//#define SIMECK
#define SIMON
#define PRECISION 12// Window size (w)


#if defined(SIMON) == defined(SIMECK)
#error Please only exactly one of SIMON and SIMECK
#endif

#define ROT(x,n) (((uint16_t)(x)<<((n)&15)) | ((uint16_t)(x)>>((16-(n))&15)))


uint8_t hamming_weight(uint16_t a)
{
	a = (a & 0x5555) + ((a >> 1) & 0x5555);
	a = (a & 0x3333) + ((a >> 2) & 0x3333);
	a = (a & 0x0f0f) + ((a >> 4) & 0x0f0f);
	return (a & 0xff) + ((a >> 8) & 0xff);
}

void compute_space(uint8_t *weight, vector<uint16_t> *space)
{
#ifdef SIMECK
	int a = 5, b = 0, c = 1; // Simeck
#endif
#ifdef SIMON
	int a = 8, b = 1, c = 2;; // Simon
#endif
	uint16_t alpha, beta, gamma, varibits, doublebits;
//#pragma omp parallel for schedule(dynamic)
	for (alpha = 0; alpha < (1 << PRECISION)-1; alpha++)
	{
		varibits = ROT(alpha, a) | ROT(alpha, b);
		doublebits = ROT(alpha, b) & (~ROT(alpha, a)) & ROT(alpha, 2 * a - b);
		weight[alpha] = hamming_weight(varibits ^ doublebits);
		for (beta = 0; beta < (1 << PRECISION); beta++)
		{
			gamma = beta ^ ROT(alpha, c);
			if (((gamma & (~varibits)) == 0) && (((gamma ^ ROT(gamma, a - b)) & doublebits) == 0))
				space[alpha].push_back(beta);
				
		}
	}
	if (PRECISION == 16) {
		weight[alpha] = 15;
		for (beta = 0; beta < (1 << PRECISION); beta++)
		{
			gamma = beta ^ ROT(alpha, c);
			if ((hamming_weight(gamma)%2 == 0))
				space[alpha].push_back(beta);
		}
	}
	else {
		varibits = ROT(alpha, a) | ROT(alpha, b);
		doublebits = ROT(alpha, b) & (~ROT(alpha, a)) & ROT(alpha, 2 * a - b);
		weight[alpha] = hamming_weight(varibits ^ doublebits);
		for (beta = 0; beta < (1 << PRECISION); beta++)
		{
			gamma = beta ^ ROT(alpha, c);
			if (((gamma & (~varibits)) == 0) && (((gamma ^ ROT(gamma, a - b)) & doublebits) == 0))
				space[alpha].push_back(beta);

		}
	}
}
double compute_diff_square_sum_pro(double** X)
{
	double max = 0;
	for (int i = 0; i < 1 << PRECISION; i++)
	{
		for (int j = 0; j < 1 << PRECISION; j++)
		{
			if (X[i][j] > 0)
			{
				max += pow(X[i][j], 2);
			}
		}
	}
	return max;
}

double compute_max_diff_pro(double** X, uint16_t out[2])
{

	double max = 0;
	for (int i = 0; i < 1 << PRECISION; i++)
	{
		for (int j = 0; j < 1 << PRECISION; j++)
		{
			if (max < X[i][j])
			{
				max = X[i][j];
				out[0] = i;
				out[1] = j;
			}
		}
	}
	return max;
}

double **one_round_diff_transition_matrix(double **X, uint8_t* weight, vector<uint16_t>* space)
{			
	double** Y = new double* [1 << PRECISION];
	for (int i = 0; i < (1 << PRECISION); i++)
	{
		Y[i] = new double[1 << PRECISION];
		memset(Y[i], 0, sizeof(double) * (1 << PRECISION));
	}
	for (int i = 0; i < (1 << PRECISION); i++)
	{
		for (int j = 0; j < (1 << PRECISION); j++)
		{
			if (X[i][j]>0)
			{
				int nb = space[i].size();
				if (nb>0)
				{
					for (int k = 0; k < nb; k++)
					{
						Y[j ^ space[i][k]][i] += pow(2, -weight[i]) * X[i][j];
					}
				}
			}
		}
	}
	for (int i = 0; i < (1 << PRECISION); i++)
	{
		memcpy(X[i], Y[i], sizeof(double) * (1 << PRECISION));
		delete[] Y[i];
	}
	delete[] Y;
	return X;
}

double** prop_square(double** X)
{
	for (int i = 0; i < 1 << PRECISION; i++)
	{
		for (int j = 0; j < 1 << PRECISION; j++)
		{
			if(X[i][j] > 0)
				X[i][j] *= X[i][j];
		}
	}
	return X;
}

void R_round_diff_square_sum(uint16_t inL, uint16_t inR, uint8_t* weight, vector<uint16_t>* space, int R)
{
	double** X = new double* [1 << PRECISION];
	for (int i = 0; i < (1 << PRECISION); i++)
	{
		X[i] = new double[1 << PRECISION];
		memset(X[i], 0, sizeof(double) * (1 << PRECISION));
	}
	X[inL][inR] = 1;
	double max_pro;
	uint16_t out[2];
	for (int r = 1; r <= R; r++)
	{

		X = one_round_diff_transition_matrix(X, weight, space);
		max_pro = compute_max_diff_pro(X, out);
		double pro_sum = compute_diff_square_sum_pro(X);
		printf("%d-round, (0x%x, 0x%x)->(0x%x, 0x%x), max_pro = 2^(%f), sum_sqaure_pro = 2^(%f)\n", r, inL, inR, out[0], out[1], log2(max_pro), log2(pro_sum));
	}
	
	for (int i = 0; i < (1 << PRECISION); i++)
	{
		delete[] X[i];
	}
	delete[] X;
}

int main()
{
	printf("w = %2d\n", PRECISION);
	clock_t s, e;
	s = clock();
	uint8_t* weight = new uint8_t[1 << PRECISION];
	vector<uint16_t>* space = new vector<uint16_t>[1 << PRECISION];
	compute_space(weight, space);
	e = clock();
	printf("Computing Space Time: %.2f s\n", (double)(e - s) / CLOCKS_PER_SEC);
	uint16_t inL = 0x8;
	uint16_t inR = 0x2022;

	R_round_diff_square_sum(inL, inR, weight, space, 4);
	delete[] weight;
	delete[] space;
	e = clock();
	printf("Running Time: %.2f s", (double)(e - s) / CLOCKS_PER_SEC);
	return 0;
}