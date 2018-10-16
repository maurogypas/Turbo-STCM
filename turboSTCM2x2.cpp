#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define termi 2     // After encoding, 2 termination symbols are sent from each constituent encoder
#define N 2050      // interleaver length (in symbols) + termi

double variance, **alpha, **beta, **gamma, **H11, **H12, **H21, **H22,
       nSymbol0[8], nSymbol1[8], nSymbol2[8], nSymbol3[8], bSymbol0[8], bSymbol1[8], bSymbol2[8], bSymbol3[8];
char nextState0[4], nextState1[4], nextState2[4], nextState3[4], preState0[4], preState1[4], preState2[4], preState3[4];

/* ***************************************************************************** */
/* Copyright:      Francois Panneton and Pierre L'Ecuyer, University of Montreal */
/*                 Makoto Matsumoto, Hiroshima University                        */
/* Notice:         This code can be used freely for personal, academic,          */
/*                 or non-commercial purposes. For commercial purposes,          */
/*                 please contact P. L'Ecuyer at: lecuyer@iro.UMontreal.ca       */
/* ***************************************************************************** */

#define W 32
#define R 32
#define M1 3
#define M2 24
#define M3 10

#define MAT0POS(t,v) (v^(v>>t))
#define MAT0NEG(t,v) (v^(v<<(-(t))))
#define Identity(v) (v)

#define V0            STATE[state_i                   ]
#define VM1           STATE[(state_i+M1) & 0x0000001fU]
#define VM2           STATE[(state_i+M2) & 0x0000001fU]
#define VM3           STATE[(state_i+M3) & 0x0000001fU]
#define VRm1          STATE[(state_i+31) & 0x0000001fU]
#define newV0         STATE[(state_i+31) & 0x0000001fU]
#define newV1         STATE[state_i                   ]

#define FACT 2.32830643653869628906e-10  // 2^-32

static unsigned state_i = 0;
static unsigned STATE[R];
static unsigned z0, z1, z2;

void InitWELLRNG1024a(unsigned *init)
{
	int j;
	state_i = 0;
	for (j = 0; j<R; j++)
		STATE[j] = init[j];
}

int WELLRNG1024a(void)
{
	// the period of the random number generator is 2^1024-1
	z0 = VRm1;
	z1 = Identity(V0) ^ MAT0POS(8, VM1);
	z2 = MAT0NEG(-19, VM2) ^ MAT0NEG(-14, VM3);
	newV1 = z1^z2;
	newV0 = MAT0NEG(-11, z0) ^ MAT0NEG(-7, z1) ^ MAT0NEG(-13, z2);
	state_i = (state_i + 31) & 0x0000001fU;
	return STATE[state_i] - 2147483648;   // returns a uniformly distributed integer in [-2^31,2^31-1]
}

void zigSet(int *k, double *w, double *f)
{
	double m = 2147483648.0, d = 3.6541528853610088, t = d, v = .00492867323399, q;
	int i;

	q = v / exp(-.5*d*d);
	k[0] = (int)(d / q*m);     k[1] = 0;
	w[0] = q / m;              w[255] = d / m;
	f[0] = 1.;               f[255] = exp(-.5*d*d);

	for (i = 254; i >= 1; i--)
	{
		d = sqrt(-2.*log(v / d + exp(-.5*d*d)));
		k[i + 1] = (int)(d / t*m);
		t = d;
		f[i] = exp(-.5*d*d);
		w[i] = d / m;
	}
}

void mapper(char *help, double *symbol)
{
	for (int i = 0; i < 4; i++)
		switch (help[i])
		{    // natural mapping
			case 0:
				symbol[i << 1] = 1.; symbol[i << 1 | 1] = 0.;
				break;
			case 1:
				symbol[i << 1] = 0.; symbol[i << 1 | 1] = 1.;
				break;
			case 2:
				symbol[i << 1] = -1.; symbol[i << 1 | 1] = 0.;
				break;
			default:
				symbol[i << 1] = 0.; symbol[i << 1 | 1] = -1.;
		}
}

void trellis()
{
	char c, i, help[4];

	FILE *fp;
	fp = fopen("trellis.txt", "r");
	if (fp == NULL)
	{
		printf("Can't open trellis.txt.\n");
		system("pause");
		exit(0);
	}

	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 4; i++)
		fscanf(fp, "%hhd", nextState0 + i);

	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 4; i++)
		fscanf(fp, "%hhd", nextState1 + i);

	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 4; i++)
		fscanf(fp, "%hhd", nextState2 + i);

	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 4; i++)
		fscanf(fp, "%hhd", nextState3 + i);


	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 4; i++)
		fscanf(fp, "%hhd", help + i);
	mapper(help, nSymbol0);

	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 4; i++)
		fscanf(fp, "%hhd", help + i);
	mapper(help, nSymbol1);

	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 4; i++)
		fscanf(fp, "%hhd", help + i);
	mapper(help, nSymbol2);

	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 4; i++)
		fscanf(fp, "%hhd", help + i);
	mapper(help, nSymbol3);


	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 4; i++)
		fscanf(fp, "%hhd", preState0 + i);

	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 4; i++)
		fscanf(fp, "%hhd", preState1 + i);

	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 4; i++)
		fscanf(fp, "%hhd", preState2 + i);

	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 4; i++)
		fscanf(fp, "%hhd", preState3 + i);


	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 4; i++)
		fscanf(fp, "%hhd", help + i);
	mapper(help, bSymbol0);

	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 4; i++)
		fscanf(fp, "%hhd", help + i);
	mapper(help, bSymbol1);

	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 4; i++)
		fscanf(fp, "%hhd", help + i);
	mapper(help, bSymbol2);

	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 4; i++)
		fscanf(fp, "%hhd", help + i);
	mapper(help, bSymbol3);

	fclose(fp);
}

void createInt(int *Int, int *deInt)
{
	// random interleaver
	int i, index, *A;

	A = new int[N - termi];
	for (i = 0; i < N - termi; i++)
		A[i] = i;

	for (i = 0; i < N - termi; i++)
	{
	here:
		index = abs(WELLRNG1024a()) & (N - termi - 1);
		if (A[index] == -1)
			goto here;
		//if (i & 1 && !(index & 1) || !(i & 1) && index & 1) // odd-even constraint
		//	goto here;

		Int[i] = index;
		deInt[index] = i;
		A[index] = -1;
	}
	delete[] A;
}

void source(char *bitStream)
{
	for (int i = 0; i<N - termi; i++)
	{
		bitStream[i << 1] = WELLRNG1024a() > 0 ? 1 : 0;
		bitStream[i << 1 | 1] = WELLRNG1024a() > 0 ? 1 : 0;
	}
}

void encode(char *bitStream, char *systematicz, char *paritiez, int *Int, int *deInt)
{
	char Upbit, Downbit, input0, input1, D1[2] = { 0 }, D2[2] = { 0 };
	static char codeword2[N - termi];
	int i;

	for (i = 0; i < N - termi; i++)
	{
		Upbit = bitStream[i << 1];
		Downbit = bitStream[i << 1 | 1];

		systematicz[i] = Upbit << 1 | Downbit;

		if ((i & 1) == 0)
			paritiez[i] = D1[0] << 1 | D1[1]; // upper encoder

		input0 = Upbit ^ D1[1];
		input1 = Downbit ^ D1[0];
		D1[0] = input0;
		D1[1] = input1;


		Upbit = bitStream[deInt[i] << 1];
		Downbit = bitStream[deInt[i] << 1 | 1];

		codeword2[i] = D2[0] << 1 | D2[1];

		input0 = Upbit ^ D2[1];
		input1 = Downbit ^ D2[0];
		D2[0] = input0;
		D2[1] = input1;
	}
	for (i = 1; i < N - termi; i += 2)
		paritiez[i] = codeword2[Int[i]];

	// 1st encoder termination
	Upbit = D1[1];
	Downbit = D1[0];

	systematicz[N - 2] = Upbit << 1 | Downbit;
	paritiez[N - 2] = D1[0] << 1 | D1[1];

	// 2nd encoder termination
	Upbit = D2[1];
	Downbit = D2[0];

	systematicz[N - 1] = Upbit << 1 | Downbit;
	paritiez[N - 1] = D2[0] << 1 | D2[1];
}

void modulate(char *codeword, double *Ic, double *Q)
{
	for (int i = 0; i < N; i++)
		switch (codeword[i])
		{    // natural mapping
			case 0:
				Ic[i] = 1.; Q[i] = 0.;
				break;
			case 1:
				Ic[i] = 0.; Q[i] = 1.;
				break;
			case 2:
				Ic[i] = -1.; Q[i] = 0.;
				break;
			default:
				Ic[i] = 0.; Q[i] = -1.;
		}
}

void wgn(double *GN, int length, int *k, double *w, double *f)
{
	// creates WGN samples by applying the Ziggurat algorithm of Marsaglia & Tsang (2000)

	int randomInt;
	double uni, uni2, x, y, r = 3.6541528853610088;
	short int i;
	int j = 0;

	while (j<length)
	{
		randomInt = WELLRNG1024a();    // gennaei enan tyxaio int, apo ton opoio 8a prokypsei to deigma 8oryvou
		i = WELLRNG1024a() & 255;      // epilegei tyxaia ena apo ta 256 strwmata tou ziggurat, gennwntas* enan kainourio int,
		                               // *symfwna me to "An Improved Ziggurat Method to Generate Normal Random Samples" tou Doornik (2005)

		if (abs(randomInt)<k[i])    // to 99.33% twn deigmatwn proerxontai apo auto to if
		{
			GN[j++] = randomInt*w[i];
			continue;
		}

		if (i == 0)
		{
			do
			{
				uni = .5 + WELLRNG1024a()*FACT;
				if (uni == 0) uni = .5 + WELLRNG1024a()*FACT;
				x = -log(uni)*0.27366123732975827;   // o antistrofos tou r
				uni = .5 + WELLRNG1024a()*FACT;
				if (uni == 0) uni = .5 + WELLRNG1024a()*FACT;
				y = -log(uni);
			} while (y + y<x*x);

			GN[j++] = randomInt > 0 ? r + x : -r - x;
			continue;
		}

		uni = randomInt*w[i];
		uni2 = .5 + WELLRNG1024a()*FACT;
		if (f[i] + uni2*(f[i - 1] - f[i]) < exp(-.5*uni*uni))
			GN[j++] = uni;

	}
}

void sisoDec1(double *noisy1I, double *susie1Q, double *noisy2I, double *susie2Q, double *LLR1, double *LLR2, double *LLR3, double *La1, double *La2, double *La3)
{
	// max-log-APP algorithm
	char numOfStates = 1, s;
	int i;
	double x, y, LL0;

	//-------------------------- gamma branch metrics calculation ----------------------------//

	for (i = 0; i < N - termi; i += 2)
	{
		for (s = 0; s < numOfStates; s++)
		{
			x = H11[0][i] + H21[0][i] * nSymbol0[s << 1] - H21[1][i] * nSymbol0[s << 1 | 1];
			y = x * (2 * noisy1I[i] - x);
			x = H11[1][i] + H21[0][i] * nSymbol0[s << 1 | 1] + H21[1][i] * nSymbol0[s << 1];
			y += x * (2 * susie1Q[i] - x);
			gamma[i][s] = y / variance;

			x = H12[0][i] + H22[0][i] * nSymbol0[s << 1] - H22[1][i] * nSymbol0[s << 1 | 1];
			y = x * (2 * noisy2I[i] - x);
			x = H12[1][i] + H22[0][i] * nSymbol0[s << 1 | 1] + H22[1][i] * nSymbol0[s << 1];
			y += x * (2 * susie2Q[i] - x);
			gamma[i][s] += y / variance;


			x = -H11[1][i] + H21[0][i] * nSymbol1[s << 1] - H21[1][i] * nSymbol1[s << 1 | 1];
			y = x * (2 * noisy1I[i] - x);
			x = H11[0][i] + H21[0][i] * nSymbol1[s << 1 | 1] + H21[1][i] * nSymbol1[s << 1];
			y += x * (2 * susie1Q[i] - x);
			gamma[i][s + 4] = La1[i] + y / variance;

			x = -H12[1][i] + H22[0][i] * nSymbol1[s << 1] - H22[1][i] * nSymbol1[s << 1 | 1];
			y = x * (2 * noisy2I[i] - x);
			x = H12[0][i] + H22[0][i] * nSymbol1[s << 1 | 1] + H22[1][i] * nSymbol1[s << 1];
			y += x * (2 * susie2Q[i] - x);
			gamma[i][s + 4] += y / variance;


			x = -H11[0][i] + H21[0][i] * nSymbol2[s << 1] - H21[1][i] * nSymbol2[s << 1 | 1];
			y = x * (2 * noisy1I[i] - x);
			x = -H11[1][i] + H21[0][i] * nSymbol2[s << 1 | 1] + H21[1][i] * nSymbol2[s << 1];
			y += x * (2 * susie1Q[i] - x);
			gamma[i][s + 8] = La2[i] + y / variance;

			x = -H12[0][i] + H22[0][i] * nSymbol2[s << 1] - H22[1][i] * nSymbol2[s << 1 | 1];
			y = x * (2 * noisy2I[i] - x);
			x = -H12[1][i] + H22[0][i] * nSymbol2[s << 1 | 1] + H22[1][i] * nSymbol2[s << 1];
			y += x * (2 * susie2Q[i] - x);
			gamma[i][s + 8] += y / variance;


			x = H11[1][i] + H21[0][i] * nSymbol3[s << 1] - H21[1][i] * nSymbol3[s << 1 | 1];
			y = x * (2 * noisy1I[i] - x);
			x = -H11[0][i] + H21[0][i] * nSymbol3[s << 1 | 1] + H21[1][i] * nSymbol3[s << 1];
			y += x * (2 * susie1Q[i] - x);
			gamma[i][s + 12] = La3[i] + y / variance;

			x = H12[1][i] + H22[0][i] * nSymbol3[s << 1] - H22[1][i] * nSymbol3[s << 1 | 1];
			y = x * (2 * noisy2I[i] - x);
			x = -H12[0][i] + H22[0][i] * nSymbol3[s << 1 | 1] + H22[1][i] * nSymbol3[s << 1];
			y += x * (2 * susie2Q[i] - x);
			gamma[i][s + 12] += y / variance;
		}
		if (i == 0) numOfStates = 4;
	}

	for (i = 1; i < N - termi; i += 2)
		for (s = 0; s < numOfStates; s++)
		{
			gamma[i][s] = 0.;
			gamma[i][s + 4] = La1[i];
			gamma[i][s + 8] = La2[i];
			gamma[i][s + 12] = La3[i];
		}

	// trellis termination
	x = H11[0][N - 2] + H21[0][N - 2];
	y = x * (2 * noisy1I[N - 2] - x);
	x = H11[1][N - 2] + H21[1][N - 2];
	y += x * (2 * susie1Q[N - 2] - x);
	gamma[N - 2][preState0[0]] = y / variance;

	x = H12[0][N - 2] + H22[0][N - 2];
	y = x * (2 * noisy2I[N - 2] - x);
	x = H12[1][N - 2] + H22[1][N - 2];
	y += x * (2 * susie2Q[N - 2] - x);
	gamma[N - 2][preState0[0]] += y / variance;


	x = -H11[1][N - 2] - H21[0][N - 2];
	y = x * (2 * noisy1I[N - 2] - x);
	x = H11[0][N - 2] - H21[1][N - 2];
	y += x * (2 * susie1Q[N - 2] - x);
	gamma[N - 2][preState1[0] + 4] = y / variance;

	x = -H12[1][N - 2] - H22[0][N - 2];
	y = x * (2 * noisy2I[N - 2] - x);
	x = H12[0][N - 2] - H22[1][N - 2];
	y += x * (2 * susie2Q[N - 2] - x);
	gamma[N - 2][preState1[0] + 4] += y / variance;


	x = -H11[0][N - 2] - H21[1][N - 2];
	y = x * (2 * noisy1I[N - 2] - x);
	x = -H11[1][N - 2] + H21[0][N - 2];
	y += x * (2 * susie1Q[N - 2] - x);
	gamma[N - 2][preState2[0] + 8] = y / variance;

	x = -H12[0][N - 2] - H22[1][N - 2];
	y = x * (2 * noisy2I[N - 2] - x);
	x = -H12[1][N - 2] + H22[0][N - 2];
	y += x * (2 * susie2Q[N - 2] - x);
	gamma[N - 2][preState2[0] + 8] += y / variance;


	x = H11[1][N - 2] + H21[1][N - 2];
	y = x * (2 * noisy1I[N - 2] - x);
	x = -H11[0][N - 2] - H21[0][N - 2];
	y += x * (2 * susie1Q[N - 2] - x);
	gamma[N - 2][preState3[0] + 12] = y / variance;

	x = H12[1][N - 2] + H22[1][N - 2];
	y = x * (2 * noisy2I[N - 2] - x);
	x = -H12[0][N - 2] - H22[0][N - 2];
	y += x * (2 * susie2Q[N - 2] - x);
	gamma[N - 2][preState3[0] + 12] += y / variance;

	//-------------------------- alpha node metrics calculation ----------------------------//

	alpha[1][0] = gamma[0][0];
	alpha[1][1] = gamma[0][4];
	alpha[1][2] = gamma[0][8];
	alpha[1][3] = gamma[0][12];

	for (i = 1; i < N - termi - 1; i++)
		for (s = 0; s < numOfStates; s++)
		{
			x = gamma[i][preState0[s]] + alpha[i][preState0[s]];
			y = gamma[i][preState1[s] + 4] + alpha[i][preState1[s]];

			x = x > y ? x : y;
			y = gamma[i][preState2[s] + 8] + alpha[i][preState2[s]];

			x = x > y ? x : y;
			y = gamma[i][preState3[s] + 12] + alpha[i][preState3[s]];

			alpha[i + 1][s] = x > y ? x : y;
		}

	//-------------------------- beta node metrics calculation -----------------------------//

	beta[N - termi - 1][preState0[0]] = gamma[N - 2][preState0[0]];
	beta[N - termi - 1][preState1[0]] = gamma[N - 2][preState1[0] + 4];
	beta[N - termi - 1][preState2[0]] = gamma[N - 2][preState2[0] + 8];
	beta[N - termi - 1][preState3[0]] = gamma[N - 2][preState3[0] + 12];

	for (i = N - termi - 1; i > 0; i--)
		for (s = 0; s<numOfStates; s++)
		{
			x = gamma[i][s] + beta[i][nextState0[s]];
			y = gamma[i][s + 4] + beta[i][nextState1[s]];

			x = x > y ? x : y;
			y = gamma[i][s + 8] + beta[i][nextState2[s]];

			x = x > y ? x : y;
			y = gamma[i][s + 12] + beta[i][nextState3[s]];

			beta[i - 1][s] = x > y ? x : y;			
		}

	//-------------------------- LLR calculation ----------------------------//

	LL0 = gamma[0][0] + beta[0][0];
	LLR1[0] = gamma[0][4] + beta[0][1] - LL0;
	LLR2[0] = gamma[0][8] + beta[0][2] - LL0;
	LLR3[0] = gamma[0][12] + beta[0][3] - LL0;

	for (i = 1; i < N - termi; i++)
	{
		x = alpha[i][0] + gamma[i][0] + beta[i][0];
		for (s = 1; s < numOfStates; s++)
		{
			y = alpha[i][s] + gamma[i][s] + beta[i][nextState0[s]];
			x = x > y ? x : y;
		}
		LL0 = x;

		x = alpha[i][0] + gamma[i][4] + beta[i][1];
		for (s = 1; s < numOfStates; s++)
		{
			y = alpha[i][s] + gamma[i][s + 4] + beta[i][nextState1[s]];
			x = x > y ? x : y;
		}
		LLR1[i] = x - LL0;

		x = alpha[i][0] + gamma[i][8] + beta[i][2];
		for (s = 1; s < numOfStates; s++)
		{
			y = alpha[i][s] + gamma[i][s + 8] + beta[i][nextState2[s]];
			x = x > y ? x : y;
		}
		LLR2[i] = x - LL0;

		x = alpha[i][0] + gamma[i][12] + beta[i][3];
		for (s = 1; s < numOfStates; s++)
		{
			y = alpha[i][s] + gamma[i][s + 12] + beta[i][nextState3[s]];
			x = x > y ? x : y;
		}
		LLR3[i] = x - LL0;
	}
}

void sisoDec2(int *deInt, double *noisy1I, double *susie1Q, double *noisy2I, double *susie2Q, double *LLR1, double *LLR2, double *LLR3, double *La1, double *La2, double *La3)
{
	// max-log-APP algorithm
	char numOfStates = 1, s;
	int i;
	double x, y, LL0;

	//-------------------------- gamma branch metrics calculation ----------------------------//

	for (i = 0; i < N - termi; i++)
	{

		if (deInt[i] & 1) // if in index i there exists a symbol from ENC #2
		{
			for (s = 0; s < numOfStates; s++)
			{
				x = H11[0][deInt[i]] + H21[0][deInt[i]] * nSymbol0[s << 1] - H21[1][deInt[i]] * nSymbol0[s << 1 | 1];
				y = x * (2 * noisy1I[i] - x);
				x = H11[1][deInt[i]] + H21[0][deInt[i]] * nSymbol0[s << 1 | 1] + H21[1][deInt[i]] * nSymbol0[s << 1];
				y += x * (2 * susie1Q[i] - x);
				gamma[i][s] = y / variance;

				x = H12[0][deInt[i]] + H22[0][deInt[i]] * nSymbol0[s << 1] - H22[1][deInt[i]] * nSymbol0[s << 1 | 1];
				y = x * (2 * noisy2I[i] - x);
				x = H12[1][deInt[i]] + H22[0][deInt[i]] * nSymbol0[s << 1 | 1] + H22[1][deInt[i]] * nSymbol0[s << 1];
				y += x * (2 * susie2Q[i] - x);
				gamma[i][s] += y / variance;


				x = -H11[1][deInt[i]] + H21[0][deInt[i]] * nSymbol1[s << 1] - H21[1][deInt[i]] * nSymbol1[s << 1 | 1];
				y = x * (2 * noisy1I[i] - x);
				x = H11[0][deInt[i]] + H21[0][deInt[i]] * nSymbol1[s << 1 | 1] + H21[1][deInt[i]] * nSymbol1[s << 1];
				y += x * (2 * susie1Q[i] - x);
				gamma[i][s + 4] = La1[i] + y / variance;

				x = -H12[1][deInt[i]] + H22[0][deInt[i]] * nSymbol1[s << 1] - H22[1][deInt[i]] * nSymbol1[s << 1 | 1];
				y = x * (2 * noisy2I[i] - x);
				x = H12[0][deInt[i]] + H22[0][deInt[i]] * nSymbol1[s << 1 | 1] + H22[1][deInt[i]] * nSymbol1[s << 1];
				y += x * (2 * susie2Q[i] - x);
				gamma[i][s + 4] += y / variance;


				x = -H11[0][deInt[i]] + H21[0][deInt[i]] * nSymbol2[s << 1] - H21[1][deInt[i]] * nSymbol2[s << 1 | 1];
				y = x * (2 * noisy1I[i] - x);
				x = -H11[1][deInt[i]] + H21[0][deInt[i]] * nSymbol2[s << 1 | 1] + H21[1][deInt[i]] * nSymbol2[s << 1];
				y += x * (2 * susie1Q[i] - x);
				gamma[i][s + 8] = La2[i] + y / variance;

				x = -H12[0][deInt[i]] + H22[0][deInt[i]] * nSymbol2[s << 1] - H22[1][deInt[i]] * nSymbol2[s << 1 | 1];
				y = x * (2 * noisy2I[i] - x);
				x = -H12[1][deInt[i]] + H22[0][deInt[i]] * nSymbol2[s << 1 | 1] + H22[1][deInt[i]] * nSymbol2[s << 1];
				y += x * (2 * susie2Q[i] - x);
				gamma[i][s + 8] += y / variance;


				x = H11[1][deInt[i]] + H21[0][deInt[i]] * nSymbol3[s << 1] - H21[1][deInt[i]] * nSymbol3[s << 1 | 1];
				y = x * (2 * noisy1I[i] - x);
				x = -H11[0][deInt[i]] + H21[0][deInt[i]] * nSymbol3[s << 1 | 1] + H21[1][deInt[i]] * nSymbol3[s << 1];
				y += x * (2 * susie1Q[i] - x);
				gamma[i][s + 12] = La3[i] + y / variance;

				x = H12[1][deInt[i]] + H22[0][deInt[i]] * nSymbol3[s << 1] - H22[1][deInt[i]] * nSymbol3[s << 1 | 1];
				y = x * (2 * noisy2I[i] - x);
				x = -H12[0][deInt[i]] + H22[0][deInt[i]] * nSymbol3[s << 1 | 1] + H22[1][deInt[i]] * nSymbol3[s << 1];
				y += x * (2 * susie2Q[i] - x);
				gamma[i][s + 12] += y / variance;
			}
		}
		else
		{
			for (s = 0; s < numOfStates; s++)
			{
				gamma[i][s] = 0.;
				gamma[i][s + 4] = La1[i];
				gamma[i][s + 8] = La2[i];
				gamma[i][s + 12] = La3[i];
			}
		}

		if (i == 0) numOfStates = 4;
	}

	// trellis termination
	x = H11[0][N - 1] + H21[0][N - 1];
	y = x * (2 * noisy1I[N - 1] - x);
	x = H11[1][N - 1] + H21[1][N - 1];
	y += x * (2 * susie1Q[N - 1] - x);
	gamma[N - 2][preState0[0]] = y / variance;

	x = H12[0][N - 1] + H22[0][N - 1];
	y = x * (2 * noisy2I[N - 1] - x);
	x = H12[1][N - 1] + H22[1][N - 1];
	y += x * (2 * susie2Q[N - 1] - x);
	gamma[N - 2][preState0[0]] += y / variance;


	x = -H11[1][N - 1] - H21[0][N - 1];
	y = x * (2 * noisy1I[N - 1] - x);
	x = H11[0][N - 1] - H21[1][N - 1];
	y += x * (2 * susie1Q[N - 1] - x);
	gamma[N - 2][preState1[0] + 4] = y / variance;

	x = -H12[1][N - 1] - H22[0][N - 1];
	y = x * (2 * noisy2I[N - 1] - x);
	x = H12[0][N - 1] - H22[1][N - 1];
	y += x * (2 * susie2Q[N - 1] - x);
	gamma[N - 2][preState1[0] + 4] += y / variance;


	x = -H11[0][N - 1] - H21[1][N - 1];
	y = x * (2 * noisy1I[N - 1] - x);
	x = -H11[1][N - 1] + H21[0][N - 1];
	y += x * (2 * susie1Q[N - 1] - x);
	gamma[N - 2][preState2[0] + 8] = y / variance;

	x = -H12[0][N - 1] - H22[1][N - 1];
	y = x * (2 * noisy2I[N - 1] - x);
	x = -H12[1][N - 1] + H22[0][N - 1];
	y += x * (2 * susie2Q[N - 1] - x);
	gamma[N - 2][preState2[0] + 8] += y / variance;


	x = H11[1][N - 1] + H21[1][N - 1];
	y = x * (2 * noisy1I[N - 1] - x);
	x = -H11[0][N - 1] - H21[0][N - 1];
	y += x * (2 * susie1Q[N - 1] - x);
	gamma[N - 2][preState3[0] + 12] = y / variance;

	x = H12[1][N - 1] + H22[1][N - 1];
	y = x * (2 * noisy2I[N - 1] - x);
	x = -H12[0][N - 1] - H22[0][N - 1];
	y += x * (2 * susie2Q[N - 1] - x);
	gamma[N - 2][preState3[0] + 12] += y / variance;


	//-------------------------- alpha node metrics calculation ----------------------------//

	alpha[1][0] = gamma[0][0];
	alpha[1][1] = gamma[0][4];
	alpha[1][2] = gamma[0][8];
	alpha[1][3] = gamma[0][12];

	for (i = 1; i < N - termi - 1; i++)
		for (s = 0; s < numOfStates; s++)
		{
			x = gamma[i][preState0[s]] + alpha[i][preState0[s]];
			y = gamma[i][preState1[s] + 4] + alpha[i][preState1[s]];

			x = x > y ? x : y;
			y = gamma[i][preState2[s] + 8] + alpha[i][preState2[s]];

			x = x > y ? x : y;
			y = gamma[i][preState3[s] + 12] + alpha[i][preState3[s]];

			alpha[i + 1][s] = x > y ? x : y;
		}


	//-------------------------- beta node metrics calculation -----------------------------//

	beta[N - termi - 1][preState0[0]] = gamma[N - 2][preState0[0]];
	beta[N - termi - 1][preState1[0]] = gamma[N - 2][preState1[0] + 4];
	beta[N - termi - 1][preState2[0]] = gamma[N - 2][preState2[0] + 8];
	beta[N - termi - 1][preState3[0]] = gamma[N - 2][preState3[0] + 12];

	for (i = N - termi - 1; i > 0; i--)
		for (s = 0; s < numOfStates; s++)
		{
			x = gamma[i][s] + beta[i][nextState0[s]];
			y = gamma[i][s + 4] + beta[i][nextState1[s]];

			x = x > y ? x : y;
			y = gamma[i][s + 8] + beta[i][nextState2[s]];

			x = x > y ? x : y;
			y = gamma[i][s + 12] + beta[i][nextState3[s]];

			beta[i - 1][s] = x > y ? x : y;
		}


	//-------------------------- LLR calculation ----------------------------//

	LL0 = gamma[0][0] + beta[0][0];
	LLR1[0] = gamma[0][4] + beta[0][1] - LL0;
	LLR2[0] = gamma[0][8] + beta[0][2] - LL0;
	LLR3[0] = gamma[0][12] + beta[0][3] - LL0;

	for (i = 1; i < N - termi; i++)
	{
		x = alpha[i][0] + gamma[i][0] + beta[i][0];
		for (s = 1; s < numOfStates; s++)
		{
			y = alpha[i][s] + gamma[i][s] + beta[i][nextState0[s]];
			x = x > y ? x : y;
		}
		LL0 = x;

		x = alpha[i][0] + gamma[i][4] + beta[i][1];
		for (s = 1; s < numOfStates; s++)
		{
			y = alpha[i][s] + gamma[i][s + 4] + beta[i][nextState1[s]];
			x = x > y ? x : y;
		}
		LLR1[i] = x - LL0;

		x = alpha[i][0] + gamma[i][8] + beta[i][2];
		for (s = 1; s < numOfStates; s++)
		{
			y = alpha[i][s] + gamma[i][s + 8] + beta[i][nextState2[s]];
			x = x > y ? x : y;
		}
		LLR2[i] = x - LL0;

		x = alpha[i][0] + gamma[i][12] + beta[i][3];
		for (s = 1; s < numOfStates; s++)
		{
			y = alpha[i][s] + gamma[i][s + 12] + beta[i][nextState3[s]];
			x = x > y ? x : y;
		}
		LLR3[i] = x - LL0;
	}
}

int main()
{
	bool test = true, genie = true;
	double *GN, *Isys, *Ipar, *noisy1I, *noisy2I, *Qsys, *Qpar, *susie1Q, *susie2Q, *IntNoisy1I, *IntSusie1Q, *IntNoisy2I, *IntSusie2Q,
	       LLR0, *LLR1, *LLR2, *LLR3,
	       *La1DEC1, *La2DEC1, *La3DEC1, *La1DEC2, *La2DEC2, *La3DEC2,
	       sigma, EbNo, EsNo,
	       wNor[256], fNor[256];
	unsigned seed[R];
	char iters, *bitStream, *systematicz, *paritiez, hd;
	int kNor[256], i, maxIters = 8, trials = 0, bitErrors = 0, frameErrors = 0, *Int, *deInt;

	//----------------------------------------------------------------------------//

	alpha = new double*[N - termi];
	for (i = 0; i < N - termi; i++)
		alpha[i] = new double[4];
	alpha[0][0] = 0.;

	beta = new double*[N - termi + 1];
	for (i = 0; i < N - termi + 1; i++)
		beta[i] = new double[4];
	beta[N - termi][0] = 0.;

	gamma = new double*[N - 1];
	for (i = 0; i < N - 1; i++)
		gamma[i] = new double[16];

	// 1st channel
	H11 = new double*[2];   // two lines
	H11[0] = new double[N]; // N columns
	H11[1] = new double[N]; // H11[0] holds the real parts, H11[1] the imaginary ones

	// 2nd channel
	H12 = new double*[2];   
	H12[0] = new double[N]; 
	H12[1] = new double[N]; 

	// 3rd channel
	H21 = new double*[2];   
	H21[0] = new double[N]; 
	H21[1] = new double[N]; 

	// 4th channel
	H22 = new double*[2];
	H22[0] = new double[N];
	H22[1] = new double[N];

	bitStream = new char[(N - termi) << 1];
	systematicz = new char[N];
	paritiez = new char[N];
	Int = new int[N - termi];
	deInt = new int[N - termi];

	LLR1 = new double[N - termi];
	LLR2 = new double[N - termi];
	LLR3 = new double[N - termi];

	La1DEC1 = new double[N - termi];
	La2DEC1 = new double[N - termi];
	La3DEC1 = new double[N - termi];
	La1DEC2 = new double[N - termi];
	La2DEC2 = new double[N - termi];
	La3DEC2 = new double[N - termi];

	GN = new double[N << 2];
	Isys = new double[N]; Ipar = new double[N];
	noisy1I = new double[N]; // real part of first receive antenna
	noisy2I = new double[N]; // real part of second receive antenna
	Qsys = new double[N];
	Qpar = new double[N];
	susie1Q = new double[N]; // imaginary part of first receive antenna
	susie2Q = new double[N]; // imaginary part of second receive antenna
	IntNoisy1I = new double[N];
	IntSusie1Q = new double[N];
	IntNoisy2I = new double[N];
	IntSusie2Q = new double[N];

	//----------------------------------------------------------------------------//

	// random number generator initialization
	for (i = 0; i < R; i++)
		seed[i] = 2345 + i;
	InitWELLRNG1024a(seed);
	zigSet(kNor, wNor, fNor);

	EbNo = 3.; // Eb/No in dB
	EsNo = EbNo;
	variance = pow(10, -EsNo / 10) / 2.;
	sigma = sqrt(variance);
	trellis();	

	//----------------------------------------------------------------------------//

	while (bitErrors < 1000)
	{
		if ((trials & 31) == 0)
			createInt(Int, deInt); // creates a different random interleaver every 32 trials

		source(bitStream);
		encode(bitStream, systematicz, paritiez, Int, deInt);

		modulate(systematicz, Isys, Qsys);
		modulate(paritiez, Ipar, Qpar);

		// iid fading
		// create 1st channel
		wgn(H11[0], N, kNor, wNor, fNor); // real part
		wgn(H11[1], N, kNor, wNor, fNor); // imaginary part

		// create 2nd channel
		wgn(H12[0], N, kNor, wNor, fNor); // real part
		wgn(H12[1], N, kNor, wNor, fNor); // imaginary part

		// create 3rd channel
		wgn(H21[0], N, kNor, wNor, fNor); // real part
		wgn(H21[1], N, kNor, wNor, fNor); // imaginary part

		// create 4th channel
		wgn(H22[0], N, kNor, wNor, fNor); // real part
		wgn(H22[1], N, kNor, wNor, fNor); // imaginary part
		
		// normalize, i.e. E[|h|^2] = 1
		for (i = 0; i < N; i++) { H11[0][i] /= sqrt(2.); H11[1][i] /= sqrt(2.); H12[0][i] /= sqrt(2.); H12[1][i] /= sqrt(2.); }
		for (i = 0; i < N; i++) { H21[0][i] /= sqrt(2.); H21[1][i] /= sqrt(2.); H22[0][i] /= sqrt(2.); H22[1][i] /= sqrt(2.); }
		
		/*
		// block fading		
		wgn(H11[0], 1, kNor, wNor, fNor); // real part
		wgn(H11[1], 1, kNor, wNor, fNor); // imaginary part

										  // create 2nd channel
		wgn(H12[0], 1, kNor, wNor, fNor); // real part
		wgn(H12[1], 1, kNor, wNor, fNor); // imaginary part

										  // create 3rd channel
		wgn(H21[0], 1, kNor, wNor, fNor); // real part
		wgn(H21[1], 1, kNor, wNor, fNor); // imaginary part

										  // create 4th channel
		wgn(H22[0], 1, kNor, wNor, fNor); // real part
		wgn(H22[1], 1, kNor, wNor, fNor); // imaginary part

		// normalize, i.e. E[|h|^2] = 1		
		H11[0][0] /= sqrt(2.); H11[1][0] /= sqrt(2.);
		for (i = 1; i < N; i++) { H11[0][i] = H11[0][0]; H11[1][i] = H11[1][0]; }

		H12[0][0] /= sqrt(2.); H12[1][0] /= sqrt(2.);
		for (i = 1; i < N; i++) { H12[0][i] = H12[0][0]; H12[1][i] = H12[1][0]; }

		H21[0][0] /= sqrt(2.); H21[1][0] /= sqrt(2.);
		for (i = 1; i < N; i++) { H21[0][i] = H21[0][0]; H21[1][i] = H21[1][0]; }

		H22[0][0] /= sqrt(2.); H22[1][0] /= sqrt(2.);
		for (i = 1; i < N; i++) { H22[0][i] = H22[0][0]; H22[1][i] = H22[1][0]; }
		*/
		
		wgn(GN, N << 2, kNor, wNor, fNor);

		for (i = 0; i < N; i++)
		{
			// 1st receive antenna
			// Re{.}
			noisy1I[i] = Isys[i] * H11[0][i] - Qsys[i] * H11[1][i];  // due to 1st path
			noisy1I[i] += Ipar[i] * H21[0][i] - Qpar[i] * H21[1][i]; // due to 2nd path
			noisy1I[i] += sigma*GN[i];                               // add noise

			// Im{.}
			susie1Q[i] = Isys[i] * H11[1][i] + Qsys[i] * H11[0][i];
			susie1Q[i] += Ipar[i] * H21[1][i] + Qpar[i] * H21[0][i];
			susie1Q[i] += sigma*GN[N + i];


			// 2nd receive antenna
			// Re{.}
			noisy2I[i] = Isys[i] * H12[0][i] - Qsys[i] * H12[1][i];  // due to 1st path
			noisy2I[i] += Ipar[i] * H22[0][i] - Qpar[i] * H22[1][i]; // due to 2nd path
			noisy2I[i] += sigma*GN[2*N + i];                         // add noise

			// Im{.}
			susie2Q[i] = Isys[i] * H12[1][i] + Qsys[i] * H12[0][i];
			susie2Q[i] += Ipar[i] * H22[1][i] + Qpar[i] * H22[0][i];
			susie2Q[i] += sigma*GN[3*N + i];
		}

		for (i = 0; i < N - termi; i++)
		{
			IntNoisy1I[Int[i]] = noisy1I[i];
			IntSusie1Q[Int[i]] = susie1Q[i];

			IntNoisy2I[Int[i]] = noisy2I[i];
			IntSusie2Q[Int[i]] = susie2Q[i];
		}
		for (; i < N; i++)
		{
			IntNoisy1I[i] = noisy1I[i];
			IntSusie1Q[i] = susie1Q[i];

			IntNoisy2I[i] = noisy2I[i];
			IntSusie2Q[i] = susie2Q[i];
		}

		// turbo decoder
		for (i = 0; i < N - termi; i++) { La1DEC1[i] = La2DEC1[i] = La3DEC1[i] = 0.; }
		for (iters = 0; iters < maxIters; iters++)
		{
			sisoDec1(noisy1I, susie1Q, noisy2I, susie2Q, LLR1, LLR2, LLR3, La1DEC1, La2DEC1, La3DEC1);

			for (i = 0; i < N - termi; i++)
			{
				La1DEC2[Int[i]] = LLR1[i] - La1DEC1[i];
				La2DEC2[Int[i]] = LLR2[i] - La2DEC1[i];
				La3DEC2[Int[i]] = LLR3[i] - La3DEC1[i];
			}

			sisoDec2(deInt, IntNoisy1I, IntSusie1Q, IntNoisy2I, IntSusie2Q, LLR1, LLR2, LLR3, La1DEC2, La2DEC2, La3DEC2);

			if (genie)
			{
				if (iters >= 3 && iters < maxIters - 1)
				{
					for (i = 0; i < N - termi; i++)
					{
						LLR0 = 0.;
						hd = 0;
						if (LLR1[i] > LLR0) { LLR0 = LLR1[i]; hd = 1; }
						if (LLR2[i] > LLR0) { LLR0 = LLR2[i]; hd = 2; }
						if (LLR3[i] > LLR0) hd = 3;
						if (hd >> 1 ^ bitStream[deInt[i] << 1]) break;
						if (hd & 1 ^ bitStream[deInt[i] << 1 | 1]) break;
					}
					if (i == N - termi) goto there;
				}
			}

			if (iters < maxIters - 1)
				for (i = 0; i < N - termi; i++)
				{
					La1DEC1[deInt[i]] = LLR1[i] - La1DEC2[i];
					La2DEC1[deInt[i]] = LLR2[i] - La2DEC2[i];
					La3DEC1[deInt[i]] = LLR3[i] - La3DEC2[i];
				}
		}

		// error counting
		for (i = 0; i < N - termi; i++)
		{
			LLR0 = 0.;
			hd = 0;
			if (LLR1[i] > LLR0) { LLR0 = LLR1[i]; hd = 1; }
			if (LLR2[i] > LLR0) { LLR0 = LLR2[i]; hd = 2; }
			if (LLR3[i] > LLR0) hd = 3;

			if (hd >> 1 ^ bitStream[deInt[i] << 1])
			{
				bitErrors++;
				if (test) {	frameErrors++; test = false; }
			}
			if (hd & 1 ^ bitStream[deInt[i] << 1 | 1])
			{
				bitErrors++;
				if (test) {	frameErrors++; test = false; }
			}
		}
		if (test == false)
			test = true;

	there:
		trials++;
		printf("%2d %d %d %d %.2e %.2e\n", iters, trials, bitErrors, frameErrors, double(bitErrors) / double(N - termi) / double(trials) / 2., double(frameErrors) / double(trials));

	}

	printf("\nEbNo=%.2f bitErrors=%d BER=%.2e FER=%.2e\n", EbNo, bitErrors, double(bitErrors) / double(N - termi) / double(trials) / 2., double(frameErrors) / double(trials));

	system("pause");
	return 0;
}
