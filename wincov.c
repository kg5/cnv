#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#define SN 111
#define EL 2500000
#define L 25000000

char sl[10000][256];
int exon[25000][2], winl, s, Len, cumsum[25000];
int *cov, dep[L], sn, cn;
float winMeanDep(int *, int, int);

int main(int argc, char *argv[])
{
	if(argc != 5)
	{
		printf("Usage: %s chr[chr*] sample_number window_length shift\n", argv[0]);
		exit(1);
	}
	sn = atoi(argv[2]);
	winl = atoi(argv[3]);
	s = atoi(argv[4]);
	int en = 0, len = 0;
	int i = 0, j = 0;
	FILE *rf, *covf, *outf, *locf, *mdf;
	char regionf[256], locn[256], outn[256];
	sprintf(regionf, "%s.loc", argv[1]);
	if((rf = fopen(regionf, "r")) == NULL)
	{
		printf("Can't open file: %s\n", regionf);
		exit(1);
	}
	while(fscanf(rf, "%d", &exon[en][0]) == 1)
		fscanf(rf, "%d", &exon[en++][1]);
	if(fclose(rf) != 0)
	{
		printf("Error closing file:%s\n", regionf);
		exit(1);
	}
	for(i = 0; i < en; i++)
	{
		Len += exon[i][1] - exon[i][0] + 1;
		cumsum[i+1] = Len;
	}
	sprintf(outn, "%s_W%dS%d.dep", argv[1], winl, s);
	if((outf = fopen(outn, "w")) == NULL)
	{
		printf("Can't open file: %s\n", outn);
		exit(1);
	}
	sprintf(locn, "%s.loc", outn);
	if((locf = fopen(locn, "w")) == NULL)
	{
		printf("Can't open file: %s\n", locn);
		exit(1);
	}
	char covn[256], mdfn[256];
	sprintf(covn, "%s.dep", argv[1]);
	if((covf = fopen(covn, "r")) == NULL)
	{
		printf("Can't open file %s\n", covn);
		exit(1);
	}
	sprintf(mdfn, "%s_W%dS%d.mdep", argv[1], winl,s);
	if((mdf = fopen(mdfn, "w")) == NULL)
	{
		printf("Can't write file: %s\n", mdfn);
		exit(1);
	}
	while(fscanf(covf, "%d", dep+len) == 1)
		len++;
	if(fclose(covf) != 0)
	{
		printf("Error closing file: %s\n", covn);
		exit(1);
	}
	if(len > L)
	{
		printf("length %d exceeds preset:%d\n", len, L);
		exit(1);
	}
	cn = len / sn;
	if(cn != Len)
	{
		printf("cov length not equal capture region\n");
		exit(1);
	}
	for(i = 0; i < sn; i++)
	{
		cov = &dep[i*cn];
		float md = 0;
		int wn = 0;
		for(j = 0; j < en; j++)
		{
			int k;
			for(k = exon[j][0]; k + winl < exon[j][1]; k += s)
			{
			  int winb = cumsum[j] + k - exon[j][0] + 1;
				int kend = k + winl - 1;
				if(k+s+winl > exon[j][1])
					kend = exon[j][1];
				float meanD = winMeanDep(cov, winb, kend - k + 1);
				md += meanD;
				wn++;
				fprintf(outf, "%.2f\t", meanD);
				if(i == 0)
					fprintf(locf, "%d\t%d\n", k, kend);
			}
		}
		fprintf(outf, "\n");
		md /= (wn+0.0);
		fprintf(mdf, "%f\n", md);
	}
	if(fclose(outf) != 0 || fclose(locf) != 0)
	{
		printf("Error closing file\n");
		exit(1);
	}
	return 0;
}

float winMeanDep(int cov[], int loc, int winl)
{
	float depth = 0;
	int k;
	for(k = loc; k < loc + winl; k++)
		depth += cov[k];
	depth /= (winl + 0.0);
	return(depth);
}
