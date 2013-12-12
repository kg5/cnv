#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#define SN 111
#define EL 2500000
#define L 25000000

char lib[SN][256];
int exon[25000][2], winl, s, Len, cumsum[25000];
int cov[L], sn, cn;
int gc[250000000];
float winMeanDep(int *, int, int);
float gcPr(int gc[], int, int);

int main(int argc, char *argv[])
{
	if(argc != 4)
	{
		printf("Usage: %s chr[chr*] window_length shift\n", argv[0]);
		exit(1);
	}
	winl = atoi(argv[2]);
	s = atoi(argv[3]);
	int en = 0;
	int i = 0, j = 0;
	FILE *rf, *outf, *locf, *gcf, *sl, *mdf;
	char regionf[256], locn[256], gcfn[256], outn[256];
	sprintf(regionf, "%s.loc", argv[1]);
	sprintf(gcfn, "/vol6/home/bgi_caofei/program/cnv/hg19gc/%s.gc",argv[1]);
	if((gcf = fopen(gcfn, "r")) == NULL)
	{
		printf("Can't open file:%s\n", gcfn);
		exit(1);
	}
	while(fscanf(gcf,"%1d",gc+i) == 1)
		i++;
	if((rf = fopen(regionf, "r")) == NULL)
	{
		printf("Can't open file: %s\n", regionf);
		exit(1);
	}
	while(fscanf(rf, "%d", &exon[en][0]) == 1)
		fscanf(rf, "%d", &exon[en++][1]);
	if((sl = fopen("sample.list", "r")) == NULL)
	{
		printf("Can't open sample.list");
		exit(1);
	}
	char gender[5], bamp[1024];
	while(fscanf(sl, "%s\t%s\t%s", gender,lib[sn],bamp) == 1)
		sn++;
	if(fclose(rf) != 0 || fclose(gcf) != 0 || fclose(sl) != 0)
	{
		printf("Error closing file:%s\n", regionf);
		exit(1);
	}
	for(i = 0; i < en; i++)
	{
		Len += exon[i][1] - exon[i][0] + 1;
		cumsum[i+1] = Len;
	}
	sprintf(outn, "%s_W%dS%d.cov", argv[1], winl, s);
	if((outf = fopen(outn, "w")) == NULL)
	{
		printf("Can't open file: %s\n", outn);
		exit(1);
	}
	sprintf(locn, "%s_W%dS%d.loc", argv[1], winl, s);
	if((locf = fopen(locn, "w")) == NULL)
	{
		printf("Can't open file: %s\n", locn);
		exit(1);
	}
	char mdfn[256];
	sprintf(mdfn, "%s_W%dS%d.mdep", argv[1], winl,s);
	if((mdf = fopen(mdfn, "w")) == NULL)
	{
		printf("Can't write file: %s\n", mdfn);
		exit(1);
	}
	if(Len > L)
	{
		printf("length %d exceeds preset:%d\n", Len, L);
		exit(1);
	}
	for(i = 0; i < sn; i++)
	{
		FILE *covf, *gcdf;
		char covn[256], gcdn[256];
		int len = 0;
	  sprintf(covn, "COV/%s/%s.cov", lib[i], argv[1]);
	  sprintf(gcdn, "COV/%s/%s_W%dS%d.gcd", lib[i], argv[1],winl,s);
		if((gcdf = fopen(gcdn, "w")) == NULL)
		{
			printf("Can't write file:%s\n", gcdn);
			exit(1);
		}
		if((covf = fopen(covn, "r")) == NULL)
		{
			printf("Can't open file %s\n", covn);
			exit(1);
		}
		while(fscanf(covf, "%d", cov+len) == 1)
			len++;
		if(len != Len)
		{
			printf("sample %s, %s: cov length not equal capture region\n", lib[i], argv[1]);
			exit(1);
		}
		float md = 0;
		int wn = 0;
		for(j = 0; j < en; j++)
		{
			int k;
			for(k = exon[j][0]; k + winl < exon[j][1]; k += s)
			{
			  int winb = cumsum[j] + k - exon[j][0];
				int wend = winb + winl - 1;
				if(k+s+winl-1 >= exon[j][1])
					wend = exon[j][1];
				float meanD = winMeanDep(cov, winb, wend - winb + 1);
				float gcP = gcPr(gc, k, wend - winb + 1);
				md += meanD;
				wn++;
				fprintf(outf, "%.2f\t", meanD);
				fprintf(gcdf, "%f\t%f\n", gcP, meanD);
				if(i == 0)
					fprintf(locf, "%d\t%d\n", k, k+wend-winb);
			}
		}
		fprintf(outf, "\n");
		md /= (wn+0.0);
		fprintf(mdf, "%f\n", md);
		if(fclose(covf) != 0 || fclose(gcdf) != 0)
		{
			printf("Error closing file\n");
			exit(1);
		}
	}
	if(fclose(outf) != 0 || fclose(locf) != 0)
	{
		printf("Error closing file\n");
		exit(1);
	}
	return 0;
}

float winMeanDep(int cov[], int loc, int wlen)
{
	float mdepth = 0;
	int k;
	for(k = loc; k < loc + wlen; k++)
		mdepth += cov[k];
	mdepth /= (wlen + 0.0);
	return(mdepth);
}

float gcPr(int gc[], int begin, int wlen)
{
	float gcp = 0;
	int k, Ncount = 0;
	for(k = begin; k < begin + wlen; k++)
	{
		if(gc[k] < 5)
			gcp += gc[k];
		else
			Ncount++;
	}
	gcp /= (wlen - Ncount + 0.0);
	return(gcp);
}
