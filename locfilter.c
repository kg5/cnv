#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#define SN 111 
#define EL 1000000
#define DT 5
#define RN 50000
#define MINRL 20

char sl[10000][256];
int exon[RN][2], winl, s, Len, cumsum[25000];//nexon[RN][2];
int cov[SN][EL], reg[EL];
int isEmpty(int cov[SN][EL], int sn, int coln);

int main(int argc, char *argv[])
{
	if(argc != 6)
	{
		printf("Usage: %s sample_list cov_dir bed_dir chr[chr*] outdir\n", argv[0]);
		exit(1);
	}
	FILE *slf;
	if((slf = fopen(argv[1], "r")) == NULL)
	{
		printf("Can't open file: %s\n", argv[1]);
		exit(1);
	}
	int sn = 0, en = 0; // nen = 0;
	char sampleID[256], bamPath[512];
	while(fscanf(slf, "%s", sampleID) == 1)
	{
		fscanf(slf, "%s", sl[sn]);
		fscanf(slf, "%s", bamPath);
		sn++;
	}
	if(sn > SN)
	{
		printf("sample number exceeds preset Maximum!\n");
		exit(1);
	}
	if(fclose(slf) != 0)
	{
		printf("Error closing file\n");
		exit(1);
	}
	int i = 0, j = 0;
	FILE *rf, *covf, *outf, *locf, *xregf;
	char regionf[256], chrn[128], locn[256], outn[256], xregfn[256];
	sprintf(regionf, "%s/%s", argv[3], argv[4]);
	if((rf = fopen(regionf, "r")) == NULL)
	{
		printf("Can't open file: %s\n", regionf);
		exit(1);
	}
	while(fscanf(rf, "%d", &exon[en][0]) == 1)
		fscanf(rf, "%d", &exon[en++][1]);
	for(i = 0; i < en; i++)
	{
		Len += exon[i][1] - exon[i][0] + 1;
		cumsum[i+1] = Len;
	}
	if(Len > EL)
	{
		printf("capture region length > preset length.\nDo chage EL\n");
//		exit(1);
	}
	sprintf(xregfn, "%s/%s.loc", argv[5], argv[4]);
	sprintf(outn, "%s/%s_deleted.cov", argv[5],argv[4]);
	if((xregf = fopen(xregfn, "w")) == NULL)
	{
		printf("Can't open file:%s\n", xregfn);
		exit(1);
	}
	if((outf = fopen(outn, "w")) == NULL)
	{
		printf("Can't open file: %s\n", outn);
		exit(1);
	}
	if(fclose(rf) != 0)
	{
		printf("Error closing file:%s\n", regionf);
		exit(1);
	}
	for(i = 0; i < sn; i++)
	{
		char covn[256];
		sprintf(covn, "%s/%s/%s.cov", argv[2], sl[i], argv[4]);
		if((covf = fopen(covn, "r")) == NULL)
		{
			printf("Can't open file %s\n", covn);
			exit(1);
		}
		int rc, len = 0, exn = 0;
		while(fscanf(covf, "%d", &cov[i][len]) == 1)
			len++;
		if(len != Len)
		{
			printf("cov length not equal capture region: sample %d", i+1);
			exit(1);
		}
		if(fclose(covf) != 0)
		{
			printf("Error closing file: %s\n", covn);
			exit(1);
		}
	}
//	printf("check mem consumption:\nplease input something:\n");
//	char inp[256];
//	while(fscanf(stdin, "%s", inp) == 1)
//		printf("%s\n", inp);
//	exit(1);
	for(i = 0; i < en; i++)
	{
		int regLen = exon[i][1] - exon[i][0] + 1;
		for(j = 0; j < regLen; j++)
		{
			int loca = cumsum[i] + j;
		  if(isEmpty(cov, sn, loca))
			{
				reg[j] = 0;
				int k;
				fprintf(outf, "loc_%d:\t", exon[i][0]+j);
				for(k = 0; k < sn; k++)
				{
					fprintf(outf, "%d ", cov[k][loca]);
					cov[k][loca] = -1;
				}
				fprintf(outf,"\n");
			}
			else
				reg[j] = 1;
		}
		int subreg[1000][2], subn = 0, isin = 0;
		for(j = 0; j < regLen; j++)
		{
			if(reg[j] == 1) 
			{
				if(isin == 0)
				{
					subreg[subn][0] = exon[i][0] + j;
					isin = 1;
				}
			}
			else
			{
				if(isin == 1)
				{
					subreg[subn][1] = exon[i][0] + j - 1;
					subn++;
					isin = 0;
				}
			}
		}
		if(reg[regLen-1] == 1 && isin == 1)
		{
			subreg[subn][1] = exon[i][1];
			subn++;
		}
		for(j = 0; j < subn; j++)
		{
			int subregLen = subreg[j][1] - subreg[j][0] + 1;
			if(subregLen > MINRL)
			{
				fprintf(xregf, "%d\t%d\n", subreg[j][0], subreg[j][1]);
//				nexon[nen][0] = subreg[j][0];
//				nexon[nen][1] = subreg[j][1];
//				nen++;
			}
			else
			{
				int k;
				for(k = subreg[j][0]; k <= subreg[j][1]; k++)
				{
					int location = cumsum[i] + k - exon[i][0];
					int m;
					fprintf(outf, "loc_%d:\t", k);
					for(m = 0; m < sn; m++)
					{
						fprintf(outf, "%d ", cov[m][location]);
						cov[m][location] = -1;
					}
					fprintf(outf, "\n");
				}
			}
		}
	}
  for(i = 0; i < sn; i++)
	{
		char depfn[256];
		sprintf(depfn, "%s/%s/%s.dep", argv[2], sl[i], argv[4]);
		FILE *depf;
		if((depf = fopen(depfn, "w")) == NULL)
		{
			printf("Can't write file:%s\n", depfn);
			exit(1);
		}	
		for(j = 0; j < Len; j++)
			if(cov[i][j] >= 0)
				fprintf(depf, "%d\n", cov[i][j]);
		if(fclose(depf) != 0)
		{
			printf("Error closing file:%s\n", depfn);
			exit(1);
		}
	}
	if(fclose(xregf) != 0)
	{
		printf("Error closing file: %s\n", xregfn);
		exit(1);
	}
	if(fclose(outf) != 0)
	{
		printf("Error closing file: %s\n", outn);
		exit(1);
	}
	return 0;
}

int isEmpty(int cov[SN][EL], int sn, int coln)
{
	int count = 0, secount = 0, k;
	float mean = 0, var = 0;
	for(k = 0; k < sn; k++)
	{
		mean += cov[k][coln];
		if(cov[k][coln] <= DT)
			count++;
		else if(cov[k][coln] <= 2 * DT)
			secount++;
	}
	mean /= sn;
	for(k = 0; k < sn; k++)
		var += (cov[k][coln] - mean)*(cov[k][coln] - mean);
	if(var > sn * mean * mean)
		return(1);
	if(count > sn - 5)
		return(1);
	else
	{
		if(count >= (sn * 2 / 3) && count + secount == sn)
			return(1);
		else 
			return(0);
	}
}
