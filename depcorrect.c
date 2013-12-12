#include<stdlib.h>
#include<stdio.h>
#include<string.h>

#define L 250000

int sn, coln, wlen, slen, chrN;
char chrl[24][16];

int main(int argc, char *argv[])
{
	if(argc != 4)
	{
		printf("Usage: %s sample_number win_len shift_len\n", argv[0]);
		exit(1);
	}
  float *cov, *mdp, *loed;
	sn = atoi(argv[1]);
	wlen = atoi(argv[2]);
	slen = atoi(argv[3]);
	cov = (float *)malloc(sn * L * sizeof(float));
	loed = (float *)malloc(sn * L * sizeof(float));
	mdp = (float *)malloc(sn * sizeof(float));
	if(cov == NULL || mdp == NULL || loed == NULL)
	{
		printf("Memory allocation failed. Plz check\n");
		exit(1);
	}
	char gender[16], bamp[512];
	FILE *clf, *mdf;
	if((clf = fopen("bed/chr.list","r")) == NULL)
	{
		printf("Can't open bed/chr.list\n");
		exit(1);
	}
	while(fscanf(clf, "%s", chrl[chrN]) == 1)
		chrN++;
	if((mdf = fopen("sample.mdep","r")) == NULL)
	{
		printf("Can't open sample.mdep\n");
		exit(1);
	}
	int i, j = 0, k;
	while(fscanf(mdf, "%f", &mdp[j]) == 1)
		j++;
	if(j != sn)
	{
		printf("sample.mdep lines not equal sample number!\n");
		exit(1);
	}
	if(fclose(clf) != 0 || fclose(mdf) != 0)
	{
		printf("Error closing file\n");
		exit(1);
	}
	for(i = 0; i < chrN; i++) 
	{
		FILE *covf, *depf, *ldpf;
		char covfn[256], depfn[256], ldpfn[256];
		sprintf(covfn, "%s_W%dS%d.cov", chrl[i], wlen, slen);
		sprintf(depfn, "%s_W%dS%d.dep", chrl[i], wlen, slen);
		sprintf(covfn, "%s_W%dS%d.ldep", chrl[i], wlen, slen);
		if((covf = fopen(covfn,"r")) == NULL)
		{
			printf("Can't open file %s\n", covfn);
			exit(1);
		}
		int covl = 0, loel = 0;
		while(fscanf(covf, "%f", &cov[covl]) == 1)
			covl++;
		coln = covl / sn;
		if((ldpf = fopen(ldpfn,"r")) == NULL)
		{
			printf("Can't open file: %s\n", ldpfn);
			exit(1);
		}
		while(fscanf(ldpf, "%f", &loed[loel]) == 1)
			loel++;
		if(loel != covl)
		{
			printf("cov file not the same scale of ldep file\n");
			exit(1);
		}
		if((depf = fopen(depfn,"w")) == NULL)
		{
			printf("Can't write file:%s\n", depfn);
			exit(1);
		}
		int cord;
		for(j = 0; j < sn; j++)
		{
			for(k = 0; k < coln; k++)
			{
				if(loed[j*coln+k] > 0)
					cord = (int)(cov[j*coln+k] * mdp[j] / loed[j*coln+k] + 0.5);
				else
				{
					printf("loess prediction = 0, %s, row:%d, col:%d\n",chrl[i],j+1,k+1);
					exit(1);
				}
				fprintf(depf, "%d\t", cord);
			}
			fprintf(depf,"\n");
		}
		if(fclose(ldpf) != 0 || fclose(depf) != 0)
		{
			printf("Error closing file\n");
			exit(1);
		}
	}
	free(cov);
	free(loed);
	free(mdp);
	return 0;
}
