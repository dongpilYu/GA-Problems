#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int power(int num)
{
	int result = 1;

	for(int i=0;i<num;i++)
		result *= 2;

	return result;
}

double logB(double x, double base)
{
	if(x==0)
		return 0;
	return log(x) / log(base);
}

int nchoosek(int a,int b)
{
	int i,j,n=1,r=1;

	for(i=a;i>=a-b+1;i--)
		n=n*i;
	for(j=b;j>=1;j--)
		r=r*j;

	return n/r;
}


int main(int argc,char**argv)
{
	int Genes=atoi(argv[2]);
	int NK_K;
	char* type=(char*)malloc(sizeof(char)*20);
	type=argv[1];
	char c_genes[20];
	char c_k[20];
	FILE* rFile;
	FILE* wFile;
	char result[40]="result_";
	strcat(result,type);

	char** solTable = (char**)malloc(sizeof(char*)*power(Genes));
	int* fitTable = (int*)malloc(sizeof(int)*power(Genes));
	double* pSumTable = (double*)malloc(sizeof(double)*(Genes+1));

	//memset(fitTable,0,sizeof(int)*power(Genes));
	//memset(pSumTable,0,sizeof(double)*(Genes+1));

	int Hx = 1;
	double Hy = 0;
	double* Hxy = (double*)malloc(sizeof(double)*(Genes));
	double* Ixy = (double*)malloc(sizeof(double)*(Genes));
	memset(Hxy,0,sizeof(double)*(Genes));
	memset(Ixy,0,sizeof(double)*(Genes));
	int Hxx = 2;
	double** Hxxy = (double**)malloc(sizeof(double*)*(Genes-1));
	double** Ixxy = (double**)malloc(sizeof(double*)*(Genes-1));
	double** epsilon = (double**)malloc(sizeof(double*)*(Genes-1));
	double eta = 0;

	double percent = (double)1 / power(Genes);

	for(int i=0;i<Genes-1;i++)
	{
		Hxxy[i]=(double*)malloc(sizeof(double)*(Genes-1));
		Ixxy[i]=(double*)malloc(sizeof(double)*(Genes-1));
		epsilon[i]=(double*)malloc(sizeof(double)*(Genes-1));
	}

	for(int i=0;i<power(Genes);i++)
	{
		solTable[i] = (char*)malloc(sizeof(char)*(Genes));
	//	memset(solTable[i],0,sizeof(solTable[i]));
	}
	
	sprintf(c_genes,"_%d",Genes);
	strcat(result,c_genes);
	
	if(strcmp(type,"nk")==0)
	{
		NK_K = atoi(argv[3]);
		sprintf(c_k,"_%d",NK_K);
		strcat(result,c_k);
	}
	

	if(strcmp(type,"royal")==0 || strcmp(type,"onemax")==0 || strcmp(type,"random")==0 || strcmp(type,"nk")==0 || strcmp(type,"deception")==0 )
	{
		if(strcmp(type,"nk")==0)
		{
			strcat(type,c_genes);
			strcat(type,c_k);
		}
		else
			strcat(type,c_genes);
		
		rFile = fopen(type,"r");
		wFile = fopen(result,"a+");	

		for(int i=0;i<power(Genes);i++)
		{
			double dFit=0;
			fscanf(rFile,"%s %lf",solTable[i],&dFit);
			fitTable[i]=(int)dFit;
			//printf("%d\n",fitTable[i]);
		}
		for(int i=0;i<=Genes;i++)
		{
			for(int j=0;j<power(Genes);j++)
			{
				if(fitTable[j]==i)
					pSumTable[i] += percent;
			}

		}

		for(int i=0;i<=Genes;i++)
		{
			double p = pSumTable[i];
			Hy = Hy - p * logB(p,2.0);
		}
		
		for(int i=0;i<Genes;i++)
		{
			double p;
			Hxy[i] = 0;
			Ixy[i] = 0;

			for(int j=0;j<2;j++)
			{
				for(int k=0;k<=Genes;k++)
				{
					p=0;
		
					if(j==0)
					{
						
						for(int l=0;l<power(Genes);l++)
						{
							if(solTable[l][i] == '0' && fitTable[l] == k)
								p += percent;
						}
						Hxy[i] = Hxy[i] - p * logB(p,2.0);
					}
					else if(j==1)
					{
						for(int l=0;l<power(Genes);l++)
						{
							if(solTable[l][i] == '1' && fitTable[l] == k)
								p += percent;
						}
						Hxy[i] = Hxy[i] - p * logB(p,2.0);
					}
				}
			}
			Ixy[i] = Hx + Hy - Hxy[i];

		}
		
		double HXY=0;
		double IXY=0;
		for(int i=0;i<Genes;i++)
		{
			HXY += Hxy[i];
			IXY += Ixy[i];
		}
		HXY = HXY/Genes;
		IXY = IXY/Genes;
		//fprintf(wFile,"Hxy : %lf\n",HXY);
		//fprintf(wFile,"Ixy : %lf\n",IXY);

		for(int i=0;i<Genes-1;i++)
		{
			for(int j=i;j<Genes-1;j++)
			{
				double p=0;
				for(int l=0;l<=Genes;l++)
				{
					
					for(int k=0;k<power(Genes);k++)
					{
						if(solTable[k][i] == '0' && solTable[k][j+1] == '0' && fitTable[k] == l)
						{
							p += percent;
						}
					}
					
					Hxxy[i][j] = Hxxy[i][j] - p * logB(p,2.0);
					p=0;
					for(int k=0;k<power(Genes);k++)
					{
						if(solTable[k][i] == '0' && solTable[k][j+1] == '1' && fitTable[k] == l)
						{
							p += percent;
						}
					}
					
					Hxxy[i][j] = Hxxy[i][j] - p * logB(p,2.0);
					p=0;
					for(int k=0;k<power(Genes);k++)
					{
						if(solTable[k][i] == '1' && solTable[k][j+1] == '0' && fitTable[k] == l)
						{
							p += percent;
						}
					}
		
					Hxxy[i][j] = Hxxy[i][j] - p * logB(p,2.0);
					p=0;
					for(int k=0;k<power(Genes);k++)
					{
						if(solTable[k][i] == '1' && solTable[k][j+1] == '1' && fitTable[k] == l)
						{
							p += percent;
						}
					}
	
					Hxxy[i][j] = Hxxy[i][j] - p * logB(p,2.0);
					p=0;
				}
				
				Ixxy[i][j] = Hxx + Hy - Hxxy[i][j];
				
			}
		}

		double IXXY=0;
		double HXXY=0;
		for(int i=0;i<Genes-1;i++)
		{
			for(int j=i;j<Genes-1;j++)
			{
				HXXY += Hxxy[i][j];
				IXXY += Ixxy[i][j];
				if(Ixxy[i][j] == 0)
					epsilon[i][j] = 0;
				else
					epsilon[i][j] = 1 - (Ixy[i] + Ixy[j+1]) / Ixxy[i][j];
				eta += fabs(epsilon[i][j]);
			}
		}
		HXXY = HXXY / nchoosek(Genes,2);
		IXXY = IXXY / nchoosek(Genes,2);
		eta = eta / nchoosek(Genes,2);
	
		//fprintf(wFile,"Hxxy : %lf\n",HXXY);
		//fprintf(wFile,"Ixxy : %lf\n",IXXY);
		fprintf(wFile,"%lf\n",eta);

	}

	

	fclose(wFile);
	fclose(rFile);
}

