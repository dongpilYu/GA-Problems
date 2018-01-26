#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <sys/time.h>

int power(int num)
{
	int result=1;

	for(int i=0;i<num;i++)
		result*=2;

	return result;
}
double landscape(int NK_N,int NK_K,char* binary)
{
	double** land=(double**)malloc(sizeof(double*)*NK_N);
	int** rand=(int**)malloc(sizeof(int*)*NK_N);
	char input[30];
	sprintf(input,"nk%d_%d",NK_N,NK_K);
	FILE* fpp=fopen(input,"r");
	
	if(fpp==NULL)
		printf("strerror(errno) : %s\n",strerror(errno));	

	for(int i=0;i<NK_N;i++)
	{
		land[i]=(double*)malloc(sizeof(double)*power(NK_K)*2);
		rand[i]=(int*)malloc(sizeof(int)*NK_K);
	}

	for(int i=0;i<NK_N;i++)
		for(int j=0;j<power(NK_K)*2;j++)
			fscanf(fpp,"%lf",&land[i][j]);
		
	for(int i=0;i<NK_N;i++)
		for(int j=0;j<NK_K;j++)
			fscanf(fpp,"%d",&rand[i][j]);
		
	double sum=0;
	int num=0;

	
	for(int i=0;i<NK_N;i++)
	{
		int bin=0;
		int k=NK_K;
		int i_binary=binary[i]-'0';
		
		bin += i_binary * power(k);
		for(int j=0;j<NK_K;j++)
		{
			k--;
			num=rand[i][j];
			i_binary=binary[num]-'0';
			bin += i_binary * power(k);
			
		}
		
		sum += land[i][bin];
	}
	fclose(fpp);

	for(int i=0;i<NK_N;i++)
	{
		free(land[i]);
		free(rand[i]);
	}
	
	free(land);
	free(rand);

	return sum;
}
int royalRoad(int genes,char* binary)
{
	int sum=0;
	int correct=0;
	int num=genes/2;

	for(int i=0;i<num;i++)
	{
		for(int j=0;j<2;j++)
		{
			int num=2*i+j;

			if(binary[num]=='0')
			{
				correct=1;
				break;
			}
		}
		if(correct==0)
			sum += 2;

		correct=0;
	}

	return sum;
	
}
int main(int argc,char** argv)
{
	int genes;
	int NK_K;
	char* type=(char*)malloc(sizeof(char)*20);
	char* binary=(char*)malloc(sizeof(char)*(genes+1));
	char c_genes[10];
	char c_k[10];
	FILE* fp;

	struct timeval t;
	gettimeofday(&t,NULL);
	srand(t.tv_usec*t.tv_sec*getpid());

	genes=atoi(argv[1]);
	type=argv[2];
	if(strcmp(type,"nk")==0)
	{	
		NK_K=atoi(argv[3]);
		sprintf(c_k,"_%d",NK_K);
	}
	sprintf(c_genes,"_%d",genes);
		
	if(strcmp(type,"royal")==0)
	{
		strcat(type,c_genes);
		fp = fopen(type,"a+");

		for(int i=0;i<power(genes);i++)
		{
			int ge=i;
			int divider = genes-1;
			int fitness = 0;

			for(int j=0;j<genes;j++)
			{
				if(ge / power(divider) == 1)
				{
					binary[j]='1';
					ge -= power(divider);	
				}
				else
				{
					binary[j]='0';
				}
				divider--;
				
			}
			binary[genes]='\0';
			fitness=royalRoad(genes,binary);
			fprintf(fp,"%s %d\n",binary,fitness);			
		}

	}

	else if(strcmp(type,"nk")==0)
	{
		strcat(type,c_genes);
		strcat(type,c_k);
		fp = fopen(type,"a+");	
		printf("123");	
		for(int i=0;i<power(genes);i++)
		{
			int ge=i;
		
			int divider = genes-1;
			double fitness = 0;
		
			for(int j=0;j<genes;j++)
			{
				
				if(ge / power(divider) == 1)
				{
					binary[j]='1';
					ge -= power(divider);	
				}
				else
					binary[j]='0';
				
				divider--;
			}
			binary[genes]='\0';
			fitness=landscape(genes,NK_K,binary);
			fprintf(fp,"%s %lf\n",binary,fitness);			
		}
		
	}

	else if(strcmp(type,"onemax")==0)
	{
		strcat(type,c_genes);
		fp = fopen(type,"a+");
	
		for(int i=0;i<power(genes);i++)
		{
			int ge=i;
			int divider = genes-1;
			int fitness = 0;

			for(int j=0;j<genes;j++)
			{
				if(ge / power(divider) == 1)
				{
					binary[j]='1';
					ge -= power(divider);	
					fitness += 1;
				}
				else
				{
					binary[j]='0';
				}
				divider--;
				
			}
			binary[genes]='\0';
			fprintf(fp,"%s %d\n",binary,fitness);			
		}
	}
	else if(strcmp(type,"random")==0)
	{
		strcat(type,c_genes);
		fp = fopen(type,"a+");

		for(int i=0;i<power(genes);i++)
		{
			int ge=i;
			int divider = genes-1;
			int fitness = 0;

			for(int j=0;j<genes;j++)
			{
				if(ge / power(divider) == 1)
				{
					binary[j]='1';
					ge -= power(divider);	
				}
				else
				{
					binary[j]='0';
				}
				divider--;
				
			}
			binary[genes]='\0';
			fitness = rand() % genes;
			fprintf(fp,"%s %d\n",binary,fitness);			
		}
	
	}
	else if(strcmp(type,"deception")==0)
	{
		strcat(type,c_genes);
		fp = fopen(type,"a+");
	
		for(int i=0;i<power(genes);i++)
		{
			int ge=i;
			int divider = genes-1;
			int fitness = 0;

			for(int j=0;j<genes;j++)
			{
				if(ge / power(divider) == 1)
				{
					binary[j]='1';
					ge -= power(divider);	
					fitness += 1;
				}
				else
				{
					binary[j]='0';
				}
				divider--;
				
			}

			if(i == 0)
				fitness = genes;
			else
				fitness -= 1;


			binary[genes]='\0';
			fprintf(fp,"%s %d\n",binary,fitness);			
		}
	
	}
		
	fclose(fp);
}

