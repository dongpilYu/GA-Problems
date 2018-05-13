#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>

#define INDIVIDUAL 200 
#define GENERATION 1000
#define PROBABILITY_CROSSOVER 0.7
#define PROBABILITY_MUTATION 0.1
#define PROBABILITY_TORNAMENT 0.7
#define NUM_OF_ELITE 4
#define ROYAL_CONSTANT 1

int GENES=0;
int NK_N=0;
int NK_K=0;
int RANDOM=0;
int NEXTDOOR=1;
char* TYPE;
char* LOCATION;
/* NK_K는 반드시 0 <= NK_K <= NK_N-1 을 만족해야 한다. */

/* 유전자의 개수를 바꿀 때 조심해야 합니다. 
royalRoad 함수에서 8로 설정했기 때문에 8의 배수가 되어야 합니다. */

int **population; 
int **next_population; 
int *randomArr;

typedef struct fitness
{
	double ideal;
	double average;
	int indexOfIdeal;
	double fit[INDIVIDUAL];	
}Fitness;
/* 세대 별 최적, 평균, 모든 fit을 저장하는 구조체 */ 

typedef struct nk
{
	int epi;
	int **Rand;
	double **land;

}NK;

int power(int num);
void nk_create(NK* landscape);
void nk_free(NK* landscape);
double nk_fitness(int index,NK* landscape);
void initialize();
void crossover(int* chromosome1,int* chromosome2,int index1,int index2);
void mutation(int* chromosome);
int royalRoad(int index);
int sumOfGene(int index);
int tornament(Fitness* generation);
void elitism(int** eliteGenes);
void fitnessCheck(Fitness* generation);
void fitnessCheckForNK(Fitness* generation,NK* landscape);
void resultToFile(Fitness gener[GENERATION],int num);

int power(int num)
{
	int result=1;

	for(int i=0;i<num;i++)
		result*=2;

	return result;
}

void nk_create(NK* landscape)
{
	int i,j;
	int pow=power(NK_K);
	int exist=0;

	if(strcmp(TYPE,"RANDOM")==0)
		landscape->epi=RANDOM;

	else if(strcmp(TYPE,"NEXTDOOR")==0)
		landscape->epi=NEXTDOOR;


	landscape->Rand=(int**)malloc(sizeof(int*)*NK_N);
	landscape->land=(double**)malloc(sizeof(double*)*NK_N);

	for(i=0;i<NK_N;i++)
	{
		landscape->land[i]=(double*)malloc(sizeof(double)*pow*2);
		landscape->Rand[i]=(int*)malloc(sizeof(int)*NK_K);
	}
	for(i=0;i<NK_N;i++)
	{
		for(j=0;j<pow*2;j++)
			landscape->land[i][j] = ( (double)rand() / RAND_MAX);
		
		for(j=0;j<NK_K;j++)
		{
			landscape->Rand[i][j] = rand() % GENES;
			if(landscape->Rand[i][j]==i)
			{
				j--;
				continue;
			}

			for(int idx=0;idx<j;idx++)
				if(landscape->Rand[i][j] == landscape->Rand[i][idx])
					exist=1;
			if(exist)
			{
				j--;
				exist=0;
				continue;
			}
		}
	}
	
}

double nk_fitness(int index, NK* landscape)
{
	
	int i,j;
	double sum=0;
	int num=0;
	int binary=0;

	if(landscape->epi == RANDOM)
	{
		for(i=0;i<NK_N;i++)
		{
			binary=0;
			int k=NK_K;

			binary += population[index][i] * power(k);
			for(j=0;j<NK_K;j++)
			{
				k--;
				num=landscape->Rand[i][j];
				binary += (population[index][num] * power(k));	
			}
			sum += landscape->land[i][binary];
		}
			
		return sum;
	}

	else if(landscape->epi == NEXTDOOR)
	{			
		for(i=0;i<NK_N;i++)
		{
			binary=0;
			int k=NK_K;

			for(j=0;j<NK_K+1;j++)
			{
				if( i+j >= GENES)
					num = (i+j) % GENES;
				else
					num = i+j;

				binary += (population[index][num] * power(k));
				k--;
			}
			sum += (landscape->land[i][binary]);	
		}
		return sum;	
	}
}

void nk_free(NK* landscape)
{
	int i;
	for(i=0;i<NK_N;i++)
	{
		free(landscape->land[i]);
		free(landscape->Rand[i]);
	}
	free(landscape);
}

void initialize()
{
    int i, j;
    for (i = 0; i < INDIVIDUAL; i++)
	{
        for (j = 0; j < GENES; j++)
		{
            population[i][j] = rand() % 2;
        }
    }
}
 
void crossover(int* chromosome1, int* chromosome2, int index1, int index2)
{
    int j;
    int gene = rand() % GENES;
   
    for (j = 0; j < gene; j++)
	{
        chromosome1[j] = population[index1][j];
    }
   
    for (j = gene; j <GENES; j++)
	{
        chromosome1[j] = population[index2][j];
    }
   
    for (j = 0; j < gene; j++)
	{
        chromosome2[j] =population[index2][j];
    }
   
    for (j = gene; j < GENES; j++)
	{
        chromosome2[j] = population[index1][j];
    }
}

void mutation(int* chromosome)
{
    int gene = rand() % GENES;
   
	if (chromosome[gene] == 0)
        chromosome[gene] = 1;
     
	else 
        chromosome[gene] = 0;    
}

int royalRoad(int index)
{
	int i,j;
	int sum=0;
	int correct=0;
	int num=GENES/ROYAL_CONSTANT;

	for(i=0;i<num;i++)
	{	
		for(j=0;j<ROYAL_CONSTANT;j++)
		{
			int num=ROYAL_CONSTANT*i+j;

			if(population[index][num]==0)
			{
				correct=1;
				break;
			}
		}
	
		if(correct==0)
			sum+=ROYAL_CONSTANT;

		correct=0;
	}
	return sum;
}
 
int sumOfGene(int index)
{
    int j;
    int sum = 0;
    for (j = 0; j < GENES; j++)
        sum += population[index][j];
    
    return sum;
}

int randomFitness(int index)
{
	int seed=0;

	for(int i=0;i<GENES;i++)
		seed = power(i) * population[index][i];
	
	
	return randomArr[seed];
}

int tornament(Fitness* generation)
{
	int rand1,rand2;
 	double prob;

	rand1=rand()%INDIVIDUAL;
	rand2=rand()%INDIVIDUAL;
	prob=((double)rand()/(RAND_MAX));
	
	if(prob < PROBABILITY_TORNAMENT)
	{
		if(generation->fit[rand1]>generation->fit[rand2])
			return rand1;
		else
			return rand2;
	}
	else
	{
		if(generation->fit[rand1]>generation->fit[rand2])
			return rand2;
		else
			return rand1;
	}
}
/* 랜덤으로 두개를 골라서 두 개중의 더 적합한 것의 인덱스를 반환 */

void elitism(int** eliteGenes)
{
	int i,j;
	for(i=0;i<NUM_OF_ELITE;i++)
	{
		for(j=0;j<GENES;j++)
			population[i][j]=eliteGenes[i][j];
	}	
}

void fitnessCheck(Fitness* generation)
{
	int ideal=0;
	int indexOfIdeal=0;
	int sum=0;	
	int i=0;
	
	for(i=0;i<INDIVIDUAL;i++)
	{
		//generation->fit[i]=royalRoad(i);
		//generation->fit[i]=sumOfGene(i);
		generation->fit[i]=randomFitness(i);
		sum+=generation->fit[i];

		if(generation->fit[i]>ideal)
		{
			ideal=generation->fit[i];
			indexOfIdeal=i;	
		}	
	}
	generation->ideal=ideal;
	generation->average=(double)sum/INDIVIDUAL;
	generation->indexOfIdeal=indexOfIdeal;	 
}


void fitnessCheckForNK(Fitness* generation,NK* landscape)
{
	double ideal=0;
	int indexOfIdeal=0;
	double sum=0;
	int i=0;

	for(i=0;i<INDIVIDUAL;i++)
	{
		generation->fit[i]=nk_fitness(i,landscape);
		sum+=generation->fit[i];
	
		if(generation->fit[i]>ideal)
		{
			ideal=generation->fit[i];
			indexOfIdeal=i;
		}
	}
	generation->ideal=ideal;
	generation->average=(double)sum/INDIVIDUAL;
	generation->indexOfIdeal=indexOfIdeal;
}

void resultToFile(Fitness gener[GENERATION],int num)
{
    FILE *fp;
    int index;
        	
    fp = fopen(LOCATION, "a+");

	//if(num != GENERATION)
	fprintf(fp,"%lf\n",gener[num-1].ideal);
	/*
	fprintf(fp,"#fitness1.dat\n");
	fprintf(fp,"#First data block (index ideal)\n#X Y\n");
	
	for(index=0;index<num;index++)
		fprintf(fp," %d %lf\n",index+1,gener[index].ideal);	

	fprintf(fp,"\n\n#Second data block (index average)\n#X Y\n");
	for(index=0;index<num;index++)
		fprintf(fp," %d %lf\n",index+1,gener[index].average);
	*/
	fclose(fp);
}
/* 이 때의 num은 세대를 다 거치기 전 끝난 경우를 위한 변수 */

int main(int argc,char** argv) 
{
	struct timeval t;
	gettimeofday(&t,NULL);
	srand(t.tv_usec * t.tv_sec * getpid());

	TYPE=(char*)malloc(sizeof(char)*30);
	TYPE=argv[1];	
	LOCATION=(char*)malloc(sizeof(char)*30);
	LOCATION=argv[2];
	NK_K=atoi(argv[3]);
	GENES=atoi(argv[4]);
	NK_N=GENES;

   	int i, j, n;	
	population=(int**)malloc(sizeof(int*)*INDIVIDUAL);
	next_population=(int**)malloc(sizeof(int*)*INDIVIDUAL);
	randomArr=(int*)malloc(sizeof(int)*power(GENES));
	int randOpt = rand() % GENES + 1;

	for(i=0;i<INDIVIDUAL;i++)
	{
		population[i]=(int*)malloc(sizeof(int)*GENES);
		next_population[i]=(int*)malloc(sizeof(int)*GENES);
	}	

   	initialize(); 
 
   	Fitness* generation=(Fitness*)malloc(sizeof(Fitness)*GENERATION);

	NK* landscape=(NK*)malloc(sizeof(NK));
	nk_create(landscape);
	
    for (i = 0; i < GENERATION; i++)
    {
    	// fitnessCheck(&generation[i]);
		
    	/* 적합도 계산을 위한 함수 호출 */
		 fitnessCheckForNK(&generation[i],landscape);
		/* NK_LANDSCAPE를 이용한 적합도 계산 함수 */
	
			
		if(i==0)
		{
			FILE* fpp;
			char resultFile[200];
			sprintf(resultFile,"land%d_%d",NK_N,NK_K);
			fpp = fopen(resultFile,"a+");
			for(int a=0;a<NK_N;a++)
			{
				for(int b=0;b<power(NK_K)*2;b++)	
					fprintf(fpp,"%lf ",landscape->land[a][b]);
				fprintf(fpp,"\n");
			}

			for(int a=0;a<NK_N;a++)
			{
				for(int b=0;b<NK_K;b++)
					fprintf(fpp,"%d ",landscape->Rand[a][b]);
				fprintf(fpp,"\n");
			}

			for(int a=0;a<INDIVIDUAL;a++)
			{
			
				if(a%10==0)
					fprintf(fpp,"\n");
				fprintf(fpp,"%lf ",generation[0].fit[a]);
			}
		}
		
		int ideal_individuo = generation[i].indexOfIdeal;
		double ideal_num = generation[i].ideal;     

		int k, l;     
			
		/*
		printf("%d Generation: \n", i + 1);
		k=ideal_individuo; 
			 
		for (l = 0; l < GENES; l++)
			printf("%d ", population[k][l]);    
	
		printf("\n");
	   	printf("%lf", ideal_num);
	   	printf("\n");
		*/
		/* 출력하는 부분 */

		
		if (ideal_num == GENES)    
		    break;	
        
		int indices_selec[INDIVIDUAL];

        for(j = 0; j < INDIVIDUAL; j++)
            indices_selec[j] = tornament(&generation[i]);
        
		/* tornament 방식으로 selection */

		for(j = 0; j < INDIVIDUAL; j++)
			for(int idx=0;idx<GENES;idx++)	
				next_population[j][idx]=population[j][idx];

        for(j = 0; j < INDIVIDUAL; j+=2)
		{
            double prob = ((double) rand() / (RAND_MAX));
	
            if(prob > PROBABILITY_CROSSOVER)
                continue;
               
            int* chromosome1=(int*)malloc(sizeof(int)*GENES);
            int* chromosome2=(int*)malloc(sizeof(int)*GENES);
           
            crossover(chromosome1,chromosome2,indices_selec[j],indices_selec[j+1]);
            for(n = 0; n < GENES; n++)
	   		{
               next_population[j][n] = chromosome1[n];
               next_population[j+1][n] = chromosome2[n];
            }
        }
		/* crossover 수행  */
       
		for(j = 0; j < INDIVIDUAL; j+=2)
		{
			double prob = ((double) rand() / (RAND_MAX));
       	    if(prob > PROBABILITY_MUTATION)
				continue;
           
          int* res=(int*)malloc(sizeof(int)*GENES);
            for(n = 0; n < GENES; n++)
                res[n]=next_population[j][n];
           
            mutation(res);
           
            for(n = 0; n <GENES; n++)
				next_population[j][n]=res[n];
		}
		/* mutation 수행 */
		
		int tmp[INDIVIDUAL];
		int** eliteGenes=(int**)malloc(sizeof(int*)*NUM_OF_ELITE);		

		for(int m=0;m<NUM_OF_ELITE;m++)
			eliteGenes[m]=(int*)malloc(sizeof(int)*GENES);

		int elite[NUM_OF_ELITE];
		int max=0;	
		int index=0;

		for(int m=0;m<INDIVIDUAL;m++)
			tmp[m]=generation[i].fit[m];
	
		for(int idx=0;idx<NUM_OF_ELITE;idx++)
		{
			for(int m=0;m<INDIVIDUAL;m++)
			{
				if(tmp[m]>max)
				{
					max=tmp[m];
					index=m;
				}
			}	
			tmp[index]=0;

			for(int w=0;w<GENES;w++)
				eliteGenes[idx][w]=population[index][w];
		}
		
	    for(j = 0; j < INDIVIDUAL; j++)
	    {	    
			for(n = 0; n < GENES; n++)
			     population[j][n] = next_population[j][n];
	    }
		/* replacement 수행 */

	    elitism(eliteGenes);
		/* elitism 수행 */		
	
	} 
	if(i==GENERATION)		
		resultToFile(generation,i);
	else 
		resultToFile(generation,i+1);
	/* 파일에 결과를 출력 */
	nk_free(landscape);
	


}

