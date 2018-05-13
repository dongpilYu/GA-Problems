# One Max-Royal Road-NK Landscape

I implement One-Max, Royal Road, NK-Landscape algorithm and experimentally compute problem epistasis of problem respectively. Problem epistasis is measured by a framework which is based on Shannon's information theory 

## One-Max 

One-Max problem is a simple problem consisting in maximizing the number of ones of a bitstring. 

## Royal Road 

Royal Road function is a function propsed by Forrest and Mitchell to investigate precisely and quantiatively how schema processing actually takes place during the typical evolution of a genetic algorithm.

## NK-Landscape 

NK-Landscape model is a model proposed by Kauffman to define a family of fitness functions that have various dimensions of search space and degrees of epistasis. The functions are tuned by two parameters: N and K. The parameters N and K determines the dimensions of the problem space and the degree of epistasis between genes constituting a chromosome, respectively. 

## Epistasis 
Epistasis means the interaction between genes. It is observed in most GA-hard problems.  

## How to use code

```markdown
gcc -o GA_Problems GA_Problems.c
./ GA_Problems RANDOM temp 1 4 
# RANDOM, NEXTDOOR - NK Landscape
# temp - File name
# 1 - NK_K
# 4 - NK_N
gcc -o Epistasis Epistasis.c
./Epistasis 4 nk 1
# 4 - NK_K
# nk - Type of Genetic Algorithm
# 1 - NK_N
