#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <assert.h>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <random>   
#include <chrono>   
#include <cfloat>
#include <limits>
#include "geneticAlgorithm.h"
#include "database.h"

using namespace std;


#define dbg 0
#define time_program 0
#define time_ga 0
#define saveFile 0


/*
 * Sanity check for genetic_algorithm().
 */
#define SANITY_CHECK()         \
	assert(popsize > 0);       \
	assert(ngen > 0);          \
	assert(g->m_rate >= 0.0);  \
	assert(g->m_rate <= 1.0);  \
	assert(g->c_rate >= 0.0);  \
	assert(g->c_rate <= 1.0);  \
	assert(g->e_rate > 0.0);   \
	assert(g->e_rate < 1.0);   \


bool ordenationDistance(const struct Individual *one, 
						const struct Individual *two)
{return (one->distance < two->distance);}

bool ordenationFitness(const struct Individual *one, 
					   const struct Individual *two)
{return (one->fitness < two->fitness);}

bool ordenationDensity(const struct Individual *one, 
					   const struct Individual *two)
{return (one->density < two->density);}

bool ordenationWeighted_sum(const struct Individual *one, 
							const struct Individual *two)
{return (one->weighted_sum < two->weighted_sum);}



/*---------------------------------------------------------------------------*/
/**
 * @brief  Creates a new individual.
 *
 * @returns   A new individual allocates in primary memory.
 */
/*---------------------------------------------------------------------------*/
Individual::Individual()
{ 
	aminoacids_chain = new int**[nproteins];

	for(int i = 0; i < nproteins; i++)
		aminoacids_chain[i] = new int*[database.maxaminoacids]
	
	objectives.resize(numberOfObjectives);
	
	density = 0;
	raw_fitness = 0;
	fitness = 0;
	weighted_sum = 0;	
	distance = 0;	
	matching = 0;
}

 
/*---------------------------------------------------------------------------*/
/**<
 * @brief  Destroys an individual deallocating the memory previously allocated.
 *
 * @param individual The individual that will be destroy.
 */
/*---------------------------------------------------------------------------*/
void individual_destroy(Individual *individual)
{
	for(int i = 0; i < nproteins;i++ )
		free(individual->aminoacids_chain[i]);

	free(individual->aminoacids_chain);
	
	individual->objectives.resize(0);
	delete individual;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief  Creates a random individual.
 *
 * @returns   
 */
/*---------------------------------------------------------------------------*/
Individual *individual_random()
{	
	int i = 0, k = 0;
	vector <int> aminoacids_list;
	
	#if(dbg>0)
	fprintf(stderr, "info: Individual random\n");
	#endif
	
	Individual *g = new Individual();
	
    /* create Individual.*/
    for (i = 0; i < nproteins; i++)
    {	
		for( k = 0; k < database.naminoacids[i]; k++)
			aminoacids_list.push_back(k);	
			
		for(int j = k; j < database.maxaminoacids;j++)
			aminoacids_list.push_back(-1);
		
		 std::random_shuffle(aminoacids_list.begin(), 
				 aminoacids_list.end());
			
		for( k = 0; k < database.maxaminoacids;k++)
		{	
			g->aminoacids_chain[i][k] = aminoacids_list.back();
			aminoacids_list.pop_back();
		}
		aminoacids_list.clear();		
	} 

	return g;
}

/* --------------------------------------------------------------------------*/
/**
 * @brief  Asserts if a individual has a feature.
 *
 * @param g 
 * @param feature
 * @param n
 *
 * @returns 1 if the individual has the feature, and 0 otherwise.  
 */
/*----------------------------------------------------------------------------*/
int has_feature(int *g, int feature, int n)
{	
	int resp = 0;

	for (int i = 0; i < n; i++)
	{
		if (g[i] == feature)
		{
			resp = 1;
			i = n;
		}
	}
		return resp;
}


/* --------------------------------------------------------------------------*/
/**
 * @brief  
 *
 * @param individual
 *
 * @returns   
 */
/*---------------------------------------------------------------------------*/
Individual *copy_individual(Individual *individual)
{
	Individual *g = new Individual();

	for(int wprotein = 0; wprotein < nproteins;wprotein++)
		memcpy(g->aminoacids_chain[wprotein],
				individual->aminoacids_chain[wprotein], 
				database.maxaminoacids*sizeof(int));
	
	for(int i = 0; i < numberOfObjectives;i++)
		g->objectives[i] = individual->objectives[i];
		
	g->density = individual->density;
	g->raw_fitness = individual->raw_fitness;
	g->fitness = individual->fitness;
	g->weighted_sum = individual->weighted_sum;	
	g->distance = individual->distance;
	g->matching = individual->matching;
		
	return g;
}


/* --------------------------------------------------------------------------*/
/**
 * @brief  
 *
 * @param offspring
 * @param individual1
 * @param individual2
 */
/*---------------------------------------------------------------------------*/
void aminoacids_crossover(Individual *offspring, 
						  Individual *individual1, 
						  Individual *individual2)
{	
	int point1 = 0;    		 		/* First crossover point.  */
	int point2 = 0;     	 		/* Second crossover point. */
	int nbegin, nmiddle, nend; 	 	/* Size of gene parts.     */
	int i, j; 						/* Loop index.			   */
	int *begin, *middle, *end;		/* Gene parts.             */
	int **map = new int**[2];
	bool search;

	/* Sanity check. */
	assert(n >= 0);
		
	/* Generate crossover points. */
	int tmp;
	do
	{
		point1 = rand()%database.maxaminoacids;
		point2 = rand()%database.maxaminoacids;
		
	} while (((point1<=1) || (abs(point1-point2)<=1) || 
			  (point2<=1)|| ((database.maxaminoacids-point2)<=1)));

	if (point1 > point2)
	{
		tmp = point1;
		point1 = point2;
		point2 = tmp;
	}
	/* Size of gene parts. */
	nbegin = point1 - 1 ;
	nmiddle = (point2 - 1) - nbegin;
	nend = database.maxaminoacids - (nbegin+nmiddle);

	for(int wprotein = 0; wprotein < nproteins; wprotein++)
	{		
		/* Gene parts. */
		begin  = new int *[nbegin];
		middle = new int *[nmiddle];
		end    = new int *[nend];
		
		map[0] = new int *[nmiddle];
		map[1] = new int *[nmiddle];
		
		memcpy(map[0],individual1->aminoacids_chain[wprotein]+ nbegin,
			   	nmiddle*sizeof(int));
		memcpy(map[1],individual2->aminoacids_chain[wprotein] + nbegin,
			   	nmiddle*sizeof(int));
	
		
		memcpy(begin, individual1->aminoacids_chain[wprotein], 
				nbegin *sizeof(int));
		memcpy(middle,individual2->aminoacids_chain[wprotein] + nbegin,
			   	nmiddle*sizeof(int));
		memcpy(end, 
				individual1->aminoacids_chain[wprotein]+nbegin+nmiddle,
				nend*sizeof(int));
		
		do
		{
			search = false;
			
			for ( i = 0; i < nmiddle; i++)
			{
				for ( j = 0; j < nbegin; j++)
				{
					if (map[1][i] == begin[j] and begin[j] != -1)
					{					
						begin[j] = map[0][i];
						search = true;
						break;				
					}
				}			
			
				for ( j = 0; j < nend; j++)
				{
					if (map[1][i] == end[j] and  end[j] != -1)
					{		
						end[j] = map[0][i];				
						search = true;
						break;
					}
				}
			}				
		}while(search);			
		
				
		memcpy(offspring->aminoacids_chain[wprotein], begin, nbegin*sizeof(int)); 
		memcpy(offspring->aminoacids_chain[wprotein] + nbegin, middle, nmiddle*sizeof(int));
		memcpy(offspring->aminoacids_chain[wprotein] + nbegin + nmiddle, end, nend*sizeof(int));

		/* House keeping. */
		delete begin;
		delete middle;
		delete end;
		delete map[0];
		delete map[1];
	}
	delete map;	
}


/* --------------------------------------------------------------------------*/
/**
 * @brief  
 *
 * @param offspring
 * @param mom
 * @param dad
 */
/*---------------------------------------------------------------------------*/
void individual_crossover(Individual *offspring, Individual *mom,
	   					  Individual *dad)
{
	#if(dbg>0)
	fprintf(stderr, "info: Individual crossover...\n");
	#endif
	
	aminoacids_crossover(offspring, mom,  dad);
}


/*---------------------------------------------------------------------------*/
/**
 * @brief  
 *
 * @param individual
 */
/*---------------------------------------------------------------------------*/
void individual_mutation(Individual *individual)
{	
	#if(dbg>0)
	fprintf(stderr, "info: Individual mutation...\n");
	#endif	
	
	int i, j;     		/* Mutation point. */
	int aminoacid;  	/* Feature.        */
	int temp;
	for(int wprotein = 0; wprotein < nproteins; wprotein++)
	{
		i = rand()%database.maxaminoacids;
		aminoacid = rand()%database.naminoacids[wprotein];
		
		if(not has_feature(individual->aminoacids_chain[wprotein], 
					aminoacid, database.maxaminoacids))
			individual->aminoacids_chain[wprotein][i] = aminoacid;
		else
		{
			j = rand()%database.maxaminoacids;
			temp = individual->aminoacids_chain[wprotein][i];
			individual->aminoacids_chain[wprotein][i] = 
							individual->aminoacids_chain[wprotein][j];
			individual->aminoacids_chain[wprotein][j] = temp;
		}
	}
}


/*---------------------------------------------------------------------------*/
/**
 * @brief  
 *
 * @param ind
 */
/*---------------------------------------------------------------------------*/
void objectiveOne(Individual *ind)
{
	#if(dbg>0)
	fprintf(stderr, "info: objective one...\n");
	#endif	
	
	double matches_aminoacids = 0.0000001;
	
	for(int wProtein1 = 0; wProtein1 < nproteins; wProtein1++)
	{
		for(int wProtein2 = 0; wProtein2 < wProtein1; wProtein2++)
		{	
			for(int wAminoacid = 0; wAminoacid < database.maxaminoacids; wAminoacid++)
			{
				int first_aminoacid  = ind->aminoacids_chain[ wProtein1 ][ wAminoacid ];
				int second_aminoacid = ind->aminoacids_chain[ wProtein2 ][ wAminoacid ];
					
				if(first_aminoacid != -1 and second_aminoacid != -1)
				{					
					if (strcmp( (database.data_aminoacids[wProtein1][ first_aminoacid ]->aminoacids_aa).c_str(),
					   (database.data_aminoacids[wProtein2][ second_aminoacid ]->aminoacids_aa).c_str()) == 0) 
						
					{
						matches_aminoacids++;
					}
				}	
			}
		}
	}
			
	ind->objectives[0] = (1/matches_aminoacids);
	ind->matching = matches_aminoacids;
}


/*---------------------------------------------------------------------------*/
/**
 * @brief  
 *
 * @param i1
 */
/*---------------------------------------------------------------------------*/
void objectiveTwo(Individual *i1)
{
	#if(dbg>0)
	fprintf(stderr, "info: objective two...\n");
	#endif	

	double x, y, z, dist = 0;
	
	for(int wAminoacid = 0; wAminoacid < database.maxaminoacids; wAminoacid++)
	{
		
		int first_aminoacid  = i1->aminoacids_chain[ 0 ][ wAminoacid ];
		int second_aminoacid = i1->aminoacids_chain[ 1 ][ wAminoacid ];
		
		if(first_aminoacid != -1 and second_aminoacid != -1)
		{
			x = database.data_aminoacids[0][first_aminoacid]->listAtoms[1]->x - database.data_aminoacids[1][second_aminoacid]->listAtoms[1]->x;
			y = database.data_aminoacids[0][first_aminoacid]->listAtoms[1]->y - database.data_aminoacids[1][second_aminoacid]->listAtoms[1]->y;
			z = database.data_aminoacids[0][first_aminoacid]->listAtoms[1]->z - database.data_aminoacids[1][second_aminoacid]->listAtoms[1]->z;

			dist += pow(x,2)+ pow(y,2) + pow(z,2); /* Calculating distance by euclidean formula. */
		}
	}          
			
	// Initializing objectives
	i1->objectives[1] = (1/(float)
	(1 + (dist/(float)(pow(1.24*pow(((database.maxaminoacids -15) -1.8), 1.0/3.0), 2)))));

}


/*---------------------------------------------------------------------------*/
/**
 * @brief  
 *
 * @param ind
 */
/*---------------------------------------------------------------------------*/
void objectiveThree(Individual *ind)
{
	double gaps = 0;
	
	for(int wProtein1 = 0; wProtein1 < nproteins; wProtein1++)
	{
		for(int wProtein2 = 0; wProtein2 < wProtein1; wProtein2++)
		{	
			for(int wAminoacid = 0; wAminoacid < database.maxaminoacids; wAminoacid++)
			{
				int first_aminoacid  = ind->aminoacids_chain[ wProtein1 ][ wAminoacid ];
				int second_aminoacid = ind->aminoacids_chain[ wProtein2 ][ wAminoacid ];
				
				if(first_aminoacid == -1 or second_aminoacid == -1)
				{
					gaps++;
				}
			}
		}
	}
	
	ind->objective[2] = 1/gaps;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief  
 *
 * @param popsize
 * @param ngen
 * @param crossoverRate
 * @param mutationRate
 */
/*---------------------------------------------------------------------------*/
void protein_matches(int popsize, int ngen, double crossoverRate, double mutationRate)
{
	/*
	*  Configuration genome.
	*/
	struct genome problem = 
	{
			mutationRate,       /* Mutation rate.    */
			crossoverRate,      /* Crossover rate.   */
			1.00            	/* Replacement rate. */
	};
	
	
	geneticAlgorithm(&problem, popsize, ngen);
} 

/*============================================================================*
 *                              genetic utilities                             *
 *============================================================================*/

/*---------------------------------------------------------------------------*/
/**
 * @brief  
 *
 * @param individual
 */
/*---------------------------------------------------------------------------*/
void print_individual(Individual *individual)
{
	#if(dbg>0)
	fprintf(stderr, "info: print individual\n");
	#endif
	
	fprintf(stderr,"\n");
		
		for(int j = 0; j < nproteins; j++)
		{		
			for(int k = 0; k < database.maxaminoacids;k++)
			{
				//fprintf(stderr,"%d ",individual->aminoacids_chain[j][k]);
				if(individual->aminoacids_chain[j][k] != -1)
					fprintf(stderr, "%s ",((database.data_aminoacids[j][individual->aminoacids_chain[j][k]])->aminoacids_aa).c_str());
				else
					fprintf(stderr, "- ");
			}	
			fprintf(stderr,"\n");
		}
		
		fprintf(stderr,"\n");
		for(int i = 0; i < numberOfObjectives; i++)
			fprintf(stderr,"Objective %d: %f \n",i,individual->objectives[i]);
		
		fprintf(stderr,"Density: %f \n",individual->density);
		fprintf(stderr,"Raw_fitness: %f \n",individual->raw_fitness);
		fprintf(stderr,"Fitness: %f \n",individual->fitness);
		fprintf(stderr,"Weight_sum: %f \n",individual->weighted_sum);
		fprintf(stderr,"Distance: %f \n",individual->distance);	
		fprintf(stderr, "Matching: %d\n ", individual->matching);
		fprintf(stderr,"\n");
	
}


/*---------------------------------------------------------------------------*/
/**
 * @brief  
 *
 * @param g
 * @param population
 * @param selected
 */
/*---------------------------------------------------------------------------*/
void reproduce(genome *g, vector<Individual *> &population, vector<Individual *> &selected)
{
	#if(dbg>0)
	fprintf(stderr, "info: reproduce\n");
	#endif
	
	/* Sanity check */
	assert(!selected.empty());
	
	for(Individual *i: population)		
		individual_destroy(i);
		
	population.resize(0);
		
	/* Make new population to construct offsprings of P. */
	while (population.size() < popsize) 
	{
		if (((double) rand() / RAND_MAX) < g->c_rate)
		{
			Individual *childSolution = new Individual();
			int r1 = rand()%popsize;
			int r2 = rand()%popsize;
			
			individual_crossover(childSolution, selected[r1], selected[r2]);
			
			if (((double) rand() / RAND_MAX) < g->m_rate) 
				individual_mutation(childSolution);
		
			population.push_back(childSolution);
		}
	}
	
	for(Individual *i: selected)		
		individual_destroy(i);
		
	selected.resize(0);
}

/*---------------------------------------------------------------------------*/
/**
 * @brief  
 *
 * @param population
 */
/*---------------------------------------------------------------------------*/
void calculate_objectives(vector<Individual *> &population)
{
	#if(dbg>0)
	fprintf(stderr, "info: calculate objectives\n");
	#endif
	
	for(int i = 0; i < population.size();i++)
	{
		objectiveOne(population[i]);
		objectiveTwo(population[i]);
	}
}

/*---------------------------------------------------------------------------*/
/**
 * @brief  
 *
 * @param individual
 *
 * @returns   
 */
/*---------------------------------------------------------------------------*/
bool Individual::dominates(Individual *individual)
{
	#if(dbg>0)
	fprintf(stderr, "info: dominates \n");
	#endif
	
	for (int i = 0; i < numberOfObjectives; i++) 
		if (this->objectives[i] > individual->objectives[i])
			return false;
	
	return true;
}


/*---------------------------------------------------------------------------*/
/**
 * @brief  
 *
 * @param Union
 */
/*---------------------------------------------------------------------------*/
void calculate_dominated(vector<Individual *> &Union)
{
	#if(dbg>0)
	fprintf(stderr, "info: calculate dominated \n");
	#endif				
	
	
	for (Individual *p1 : Union)
	{
		if(!p1->dom_set.empty())
			p1->dom_set.clear();

		for (Individual *p2 : Union)
		{			
			if (p1 == p2) 
				continue;

			if (p1->dominates(p2))
				p1->dom_set.push_back(p2);				
		}

	}	
	 
	
}


/*---------------------------------------------------------------------------*/
/**
 * @brief  
 *
 * @param p1
 * @param Union
 *
 * @returns   
 */
/*---------------------------------------------------------------------------*/
int calculate_raw_fitness(Individual *p1, vector<Individual *> &Union)
{
	#if(dbg>0)
	fprintf(stderr, "info: raw_fitness\n");
	#endif
	
	int sum = 0;
	
	for (Individual *p2 : Union)
	{
		if (p1 == p2) 
			continue;

		if (p2->dominates(p1))
			sum = sum + p2->dom_set.size();		
	}
	
	return sum;
}


/*---------------------------------------------------------------------------*/
/**
 * @brief  
 *
 * @param objectivesI1
 * @param objectivesI2
 *
 * @returns   
 */
/*---------------------------------------------------------------------------*/
double euclidian_distance(vector<double> objectivesI1, vector<double> objectivesI2)
{
	#if(dbg>0)
	fprintf(stderr, "info: euclidian distance\n");
	#endif
	
	double sum = 0.0;
	
	for(int i = 0; i < objectivesI1.size(); i++)
		sum += pow((objectivesI1[i] - objectivesI2[i]), 2);
	
	return sqrt(sum);
}


/*---------------------------------------------------------------------------*/
/**
 * @brief  
 *
 * @param p1
 * @param Union
 *
 * @returns   
 */
/*---------------------------------------------------------------------------*/
double calculate_density(Individual *p1, vector<Individual *> &Union)
{
	#if(dbg>0)
	fprintf(stderr, "info: calculate density\n");
	#endif
	
	for(Individual *p2 : Union)
		p2->distance = euclidian_distance(p1->objectives, p2->objectives);
	
	vector<Individual *> list = sortVector(Union, "distance");
	
	int k = sqrt(list.size());
		
	return (1.0/(list[k]->distance + 2.0));
	
}

/*---------------------------------------------------------------------------*/
/**
 * @brief  
 *
 * @param population
 * @param archive
 */
/*---------------------------------------------------------------------------*/
void calculate_fitness(vector<Individual *> &population, vector<Individual *> &archive)
{
	#if(dbg>0)
	fprintf(stderr, "info: calculate fitness\n");
	#endif
	
	calculate_objectives(population);
	
	vector<Individual *> Union;
	Union.reserve(population.size() + archive.size());
	Union.insert(Union.end(), population.begin(), population.end());
	Union.insert(Union.end(), archive.begin(), archive.end());
	
	calculate_dominated(Union);
	
	for (Individual *p : Union)
	{
		p->raw_fitness = calculate_raw_fitness(p, Union);
		p->density = calculate_density(p, Union);
		p->fitness = p->raw_fitness + p->density;
	}
	
	Union.clear();
}


/*---------------------------------------------------------------------------*/
/**
 * @brief  
 *
 * @param population
 * @param archive
 */
/*---------------------------------------------------------------------------*/
void environmetal_selection(vector<Individual *> &population, vector<Individual *> &archive)
{
	#if(dbg>0)
	fprintf(stderr, "info: environmental selection\n");
	#endif
	
	vector<Individual *> Union;
	Union.reserve(population.size() + archive.size());
	Union.insert(Union.end(), population.begin(), population.end());
	Union.insert(Union.end(), archive.begin(), archive.end());
	
	vector<Individual *> enviroment;
	
	for (Individual *p : Union)
		if(p->fitness < 1.0)
			enviroment.push_back(copy_individual(p));
	
	if(enviroment.size() < archiveSize)
	{
		#if(dbg>0)
		fprintf(stderr, "info: archive shorter than environment\n");
		#endif
		/* TEM QUE SER EM ORDEM DECRESCENTE*/
		std::sort(Union.rbegin(), Union.rend(), ordenationFitness);
		
		int count = 0;
		while(enviroment.size() < archiveSize)
		{
			if(Union[count]->fitness >= 1.0)
				enviroment.push_back(copy_individual(Union[count]));
				
			count++;
		}		
	}
	else if(enviroment.size() > archiveSize)
	{
		#if(dbg>0)
		fprintf(stderr, "info: environment greater than archive \n");
		#endif
		do
		{
			int k = sqrt(enviroment.size());
						
			for(Individual *p1: enviroment)
			{
				p1->density = 0;
				
				for(Individual *p2: enviroment)
				{
					if (p1 == p2) 
						continue;
					
					p2->distance = euclidian_distance(p1->objectives, p2->objectives);
				}
				vector<Individual *> list = sortVector(enviroment, "distance");
				
				//for(int i = 0; i < k; i++)
				p1->density = list[k]->distance;			
			}	
			std::sort(enviroment.rbegin(), enviroment.rend(), ordenationDensity);  
			enviroment.erase(enviroment.begin(),enviroment.begin()+1);
			
		}while(enviroment.size() <= archiveSize);
	}
	
	/* House Keeping. */
	for(Individual *i: archive)
		individual_destroy(i);
	
	archive.clear();
	
	/* New archive generate from the environment. */
	archive = enviroment;	
}


/*---------------------------------------------------------------------------*/
/**
 * @brief  
 *
 * @param oldVector
 * @param type
 *
 * @returns   
 */
/*---------------------------------------------------------------------------*/
vector<Individual *> sortVector(vector<Individual *> oldVector, string type)
{
	vector<Individual *> copy;
	
	
	if(type.compare("distance") == 0)
		std::sort(oldVector.begin(), oldVector.end(), ordenationDistance);
	else if(type.compare("fitness") == 0)
		std::sort(oldVector.begin(), oldVector.end(), ordenationFitness);
	else if(type.compare("density") == 0)
		std::sort(oldVector.begin(), oldVector.end(), ordenationDensity);
	
	copy = oldVector;
	
	return copy;
	
}


/*---------------------------------------------------------------------------*/
/**
 * @brief  
 *
 * @param archive
 */
/*---------------------------------------------------------------------------*/
void calculate_weighted_sum(vector<Individual *> &archive)
{
	#if(dbg>0)
	fprintf(stderr, "info: calculate weighted sum\n");
	#endif
	
	for(int i = 0; i < archive.size();i++)
	{
		double sum = 0;
		
		for(int j = 0; j < numberOfObjectives;j++)
			sum +=archive[i]->objectives[j]; 
		
		archive[i]->weighted_sum = sum;
	}
}


/*---------------------------------------------------------------------------*/
/**
 * @brief  
 *
 * @param selected
 * @param population
 * @param tournament
 */
/*---------------------------------------------------------------------------*/
void binary_tournament(vector<Individual *> &selected, 
					   vector< Individual *> &population, unsigned tournament)
{
	#if(dbg>0)
	fprintf(stderr, "info: tournament\n");
	#endif
	
	/* Sanity check */
	assert(!population.empty());
		
	unsigned i = 0, j = 0; /* Loop index*/
	vector<Individual *> ranking;
			
	while(i < popsize)
	{
		j = 0;
		
		while(j < tournament)
		{
			unsigned pos  = rand()%popsize;
			ranking.push_back(population[pos]);		
			j++;
		}			
				
		std::sort(ranking.begin(),ranking.end(),ordenationFitness);
		selected.push_back(copy_individual(ranking[0]));
				
		ranking.resize(0);
									
		i++;
	}
}

/*============================================================================*
 *                              genetic Algoritmh                             *
 *============================================================================*/

/*---------------------------------------------------------------------------*/
/**
 * @brief  
 *
 * @param g
 * @param popsize
 * @param ngen
 */
/*---------------------------------------------------------------------------*/
void geneticAlgorithm(genome *g, int popsize, int ngen)
{

	//srand (time(NULL));
	fprintf(stderr,"popsize: %d\n", popsize);
	fprintf(stderr,"ngen: %d\n", ngen);
	fprintf(stderr,"Crossover rate: %.2f\n", g->c_rate);
	fprintf(stderr,"Repalcement rate: %.2f\n", g->r_rate);
	fprintf(stderr,"Mutation rate: %.2f\n", g->m_rate);
	fprintf(stderr, "Number od objectives: %d\n\n", numberOfObjectives);
	
	
	vector<Individual *> population;/* Current population.       					 */
	vector<Individual *> archive;	/* Archive to storage individuals non-dominated. */
	vector<Individual *> selected;	/* Individuals selected from the population.	 */
	
	SANITY_CHECK();
	
	for(int count = 0; count < popsize; count++)
		population.push_back(individual_random());
	
	for(int count = 0; count < ngen; count++)
	{	
		/* Calculates fitness for each individual. */
		calculate_fitness(population, archive);	
		/* Enviromental selection using non-dominated individual and trucator operator. */
		environmetal_selection(population, archive);
		/*	Calculated weighted sum of objective. */	
		calculate_weighted_sum(archive);
		/* Sort population to show the best individuals from the population. */
		std::sort(archive.begin(), archive.end(), ordenationWeighted_sum);
		print_individual(archive[0]);
		/* Tournament to select the most probable. */
		binary_tournament(selected, population, 2);			
		/* Executes mutation and crossover. */
		reproduce(g, population, selected);			
	}	
	
	/* House Keeping. */
	for(Individual *i: archive)		
		individual_destroy(i);
				
	for(Individual *i: population)		
		individual_destroy(i);	
}
