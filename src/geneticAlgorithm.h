#ifndef _GENETIC_ALGORITHM_H_
#define _GENETIC_ALGORITHM_H_
#include <vector>

	/**
	 *  Gene.
	 */
	struct Individual
	{
		int **aminoacids_chain;    		 /* Chains of amino acids.  */
		double density;
		double raw_fitness;
		double fitness;
		double weighted_sum;	
		double distance;	
		int matching;

		std::vector<double> objectives;
		std::vector<Individual *> dom_set;
		
		Individual();
		bool dominates(Individual *individual);
		
	};

	/*
	 * Genome.
	 */
	struct genome
	{
		/* Attributes. */
		double m_rate; /* Mutation rate.    */
		double c_rate; /* Crossover rate.   */
		double r_rate; /* Replacement rate. */
	};
	
	/**
	 * Creates a Individual.
	 */
	Individual *individual_create();
		
	/**
	 * Destroys an individual.
	 */
	void individual_destroy(Individual *individual);
	
	/**
	 * Mates two individuals.
	 */
	void individual_crossover(Individual *offspring, Individual *mom, Individual *dad);
	
	/**
	 * Mates two amino acids chain.
	 */
	void aminoacids_crossover(Individual *offspring,  Individual *mom, Individual *dad);
	
	/**
	 * Mutates a individual.
	 */
	void individual_mutation(Individual *individual);
	
	/**
	 * Generates a random individual.
	 */
	Individual *individual_random();

	/**
	 *  One if the Individual has the feature, and zero otherwise.
	 */
	int has_feature(int *g, int feature, int n);
	
	/**
	 * Genetic algorithm.
	 */
	void geneticAlgorithm(genome *g, int popsize, int ngen);

	/**
	 * Call genetic algorithm for protein matches
	 */
	void protein_matches(int popsize, int ngen, double crossoverRate, double mutationRate);
	
	/**
	 * Prints the individual and his respective amino acids and atoms.
	 */
	void print_individual(Individual *individual);
	
	/**
	 * Evaluates each individual 
	 */
	void individual_evaluate(std::vector<struct Individual> &population);
			
	/**
	 * Funtion used to make a copy of a Individual.
	 */
	Individual *copy_individual(Individual *individual);
	
	/**
	 * 
	 */
	void binary_tournament(std::vector<Individual *> &parents, std::vector< Individual *> &population, unsigned tournament);

	/**
	 * 
	 */
	void calculate_weighted_sum(std::vector<Individual *> &archive);
	
	/**
	 * 
	 */
	std::vector<Individual *> sortVector(std::vector<Individual *> oldVector, std::string type);
		
	/**
	 *
	 */
	void environmetal_selection(std::vector<Individual *> &population, std::vector<Individual *> &archive);
	
	/**
	 * 
	 */
	void calculate_fitness(std::vector<Individual *> &population, std::vector<Individual *> &archive);
	
	/**
	 * 
	 */
	double calculate_density(Individual *p1, std::vector<Individual *> &Union);
	
	/**
	 * 
	 */
	int calculate_raw_fitness(Individual *p1, std::vector<Individual *> &Union);
	
	/**
	 * 
	 */
	void calculate_dominated(std::vector<Individual *> &Union);
	
	/**
	 * 
	 */
	void calculate_objectives(std::vector<Individual *> &population);
	
	/**
	 * 
	 */
	void reproduce(genome *g, std::vector<Individual *> &population, std::vector<Individual *> &selected);

	/**
	 * 
	 */
	void objectiveOne(Individual *);
	
	/**
	 * 
	 */
	void objectiveTwo(Individual *);
	
	/**
	 * 
	 */
	void objectiveThree(Individual *);
	
	/**
	 * 
	 */
	double euclidian_distance(std::vector<double> objectivesI1, std::vector<double> objectivesI2);

	
	
	/* Global parameters. */
	extern int gen;       			 /* Number of generations elapsed. 	*/
	extern int popsize;   			 /* Population size.               	*/
	extern int archiveSize;			 /* Archive size. 					*/
	extern genome *g; 				 /* Genome.                        	*/
	extern int numberOfObjectives;
	extern double mutationRate;
	extern double crossoverRate;


#endif /* _GENETIC_ALGORITHM_H_ */
