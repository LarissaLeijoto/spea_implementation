#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <time.h>
#include <omp.h>
#include <time.h>
#include "database.h"
#include "util.h"
#include "geneticAlgorithm.h"

using namespace std;

	
#define dbg  0
#define time_program 0

/* Program parameters. */
static const char **filenames = NULL; /* Name of input files.               */
int nproteins = 0;               	  /* Number of proteins (input files).  */
int popsize = 10;		         	  /* Population size.                   */
int ngen = 50;          		 	  /* Number of generations.             */
int numberOfObjectives = 1;			  /* Number of external individuals.	*/
int archiveSize = 10;				  /* Archive size.						*/
double mutationRate = 0.1;
double crossoverRate = 0.8;
/**
 * @brief Database.
 */
struct database database;


/**
 * @brief Prints program usage and exits.
 * 
 * @details Prints program usage and exits.
 */
static void usage(void)
{
	printf("Usage: genetic algorithm <popsize><Archive size><ngen><Crosssover Rate><Mutation Rate><Number Of Objectivies><protein files>\n");
	exit(EXIT_SUCCESS);
}

/**
 * @brief Reads command line arguments.
 * 
 * @details Reads and parses command line arguments pointed to by @p argv.
 * 
 * @todo Read and parse command line arguments.
 */
static void readargs(int argc, char **argv)
{
	
	/* Missing arguments. */
	cout<<argc<<endl<<endl;
	if (argc < 8)
		usage();
	
	popsize = atoi(argv[1]);
	archiveSize = atoi(argv[2]);
	ngen = atoi(argv[3]);
	crossoverRate = stod(argv[4]);
	mutationRate = stod(argv[5]);
	numberOfObjectives = atoi(argv[6]);	
	
	
	/* Count the number of proteins. */
	for (int i = 7; argv[i] != NULL; i++)
		nproteins++;
		
	filenames = (const char **)smalloc(nproteins*sizeof(char *));
	
	/* Extract protein files. */
	for (int i = 0; i < nproteins; i++)
		filenames[i] = argv[7 + i];
	
	/* Assert program parameters. */
	 if (popsize == 0)
		error("invalid population size");
	else if (ngen == 0)
		error("invalid number of generations");
	else if(numberOfObjectives < 1)
		error("invalid number of objectives");
	else if(archiveSize == 0)
		error("invalid archive size");
}

int main(int argc, char **argv)
{
	clock_t start_time = 0;
	start_time = clock();
	srand( (unsigned)time(NULL) );
	//srand( 198745 );
	//omp_set_num_threads(4);

	readargs(argc, argv);
	
	#if(dbg>0)
	fprintf(stderr, "info: parsing database! Ok \n");
	#endif
	/* Parse database in order to determine the largest number of amino acids 
	   among all proteins.*/
	database_parse(filenames, nproteins);
	
	#if(dbg>0)
	fprintf(stderr, "info: reading database! Ok\n");
	#endif
	/* Read database. */
	database_read(filenames, nproteins);
	
	//readMatrix();
	
	#if(dbg>0)
	fprintf(stderr, "info: proteins matches...\n");
	#endif	
	protein_matches(popsize, ngen, crossoverRate, mutationRate);
	
	//print_base();
	
	/* House keeping. */
	database_destroy();
	
	double time_in_seconds = (clock() - start_time) / (double)CLOCKS_PER_SEC;
	fprintf(stderr, "program time main: %.2f\n",time_in_seconds );
	
	
	free(filenames);
	
	return (EXIT_SUCCESS);
}
