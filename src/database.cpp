#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <cstdio>
#include <cstdlib>
#include <limits.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "database.h"

/* program's debug*/
#define dbg 1
#define time_program 0

using namespace std;

/**
 * @brief Creates a atom.
 * @returns A new atom.
 */
	struct Atom *atom_create(void)
	{
		struct Atom *at;
		
		at = new Atom;
		at->atoms = "\0";
		at->x = 0; 
		at->y = 0;
		at->z = 0;
		at->id = 0;
		
		return (at);
	}
/**
 * @brief Creates a amino acid.
 * @returns A new amino acid.
 */	
	struct Aminoacid *aminoacid_create(void)
	{
		struct Aminoacid *aa ;
		aa = new Aminoacid;
		aa->aminoacids_aa = "\0";
		
		return (aa);
	}

/**
 * @brief Parses the database.
 * @details Parses the database so that we can determine the largest number of
 *          amino acids among all proteins.
 *
 * @param filenames Name of input files.
 * @param nproteins Number of proteins (number of input files).
 * @returns The largest number of amino acids among all proteins.
 * 
 * @todo Check for bad file format.
 */
void database_parse(const char **filenames, int nproteins)
{
	/* Sanity checks. */
	assert(filenames != NULL);
	assert(nproteins != 0);
	
	/** Amino acids initialization */
	database.maxaminoacids = 0;
        database.minaminoacids = INT_MAX;
	database.naminoacids = new int*[nproteins]; /* Number of amino acids of each protein */

	/** Atoms initialization */
	database.maxatoms = 0;
	database.minatoms = INT_MAX;
	database.natoms =  new int*[nproteins];   /* Number of atoms of each protein */

	/* Find largest and smallest number of amino acids among all proteins. */
	for (int wprotein = 0; wprotein < nproteins; wprotein++)
	{
		FILE *wfile; 		/* Working file. 				    */
		char *line;         /* Working line.    				*/
		int cToken;			/* Count token.     				*/
		string lable;       /* Define if is a atom or amino acid*/
		
		/* Open working file. */
		wfile = fopen(filenames[wprotein], "r");
		if (wfile == NULL)
			error ("cannot open input file");
			
		/* Get number of amino acids and atoms. */
		database.natoms[wprotein] = 0;
		
		while (!feof(wfile))
		{
			char *token;       	/* Working token.  */

			line = readline(wfile);
			cToken = 0;
			
			/* Read line. */
			token = strtok(line," ");
			while (token != NULL)
			{
				if(cToken == 0)
					lable = string(token);
				
				token = strtok(NULL," ");
				cToken++;
			}
				if(strcmp(lable.c_str(),"ATOM") == 0)
					database.natoms[wprotein]++;

			free(line);
			free(token);
		}
		
		database.natoms[wprotein]--;

		/* Largest number of atoms found. */
		if (database.maxatoms < database.natoms[wprotein])
			database.maxatoms = database.natoms[wprotein];
		/* smallest  number of atoms found. */
	    if (database.minatoms > database.natoms[wprotein])
			database.minatoms = database.natoms[wprotein];
		
		fclose(wfile);
	}
		/* Invalid number of amino acids. */
		if (database.maxatoms == 0)
			error("invalid number of atoms");
		/* Invalid number of amino acids. */
		if (database.minatoms == 0)
			error("invalid number of atoms");
}

/**
 * @brief Reads the database into memory.
 * @details Reads the database into memory so we can speedup 
 * computation.
 * @param filenames   Name of input files.
 * @param nproteins   Number of proteins (number of input files).
 *
 * @returns The database.
 * 
 * @todo Check for bad file format.
 */
void database_read (const char **filenames, int nproteins)
{
	/* Sanity check. */
	assert(filenames != NULL);
	assert(nproteins != 0);
	
	
	#if(dbg>0)
	fprintf(stderr,"Max number of atoms: %d \n", database.maxatoms);
	fprintf(stderr,"Number of protein: %d \n", nproteins);
    #endif		
	
	/* Allocate database. */
	database.data_atoms = new Atom**[nproteins];
	
	for (int i = 0; i < nproteins; i++)
		database.data_atoms[i] = new Atom*[database.maxatoms];
		
	/* Read database. */
	for (int wprotein = 0; wprotein < nproteins; wprotein++)
	{
		int wline;		/* Working line.	*/
		int watom;		/* Current atom.    */
   		int cToken;		/* Count token.     */
		FILE *wfile;	/* Working file.    */
		char *line;     /* Working line.    */
		char *token;    /* Working token.   */

		/* Open working file. */
		wfile = fopen(filenames[wprotein], "r");
		
		if (wfile == NULL)
			error ("cannot open input file");
		
		/* Read amino acid and atom. */
		wline = 0;
		watom = 0;
		
		while (!feof(wfile))
		{
			line = readline(wfile);
			/* Read line. */
			token = strtok(line, " ");
			if(token != NULL)
			{
				if(strcmp(token,"ATOM") == 0)
				{
					Atom *at = atom_create();

					cToken = 0;
					while (token != NULL)
					{
						if(token!=NULL)
						{							

							token = strtok(NULL, " ");
							if(cToken == 0)
							{
								//
								//fprintf(stderr, "%s",token);
							}
							else if(cToken == 1)
							{
								//fprintf(stderr, "%s",token);
								at->atoms = token;
							}
							else if(cToken == 2)
							{
								//fprintf(stderr, "%s",token);
								at->aminoacids_at = token;
								//at->chain_at = token;
							}
							else if(cToken == 3)
							{
								//fprintf(stderr, "%s",token);
								at->id = atoi(token);
							}
							else if(cToken == 4)
							{
								//fprintf(stderr, "%s",token);
								at->x = (float)atof(token);
							}
							else if(cToken == 5)
							{
								//fprintf(stderr, "%s",token);
								at->y = (float)atof(token);
							}
							else if(cToken == 6)
							{
								//fprintf(stderr, "%s",token);
								at->z = (float)atof(token);
							}
							//cerr<<endl;
								
							cToken++;
						}	
					}
				 database.data_atoms[wprotein][watom] = at	;
				 watom++;
				}
			}		
			wline++;
			free(line);
			free(token);
		}
		fclose(wfile);
	}
	
	int previousAminoacid = INT_MAX;
		
	for (int wprotein = 0; wprotein < nproteins; wprotein++)
	{
		database.naminoacids[wprotein] = 0;
		
		for (int watom = 0; watom < database.natoms[wprotein]; watom++)
		{
			if(previousAminoacid != database.data_atoms[wprotein][watom]->id) 
			{
				previousAminoacid = database.data_atoms[wprotein][watom]->id;
				database.naminoacids[wprotein]++;
			}
		}
		
		/* Largest number of amino acids found. */
		if (database.maxaminoacids < database.naminoacids[wprotein])
			database.maxaminoacids = database.naminoacids[wprotein];
		/* smallest  number of amino acids found. */
	    if (database.minaminoacids > database.naminoacids[wprotein])
			database.minaminoacids = database.naminoacids[wprotein];
		cout<<database.naminoacids[wprotein]<<endl;
	}	
	
	
	/* Allocate database. */
	database.data_aminoacids = new Aminoacid**[nproteins];
	
	for (int i = 0; i < nproteins; i++)
		database.data_aminoacids[i] = new Aminoacid*[database.maxaminoacids];	
	
	
	for (int wprotein = 0; wprotein < nproteins; wprotein++)
	{		
		int waminoacid = -1;
		previousAminoacid = INT_MAX;

		for (int watom = 0; watom < database.natoms[wprotein]; watom++)
		{
			if(previousAminoacid != database.data_atoms[wprotein][watom]->id) 
			{
				waminoacid++;
				previousAminoacid = database.data_atoms[wprotein][watom]->id;
				database.data_aminoacids[wprotein][waminoacid] = aminoacid_create();
				database.data_aminoacids[wprotein][waminoacid]->aminoacids_aa = database.data_atoms[wprotein][watom]->aminoacids_at;
				database.data_aminoacids[wprotein][waminoacid]->listAtoms.push_back(database.data_atoms[wprotein][watom]);
			}
			else
				database.data_aminoacids[wprotein][waminoacid]->listAtoms.push_back(database.data_atoms[wprotein][watom]);
		}
		
	}	
	
}

void readMatrix()
{
	ifstream fileBlosum;
	fileBlosum.open("../blosum.txt");
	int number;
	
	int **Blosum = (int **)smalloc(23*sizeof(int*));
	
	for(int count = 0; count < 23 ;count++)
		Blosum[count] = (int *)smalloc(23*sizeof(int));
	
	for(int line = 0; line < 23 ;line++)
	{
		for(int column = 0; column < 23;column++)
		{
			fileBlosum >> number;
			Blosum[line][column] = number;
			cout<< Blosum[line][column] <<" ";
		}
		cout<<endl;
	}	
	
	fileBlosum.close();
	
	cout<<endl<<endl<<endl;
	
	ifstream fileClesum;
	fileClesum.open("../clesum.txt");
	
    int **Clesum = (int **)smalloc(17*sizeof(int*));
	
	for(int count = 0; count < 17 ;count++)
		Clesum[count] = (int *)smalloc(17*sizeof(int));
		
    for(int line = 0; line < 17 ;line++)
	{
		for(int column = 0; column < 17;column++)
		{
			fileClesum >> number;
			Clesum[line][column] = number;
			cout<<Clesum[line][column]<<" ";
		}
		cout<<endl;
	}	
	fileClesum.close();
}


/**
 * @brief Destroys the database.
 * 
 * @details Destroys the database, freeing underlying resources.
 */
void database_destroy()
{		
	int i = 0, j = 0; /* Loop index*/
	
	for (i = 0; i < nproteins; i++)
	{
		for(j = 0; j < database.naminoacids[i];j++ )
			delete database.data_aminoacids[i][j];
			
		for(j = 0; j < database.natoms[i];j++ )
			delete  database.data_atoms[i][j];
			
		delete [] database.data_aminoacids[i];
		delete [] database.data_atoms[i];
	}	
	
	delete [] database.data_aminoacids;
	delete [] database.data_atoms;
	
	free(database.naminoacids);
	free(database.natoms);	
}

/**
 * @brief Print the database.
 * 
 * @details Prints the database, when is necessary.
 */

void print_base(void)
{
	for(int i = 0; i < nproteins;i++)
	{
		for(int k = 0;k < database.naminoacids[i];k++)
		{
			//fprintf(stderr,"%d %s %s\n",k,((database.data_aminoacids[i][k])->chain_aa).c_str(), ((database.data_aminoacids[i][k])->aminoacids_aa).c_str());
		}
	}
	
	for(int i = 0; i < nproteins;i++)
	{
		for(int k = 0;k < database.natoms[i];k++)
		{
			fprintf(stderr,"%d %s %s %u %.3f %.3f %.3f\n",k,(database.data_atoms[i][k])->aminoacids_at.c_str(),
			(database.data_atoms[i][k])->atoms.c_str(), 
			//(database.data_atoms[i][k])->chain_at.c_str(), 
			(database.data_atoms[i][k])->id, 
			(database.data_atoms[i][k])->x, 
			(database.data_atoms[i][k])->y, 
			(database.data_atoms[i][k])->z);
		}
	}	
	
	
}
