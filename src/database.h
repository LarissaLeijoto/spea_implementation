#include <string.h>
#include <string>
#include <vector>


#ifndef DATABASE_H_
#define DATABASE_H_


	/**
	* @brief Atoms 
	*/
	struct Atom
	{
		std::string aminoacids_at;
		std::string atoms;
		float x,y,z;
		std::string chain_at;
		int id;
	};
	
	/**
	* @brief Amino acids 
	*/
	struct Aminoacid
	{
		std::string aminoacids_aa;
		std::vector<Atom*> listAtoms;
	};

	/**
	 * @brief Database.
	 */
	struct database
	{
	    /** Data of amino acids and atoms*/
		Aminoacid ***data_aminoacids;/* amino acids data.			 */
		Atom ***data_atoms;      	 /* atoms data.					 */
		
	    int maxaminoacids;		 	/*largest number of amino acids. */
	    int minaminoacids; 	 		/*smallest number of amino acids.*/        
		int *naminoacids;  	 		/*number of amino acids 
									   for each protein.  */
	 	
		int maxatoms;			 	/*largest number of atoms.       */
	    int minatoms; 			 	/*smallest number of atoms  	 */
	    int *natoms; 			 	/*number of atoms for each 
									   protein. 		*/
	};
	
	
	/* Forward definitions. */
	extern void database_read(const char **, int);
	extern void database_parse(const char **, int);
	extern void protein_matches(int popsize, int ngen);
	extern void database_destroy(void);
	extern void print_base(void);
	extern void readMatrix();
	
	/* Forward definitions. */
	extern int nproteins;
	extern int **Blosum;
	extern int **Clesum;
	extern struct database database;
	

	
#endif /* DATABASE_H_ */
