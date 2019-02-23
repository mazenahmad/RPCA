
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "arrays.h"
#include "strTools.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////========================================================================================/////

/////////////////////////////////////////////////////////////////////////////////////////////////
 void writepdblines2file(FILE *fout,
						char * title,
						int model,
						int natoms,
						t_atoms  *atoms,
						int *index,
						rvec *x,
						int *x_index,
						double *bfactors)
{

	if (fout  == NULL)
	{
		fprintf(stderr, "Cannot write to the pdb file, exit\n");
		exit(1);
	}
	/* THE pdb format (for ATOM/HETATOM lines) */
	static const char *pdbformat ="%-6s%5u %-4.4s%c%3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           %-2s \n"  ;
	int i,ndx, x_id;
	const char *Record= "ATOM  ";
	char altLoc = ' ';
	float occup =1 , bfac;
	if(model) fprintf(fout,"MODEL %d\n",model);
	fprintf(fout,"TITLE %s\n",title);
	for (i = 0; i < natoms; i++)
	{
		if(index) ndx=index[i];
		else      ndx=i;
		if(x_index) x_id= x_index[i];
		else        x_id= i;
		if(bfactors) bfac= bfactors[i];
		else bfac=0;


		fprintf(fout, pdbformat,
                                    		Record,
											ndx+1,//atoms->atom[ndx].atomnumber,
											*(atoms->atomname[ndx]),
		                                    altLoc,
											*(atoms->resinfo[atoms->atom[ndx].resind].name),
											atoms->resinfo[atoms->atom[ndx].resind].chainid ,
											atoms->resinfo[atoms->atom[ndx].resind].nr,
		                                    x[x_id][0]*10,
											x[x_id][1]*10,
											x[x_id][2]*10,
		                                    occup,
		                                    bfac,
											atoms->atom[ndx].elem );
	}
	if(model) fprintf(fout,"ENDMDL\n");

}

/////////////////////////////////////////////////////////////////////////////////////////////////

void writepdb(const char *fname,
						char * title,
						int model,
						int natoms,
						t_atoms  *atoms,
						int *index,
						rvec *x,
						int *x_index,
						double *bfactors,
						int append)
{

	FILE *fout;
	if(append) fout = fopen(fname, "a");
	else fout = fopen(fname, "w");
	if (fout  == NULL)
	{
		fprintf(stderr, "Cannot open file %s for writing, exit\n",	fname);
		exit(1);
	}
	writepdblines2file(fout, title, model,natoms,atoms,index,x,x_index,bfactors);
	//fprintf(fout,"END   \n");
	fclose(fout);

}

/////////////////////////////////////////////////////////////////////////////////////////////////

////========================================================================================/////

/////////////////////////////////////////////////////////////////////////////////////////////////
/************************************************************************************************
 ********************  struct for the info of one atom in pdb file (one line ) ******************
 ***********************************************************************************************/
/////////////////////////////////////////////////////////////////////////////////////////////////

typedef struct
{
	char  Record[7];              /* ATOM or HETATEM                                           */
	int  atomnr;                  /* PDB atom number                                           */
	char name[5];                 /* atom name                                                 */
	char altLoc;                  /* Alternate location                                        */
	char resname[5];              /* residue name                                              */
	char chid ;                   /* chainID                                                   */
	int  resi;                    /* res number                                                */
	char insert[2];               /* Code for insertion of residues.                           */
	real x;
	real y;
	real z ;
	real occup;                   /* Occupancy                                                 */
	real bfac;                    /* B-factor                                                  */
	char element[3];              /* Element symbol, right-justified                           */
} pdbAtom_line_t;

typedef struct
{
		pdbAtom_line_t *array;
		int used;
		int size;
} pdbAtom_line_t_cap;
////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/************************************************************************************************
 *							   read the information in an atom line 							*
 *																								*
*************************************************************************************************/



pdbAtom_line_t read_pdbatomline(char line[] )
{
	pdbAtom_line_t pdbatom;
	char tmp[12], tmp2[7];
	int  l,i;
	l=0;
	for(i=0; (i<6); i++,l++) pdbatom.Record[i]=line[l]; /*atom or hetatm*/
	l=6;
    sscanf(line + 6, "%d", &(pdbatom.atomnr)); /* atom number */
    l=12;
    for(i=0; (i<4); i++,l++) pdbatom.name[i]=line[l];/* atom name */
    l=16;
    //sscanf(line +l, "%c", &(pdbatom.altLoc));
    pdbatom.altLoc =line[l];
    l=17;
    for(i=0; (i<3); i++,l++) pdbatom.resname[i]=line[l];
    l=21;
    sscanf(line + l, "%c", &(pdbatom.chid)); /* chain id number */
    l= 22;
    sscanf(line + l, "%d", &(pdbatom.resi)); /* residue number */
    /* X,Y,Z Coordinate */
    l=30;
    for(i=0; (i<8); i++,l++) tmp[i]=line[l];/* x */
    pdbatom.x= strtod(tmp,NULL)*0.1;
    l=38;
    for(i=0; (i<8); i++,l++) tmp[i]=line[l];/* y */
	pdbatom.y= strtod(tmp,NULL)*0.1;
	l=46;
	for(i=0; (i<8); i++,l++) tmp[i]=line[l];/* z */
	pdbatom.z= strtod(tmp,NULL)*0.1;
	l=54;
	for(i=0; (i<6); i++,l++) tmp2[i]=line[l];/* occup */
	pdbatom.occup= strtod(tmp2,NULL);
	l=60;
	for(i=0; (i<6); i++,l++) tmp2[i]=line[l];/* bfac */
	pdbatom.bfac= strtod(tmp2,NULL);
	l=76;
	for(i=0; (i<2); i++,l++) pdbatom.element[i]=line[l]; /* element*/
	 return pdbatom;
}
/////////////////////////////////////////////////////////////////////////////////////////////////

////========================================================================================/////
////========================================================================================/////

/////////////////////////////////////////////////////////////////////////////////////////////////
/************************************************************************************************
* function to read a pdb file  and putting the information of the atoms                         *
* in a structure:
* altLoc : which altLoc to read
* pdbatoms: structure for pdbAtoms_t
* pdbrecords: structure for pdbrecord_t
* pdbAtoms_t	 *pdbatoms;
* INIT_CAP_STRUC(pdbatoms ,pdbAtom_line_t,100) ;
* INIT_PDBINF(pdbrecords);
*************************************************************************************************/
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

 void Readpdb2atoms(	char *fname,
		 	 	 	 	 pdbAtom_line_t_cap *pdbAtoms,
		 	 	 	 	char altLoc  )
{

	 //if (altLoc == NULL) altLoc= ' ';
	 FILE   *fpin;
	#define STRLEN 4096
	if ((fpin = fopen(fname, "r")) == NULL)  PRINT_ERR_OPENFILE(fname);
	char line[STRLEN];
	pdbAtom_line_t tmpatom;
	for (;;)
	{
		if (fgets(line, STRLEN, fpin) == NULL) break;
		if(( ! strncmp("ATOM  ", line, 6) || ! strncmp("HETATM", line, 6))  && line[16]== altLoc )
		{
			tmpatom= read_pdbatomline(line );
			INSERT_IN_CAP_STRUC(pdbAtoms,tmpatom,+100);
		}
		if ( ! strncmp("END   ", line, 6) ||  ! strncmp("ENDMDL", line, 6) )
			break;

	}
	fclose(fpin);
	SHRINK_CAP_STRUC(pdbAtoms);

}
