///////////////////////////////////////////////////////////////////////////////////////////////////////////
/**********************************************************************************************************
 * arrays.h
 *
 *  Created on: Mar 29, 2013
 *      Author: mazen
 *********************************************************************************************************/
///////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ARRAYS_H_
#define ARRAYS_H_

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>





#ifndef FALSE
#define  FALSE   (0)
#endif
#ifndef TRUE
#define  TRUE    (1)
#endif


/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////





#define PRINT_ERR_LOC fprintf(stderr,"\nThe error is in function %s(): file %s:line %d \n" \
										, __func__,	__FILE__,__LINE__)

#define PRINT_MEM_ERR  do										\
						{										\
							fprintf(stderr,"Out of memory\n");	\
							PRINT_ERR_LOC;						\
							exit(1);							\
						} while (0)

static void printmessge(char *message)
{
	fprintf(stderr,"Error: %s\n",message);
}

#define PRINT_FATAL_ERR_MES(message)							\
						do								\
						{								\
							printmessge(message);		\
							PRINT_ERR_LOC;				\
							exit(1);					\
						} while (0)

#define PRINT_ERR_OPENFILE(file)						\
						do								\
						{								\
							fprintf(stderr,"Can not open file %s\n",file);\
							PRINT_ERR_LOC;				\
							exit(1);					\
						} while (0)

#define FUNCTION_CALLED() 		printf( " %s called\n " , __func__)
#define FUNCTION_RETURNS()   	printf( " %s returns\n " , __func__)

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////   1D array   ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
////// declare and allocate 1D array of type type 
#define ONE_DIM_ARRAY(ptr,type,n)  					\
	type * ptr ;									\
	do												\
	{												\
		if((ptr= (type *) malloc( (n) *sizeof(type)) )==NULL) PRINT_MEM_ERR;\
	} while(0)
		
/////////////////////////////////////////////////////////////////////////////////////
/////// allocate (only) predeclared 1D array of type type
#define ALLO_ONE_DIM_ARRAY(ptr,type,n)  		   \
	do												\
	{												\
		if((ptr= (type *) malloc( (n) *sizeof(type)) )==NULL) PRINT_MEM_ERR;\
	} while(0)
/////////////////////////////////////////////////////////////////////////////////////
/////// zero 1D array of type type
#define ZERO_ONE_DIM_ARRAY(ptr,n)  		            \
	do												\
	{												\
		int ii;									    \
	for (ii=0; ii< n; ii++) ptr[ii]= 0;             \
	} while(0)

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////   2D array   ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

////// declare and allocate 2D array of type type 

#define TWO_DIM_ARRAY(ptr,type,row, col)  			\
	type ** ptr ;									\
	do												\
	{												\
		int ii ;								    \
		if((ptr= (type **) malloc( (row) *sizeof(type*)) )==NULL) PRINT_MEM_ERR;\
		for(ii=0; ii< row; ii++) if((ptr[ii]= (type *) malloc( (col) *sizeof(type)) )==NULL) PRINT_MEM_ERR;\
	} while(0)

/////////////////////////////////////////////////////////////////////////////////////
////// allocate (only) predeclared 2D array of type type 

#define ALLO_TWO_DIM_ARRAY(ptr,type,row, col)  		\
	do												\
	{												\
		int ii ;								    \
		if((ptr= (type **) malloc( (row) *sizeof(type*)) )==NULL) PRINT_MEM_ERR;\
		for(ii=0; ii< row; ii++) if((ptr[ii]= (type *) malloc( (col) *sizeof(type)) )==NULL) PRINT_MEM_ERR;\
	} while(0)
/////////////////////////////////////////////////////////////////////////////////////
////////////////
#define FREE_TWO_DIM_ARRAY(ptr,n)  			        \
	do												\
	{												\
		int ii ;								    \
		for(ii = 0; ii < n; ii++)	free(ptr[ii]);  \
		free(ptr);									\
	} while(0)
		
/////////////////////////////////////////////////////////////////////////////////////
/////// zero 2D array of type type
#define ZERO_TWO_DIM_ARRAY(ptr,n, m)  		         				    \
	do																	\
	{																	\
		int ii, jj;									   					\
		for (ii=0; ii< n; ii++) for(jj=0; jj< m; jj++) ptr[ii][jj]= 0;  \
	} while(0)

/////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////

/*********************************************************************************************************
 **	*	*	*	*	*	*	*			dynamically growing  arrays 		*	*	*	*	*	*	*	**
 *********************************************************************************************************/

/////////////////////////////////////////////////////////////////////////////////////
/************************************************************************************
 *  --------------------------Automatic type defining ------------------------------*
 *  A new capsule structure (type_cap) for any type (including structures)      	*
 *  of arrays. This capsule structure contains dynamically growing  array of the 	*
 *  type (type). 																	*
 *  The following functions can be used to handle this array.  						*
 *         																			*
 *  type : the type of the array that we want to make it dynamically growing        *
 *  the macro will return a new structure type with the pasted name type_cap		*
 ************************************************************************************/
/////////////////////////////////////////////////////////////////////////////////////

#define DEF_CAP_STRUC_TYPE(type)    \
	typedef struct					\
	{								\
		type *array; 				\
		int used;					\
		int size; 					\
	} type##_cap


/////// for array of real type
/// DEF_CAP_STRUC_TYPE(real);
/***
typedef struct
	{
		real *array;
		int used;
		int size;
	} real_cap;
***/
/////////////////////////////////////////////////////////////////////////////////////

/************************************************************************************
 *  -------------Declaration and initialization of the capsule structure------------*
 * declare and initialize the capsule structure that contain the dynamically		*
 * growing array of type (type).         											*
 *  type : the type of the array that we want to make it dynamically growing        *
 *  																				*
 *  ptr: pointer to the capsule structure that contains the array  of type "type" 	*
 *  																				*
 * initialSize: the initial size of the array to allocate( Nr. of elements).  		*
 ************************************************************************************/
/////////////////////////////////////////////////////////////////////////////////////

#define INIT_CAP_STRUC(ptr,type,initialSize)  												\
			if((ptr= (type##_cap *)	malloc( sizeof(type##_cap)) )==NULL) PRINT_MEM_ERR;		\
			do {																			\
				if ((ptr->array =(type *) malloc( (initialSize) * sizeof(type)))==NULL)		\
					PRINT_MEM_ERR;															\
				ptr->used = 0;																\
				ptr->size = initialSize;													\
			} while(0)
/////////////////////////////////////////////////////////////////////////////////////


/************************************************************************************
 *  ------------------ Insert and element in the capsule array ---------------------*
 *  																				*
 *  ptr: pointer to the capsule structure that contains the array  of type "type" 	*
 *  																				*
 * initialSize: the initial size of the array to allocate( Nr. of elements).  		*
 * 																					*
 * Increment: the way how to increment the array in case it is full ex. *2 or +100  *
 ************************************************************************************/

#define INSERT_IN_CAP_STRUC(ptr,element,incremnet)			    \
do																\
{																\
	if (ptr->used == ptr->size)									\
	{															\
		ptr->size = ptr->size incremnet;						\
		if ((ptr->array = realloc(ptr->array, ptr->size * sizeof(ptr->array[0])))==NULL)\
		{														\
				PRINT_MEM_ERR;									\
		}														\
	}															\
	ptr->array[ptr->used++] = (element);						\
} while(0)

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


/************************************************************************************
 *  --------------- Shrink the capsule array to it is used size---------------------*
 *  																				*
 *  ptr: pointer to the capsule structure that contains the array  of type "type" 	*
 *  																				*
 ************************************************************************************/

#define SHRINK_CAP_STRUC(ptr)								    \
do																\
{																\
	if (ptr->used < ptr->size)									\
	{															\
		(ptr)->size = ptr->used ;         						\
		if(ptr->used == 0)	(ptr)->size =1;						\
		ptr->array = realloc(ptr->array, ptr->size * sizeof(ptr->array[0]));\
	}															\
} while(0)

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/************************************************************************************
 *  --------------- resize the capsule array to it is used size---------------------*
 *  																				*
 *  ptr: pointer to the capsule structure that contains the array  of type "type" 	*
 *  																				*
 * Increment: the way how to increment the array in case it is full ex. *2 or +100  *
 ************************************************************************************/

#define EXTEND_CAP_STRUC(ptr,increment)							\
do																\
{																\
	ptr->size = ptr->size increment	;							\
	if ((ptr->array = realloc(ptr->array, ptr->size * sizeof(ptr->array[0])))==NULL)\
					PRINT_MEM_ERR;								\
} while(0)

	/////////////////////////////////////////////////////////////////////////////////////

/************************************************************************************
 *  ------------------------ Delete capsule structure ------------------------------*
 *  ptr: pointer to the capsule structure that contains the array  of type "type" 	*
 ************************************************************************************/
#define DELETE_CAP_STRUC(ptr)		\
do									\
{									\
  free(ptr->array);					\
  ptr->used= ptr->size=0;			\
  free(ptr);						\
}while(0)
//////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------//

#endif /* ARRAYS_H_ */
