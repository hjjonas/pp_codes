// storage.h

// dynamics structs 

#ifndef STORAGE_H
#define STORAGE_H

/*------------------GLOBAL VARIABLES-----------------------------------------*/
/*---------------------------------------------------------------------------*/

/*------------------STRUCTURE DEFINITIONS------------------------------------*/
typedef struct IntArraytype {
  // a list type of intergers
  int *ints;
  size_t used;
  size_t size;

} IntArray ;
/*---------------------------------------------------------------------------*/

/*------------------GLOBAL FUNCTIONS-----------------------------------------*/
extern void initIntArray(IntArray *, size_t );
extern void insertIntArray(IntArray *, int ); 
extern void freeIntArray(IntArray *);
extern void removeIthElementIntArray(IntArray *, int  );
extern void removeElementXIntArray(IntArray *, int  );
extern int checkElementXIntArray(IntArray *, int );
extern void memory_allocation(void);
extern void printIntArray(IntArray * );
/*---------------------------------------------------------------------------*/

/*------------------LOCAL FUNCTIONS------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*------------------GLOBAL STRUCTURES----------------------------------------*/
/*---------------------------------------------------------------------------*/




#endif
