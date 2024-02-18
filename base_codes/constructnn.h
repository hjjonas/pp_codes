// constructnn.h

#ifndef NN_H
#define NN_H

/*------------------GLOBAL VARIABLES-----------------------------------------*/
// definitions for neighborlist
#define MAXPARTS            NPART
#define MAX_IMAGES          (MAXPARTS *8) /* is MAXPARTS NPART?*/
#define MAX_NEIGHBORS       (MAXPARTS *300)
#define MAX_IM_NEIGHBORS    300
#define STACKDEPTH          30
#define NSHIFTS             27

#define GRIDOFFSET          3
#define LJ_RZERO            1.5
#define MAXBOXSIZE          20
#define MINCELLSIZE         (LJ_RZERO/GRIDOFFSET)
// #define NCELLX              ((int) (MAXBOXSIZE / MINCELLSIZE)) gives errors?? Just hardcoded it:
#define NCELLX              40
#define NCELLY              (NCELLX)
#define NCELLZ              (NCELLX)
/*---------------------------------------------------------------------------*/

/*------------------STRUCTURE DEFINITIONS------------------------------------*/
typedef struct cell2_type {
 // for nearest neighbor
  int                n                 ,
                     stack[STACKDEPTH] ;
} Cell;

typedef struct nlist_item_type {
 // for nearest neighbor
  int                ipart             ,
                     image_id          ,
                     *first           ;
} List_item;

typedef struct list_type {
 // for nearest neighbor
  int                cellx             ,
                     celly             ,
                     cellz             ;  

  vector             cellsize          ,
                     inv_cellsize      ,
                     trans[NSHIFTS]    ;

  List_item          nli[MAX_IMAGES]    ;

  int                nlj[MAX_NEIGHBORS] ,
                     nneighbors         ,
                     nimages           ,
                     count             ;

  double             cutoff            ,
                     cutoff2           ,
                     cutoff_buffer2    ,
                     bigpart_cutoff2   ;

  Cell               cell[NCELLX][NCELLY][NCELLZ];

} List;
/*---------------------------------------------------------------------------*/

/*------------------GLOBAL FUNCTIONS-----------------------------------------*/
extern void setup_nnlist();
extern void update_nnlist(Slice *);
/*---------------------------------------------------------------------------*/

/*------------------GLOBAL STRUCTURES-----------------------------------------*/
List list;
/*---------------------------------------------------------------------------*/


#endif
