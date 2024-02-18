#include "path.h"



void initIntArray(IntArray *a, size_t initialSize) {
  //https://stackoverflow.com/questions/3536153/c-dynamically-growing-array
  //size_t initialSize is the intial size of the array. for example "5" intergers, do not keep it too small, as the array is doubled if it is too short 
  // a->used is the number of used entries, upon initialization it is zero.
  
  if(((a->ints = (int *) malloc(initialSize * sizeof(int) ))==NULL)) {;
       error("could not malloc initialSize * sizeof(int) ");
    }
  
  a->used = 0;
  a->size = initialSize;
  return;

}

void insertIntArray(IntArray *a, int element) {
  // the original array a is elongated if you want to add more elements that is currently has. Just for easiness it is doubled. 
  // NOT handy when the array becomes big!
  // a->used is the number of used entries, because a->ints[a->used++] updates a->used only *after* the array has been accessed.

  // Therefore a->used can go up to a->size 
  if (a->used == a->size-1) {
    a->size *= 2;
    if( (a->ints = realloc(a->ints, a->size * sizeof(int)))==NULL) {;
       error("could not realloc IntArray->ints ");
    }
    // a->ints = realloc(a->ints, a->size * sizeof(int)); // reallocate the elongated array with the new size.
  }
  a->ints[a->used] = element;
  a->used++;
  return;

}

void freeIntArray(IntArray *a) {
  free(a->ints); //free the memory
  a->ints = NULL;// make pointer zero
  a->used = a->size = 0;   
  return;

}

void removeIthElementIntArray(IntArray *a, int elemi ){
  // elemi stands for the ith element 
  // the array a has a->used  elements currently
  int start =0,i;

  if (a->used < elemi ){
    printf("the max length of the int arrat is %zu, but you want to remove element %d\n",a->used,elemi );
    error("Not possible");
  }

  a->used--;
  for (i=elemi;i<a->used;i++){
      a->ints[i]=a->ints[i+1];
  }
  return;
}

void removeElementXIntArray(IntArray *a, int elemX ){
  // elemX stands for element number X
  // the array a has a->used  elements currently
  int start =0,i;
  int last_int=(int)a->used-1;
  ///because you remove here one on forehand,  a[i]=a[i+1]; will not give problems
  // dprint((int)a->used);
  
  if(a->ints[last_int]==elemX){
      a->used--;
      return;
  }
  else {
    a->used--;
    for (i=0;i<a->used;i++){
      if(a->ints[i]==elemX){
        start=1;
      }
      if (start){
          a->ints[i]=a->ints[i+1];
      }
    }
  }
  
  if (start==0   ){
    // dprint((int)a->used+1);
    dprint(a->ints[last_int]);
    printf("you want to remove elemnt %d \n",elemX);
    printf("the remaining elements are :\n");

    for (i=0;i<a->used;i++){
      printf(" i = %d  ",i);
      dprint(a->ints[i]);
    }
    error("you are trying to remove an element that is not in the list.");
  }

  return;

}

int checkElementXIntArray(IntArray *a, int elemX ){
  // check if elementX is in the array
  // the array a has a->used  elements currently
  int check =0,i;

  for (i=0;i<a->used;i++){
      if(a->ints[i]==elemX){
        check=1;
        break;
      }      
  }

  return check;

}

void printIntArray(IntArray *a ){
  // check if elementX is in the array
  // the array a has a->used  elements currently
  int i;
  for (i=0;i<a->used;i++){
    printf(" item %d has value %d\n", i,a->ints[i]);
     
  }

  return ;

}




//    #####  #####  #####  #####  ##### BondTimeInfoArray ##### #####  #####  #####  ##### 
//    #####  #####  #####  #####  ##### BondTimeInfoArray ##### #####  #####  #####  ##### 
//    #####  #####  #####  #####  ##### BondTimeInfoArray ##### #####  #####  #####  ##### 
//    #####  #####  #####  #####  ##### BondTimeInfoArray ##### #####  #####  #####  ##### 
//    #####  #####  #####  #####  ##### BondTimeInfoArray ##### #####  #####  #####  ##### 
//    #####  #####  #####  #####  ##### BondTimeInfoArray ##### #####  #####  #####  ##### 


// void initBondTimeInfoArray(BondTimeInfoArray *a, size_t initialSize) {
//   //https://stackoverflow.com/questions/3536153/c-dynamically-growing-array
//   //size_t initialSize is the intial size of the array. for example "5" intergers, do not keep it too small, as the array is doubled if it is too short 
//   // a->used is the number of used entries, upon initialization it is zero.
  
//   // if(((a->ints = (int *) malloc(initialSize * sizeof(int) ))==NULL)) {;
//   if(((a->bondinfos = (BondTimeInfo *) malloc(initialSize * sizeof(BondTimeInfo) ))==NULL)) {;
//        error("could not malloc initialSize * sizeof(BondTimeInfoArray) ");
//     }
//   // printf("initBondTimeInfoArray done \n");
//   a->used = 0;
//   a->size = initialSize;
//   return;

// }

// void insertBondTimeInfoArray(BondTimeInfoArray *a, BondTimeInfo element) {
//   // the original array a is elongated if you want to add more elements that is currently has. Just for easiness it is doubled. 
//   // NOT handy when the array becomes big!
//   // a->used is the number of used entries, because a->ints[a->used++] updates a->used only *after* the array has been accessed.

//   // Therefore a->used can go up to a->size 
//   if (a->used == a->size-1) {
//     a->size *= 2;
//     if( (a->bondinfos = realloc(a->bondinfos, a->size * sizeof(BondTimeInfo)))==NULL) {;
//        error("could not realloc a->bondinfos ");
//     }
//     // a->bondinfos = realloc(a->bondinfos, a->size * sizeof(BondTimeInfoArray)); // reallocate the elongated array with the new size.
//   }
//   a->bondinfos[a->used] = element;
//   a->used++;
//   return ;

// }

// void insertBondTimeInfoArray_atIth(BondTimeInfoArray *a, BondTimeInfo element, int ith_place) {
//   // the original array a is elongated if you want to add more elements that is currently has. Just for easiness it is doubled. 
//   // NOT handy when the array becomes big!
//   // a->used is the number of used entries, because a->ints[a->used++] updates a->used only *after* the array has been accessed.

//   // Therefore a->used can go up to a->size 
//   if (a->used == a->size-1) {
//     a->size *= 2;
//     if( (a->bondinfos = realloc(a->bondinfos, a->size * sizeof(BondTimeInfo)))==NULL) {;
//        error("could not realloc a->bondinfos ");
//     }
//   }

//   for (int i=a->used;i>ith_place;i--){
//       a->bondinfos[i+1]=a->bondinfos[i];
//   }
//   a->bondinfos[ith_place]=element;
//   a->used++;

//   return ;

// }

// void freeBondTimeInfoArray(BondTimeInfoArray *a) {
//   free(a->bondinfos); //free the memory
//   a->bondinfos = NULL;// make pointer zero
//   a->used = a->size = 0;   
//   return;

// }

// void removeIthElementBondTimeInfoArray(BondTimeInfoArray *a, int elemi ){
//   // elemi stands for the ith element 
//   // the array a has a->used  elements currently
//   int i;

//   if (a->used < elemi ){
//     printf("the max length of the  array is %zu, but you want to remove element %d\n",a->used,elemi );
//     free(a);
//     error("Not possible");
//   }

//   a->used--;
//   for (i=elemi;i<a->used;i++){
//       a->bondinfos[i]=a->bondinfos[i+1];
//   }
//   return;
// }

// int checkElementXBondTimeInfoArray_atIth(BondTimeInfoArray *a, int ipart , int jpart , int start_at ){
//   // check if elementX is in the array
//   // the array a has a->used  elements currently
//   int check =-1,i;

//   if (start_at>a->used){
//     printf("the max length of the  array is %zu, but you want to check ith element %d\n",a->used,start_at );
//     free(a);
//     error("Not possible");
//   }
//   if((a->used>0) && (a->bondinfos[start_at].ipart==ipart) && (a->bondinfos[start_at].jpart==jpart) ){
//       check=start_at;
//   }     
  

//   return check;

// }


// int checkElementXBondTimeInfoArray(BondTimeInfoArray *a, int ipart , int jpart  ){
//   // check if elementX is in the array
//   // the array a has a->used  elements currently
//   int check =-1,i;

//   for (i=0;i<a->used;i++){
//     if(a->bondinfos[i].ipart==ipart ){
//       if (a->bondinfos[i].jpart==jpart){
//         check=i;
//         break;
//       }
//     }     
//   }

//   return check;

// }

// void printBondTimeInfoArray(BondTimeInfoArray *a ){
//   // check if elementX is in the array
//   // the array a has a->used  elements currently
//   int i;
//   printf(" number of tracked bonds  is %zu \n",a->used);
//   for (i=0;i<a->used;i++){
//     printf(" bond %d connected at %3.2lf has connection between ipart %d and jpart %d \n", i,a->bondinfos[i].start_time,a->bondinfos[i].ipart,a->bondinfos[i].jpart);
     
//   }

//   return ;

// }



