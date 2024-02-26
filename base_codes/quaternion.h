// quaternion.h

#ifndef QUATERNION_H
#define QUATERNION_H

/*------------------GLOBAL VARIABLES-----------------------------------------*/
/*---------------------------------------------------------------------------*/

/*------------------STRUCTURE DEFINITIONS------------------------------------*/


typedef struct quattensor_type {

    vector        w,x,y,z;

} quattensor;
/*---------------------------------------------------------------------------*/

/*------------------GLOBAL FUNCTIONS-----------------------------------------*/
extern double langrange_multiplier_quat(quaternion, quaternion);
extern quaternion quatVecToVec(vector, vector);
extern tensor getrotmatrix(quaternion);
extern quattensor getquatmatrix(quaternion);
extern quaternion RandomQuaternion();
extern quaternion RandomQuaternionRange(double);
extern quaternion QuaternionXaxis( double );
extern quaternion QuaternionYaxis( double );
extern quaternion QuaternionZaxis( double );
extern void rotate_quaternion_ipart( Slice *, int , quaternion  );
extern quaternion RotateQuaternion(vector  , double ) ;
extern quaternion normalize_quat(quaternion );
/*---------------------------------------------------------------------------*/

/*------------------LOCAL FUNCTIONS------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*------------------GLOBAL STRUCTURES-----------------------------------------*/
/*---------------------------------------------------------------------------*/




#endif
