
// random.h

#ifndef RANDOM_H
#define RANDOM_H

/*------------------GLOBAL FUNCTIONS------------------------------------*/
extern double RandomNumber(void);
extern void InitializeRandomNumberGenerator(double);
extern double RandomGaussianNumber();
extern vector RandomBrownianVector(double);
extern vector RandomVector(double);
extern double RandomVelocity(double);
extern vector RandomUnitVector(void);
extern int Random_urandom(void);
extern vector RandomVector3D(vector ) ;
extern double RandomNumberRange(double , double );
extern int RandomIntegerRange(int  , int );
extern vector RandomVector1(void) ;
/*---------------------------------------------------------------------------*/

#endif