#ifndef _BASEDEFINE_H_
#define _BASEDEFINE_H_

#ifdef min
#undef min
#endif
#define min(x,y)((x)>(y)?(y):(x))

#ifdef max
#undef max
#endif
#define max(x,y)((x)>(y)?(x):(y))

#ifdef NULL
#undef NULL
#endif
#define NULL 0

#ifdef FREE
#undef FREE
#endif
#define FREE(x) {if(x!=NULL){free(x);x=NULL;}}

#ifdef M_FREE
#undef M_FREE
#endif
#define M_FREE(x) if((x)!=NULL){free(x);x=NULL;}

#define GrIdx(n,m) (((n))*((n)+1)/2+(m)-3)

#define LegIdx(n,m) (((n))*((n)+1)/2+(m))

#define PI   3.1415926535897932384626

#define RADDEG    (PI/180.0)

#endif
