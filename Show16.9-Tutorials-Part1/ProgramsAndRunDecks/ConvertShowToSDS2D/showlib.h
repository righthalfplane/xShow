#include <stdio.h>
struct showdata{
    FILE *out;
    long nodeCount;
    long valueCount;
    long normalCount;
    long dxyzCount;
    long rsCount;
    long localSystem;
    long elementCount;
    long nodesMax;
    long stressCount;
    long pc;
    long *elements;
    double *x;
    double *y;
    double *z;
    double *v;
    double *nx;
    double *ny;
    double *nz;
    double *dx;
    double *dy;
    double *dz;
    double *r;
    double *s;
    double *stress;
};
typedef struct showdata *showPtr;

showPtr showStart(char *name,char *mode);
int showEnd(showPtr s);
int showDone(showPtr s);
int showRS(showPtr s,long nodeCount,double *rr,double *ss);
int showDisplacements(showPtr s,long nodeCount,double *dx,double *dy,double *dz);
int showValues(showPtr s,long nodeCount,double *v);
int showNormals(showPtr s,long nodeCount,double *nx,double *ny,double *nz);
int showNodes(showPtr s,long nodeCount,double *x,double *y,double *z);
int showElements(showPtr s,long elementCount,long nodesMax,long *elements);
int showStress(showPtr s,long nodeCount,double *stress);
