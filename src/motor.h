#include "fem.h"

#define MAXSTRING 50
typedef char string[MAXSTRING];

typedef struct {
    int *elem;
    double *X;
    double *Y;
    int nElem;
    int nNode;
    int nLocalNode;
    int nDomain;
    int *nElemDomain;
    string *nameDomain;
    int *domain;
} motorMesh;

typedef struct {
    int size; // = mesh->nNode
    motorMesh *mesh;
    double *a;
    double time;
    double theta;
    double omega;
    int *movingNodes;
    double inertia;
    double L;
    double *js;
    double *mu;
    int nonLinearFlag;
    int nHystereticCurve;
    const double *hystereticCurveH;
    const double *hystereticCurveB;
} motor;


void 	  	    motorMeshWrite(const motorMesh *theMotorMesh, const char *filename);
motorMesh 	   *motorMeshRead(const char *filename);
femMesh        *motorDomainCreate(const motorMesh *theMotorMesh, int iDomain);

motor          *motorCreate(motorMesh *theMotorMesh);
void            motorPrintInfos(const motor *theMotor);
void            motorComputeMagneticPotential(motor *myMotor);
double          motorComputeCouple(motor *myMotor);
void            motorAdaptMesh(motor *myMotor, double delta);
void            motorComputeCurrent(motor *myMotor);
