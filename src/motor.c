#include "motor.h"

typedef struct {
    int size;
    motor *motor;
    int *dirichlet;
    femSolver *solver;
    femMesh *meshEntier;
} femMotorProblem;

double perimetre(double *X, double *Y) {
    double res = 0.0;
    res += sqrt( (X[0]-X[1])*(X[0]-X[1]) + (Y[0]-Y[1])*(Y[0]-Y[1]) );
    res += sqrt( (X[2]-X[1])*(X[2]-X[1]) + (Y[2]-Y[1])*(Y[2]-Y[1]) );
    res += sqrt( (X[0]-X[2])*(X[0]-X[2]) + (Y[0]-Y[2])*(Y[0]-Y[2]) );
    return res;
}

double jaco(double *X, double *Y) {
    return (X[1]-X[0])*(Y[2]-Y[0]) - (X[2]-X[0])*(Y[1]-Y[0]);
}

void motorAdaptMesh(motor *theMotor, double delta) {
    motorMesh *theMesh = theMotor->mesh;

    double x,y;
    for(int i = 0; i < theMesh->nNode; ++i){
        if  (theMotor->movingNodes[i] == 1){
            x = theMesh->X[i]*cos(delta) - theMesh->Y[i]*sin(delta);
            y = theMesh->X[i]*sin(delta) + theMesh->Y[i]*cos(delta);
            theMesh->X[i] = x;
            theMesh->Y[i] = y; }}
    theMotor->theta += delta;

    int *elem = theMotor->mesh->elem;
    double *X = theMotor->mesh->X;
    double *Y = theMotor->mesh->Y;
    double Xloc[3], Yloc[3];

    int debutElemDomain = 0;
    for (int i = 0; i < 11; i++) {debutElemDomain += theMotor->mesh->nElemDomain[i];}

    double max, ratio; int maxIdx;
    double angle1, angle2, angle;

    for (int k = debutElemDomain; k < theMotor->mesh->nElem; k++) {
        max = 0;
        maxIdx = elem[3*debutElemDomain];

        Xloc[0] = X[elem[3*k]];
        Xloc[2] = X[elem[3*k+2]];
        Yloc[0] = Y[elem[3*k]];
        Yloc[2] = Y[elem[3*k+2]];

        for (int i = debutElemDomain; i < theMotor->mesh->nElem; i++) {
            for (int j = 0; j < 3; j++) {

                angle1 = fmin(atan2(Yloc[0], Xloc[0]), atan2(Yloc[2], Xloc[2]));
                angle2 = fmax(atan2(Yloc[0], Xloc[0]), atan2(Yloc[2], Xloc[2]));
                angle  = atan2(Y[elem[3*i + j]], X[elem[3*i + j]]);

                if (angle >= angle1-0.1 && angle <= angle2+0.1) {
                    Xloc[1] = X[elem[3*i + j]];
                    Yloc[1] = Y[elem[3*i + j]];

                    ratio = jaco(Xloc, Yloc)/(perimetre(Xloc, Yloc)*perimetre(Xloc, Yloc));
                    if (ratio > max) {
                        max = ratio;
                        maxIdx = elem[3*i + j];
                    }
                }
            }
        }

        elem[3*k + 1] = maxIdx;
    }
}

double motorComputeCouple(motor *theMotor) {
    int *elem = theMotor->mesh->elem;
    double *X = theMotor->mesh->X;
    double *Y = theMotor->mesh->Y;
    double *A = theMotor->a;

    int iDomain = 10; // domaine rotor gap
    int nElemDomain = theMotor->mesh->nElemDomain[iDomain];

    double w[] = {1.0/6.0, 1.0/6.0, 1.0/6.0};
    double xsi[] = {1.0/6.0, 1.0/6.0, 2.0/3.0}, eta[] = {1.0/6.0, 2.0/3.0, 1.0/6.0};
    double phi[3][3] = {{2.0/3.0, 1.0/6.0, 1.0/6.0},
                      {1.0/6.0, 1.0/6.0, 2.0/3.0},
                      {1.0/6.0, 2.0/3.0, 1.0/6.0}};
    double dphidxsi[] = {-1.0, 1.0, 0.0};
    double dphideta[] = {-1.0, 0.0, 1.0};
    double I = 0.0;

    double Xloc[3], Yloc[3], Aloc[3], dphidx[3], dphidy[3];

    int debutElemDomain = 0;
    for (int i = 0; i < iDomain; i++) {debutElemDomain += theMotor->mesh->nElemDomain[i];}

    // calcul de d
    double d1 = fabs(sqrt(X[elem[debutElemDomain*3]]  *X[elem[debutElemDomain*3]]   + Y[elem[debutElemDomain*3]]  *Y[elem[debutElemDomain*3  ]])
                   - sqrt(X[elem[debutElemDomain*3+1]]*X[elem[debutElemDomain*3+1]] + Y[elem[debutElemDomain*3+1]]*Y[elem[debutElemDomain*3+1]]));
    double d2 = fabs(sqrt(X[elem[debutElemDomain*3]]  *X[elem[debutElemDomain*3]]   + Y[elem[debutElemDomain*3]]  *Y[elem[debutElemDomain*3  ]])
                   - sqrt(X[elem[debutElemDomain*3+2]]*X[elem[debutElemDomain*3+2]] + Y[elem[debutElemDomain*3+2]]*Y[elem[debutElemDomain*3+2]]));
    double d3 = fabs(sqrt(X[elem[debutElemDomain*3+2]]*X[elem[debutElemDomain*3+2]] + Y[elem[debutElemDomain*3+2]]*Y[elem[debutElemDomain*3+2]])
                   - sqrt(X[elem[debutElemDomain*3+1]]*X[elem[debutElemDomain*3+1]] + Y[elem[debutElemDomain*3+1]]*Y[elem[debutElemDomain*3+1]]));
    double d = fmax(d1, d2);
    d = fmax(d, d3);

    for (int iElem = debutElemDomain; iElem < debutElemDomain+nElemDomain; iElem++) {
        for (int i = 0; i < 3; i++) {
            Xloc[i] = X[elem[iElem*3 + i]];
            Yloc[i] = Y[elem[iElem*3 + i]];
            Aloc[i] = A[elem[iElem*3 + i]];
        }

        double dxdxsi = 0;
        double dxdeta = 0;
        double dydxsi = 0;
        double dydeta = 0;
        for (int i = 0; i < 3; i++) {
            dxdxsi += Xloc[i]*dphidxsi[i];
            dxdeta += Xloc[i]*dphideta[i];
            dydxsi += Yloc[i]*dphidxsi[i];
            dydeta += Yloc[i]*dphideta[i];
        }
        double jac = dxdxsi * dydeta - dxdeta * dydxsi;

        for (int i = 0; i < 3; i++) {
            dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
            dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
        }

        double dadx = 0.0;
        double dady = 0.0;
        for (int i = 0; i < 3; i++) {
            dadx += Aloc[i] * dphidx[i];
            dady += Aloc[i] * dphidy[i];
        }

        for (int i = 0; i < 3; i++) {
            double Xint = Xloc[0]*phi[i][0] + Xloc[1]*phi[i][1] + Xloc[2]*phi[i][2];
            double Yint = Yloc[0]*phi[i][0] + Yloc[1]*phi[i][1] + Yloc[2]*phi[i][2];
            double angle = atan(Yint/Xint);
            if (Xint < 0.0) {angle += M_PI;}

            double r = sqrt(Xint*Xint + Yint*Yint);
            double dadr  =  cos(angle)*dadx + sin(angle)*dady;
            double dadth = (-sin(angle)*dadx + cos(angle)*dady)*r;

            I += dadr * dadth * jac * w[i];
        }
    }

    printf("C = %f, w = %f, theta = %f\n", -1*I*0.06/(4*M_PI*1e-7*d), theMotor->omega, theMotor->theta);
    return -1*I*0.06/(4*M_PI*10e-7*d);
}

void motorComputeCurrent(motor *theMotor) {
    double offsetCurrent = 0.3;
    double amplitude = theMotor->theta + offsetCurrent;
    while (amplitude < 0.0) {amplitude += M_PI/2;}
    while (amplitude >= M_PI/2) {amplitude -= M_PI/2;}
    double courant = 8.8464*1e5;

    for (int i = 1; i <= 6; i++) {theMotor->js[i] = 0.0;}

    printf("amplitude = %f\n", amplitude);
    /*if (amplitude < M_PI/12) {
        theMotor->js[3] = courant;
        theMotor->js[4] = -courant;
    } else if (amplitude < 2*M_PI/12) {
        theMotor->js[3] = courant;
        theMotor->js[4] = -courant;
        theMotor->js[1] = courant;
        theMotor->js[2] = -courant;
    } else if (amplitude < 3*M_PI/12) {
        theMotor->js[1] = courant;
        theMotor->js[2] = -courant;
    } else if (amplitude < 4*M_PI/12) {
        theMotor->js[1] = courant;
        theMotor->js[2] = -courant;
        theMotor->js[5] = courant;
        theMotor->js[6] = -courant;
    } else if (amplitude < 5*M_PI/12) {
        theMotor->js[5] = courant;
        theMotor->js[6] = -courant;
    } else if (amplitude < 6*M_PI/12) {
        theMotor->js[3] = courant;
        theMotor->js[4] = -courant;
        theMotor->js[5] = courant;
        theMotor->js[6] = -courant;
    }*/
    if (amplitude < M_PI/6) {
        theMotor->js[3] = courant;
        theMotor->js[4] = -courant;
    } else if (amplitude < 2*M_PI/6) {
        theMotor->js[1] = courant;
        theMotor->js[2] = -courant;
    } else if (amplitude < 3*M_PI/6) {
        theMotor->js[5] = courant;
        theMotor->js[6] = -courant;
    }
}

femMotorProblem *femMotorCreate(motor *theMotor, femSolverType solverType) {
    femMotorProblem *problem = malloc(sizeof(femMotorProblem));

    motorMesh *mesh = theMotor->mesh;
    femMesh *meshEntier = malloc(sizeof(femMesh));
    meshEntier->elem = mesh->elem;
    meshEntier->X = mesh->X;
    meshEntier->Y = mesh->Y;
    meshEntier->nNode = mesh->nNode;
    meshEntier->nElem = mesh->nElem;
    meshEntier->nLocalNode = mesh->nLocalNode;
    meshEntier->number = malloc(sizeof(int) * theMotor->size);
    femMeshRenumber(meshEntier, FEM_NO);

    switch (solverType) {
        case FEM_FULL: problem->solver = femSolverFullCreate(theMotor->size, 3); break;
        case FEM_ITER: problem->solver = femSolverIterativeCreate(theMotor->size, 3); break;
        case FEM_BAND: ;
            femMeshRenumber(meshEntier, FEM_XNUM);
            int band = femMeshComputeBand(meshEntier);
            //printf("band size = %d\n", band);
            problem->solver = femSolverBandCreate(theMotor->size, 3, band);
            break;
        default: break;
    }

    problem->dirichlet = malloc(sizeof(int) * theMotor->size);
    for (int i = 0; i < theMotor->size; i++) {problem->dirichlet[i] = 0;}
    femEdges *edges = femEdgesCreate(meshEntier);
    for (int i = 0; i < edges->nEdge; i++) {
        if (edges->edges[i].elem[1] < 0) {
            problem->dirichlet[edges->edges[i].node[0]] = 1;
            problem->dirichlet[edges->edges[i].node[1]] = 1;
        }
    }

    femEdgesFree(edges);

    problem->motor = theMotor;
    problem->meshEntier = meshEntier;
    problem->size = theMotor->size;

    return problem;
}

void femMotorFree(femMotorProblem *theProblem) {
    free(theProblem->meshEntier->number);
    femSolverFree(theProblem->solver);
    free(theProblem->dirichlet);
    free(theProblem->meshEntier);
    free(theProblem);
}

void femMotorMeshLocal(const femMotorProblem *theProblem, const int iElem, int *map, int *ctr, double *X, double *Y, double *U) {
    femMesh *theMesh = theProblem->meshEntier;

    for (int j = 0; j < 3; ++j) {
        map[j] = theMesh->elem[iElem*3+j];
        X[j]   = theMesh->X[map[j]];
        Y[j]   = theMesh->Y[map[j]];
        U[j]   = theProblem->motor->a[map[j]];
        ctr[j] = theProblem->dirichlet[map[j]];
        map[j] = theProblem->meshEntier->number[map[j]];
    }
}

void femMotorCompute(femMotorProblem *theProblem) {
    femSolver *solver = theProblem->solver;
    motorMesh *mesh = theProblem->motor->mesh;
    motor *theMotor = theProblem->motor;
    int *number = theProblem->meshEntier->number;

    double *Aloc = solver->local->A[0];
    double *Bloc = solver->local->B;
    double Uloc[3];
    int ctr[3];

    double xsi = 1.0/3.0, eta = 1.0/3.0;
    double phi[] = {1.0/3.0, 1.0/3.0, 1.0/3.0};
    double dphidxsi[] = {-1.0, 1.0, 0.0};
    double dphideta[] = {-1.0, 0.0, 1.0};

    double X[3], Y[3], dphidx[3], dphidy[3];
    int map[3];

    int nDomain = mesh->nDomain;
    int offsetIndex = 0;

    for (int iDomain = 0; iDomain < nDomain; iDomain++) {
        int nElem = mesh->nElemDomain[iDomain];

        for (int iElem = 0; iElem < nElem; iElem++) {
            femMotorMeshLocal(theProblem, iElem + offsetIndex, map, ctr, X, Y, Uloc);

            double dxdxsi = 0;
            double dxdeta = 0;
            double dydxsi = 0;
            double dydeta = 0;
            for (int i = 0; i < 3; i++) {
                dxdxsi += X[i]*dphidxsi[i];
                dxdeta += X[i]*dphideta[i];
                dydxsi += Y[i]*dphidxsi[i];
                dydeta += Y[i]*dphideta[i];
            }
            double jac = dxdxsi * dydeta - dxdeta * dydxsi;

            for (int i = 0; i < 3; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }
            for (int i = 0; i < 3; i++) {
                for(int j = 0; j < 3; j++) {
                    Aloc[i*3+j] = (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]) * jac / (2 * theMotor->mu[iDomain]);
                }
            }
            for (int i = 0; i < 3; i++) {
                Bloc[i] = phi[i] * jac * theMotor->js[iDomain] / 2;
            }
            for (int i = 0; i < 3; i++) {
                if (ctr[i] == 1) {
                    femFullSystemConstrain(solver->local, i, 0.0);
                }
            }
            femSolverAssemble(solver, Aloc, Bloc, Uloc, map, 3);
        }
        offsetIndex += nElem;
    }

    double *soluce = femSolverEliminate(solver);

    for (int i = 0; i < theProblem->size; i++) {
        theProblem->motor->a[i] += soluce[number[i]];
    }
}

void motorComputeMagneticPotential(motor *theMotor) {
    femSolverType type = FEM_BAND;
    femMotorProblem *theProblem = femMotorCreate(theMotor, type);

    for (int i = 0; i < theMotor->size; i++) {theMotor->a[i] = 0.0;}

    int testConvergence;
    do {
        femMotorCompute(theProblem);
        testConvergence = femSolverConverged(theProblem->solver);
    } while (testConvergence == 0);
    if (testConvergence == -1)  printf("Iterative solver stopped afer a maximum number of iterations\n");

    femMotorFree(theProblem);
}
