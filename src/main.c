
/*
 *
 * ==========================================================================
 *
 *  Simulation par éléments finis d'un motor à reluctance variable
 *  Projet du cours LEPL1110 : ann�e 20-21
 *
 *  Vincent Legat
 *  Michel Henry
 *  François Henrotte
 *  Benjamin Legat
 *
 * ==========================================================================
 */

#include "glfem.h"
#include "motor.h"

int main(void)
{
    motorMesh *theMotorMesh = motorMeshRead("../data/motor838.txt");
    motor *theMotor = motorCreate(theMotorMesh);
    motorPrintInfos(theMotor);

    femMesh *theRotor        = motorDomainCreate(theMotorMesh,8);
    femMesh *theStator       = motorDomainCreate(theMotorMesh,0);
    femMesh *theGap          = motorDomainCreate(theMotorMesh,11);
    femMesh *theCoilPositive = motorDomainCreate(theMotorMesh,1);
    femMesh *theCoilNegative = motorDomainCreate(theMotorMesh,2);
    femMesh *theRotorGap     = motorDomainCreate(theMotorMesh,10);
    femMesh *theStatorGap    = motorDomainCreate(theMotorMesh,7);
    femMesh *theRotorAir     = motorDomainCreate(theMotorMesh,9);

    const char theHelpMessage[] = {
    "   [esc] : Exit\n"
    "    R    : Restart and reset zoom, translations \n"
    "    S    : Show magnetic potential on rotor and stator \n"
    "    K    : Domains with colors \n"
    "    H    : Display or hide keyboard shortcuts \n"};
    printf("\n%s\n",theHelpMessage);
    glfemWindowCreate("EPL1110 : Switched Reluctance Motor",480,480,
            theMotorMesh->nNode,theMotorMesh->X,theMotorMesh->Y);
    glfemWindowSetHelpMessage(theHelpMessage);


    double 	theDiscreteTime = 0.0;
    double 	theStartingTime = 0.0;
    double  theTimeStep  = 0.002;
    double  theStop = 0;
    double  omega = 1.0;
    int     thePlotMode = 1;
    int 	theIteration = 0;
    char    theAction = 'K';
    char    theMessage[256];

    FILE *fp = fopen("../data/couples.txt", "w");

    do
    {
        glfemWindowUpdate();
        char action = glfemGetAction();

//
//  Gestion de l'animation en temps reel
//    - "theTime" est le temps courant de l'application (sauf si il a été remis a zero :-)
//    - Le facteur entre le temps réel et le temps de l'application permet de controler la vitesse d'execution
//    - Si necessaire, une nouvelle iteration discrete est calculee...
//      C'est donc ici que se trouve "virtuellement" la boucle sur toutes les iterations temporelles
//
//  Pour figer/ne pas figer le resultat a un temps, decommenter les deux lignes ci-dessous
//
        double theTime = (glfwGetTime() - theStartingTime) * 5;

        if (theTime >= theStop) theTime = theStop;


        if (action == 'K') thePlotMode = 0;
        if (action == 'S') thePlotMode = 1;

        if (action == 'R') {
            theStartingTime = glfwGetTime();
            theDiscreteTime = 0;
            theStop = 1;}

        glfemSetColorLine(GLFEM_BLACK);
        glfemSetColor(GLFEM_BACKGROUND);
        glfemSetLineWidth(0.0001);

//
// Pour faire apparaitre les sous-domaines : c'est ici qu'il est possible
// de modifier l'affichage :-)
//

        if (thePlotMode == 0) {
            glfemPlotMesh((femMesh*)theMotorMesh);
            glfemSetColor(GLFEM_BLUE);
            glfemPlotMesh(theRotor);
            glfemPlotMesh(theStator);
            glfemSetColor(GLFEM_GREEN);
            glfemPlotMesh(theGap);
            glfemSetColor(GLFEM_RED);
            glfemPlotMesh(theRotorGap);

            glfemSetColor(GLFEM_RED);
            for (int i = 1; i <= 6; i++) {
                if (theMotor->js[i] < 0.0) {
                    glfemPlotMesh(motorDomainCreate(theMotorMesh,i));
                }
            }
            glfemSetColor(GLFEM_BLACK);
            for (int i = 1; i <= 6; i++) {
                if (theMotor->js[i] > 0.0) {
                    glfemPlotMesh(motorDomainCreate(theMotorMesh,i));
                }
            }
        }

//
// Visualisation du potentiel magn�tique sur le rotor et le stator
//

        if (thePlotMode == 1) {
            glfemSetScale((femMesh*)theMotorMesh,theMotor->a);
            glfemPlotSolution(theRotor,theMotor->a);
            glfemPlotSolution(theStator,theMotor->a);
        }

//
// Calcul d'une itération temporelle
// C'est ici qu'on exécute le projet :-)
//
//   Calcul du potentiel magnétique A
//   Calcul du couple
//   Calcul de omega par l'équation de Newton
//   Rotation du rotor et remaillage de la bande glissante
//   Mise a jour des courants dans les inducteurs en fonction de l'angle
//

        if (theTime > theDiscreteTime) {
            theIteration += 1;
            theDiscreteTime += theTimeStep;
            motorComputeMagneticPotential(theMotor);
            double C = motorComputeCouple(theMotor);
            fprintf(fp, "%.10f %.10f %.10f %.10f\n", theDiscreteTime, theMotor->theta, C, theMotor->omega);
            // couple résistif
            //C -= 0.0002*theMotor->omega;
            omega += C * theTimeStep / theMotor->inertia;
            theMotor->omega = omega;
            motorAdaptMesh(theMotor,omega*theTimeStep); // changement de theta
            motorComputeCurrent(theMotor);
            printf("Iteration  %2d - %.2f : %14.7e \n",theIteration,theDiscreteTime,theMotor->theta);
        }
        sprintf(theMessage,"Time = %.2f iteration = %d",theDiscreteTime,theIteration);


        sprintf(theMessage, "Entrefer ");
        glfemDrawMessage(theMessage,(double[2]){16.0, 30.0});

    } while(!glfemWindowShouldClose());

    fclose(fp);

    exit(EXIT_SUCCESS);
    return 0;


}



// =========== Quelques fonctions fournies gracieusement par l'�quipe didactique ====== //
//
// Attention : il n'est pas permis de les modifier :-)
// Attention : la structure de données du problème est figée et vous ne pouvez pas
//             la modifier : c'est cette structure que le correcteur automatique
//             utilisera pour tester votre programme
//
// =====================================================================================//

static const double _hystereticCurveH[43] = { 0, 10, 20, 30, 40, 50,
      60, 70, 80, 90, 100, 125, 150, 175, 200, 250,
      300, 400, 500, 600,  700, 800, 900, 1000, 1250, 1500, 2000, 2500, 5000,
      7500,  10000, 15000, 20000, 59000, 174000, 514000, 1520000, 4470000,
      13200000, 38900000, 115000000, 339000000, 1000000000 };
static const double _hystereticCurveB[43] =  { 0.0,
      0.194880829963, 0.377143018857, 0.537767739762, 0.672888260835,
      0.783043000477, 0.871342430831,0.941778611986, 0.998183303557, 1.04378111223,
      1.08110469369, 1.14963767549, 1.19607212343, 1.22964695907, 1.25515221835,
      1.29162498935, 1.31678879432, 1.35015120537, 1.37220092877, 1.38859114656,
      1.4017440574, 1.41287024565, 1.42264180514, 1.43146158921, 1.45082466146,
      1.46784549989, 1.49819370601, 1.52578650709, 1.64314027719, 1.73458485332,
      1.8039068939,1.89568786291, 1.95213815187, 2.1390774927, 2.45827909293,
      3.32303272825, 5.85485500678, 13.2701832298, 35.2114648741, 99.8027446541,
      291.062951228, 854.036370229, 2515.3105707 };



motor *motorCreate(motorMesh *theMesh)
{
    motor *theMotor = malloc(sizeof(motor));
    theMotor->mesh = theMesh;
    theMotor->size = theMesh->nNode;

//
//  Identification des noeuds mobiles
//

    theMotor->movingNodes = malloc(sizeof(int)*theMesh->nNode);
    for (int i=0; i < theMotor->size; i++)
        theMotor->movingNodes[i] = 0;
    for (int i=0; i < theMesh->nElem; i++) {
        int domain = theMesh->domain[i];
        if (domain == 8 || domain == 9 || domain == 10 ) {
            int *elem = &(theMesh->elem[i*3]);
            for (int j=0; j < 3; j++) {
                theMotor->movingNodes[elem[j]] = 1; }}}
//
//  Initialisation des inconnues � la valeur X pour voir le maillage tourner
//
    theMotor->a = malloc(sizeof(double)*theMotor->size);
    for (int i=0; i < theMotor->size; i++) {
        theMotor->a[i] = theMesh->X[i];
    }

    theMotor->theta = 0;
    theMotor->omega = 0;
    theMotor->time  = 0;

//
//  Parametres materiels
//

    double mu_0 = 4*3.141592653589793*1e-7; //kg m /(A**2 s**2)
    double mu_r = 1e3;
    double js   = 8.8464*1e5;               // A / m**2
    theMotor->inertia = 5*1e-4;             // kg m**2
    theMotor->L       = 0.06;               // m
    theMotor->js = malloc(sizeof(double)*theMesh->nDomain);
    theMotor->mu = malloc(sizeof(double)*theMesh->nDomain);
    for(int i = 0; i < theMesh->nDomain; i++) {
        theMotor->js[i] = 0.0;
        theMotor->mu[i] = mu_0; }
    theMotor->mu[0] *= mu_r;
    theMotor->mu[8] *= mu_r;
    theMotor->js[1] = js;
    theMotor->js[2] = -js;

//
//  Bonus : magnetostatique non-lineaire
//

    theMotor->nonLinearFlag = 0;
    theMotor->nHystereticCurve = 43;
    theMotor->hystereticCurveH = _hystereticCurveH;
    theMotor->hystereticCurveB = _hystereticCurveB;



    return theMotor;
}

void motorPrintInfos(const motor *theMotor)
{
    int  size = theMotor->size;
    motorMesh *theMesh = theMotor->mesh;
    printf(" \n");
    printf(" ====== Switched Reluctance Motor Simulation ============================\n");
    for(int i = 0; i < theMesh->nDomain; i++)
        printf("    Domain %2d : %-16s : nElem = %6d, mu = %10.3e, js = %10.3e \n",i,theMesh->nameDomain[i],
                    theMesh->nElemDomain[i],theMotor->mu[i],theMotor->js[i]);
    printf("                                 : mu permeability [kg m s2/A2] - js current density [A/m2]\n");
    printf("    Number of elements           : %d\n",theMesh->nElem);
    printf("    Number of nodes              : %d\n",theMotor->size);
    printf("    Flag for non linearities     : %d\n",theMotor->nonLinearFlag);
    printf("    Rotor inertia                : %13.7e [kg m2]\n",theMotor->inertia);
    printf("    Motor axial length           : %13.7e [m]\n",theMotor->L);
    printf("    Time                         : %13.7e [s]\n",theMotor->time);
    printf("    Angular position             : %13.7e [rad]\n",theMotor->theta);
    printf("    Angular velocity             : %13.7e [rad/s]\n",theMotor->omega);
    printf("=========================================================================\n");


}

motorMesh *motorMeshRead(const char *filename)
{
    motorMesh *theMesh = malloc(sizeof(motorMesh));
    theMesh->nLocalNode = 3;

    FILE* file = fopen(filename,"r");
    if (file == NULL) Error("No mesh file !");
    int i,j,trash,*elem;

    ErrorScan(fscanf(file, "Number of nodes %d\n", &theMesh->nNode));
    theMesh->X = malloc(sizeof(double)*theMesh->nNode);
    theMesh->Y = malloc(sizeof(double)*theMesh->nNode);
    for (i = 0; i < theMesh->nNode; i++)
         ErrorScan(fscanf(file,"%d : %le %le\n",&trash,&theMesh->X[i],&theMesh->Y[i]));

    ErrorScan(fscanf(file, "Number of sub-domains %d \n", &theMesh->nDomain));
    theMesh->nameDomain = malloc(sizeof(string)*theMesh->nDomain);
    theMesh->nElemDomain = malloc(sizeof(int)*theMesh->nDomain);
    for(i = 0; i < theMesh->nDomain; i++)
     	  ErrorScan(fscanf(file, "%6d : %s : %d\n",&trash,theMesh->nameDomain[i],&theMesh->nElemDomain[i]));

    ErrorScan(fscanf(file, "Number of triangles %d  \n", &theMesh->nElem));
    theMesh->elem = malloc(sizeof(int)*3*theMesh->nElem);
    theMesh->domain = malloc(sizeof(int)*theMesh->nElem);
    for (i = 0; i < theMesh->nElem; i++) {
    	  elem = &(theMesh->elem[i*3]);
    	  ErrorScan(fscanf(file,"%d : %d %d %d %d\n",&trash,&elem[0],&elem[1],&elem[2],&theMesh->domain[i])); }

    fclose(file);
    return theMesh;
}

void motorMeshWrite(const motorMesh *theMesh, const char *filename)
{
    int i,j,*elem;

    FILE* file = fopen(filename,"w");
    fprintf(file, "Number of nodes %d \n", theMesh->nNode);
    for (i = 0; i < theMesh->nNode; i++) {
        fprintf(file,"%6d : %14.7e %14.7e \n",i,theMesh->X[i],theMesh->Y[i]); }

    fprintf(file, "Number of sub-domains %d \n", theMesh->nDomain);
    for (i = 0; i < theMesh->nDomain; i++) {
     	  fprintf(file, "%6d : %-16s : %d \n",i,theMesh->nameDomain[i],theMesh->nElemDomain[i]);}

    fprintf(file, "Number of triangles %d \n", theMesh->nElem);
    for (i = 0; i < theMesh->nElem; i++) {
    	  elem = &(theMesh->elem[i*3]);
       	fprintf(file,"%6d : %6d %6d %6d %6d \n",i,elem[0],elem[1],elem[2],theMesh->domain[i]); }

    fclose(file);
}

femMesh *motorDomainCreate(const motorMesh *theMotorMesh, int iDomain)
{
    femMesh *theMesh = malloc(sizeof(femMesh));
    theMesh->nLocalNode = 3;
    theMesh->nNode = theMotorMesh->nNode;
    theMesh->X = theMotorMesh->X;
    theMesh->Y = theMotorMesh->Y;

    int shift = 0;
    for (int i=0; i < iDomain; i++)
      shift += theMotorMesh->nElemDomain[i];
    theMesh->elem = &theMotorMesh->elem[3*shift];
    theMesh->nElem = theMotorMesh->nElemDomain[iDomain];
    return theMesh;
}
