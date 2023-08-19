#include <stdio.h>
#include <math.h>
#include <stdlib.h>					
#include <string.h>
#include <float.h>
#include <mpi.h>
/*---------------------------- CONSTANTS IN THE SYSTEM ------------------------------------------------*/
#define	N					10000
#define	DOUBLE_N				10000.0
#define NA					6500
#define	CUTOFF_SQRD				(2.50*2.50)
#define	INV_CUTOFF_SQRD				(1.0/CUTOFF_SQRD)
#define	DR					0.35
#define	MAX_NEBZ				512
#define	MAX_CELLS				(N/2) 
#define	MAX_BUDS_PER_CELL			1000
#define	TIME_STEP				0.005
#define	HALF_TIME_STEP				0.5*TIME_STEP
#define	LENGTH					1.0
#define	HALF_LENGTH				0.5
#define DENSITY 				1.20
#define TAU_T					5.0
#define	TWO_PI					6.28318530717958
#define LARGE                          		1
#define SMALL                           	0
#define DIVIDE_TO                               2
#define IMPURITY_DENSITY			0.0
#define NUMBER_IMPURITY				(int)(DOUBLE_N*IMPURITY_DENSITY)
/*-------------------- MINIMIZER SHIT ---------------------------------------------------------------------*/
#define	MAX_ITERATIONS				500000
#define	TOL					1e-18 		
#define	LINMIN_MAX_ITERATIONS			40
#define	LINMIN_G1				2.0 // factors for growing and shrinking the interval- don't change
#define	LINMIN_G2				1.25
#define	LINMIN_G3				0.5
#define	LAST_X_DEFAULT				0.0001
/*------------------------------------ BORING STUFF FROM NUMERICAL RECIPES -----------------------------------*/
#define	IA					16807
#define	IM					2147483647
#define	AM					(1.0/IM)
#define	IQ					127773
#define	IR					2836
#define	NTAB					32
#define	NDIV					(1+(IM-1)/NTAB)
#define	RNMX					(1.0-DBL_EPSILON)
#define	ITMAX					200
#define	EPS					3.0e-7
#define	FPMIN					1.0e-30
#define	UNI					((double)rand()/((double)RAND_MAX + 1.0)) 
#define SWAP(a,b) 				temp=(a);(a)=(b);(b)=temp;
#define M 					7
#define NSTACK 					50
/*-------------------------------------------------------------------------------------------------------*/
void calculateForces();
void initializeSystem();
void initializeSystemRestart();		// Here I have added this
void updateNebzLists();
int min();

/*-------------------------- GLOBALS VARIABLES -----------------------------------------------------------*/
double L; 			/*	length of <square> box 		*/
double invL;			/*	inverse of length		*/
double V;			/*	volume				*/
double T;			/*	temperature			*/
double inst_T;			/* 	instantaneous temperature	*/
double beta; 			/*  	inverse temperature 		*/
double t;			/*	time				*/
double u;			/*	potential energy		*/
double particleU;		/* 	energy of single particle	*/
double kinetic;			/*	kinetic energy			*/
double strain;			/*	strain of the system		*/
double rx[N];			/*	x component of position		*/
double ry[N];			/*	y component of position		*/
double rxUnFolded[N];		/*	x component of position		*/
double ryUnFolded[N];		/*	y component of position		*/
double px[N];			/*	x component of momentum		*/
double py[N];			/*	y component of momentum		*/
double fx[N];			/*	x component of force		*/
double fy[N];			/*	y component of force		*/
int type[N]; 			/*	type large or small		*/
double leDisplacement;  	/* 	LEE-EDWARDS displacement  	*/

int impurityParticleIndex[NUMBER_IMPURITY];
double rxImpurity[NUMBER_IMPURITY];
double ryImpurity[NUMBER_IMPURITY];

/*-------------- restart job options -------------------------------------------------------------------*/
int restartEq,restartPostEq;
int restartRun;
int iterFinish;
/*---------------------------- PARAMETERS FOR THE SYSTEM -----------------------------------------------*///change these
const double invSizeSqrd[2][2] = {{1.00,1.56250},{1.56250,1.29132}};
// [0][0] = large-large; [0][1] = [1][0] = large-small; [1][1] = small-small;
double eps[2][2] = {{1.0,1.50},{1.50,0.50}}; 
double eps4[2][2] ={{4.00,6.00},{6.00,2.00}} ;
double vc1[2][2] = {{-0.001949973872640, -0.003046834176000},{-0.003046834176000, -0.002518044773554}};
double vc2[2][2] = {{0.016266559488000, 0.016266559488000},{0.016266559488000, 0.016266559488000}};
double fc[2][2] = {{24.000000000000000, 56.250000000000000},{56.250000000000000, 15.495867768595042}};
double f1[2][2] = {{-0.000005368709120, -0.000005368709120},{-0.000005368709120, -0.000005368709120}};
double f2[2][2] = {{0.000655360000000, 0.000655360000000},{0.000655360000000, 0.000655360000000}};	
/*------------------------------------- FOR NEBZ LIST -----------------------------------------------------*/  //old values //change it
int nebz[N][MAX_NEBZ]; 		//maximum we give here MAX_NEBZ nebz. 
int numOfNebz[N]; 		//number of elements in nebz[N][MAX_NEBZ]
double maxD,listCutOffSqrd; 	//for updating nebz list.
int nebListCounter; 		//for counting how many steps have took place for nebz list updating.

/*---------------------------------------------------------------------------------------------------------------*/ //new one
int nebz[N][MAX_NEBZ]; //maximum we give here MAX_NEBZ nebz. 
int numOfNebz[N]; //number of elements in nebz[N][MAX_NEBZ]
double cellCutOff,maxD,listCutOffSqrd; //for updating nebz list.
int nebListCounter; //for counting how many steps have took place for nebz list updating.
int numInCell[MAX_CELLS];
int cells[MAX_CELLS][MAX_NEBZ];
const int shift[4][2] = {{-1,0},{-1,1},{0,1},{1,1}}; //{yShift,xShift};
/*--------------------- FOR LOGARITHAMIC STORAGE -----------------------------------------------------------*/
unsigned long long int *nextStepIndex;
int serial;
int for_ran1;
/*------------------------------ FOR MINIMIZER -------------------------------------------------------------*/
int restart;		/* whether to restart optimizer - fresh cg directions */
double lastx;		/* keeps track of typical step length */
double gtyp;		/* stores the rms gradient for linmin */
double p[3*N],q[3*N],g[3*N],h[3*N],mosesX[3*N],mosesY[3*N],scratch[3*N];
int NUM_TASK, TASK_ID, jobNo;
/**-------------------------------------------- MAIN PROGRAM ----------------------------------------------**/
//initialize with a negative integer. then just give ONE positive one. 
double ran1(int *idum){
	int j;
	int k;
	static int iy=0;
	static int iv[NTAB];
	double temp;

	if ( *idum <= 0 || !iy) { 
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) { 
			k = *idum/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=*idum/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k; 
	if (*idum < 0) *idum += IM; 
	j=iy/NDIV; 
	iy=iv[j]; 
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

void shellSort(int n, unsigned long long int *a){
	/*Sorts an array a[1..n] into ascending numerical order by Shell's method (diminishing increment
	  sort). n is input; a is replaced on output by its sorted rearrangement.*/
	int i,j,inc;
	unsigned long long int v;

	inc=1; //Determine the starting increment.
	do {
		inc *= 3;
		inc++;
	} while (inc <= n);
	do { //Loop over the partial sorts.
		inc /= 3;
		for (i=inc+1;i<=n;i++) { //Outer loop of straight insertion.
			v=a[i];
			j=i;
			while (a[j-inc] > v) { //Inner loop of straight insertion.
				a[j]=a[j-inc];
				j -= inc;
				if (j <= inc) break;
			}
			a[j]=v;
		}
	} while (inc > 1);

	return;
}
/*----------------------------------------------------------------------------------------------------------*/
int intSqrt(int number){
        int k=1;
        while ( k*k < number )
                k++;
        return k;
}

void initializeGrid(){
        int i,sqrtN,j,temp;
        double space;

	for (i=0;i<N;i++){
                if(i<NA) type[i] = 0;
		else type[i] = 1;
        }
        //permutate indices
        for (i=0; i<N; i++){
                j = (int)(DOUBLE_N*UNI);
                temp = type[i];
                type[i] = type[j];
                type[j] = temp;
        }

	V = DOUBLE_N/DENSITY;
        sqrtN = intSqrt(N);
        L = sqrt(V);
	printf("Simulation box size L = %.10f\n",L);
        space = 1.0/(double)sqrtN;
        strain = 0.0;
        for (i=0; i<N; i++){
                rx[i] = space*(double)(i%sqrtN);
                ry[i] = space*(double)(i/sqrtN);
                px[i] = 2.0*ran1(&for_ran1)-1.0;
                py[i] = 2.0*ran1(&for_ran1)-1.0;
        }
        initializeSystem();
        return;
}



void saveSnapShot(char *fileName){			
	int i;
	FILE *file;
	file = fopen(fileName,"wb");
	fprintf(file,"%f\t%.10f\t%.10f\t%.10f\n",L,u,T,1.1);
	for (i=0; i<N; i++){
		fprintf(file,"%.8f\t%.8f\t%d\n",rx[i],ry[i],type[i]);		// Changed
	}
	fclose(file);

	return;
}

void saveCompleteState(char *outFileName){
	int i;
	FILE *file;

	file = fopen(outFileName,"wb");
	fprintf(file,"%.15g\t%.8g\t%.10g\t%.10f\t%.10f\t1.1\t1.1\n",L,u,T,1.1,kinetic); 
	for (i=0; i<N; i++){
		fprintf(file,"%.8g\t%.8g\t%d\t%.8g\t%.8g\n",rx[i],ry[i],type[i],px[i],py[i]);	// Changed
	}
	fclose(file);

	return;
}

void saveCompleteRestartState(char *outFileName, int iter){
	int i;
	FILE *file;

	file = fopen(outFileName,"wb");
	fprintf(file,"%.15g\t%.8g\t%.10g\t%d\t%.10f\t1.1\t1.1\n",L,u,T,iter,kinetic); 
	for (i=0; i<N; i++){
		fprintf(file,"%.8g\t%.8g\t%d\t%.8g\t%.8g\n",rx[i],ry[i],type[i],px[i],py[i]);	// Changed
	}
	fclose(file);

	return;
}

void readSnapShot(char *fileName){
	int i;
	double dummy;
	FILE *readFile;

	readFile = fopen(fileName, "rb");					// Change here
	fscanf(readFile, "%lf",&(L));
	fscanf(readFile, "%lf",&(u));
	fscanf(readFile, "%lf",&(T));
	fscanf(readFile, "%lf",&(dummy));
	for (i=0;i<N;i++){
		fscanf(readFile, "%lf",&(rx[i]));
		fscanf(readFile, "%lf",&(ry[i]));
		fscanf(readFile, "%d",&(type[i]));
		px[i] = 0.0;	
		py[i] = 0.0;	
	}
	fclose(readFile);
	initializeSystem();

	return;
}

void readCompleteState(char *inFileName){
	int i;
	double dummy;
	FILE *file;

	strain = 0.0;//this we need for the minimizer as it is defined for general case
	file = fopen(inFileName,"rb");
	fscanf(file, "%lf",&(L));
	fscanf(file, "%lf",&(u));
	fscanf(file, "%lf",&(T));
	fscanf(file, "%lf",&(dummy));
	for (i=0;i<3;i++){				// Change here
		fscanf(file,"%lf",&dummy);
	}
	for (i=0;i<N;i++){
		fscanf(file, "%lf",&(rx[i]));
		fscanf(file, "%lf",&(ry[i]));
		fscanf(file, "%d",&(type[i]));
		fscanf(file, "%lf",&(px[i]));
		fscanf(file, "%lf",&(py[i]));
	}
	fclose(file);
	initializeSystemRestart();

	return;
}

void readCompleteRestartState(char *inFileName){
	int i;
	double dummy;
	FILE *file;

	strain = 0.0;//this we need for the minimizer as it is defined for general case
	file = fopen(inFileName,"rb");
	fscanf(file, "%lf",&(L));
	fscanf(file, "%lf",&(u));
	fscanf(file, "%lf",&(T));
	fscanf(file, "%d",&(iterFinish));
	printf("%g\t%g\t%g\n",L,u,T);
	for (i=0;i<3;i++){
		fscanf(file,"%lf",&dummy);
	}
	for (i=0;i<N;i++){
		fscanf(file, "%lf",&(rx[i]));
		fscanf(file, "%lf",&(ry[i]));
		fscanf(file, "%d",&(type[i]));
		fscanf(file, "%lf",&(px[i]));
		fscanf(file, "%lf",&(py[i]));
		printf("%g\t%g\t%d\t%g\t%g\n",rx[i],ry[i],type[i],px[i],py[i]);		// Changed
	}
	fclose(file);
	initializeSystem();

	return;
}

void fixDrift(){							// Changed
	double Px,Py;
	int i;
	Px = 0.0; Py = 0.0; 
	for (i=0;i<N;i++){
		Px += px[i]; Py += py[i];
	}
	Px = Px/DOUBLE_N; Py = Py/DOUBLE_N;
	for (i=0;i<N;i++){
		px[i] -= Px; py[i] -= Py; 
	}

	return;
}

void applyPBC(){
        int i;

        for(i=0;i<N;i++){
                if (rx[i] >= LENGTH){
                        while(rx[i] >= LENGTH)
                                rx[i] = rx[i] - LENGTH;
                }
                else if (rx[i] < 0.0){
                        while(rx[i] < 0.0)
                                rx[i] = rx[i] + LENGTH;
                }

                if (ry[i] >= LENGTH){
                        while(ry[i] >= LENGTH)
                                ry[i] = ry[i] - LENGTH;
                }
                else if (ry[i] < 0.0){
                        while(ry[i] < 0.0)
                                ry[i] = ry[i] + LENGTH;
                }


        }
        return;
}

void adjustStartingTemperature(){		// Changed
	int i;
	double xi_T;

	kinetic = 0.0;
	for(i=0;i<N;i++){
		kinetic += px[i]*px[i] + py[i]*py[i] ;	// Changed
	} 
	kinetic = 0.50*kinetic;
	inst_T = (2.0*kinetic)/(2.0*DOUBLE_N - 3.0);		//// Change here
	xi_T = sqrt(T/inst_T);

	for(i=0;i<N;i++){
		px[i] *= xi_T;
		py[i] *= xi_T;
	}
	return;
}
//initiallize system - stuff for cell subdivision and nebz list.
void initializeSystem(){						// Change here ...
	double rcutSq,sr2,sr6,sr12;
	int i,j;

	invL = 1.0/L;		// Change here
	V = L*L;
	cellCutOff = 1.0*sqrt(CUTOFF_SQRD) + DR;
	listCutOffSqrd = cellCutOff*cellCutOff;
	strain = 0.0; // notice that the strain set to be zero as the minimizer is for general condition
	leDisplacement = fmod(strain,1.0);
	if (leDisplacement < -0.5)
		leDisplacement += 1.0;
	else if (leDisplacement > 0.5)
		leDisplacement -= 1.0;
	
	for(i=0;i<N;i++){
		rxUnFolded[i] = rx[i];
		ryUnFolded[i] = ry[i];
	}
	applyPBC();
	updateNebzLists();
	calculateForces();
	nebListCounter = 0;
	maxD = 0.0;
	fixDrift();
	adjustStartingTemperature();
	return;
}
//initiallize system - stuff for cell subdivision and nebz list.
void initializeSystemRestart(){
	double rcutSq,sr2,sr6,sr12;
	int i,j;

	invL = 1.0/L;
	V = L*L;							// Change here ...
	cellCutOff = 1.0*sqrt(CUTOFF_SQRD) + DR;
	listCutOffSqrd = cellCutOff*cellCutOff;
	strain = 0.0; // notice that the strain set to be zero as the minimizer is for general condition
	leDisplacement = fmod(strain,1.0);
	if (leDisplacement < -0.5)
		leDisplacement += 1.0;
	else if (leDisplacement > 0.5)
		leDisplacement -= 1.0;

	//applyPBC();
	updateNebzLists();
	calculateForces();
	nebListCounter = 0;
	maxD = 0.0;
	//fixDrift();

	for(i=0;i<N;i++){
		rxUnFolded[i] = rx[i];
		ryUnFolded[i] = ry[i];
	}	
	//adjustStartingTemperature();			// Upto here its is changed 

	return;
}	
/*********************************** FILE OUTPUT FUNCTIONS ************************************************************/
void writeState(char *outCoordFileName, char *outStuffFileName){
	int k;
	FILE *outCoordFile, *outStuffFile;

	outCoordFile = fopen(outCoordFileName,"ab");//append...
	outStuffFile = fopen(outStuffFileName,"ab");//append...
	fprintf(outStuffFile,"%.14g\t%.14g\t%.14g\t%g\n",t,u,inst_T,kinetic);
	for (k=0; k<N; k++){
		//fprintf(outCoordFile,"%g\t%g\t%g\t%d\n",rx[k],ry[k],rz[k],type[k]);
		fprintf(outCoordFile,"%g\t%g\t%d\n",rxUnFolded[k],ryUnFolded[k],type[k]);
	}
	fclose(outCoordFile);
	fclose(outStuffFile);

	return;
}

void writeEnd(char *outFileName){
	FILE *outFile;
	outFile = fopen(outFileName,"ab");
	fprintf(outFile,"-1.1\t-1.1\t-1.1\n-1.1\n"); // mark end of input no momentum
	fclose(outFile);
	return;
}

void updateNebzLists(){					//replace it for 2d
	int i,j,x,y,current,m;
	int a,b,c,k,l,numHere,numThere,target,w,xIndex,yIndex;
	double dx,dy,r2,invCellSize;
	double rxi,ryi;
	int mSqrd;
	int ii,jj;

	double includeSqrd;
	nebListCounter = 0;
	maxD = 0.0;

	if ( N < 512 ){		//change here 
		//set all neb lists to zero
		for (i=0; i<N; i++)
			numOfNebz[i] = 0;

		for (i=0;i<N-1;i++){
			rxi = rx[i];
			ryi = ry[i];
			for (j=i+1; j<N; j++){
				dx = rx[j] - rxi;
				dy = ry[j] - ryi;
				if ( dx >= HALF_LENGTH )
					dx -= LENGTH;
				else if ( dx < -HALF_LENGTH )
					dx += LENGTH;
				if ( dy >= HALF_LENGTH )
					dy -= LENGTH;
				else if ( dy < -HALF_LENGTH )
					dy += LENGTH;
				// end of mess 
				r2 = L*L*( dx*dx + dy*dy );
				if (r2 < listCutOffSqrd){
					nebz[i][numOfNebz[i]] = j;
					nebz[j][numOfNebz[j]] = i;
					numOfNebz[i]++;
					numOfNebz[j]++;
				}
			}
		}
	}
	//this uses the cell subdivision.

	
	else{
	m =(int)(L/cellCutOff);
	mSqrd = m*m; 		//the length of cells and numInCell...
	invCellSize = (double)m;
	
	//set number in each cell to zero
	for (i=0; i<mSqrd; i++)
		numInCell[i] = 0;
	
	/*********cells sorting********/
	for (i=0;i<N;i++){
		current = m*(int)(invCellSize*ry[i]) + (int)(invCellSize*rx[i]);
		cells[current][numInCell[current]] = i;
		numInCell[current]++;
	}

	//set all neb lists to zero
	for (i=0; i<N; i++)
		numOfNebz[i] = 0;
	
	//here's where it can get messy when sheering. maybe deal with that later.
	for (i=0;i<m;i++){
		for (j=0;j<m;j++){
			current = m*i + j;
			numHere = numInCell[current];
			//first check interactions within a cell
			for (k=0;k<numHere-1;k++){
				for (l=k+1;l<numHere;l++){
					ii = cells[current][k];
					jj = cells[current][l];
					dx = rx[jj] - rx[ii];
					dy = ry[jj] - ry[ii];
					r2 = V*( dx*dx + dy*dy );
					if (r2 < listCutOffSqrd){
						nebz[ii][numOfNebz[ii]] = jj;
						nebz[jj][numOfNebz[jj]] = ii;
						numOfNebz[ii]++;
						numOfNebz[jj]++;
					}
				}
			}
			
			//now there are 4 nebz that need to be checked.
			for (w=0;w<4;w++){ //w is an index... sorry.
				a = i + shift[w][0];
				if (a<0) a = m-1;
				else if (a==m) a = 0;
				b = j + shift[w][1];
				if (b==m) b = 0;
				target = m*a + b;
				numThere = numInCell[target];
				for (k=0;k<numHere;k++){
					for (l=0;l<numThere;l++){
						ii = cells[current][k];
						jj = cells[target][l];
						dx = rx[jj] - rx[ii];
						dy = ry[jj] - ry[ii];
						// mess due to periodic boundary conditions 
						if ( dx >= HALF_LENGTH )
							dx -= LENGTH;
						else if ( dx < -HALF_LENGTH )
							dx += LENGTH;
						if ( dy >= HALF_LENGTH )
							dy -= LENGTH;
						else if ( dy < -HALF_LENGTH )
							dy += LENGTH;
						// end of mess 
						r2 = V*( dx*dx + dy*dy );
						if (r2 < listCutOffSqrd){
							nebz[ii][numOfNebz[ii]] = jj;
							nebz[jj][numOfNebz[jj]] = ii;
							numOfNebz[ii]++;
							numOfNebz[jj]++;
						}
					}
				}
			}
		}
	}
	}
	return;
}


/*-------------------------------------------------------------------------------------------------------------------------*/

							//we can use this
void calculateForces(){
	int i,j,m,k,typeI,typeJ;
	double g,dx,dy,dz,rxi,ryi,rzi,r2OverSigma2,invSigmaSqrd;
	double sigmaOverR12,temp;
	double sr2,sr6,sr8,rijSq;

	//first set 'em all to zero
	u = 0.0;
	for (i=0;i<N;i++){
		fx[i] = 0.0;
		fy[i] = 0.0;

	}

	for (i=0;i<N-1;i++){
		rxi = rx[i];
		ryi = ry[i];
		typeI = type[i];
		m = numOfNebz[i];
		for (k=0;k<m;k++){
			j = nebz[i][k];
			if (j > i){
				dx = rx[j] - rxi;
				dy = ry[j] - ryi;
				// mess due to LEES EDWARDS periodic boundary conditions 
	
				if ( dx >= HALF_LENGTH )
					dx -= LENGTH;
				else if ( dx < -HALF_LENGTH )
					dx += LENGTH;
				if ( dy >= HALF_LENGTH )
					dy -= LENGTH;
				else if ( dy < -HALF_LENGTH )
					dy += LENGTH;
				// end of mess 
				typeJ = type[j];
				invSigmaSqrd = invSizeSqrd[typeI][typeJ];
				rijSq = L*L*( dx*dx + dy*dy);
				r2OverSigma2 = L*L*( dx*dx + dy*dy)*invSigmaSqrd;
				if (r2OverSigma2 < CUTOFF_SQRD){
					sr2=1.0/r2OverSigma2;
                     			sr6=sr2*sr2*sr2;
                     			sr8=sr6*sr2;
                     			u += eps4[typeI][typeJ]*(sr6*(sr6-1.0)+vc1[typeI][typeJ]*rijSq
								 +vc2[typeI][typeJ]);
                     			g = fc[typeI][typeJ]*(sr8*(2.0*sr6-1.0)+f1[typeI][typeJ]
				                              +f2[typeI][typeJ]);
					temp = L*dx*g;
					fx[j] += temp;
					fx[i] -= temp;
					temp = L*dy*g;
					fy[j] += temp;
					fy[i] -= temp;
				}
			}
		}
	}
	return;
}		

/*this runs NVE and the update rules are 1. update coordinate first then update velocity by
  half time step and then calculate forces and then update velocity by another half time step*/
void advanceTimeNVE(){
	int i;
	double temp,dx,dy,sum;

	for (i=0; i<N; i++){
		dx = TIME_STEP*px[i] + 0.5*TIME_STEP*TIME_STEP*fx[i];
		rxUnFolded[i] += dx*invL;	
		temp = rx[i] + dx*invL;
		if (temp >= LENGTH)
			rx[i] = temp - LENGTH;
		else if (temp < 0.0)
			rx[i] = temp + LENGTH;
		else
			rx[i] = temp;


		dy = TIME_STEP*py[i] + 0.5*TIME_STEP*TIME_STEP*fy[i];
		ryUnFolded[i] += dy*invL;
		temp = ry[i] + dy*invL;
		if (temp >= LENGTH)
			ry[i] = temp - LENGTH;
		else if (temp < 0.0)
			ry[i] = temp + LENGTH;
		else
			ry[i] = temp;

		px[i] += HALF_TIME_STEP*fx[i];
		py[i] += HALF_TIME_STEP*fy[i];

		temp = dx*dx + dy*dy ;
		if (temp > maxD)
			maxD = temp;
	}
	nebListCounter++;
	if ( 2.0*((double)nebListCounter)*sqrt(maxD) > DR ){
		updateNebzLists();
		nebListCounter = 0;
		maxD = 0.0;
	}

	calculateForces();
	kinetic = 0.0;
	for (i=0; i<N; i++){
		px[i] += HALF_TIME_STEP*fx[i];
		py[i] += HALF_TIME_STEP*fy[i];
		kinetic += px[i]*px[i] + py[i]*py[i];
	}
	kinetic *= 0.5;
	T = 2.0*kinetic/(2.0*DOUBLE_N-3.0);

	return;
}

void advanceTimeNVT(){
	int i;
	double temp,dx,dy;
	double xi_T;

	for (i=0; i<N; i++){
		px[i] += HALF_TIME_STEP*fx[i];
		py[i] += HALF_TIME_STEP*fy[i];
		dy = TIME_STEP*py[i];
		
		ryUnFolded[i] += dy*invL;
		temp = ry[i] + dy*invL;
		if (temp >= LENGTH)
			ry[i] = temp - LENGTH;
		else if (temp < 0.0)
			ry[i] = temp + LENGTH;
		else
			ry[i] = temp;

		dx = TIME_STEP*px[i];
		rxUnFolded[i] += dx*invL;
		temp = rx[i] + dx*invL;
		if (temp >= LENGTH)
			rx[i] = temp - LENGTH;
		else if (temp < 0.0)
			rx[i] = temp + LENGTH;
		else
			rx[i] = temp;

		temp = dx*dx + dy*dy ;
		if (temp > maxD)
			maxD = temp;
	}

	nebListCounter++;
	if ( 2.0*((double)nebListCounter)*sqrt(maxD) > DR )
		updateNebzLists();

	calculateForces();
	kinetic = 0.0;
	for (i=0; i<N; i++){
		px[i] += HALF_TIME_STEP*fx[i];
		py[i] += HALF_TIME_STEP*fy[i];
		kinetic += px[i]*px[i] + py[i]*py[i];
	}
	kinetic = 0.5*kinetic;
	inst_T = kinetic/(DOUBLE_N-3.0);
	xi_T = sqrt(1.0 + TIME_STEP*(T/inst_T - 1.0)/TAU_T);
	//xi_T = sqrt(T/inst_T);
	for (i=0;i<N;i++){
		px[i] *= xi_T;
		py[i] *= xi_T;
	}
	return;
}

void advanceTimeNVTBrownClarke(){
	int i;
	double temp,dx,dy;
	double pxTemp,pyTemp;
	double xi_T;

	for (i=0; i<N; i++){
		dy = TIME_STEP*py[i];
		ryUnFolded[i] += dy*invL;
		temp = ry[i] + dy*invL;
		if (temp >= LENGTH)
			ry[i] = temp - LENGTH;
		else if (temp < 0.0)
			ry[i] = temp + LENGTH;
		else
			ry[i] = temp;

		dx = TIME_STEP*px[i];
		rxUnFolded[i] += dx*invL;
		temp = rx[i] + dx*invL;
		if (temp >= LENGTH)
			rx[i] = temp - LENGTH;
		else if (temp < 0.0)
			rx[i] = temp + LENGTH;
		else
			rx[i] = temp;

		temp = dx*dx + dy*dy;
		if (temp > maxD)
			maxD = temp;
	}

	nebListCounter++;
	if ( 2.0*((double)nebListCounter)*sqrt(maxD) > DR )
		updateNebzLists();

	calculateForces();
	kinetic = 0.0;
	for (i=0; i<N; i++){
		pxTemp = px[i] + HALF_TIME_STEP*fx[i];
		pyTemp = py[i] + HALF_TIME_STEP*fy[i];
		kinetic += pxTemp*pxTemp + pyTemp*pyTemp;
	}
	kinetic = 0.5*kinetic;
	inst_T = kinetic/(DOUBLE_N-3.0);
	xi_T = sqrt(T/inst_T);

	kinetic = 0.0;
	for (i=0;i<N;i++){
		px[i] = (2.0*xi_T-1.0)*px[i] + xi_T*TIME_STEP*fx[i];
		py[i] = (2.0*xi_T-1.0)*py[i] + xi_T*TIME_STEP*fy[i];
		kinetic += px[i]*px[i] + py[i]*py[i];
	}
	kinetic = 0.5*kinetic;
	inst_T = kinetic/(DOUBLE_N-3.0);

	return;
}

void advanceTimeNVTWithImpurity(){
	int i;
	double temp,dx,dy;
	double xi_T;

	for (i=0; i<N; i++){
		px[i] += HALF_TIME_STEP*fx[i];
		py[i] += HALF_TIME_STEP*fy[i];

		dy = TIME_STEP*py[i];
		ryUnFolded[i] += dy*invL;
		temp = ry[i] + dy*invL;
		if (temp >= LENGTH)
			ry[i] = temp - LENGTH;
		else if (temp < 0.0)
			ry[i] = temp + LENGTH;
		else
			ry[i] = temp;

		dx = TIME_STEP*px[i];
		rxUnFolded[i] += dx*invL;
		temp = rx[i] + dx*invL;
		if (temp >= LENGTH)
			rx[i] = temp - LENGTH;
		else if (temp < 0.0)
			rx[i] = temp + LENGTH;
		else
			rx[i] = temp;

		temp = dx*dx + dy*dy;
		if (temp > maxD)
			maxD = temp;
	}

	for(i=0; i<NUMBER_IMPURITY; i++){
		rx[impurityParticleIndex[i]] = rxImpurity[i];
		ry[impurityParticleIndex[i]] = ryImpurity[i];
		rxUnFolded[impurityParticleIndex[i]] = rxImpurity[i];
		ryUnFolded[impurityParticleIndex[i]] = ryImpurity[i];

	} 	

	nebListCounter++;
	if ( 2.0*((double)nebListCounter)*sqrt(maxD) > DR )
		updateNebzLists();

	calculateForces();
	for (i=0; i<N; i++){
		px[i] += HALF_TIME_STEP*fx[i];
		py[i] += HALF_TIME_STEP*fy[i];
	}

	for(i=0; i<NUMBER_IMPURITY; i++){
		px[impurityParticleIndex[i]] = 0.0;
		py[impurityParticleIndex[i]] = 0.0;
	}

	kinetic = 0.0;
	for(i=0; i<N; i++){
		kinetic += px[i]*px[i] + py[i]*py[i];
	} 	
	kinetic = 0.5*kinetic;
	inst_T = kinetic/((DOUBLE_N - NUMBER_IMPURITY)-1.0);
	xi_T = sqrt(1.0 + TIME_STEP*(T/inst_T - 1.0)/TAU_T);
	//xi_T = sqrt(T/inst_T);
	for (i=0;i<N;i++){
		px[i] *= xi_T;
		py[i] *= xi_T;
	}
	return;
}

void cool(double rate, double finalTemp){
	int i;
	FILE *outFile;
	time_t t0,t1;   /*for timing purposes...*/

	t0 = time(0);
	printf("Starting run at T = %f, u = %f\n\n",T,u);
	i = 0;
	while ( T > finalTemp ){
		/*if (!(i % 10)){
		  fprintf(stderr, "\r");
		  fprintf(stderr, "Running, T = : %f u = %.10f, kinetic = %.10f",T,u,kinetic);
		  }*/
		advanceTimeNVTBrownClarke();
		//advanceTimeNVT();
		T -= TIME_STEP*rate;
		i++;
	}
	t1 = time(0);
	printf("\n\n Time in seconds: %ld\n", t1-t0);

	return;
}

int prepareSavingArray(unsigned long long int runLength, int numOfOrigins, double factor){
	int i,j,k;
	unsigned long long int linearInterval,maximalInterval,index;
	int offset,current;
	char nextStepIndexFileName[128];
	FILE *outFile;

	linearInterval = runLength/numOfOrigins;
	//maximalInterval = runLength/DIVIDE_TO;
	maximalInterval = runLength;

	printf("no of origins = %d\n",numOfOrigins);
	printf("linearInterval = %lld\n",linearInterval);
	printf("maximalInterval = %lld\n",maximalInterval);

	current = 0;
	for (k=0;k<numOfOrigins;k++){
		nextStepIndex[current] = k*linearInterval;
		current++;
		offset = 1; //the smallest interval
		while (offset < maximalInterval){
			index = k*linearInterval + offset;
			if (index<runLength){
				nextStepIndex[current] = index;
				current++;
			}
			if ( (int)(offset*factor) == offset )
				offset++;
			else
				offset = (int)(offset*factor);
		}
	}

	shellSort(current, nextStepIndex);
	j=0; i=0;
	while (j<current){
		while ( nextStepIndex[j] == nextStepIndex[i] )
			j++;
		i++;
		nextStepIndex[i] = nextStepIndex[j];
	}
	current = i+1;

	printf("current = %d\n",current);
	sprintf(nextStepIndexFileName,"nextIndex_%d_%.2f_%.3d.dat",N,T,jobNo);	//Here it is changed
	outFile = fopen(nextStepIndexFileName,"wb");
	for (i=0;i<current;i++)
		fprintf(outFile,"%lld\n",nextStepIndex[i]);
	return current;
}

//run the system for <steps> time steps
void equilibrate(double duration){
	int i,steps;
	int t0,t1;      /*for timing purposes...*/
	char restartFileName[128];
	int start;

	t0 = time(0);
	steps = (int)(duration/TIME_STEP);
	sprintf(restartFileName,"%.3d/restartEq_%d_%.3f_%.3d.dat",jobNo,N,T,jobNo);		//Here it is changed

	if(restartEq){
		start = iterFinish;
		printf("Restarting equilibration at T = %g, u = %g\n\n",T,u);
	}
	else{
		start = 0;
		printf("Starting equilibration at T = %g, u = %g\n\n",T,u);
	}

	for (i=start; i<steps; i++){
		advanceTimeNVT();
		if(!(i%((int)(steps/10))))
			saveCompleteRestartState(restartFileName,i);
	}

	if(!((i+1)%((int)(steps/10))))
      		saveCompleteRestartState(restartFileName,i);

	// setting the unfolded coordinates for the production run	
	for(i=0;i<N;i++){
		rxUnFolded[i] = rx[i];
		ryUnFolded[i] = ry[i];
	}

	t1 = time(0);
	printf("equilibration at T = %g - Time in seconds: %d\n\n\n",T, t1-t0);
	return;
}

void equilibratePinned(double duration){
	int k,i,steps,stepsCompleted,between,current;
	int t0,t1;      /*for timing purposes...*/
	char stuffFileName[128], dataFileName[128],dataEqFileName[128];
	char tempFileName[128];
	FILE *stuffFile, *dataFile, *dataEqFile, *lambdaFile;
	int counter, start, intDummy, numOfFrames;
	double mu,lambda,lambdaExact;
	double dummy;
	char restartFileName[128],moveDataFile[128];
	FILE *tempFile;

	sprintf(restartFileName,"%.3d/restartPostEq_%d_%.3f_%.3d.dat",jobNo,N,T,jobNo);	//Here it is changed

	t0 = time(0);
	steps = (int)(duration/TIME_STEP);
	if(restartPostEq){
		start = iterFinish;
		printf("Restarting post eq at T = %g, u = %g\n\n",T,u);
	}
	else{
		start = 0;
		printf("Starting post eq at T = %g, u = %g\n\n",T,u);
	}

	for (k=start; k<steps; k++){
		advanceTimeNVTWithImpurity();
		if(!(k%((int)(steps/10))))
			saveCompleteRestartState(restartFileName,k);

		if ( !(k%655360) )
			fixDrift();
	}

	if(!((k+1)%((int)(steps/10))))
      		saveCompleteRestartState(restartFileName,k);

	// setting the unfolded coordinates for the production run	
	dummy=0;
	for(i=0;i<N;i++){
		rxUnFolded[i] = rx[i];
		ryUnFolded[i] = ry[i];
	}
	t1 = time(0);
	printf("Post-equilibration at T = %g - Time in seconds: %d\n\n\n",T, t1-t0);

	return;
}

void run(double duration){
	unsigned long long int k,i,steps,between,current,start,numOfFrames;
	int t0,t1;      /*for timing purposes...*/
	char stuffFileName[128], dataFileName[128];
	FILE *stuffFile, *dataFile;
	char restartFileName[128],tempFileName[128],moveDataFile[128];	//Here it is changed
	FILE *tempFile;
	double dummy;
	int intDummy;

	sprintf(tempFileName,"%.3d/tempData.dat",jobNo);
	sprintf(restartFileName,"%.3d/restartRun_%d_%.3f_%.3d.dat",jobNo,N,T,jobNo);
	sprintf(stuffFileName,"%.3d/stuff_%d_%.3f_%.3d.dat",jobNo,N,T,jobNo);				//Here it is changed
	sprintf(dataFileName,"%.3d/data_%d_%.3f_%.3d.dat",jobNo,N,T,jobNo);		
	

	t0 = time(0);
	steps = (unsigned long long int)(duration/TIME_STEP);
	between = (unsigned long long int)(1.0/TIME_STEP);
	printf("steps = %lld\n",steps);
	printf("Starting run at T = %g, u = %g\n\n",T,u);
	if(restartRun){
		printf("Restarting run at T = %g, u = %g\n\n",T,u);
		stuffFile = fopen(stuffFileName,"ab");
		dataFile = fopen(dataFileName,"rb"); // lets see whether the writing is done properly
		k = 0;
		while ( fscanf(dataFile,"%lf %lf %lf %d",&(dummy),&(dummy),&(dummy),&(intDummy)) != EOF )
			k++;
		numOfFrames = (int)k/DOUBLE_N;
		current = numOfFrames - 1;
		start = nextStepIndex[current] + 1;
		current = current + 1;
		rewind(dataFile);
		tempFile = fopen(tempFileName,"wb");
		for(k=0;k<numOfFrames;k++){
			for(i=0;i<N;i++){
				fscanf(dataFile,"%lf %lf %d",&(rx[i]),&(ry[i]),&(type[i]));
				fprintf(tempFile,"%.10g\t%.10g\t%d\n",rx[i],ry[i],type[i]);
			}
		}
		fclose(dataFile);
		fclose(tempFile);
		sprintf(moveDataFile,"mv %.3d/tempData.dat %.3d/data_%d_%.3f_%.3d.dat",jobNo,jobNo,N,T,jobNo);	//Here it is changed
		system(moveDataFile);
		dataFile = fopen(dataFileName,"ab"); // lets see whether the writing is done properly
		initializeSystem();	
	}
	else{
		printf("Starting run at T = %g, u = %g\n\n",T,u);
		current = 0;
		start = 0;
		stuffFile = fopen(stuffFileName,"wb");
		dataFile = fopen(dataFileName,"wb");
	}
	for (k=start; k<steps; k++){
		advanceTimeNVTWithImpurity();
		if ( k==nextStepIndex[current] ){
			fprintf(stuffFile,"%lld\t%.10g\t%.10g\t%.10g\n",k,u,kinetic,inst_T);
			//printf("%d\t%.10g\t%.10g\t%.10g\n",k,u,kinetic,inst_T);
			for (i=0;i<N;i++)
				fprintf(dataFile,"%.10g\t%.10g\t%d\n",rxUnFolded[i],ryUnFolded[i],type[i]);
			current++;
		}
		if(!(k%100000))
			saveCompleteRestartState(restartFileName,k);

		if ( !(k%65536) )
			fixDrift();
	}

	printf("Last index of this run =%lld\n",current);

	t1 = time(0);
	printf("data collection at T = %g - Time in seconds: %d\n\n\n",T, t1-t0);
	fclose(stuffFile);
	fclose(dataFile);
	return;
}

void selectRandomImpurity(){				//Thats ok
	int i,temp,k;
	int indices[N];
	FILE *impurityFile;
	char impurityFileName[128];

	for(i=0;i<N;i++)
		indices[i] = i;

	//permute indices:
	for (i=0; i<N; i++){
		k = (int)(ran1(&for_ran1)*DOUBLE_N); //choose random index to switch with
		temp = indices[i];
		indices[i] = indices[k];
		indices[k] = temp;
	}

	sprintf(impurityFileName,"%.3d/impurityFileIndex_%d_%.3d.dat",jobNo,N,jobNo);		//Here it is changed
	impurityFile = fopen(impurityFileName,"wb"); 	
	printf("number of impurity = %d\n",NUMBER_IMPURITY);
	for(i=0; i<NUMBER_IMPURITY; i++){
		impurityParticleIndex[i] = indices[i];
		fprintf(impurityFile,"%d\n",impurityParticleIndex[i]);	
		rxImpurity[i] = rx[impurityParticleIndex[i]];
		ryImpurity[i] = ry[impurityParticleIndex[i]];
	}
	fclose(impurityFile);

	return;
}

void readImpurity(char *impurityFile){
	int i;
	FILE *file;

	file = fopen(impurityFile,"rb");
	for(i=0; i<NUMBER_IMPURITY; i++){
		fscanf(file, "%d",&(impurityParticleIndex[i]));
		rxImpurity[i] = rx[impurityParticleIndex[i]];
		ryImpurity[i] = ry[impurityParticleIndex[i]];
	}
	fclose(file);
	return;
}

int main(int argc,char *argv[]){
	unsigned long long int runLength;
	int numOfOrigins,nextStepIndexSize,restartFlag;
	double factor,temp,duration,durationPostPin,prodDuration,T_org;
	char snapFileName[128],impurityFileName[128],dir[128];
	int ierr;

	ierr = MPI_Init( &argc, &argv );
	ierr = MPI_Comm_size ( MPI_COMM_WORLD, &NUM_TASK );
	ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &TASK_ID );

	if ( TASK_ID == 0 ){
		printf ( "\n" );
		printf ( "  MPI is initialized and running on %i processors\n", NUM_TASK );
	}	

	sscanf(argv[1],"%d",&serial);								
	jobNo = serial*NUM_TASK + TASK_ID; // this is only for these set of runs			//Here it is changed
	for_ran1 = -(int)(time(0)+jobNo);
	printf("random no seed = %d\n",for_ran1);
	ran1(&for_ran1); //initialize the random number
	srand( jobNo + 100 - for_ran1);
	sprintf(dir,"mkdir -p %.3d",jobNo);
	system(dir);

	duration =2*50000;
	durationPostPin = 0;
	prodDuration =2*50000;

	/*********** initialization ************/
	T =0.45;
	T_org = T;
	sprintf(impurityFileName,"%.3d/impurityFileIndex_%d_%.3d.dat",jobNo,N,jobNo);
	runLength = (unsigned long long int)(prodDuration/TIME_STEP);
	factor = 1.4;
	numOfOrigins = 200;

	restartEq=0;
	restartPostEq=0;
	restartRun=0;
	restartFlag = 0;
	/* RRR = 0 => fresh run */
	/* RRR = 1 => restart equilibration */
	/* RRR = 2 => restart post-equilibration */
	/* RRR = 3 => restart production */

	switch (restartFlag){
		case 0:
			initializeGrid();
			cool(1e-1,0.40);
			printf("energy after cooling u = %.10g\n",u);
			T = T_org;
			adjustStartingTemperature();
			if (duration>0) equilibrate(duration);
			selectRandomImpurity();//impurity particles are selected 
			if (durationPostPin>0) equilibratePinned(durationPostPin);
			break;
		case 1:
			restartEq=1;
			sprintf(snapFileName,"%.3d/restartEq_%d_%.3f_%.3d.dat",jobNo,N,T,jobNo);
			printf("restart equilibration from file %s\n",snapFileName);
			readCompleteRestartState(snapFileName);
			printf("energy after previous run u = %.10g\n",u);
			if (duration>0) equilibrate(duration);
			selectRandomImpurity();//impurity particles are selected 
			if (durationPostPin>0) equilibratePinned(durationPostPin);
			break;
		case 2:
			restartPostEq=1;
			printf("restart post-equilibration\n");
			sprintf(snapFileName,"%.3d/restartPostEq_%d_%.3f_%.3d.dat",jobNo,N,T,jobNo);
			readCompleteRestartState(snapFileName);
			readImpurity(impurityFileName);
			if (durationPostPin>0) equilibratePinned(durationPostPin);
			break;
		case 3:
			restartRun=1;
			printf("restart production\n");
			sprintf(snapFileName,"%.3d/restartRun_%d_%.3f_%.3d.dat",jobNo,N,T,jobNo);
			readCompleteRestartState(snapFileName);
			readImpurity(impurityFileName);
			break;
	}


	if (prodDuration>0){
		/**************prepare saving index array*****************/
		nextStepIndexSize = (numOfOrigins) * (10+(int)( log((double)((runLength)/DIVIDE_TO))/log(factor) ) );
		nextStepIndex = (unsigned long long int *)malloc(sizeof(unsigned long long int)*nextStepIndexSize);
		prepareSavingArray(runLength,numOfOrigins,factor);
		/**************run*****************/
		run(prodDuration);
		/**********************************/
	}

	free(nextStepIndex);
	ierr = MPI_Finalize ( );

	return 0;
}
