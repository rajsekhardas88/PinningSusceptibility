#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "/export/data/mpi/openmpi/1.4.5/intel/include/mpi.h"

// (N,NA) = (1000, 800) (2000, 1600), (4000, 3200) (28160, 22528)
/*---------------------------- CONSTANTS IN THE SYSTEM ------------------------------------------------*/
#define	N					10000
#define	DOUBLE_N				10000.0
#define NA					5000
#define	CUTOFF_SQRD				1.919383103666485
#define	INV_CUTOFF_SQRD				(1.0/CUTOFF_SQRD)
#define	C_0					-0.806140903539923
#define	C_2					0.7
#define	C_4					-0.156300219287607
#define	N_0					10.0
#define	DR					0.35
#define	MAX_NEBZ				512
#define	MAX_CELLS				(N/2) 
#define	MAX_BUDS_PER_CELL			1000
#define	TIME_STEP				0.005
#define	HALF_TIME_STEP				0.5*TIME_STEP
#define	LENGTH					1.0
#define	HALF_LENGTH				0.5
#define DENSITY 				0.81
#define TAU_T					5.0
#define	TWO_PI					6.28318530717958
#define LARGE                          		1
#define SMALL                           	0
#define DIVIDE_TO                               2
#define IMPURITY_DENSITY			0.000
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
double rz[N];			/*	z component of position 	*/
double rxUnFolded[N];		/*	x component of position		*/
double ryUnFolded[N];		/*	y component of position		*/
double rzUnFolded[N];		/*	z component of position 	*/
double px[N];			/*	x component of momentum		*/
double py[N];			/*	y component of momentum		*/
double pz[N];			/*	z component of momentum		*/
double fx[N];			/*	x component of force		*/
double fy[N];			/*	y component of force		*/
double fz[N];			/*	z component of force		*/
int type[N]; 			/*	type large or small		*/
double leDisplacement;  	/* 	LEE-EDWARDS displacement  	*/
int impurityParticleIndex[NUMBER_IMPURITY];
double rxImpurity[NUMBER_IMPURITY];
double ryImpurity[NUMBER_IMPURITY];
double rzImpurity[NUMBER_IMPURITY];

/*-------------- restart job options -------------------------------------------------------------------*/
int restartEq,restartPostEq;
int restartRun;
int iterFinish;
/*---------------------------- PARAMETERS FOR THE SYSTEM -----------------------------------------------*/
const double invSizeSqrd[2][2] = {{1.00,0.718184429761563},{0.718184429761563,0.510204081632653}};
// [0][0] = large-large; [0][1] = [1][0] = large-small; [1][1] = small-small;
/*------------------------------------- FOR NEBZ LIST -----------------------------------------------------*/
int nebz[N][MAX_NEBZ]; 		//maximum we give here MAX_NEBZ nebz. 
int numOfNebz[N]; 		//number of elements in nebz[N][MAX_NEBZ]
double maxD,listCutOffSqrd; 	//for updating nebz list.
int nebListCounter; 		//for counting how many steps have took place for nebz list updating.
/*------------------------------------FOR CELL SUBDIVISION ------------------------------------------------*/
double cellCutOff;
int numInCell[MAX_CELLS];
int cells[MAX_CELLS][MAX_BUDS_PER_CELL];
const int shift[13][3] = {{-1,-1,1},{-1,0,1},{-1,1,1},{0,-1,1},{0,0,1},{0,1,1},{1,-1,1},
	{1,0,1},{1,1,1},{-1,-1,0},{-1,0,0},{-1,1,0},{0,1,0}};
/*--------------------- FOR LOGARITHAMIC STORAGE -----------------------------------------------------------*/
int *nextStepIndex;
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

void shellSort(int n, int *a){
	/*Sorts an array a[1..n] into ascending numerical order by Shell's method (diminishing increment
	  sort). n is input; a is replaced on output by its sorted rearrangement.*/
	int i,j,inc;
	int v;

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

int getUpperThrdRoot(int number){
	int q=1;
	while (q*q*q<number)
		q++;
	return q;
}

void initializeGrid(){
	int i,thrdrtN,j,temp;
	double space,tempd;

	V = DOUBLE_N/DENSITY;
	thrdrtN = getUpperThrdRoot(N);
	printf("got %d\n",thrdrtN);
	L = pow(V,1.0/3.0);
	printf("Simulation box size L = %.10f\n",L);
	space = 1.0/(double)thrdrtN;

	for (i=0; i<N; i++){
		rx[i] = space*(double)( (i%(thrdrtN*thrdrtN)) % thrdrtN );
		ry[i] = space*(double)((i%(thrdrtN*thrdrtN)) / thrdrtN );
		rz[i] = space*(double)(i/(thrdrtN*thrdrtN));
	}

	//permutate indices
	for (i=0; i<N; i++){
		j = (int)(DOUBLE_N*UNI);

		tempd = rx[i];
		rx[i] = rx[j];
		rx[j] = tempd;

		tempd = ry[i];
		ry[i] = ry[j];
		ry[j] = tempd;

		tempd = rz[i];
		rz[i] = rz[j];
		rz[j] = tempd;

		px[i] = 2.0*ran1(&for_ran1)-1.0;
		py[i] = 2.0*ran1(&for_ran1)-1.0;
		pz[i] = 2.0*ran1(&for_ran1)-1.0;
	}

	for (i=0;i<N;i++){
		if(i<NA) type[i] = 0;
		else type[i] = 1;
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
		fprintf(file,"%.8f\t%.8f\t%.8f\t%d\n",rx[i],ry[i],rz[i],type[i]);
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
		fprintf(file,"%.8g\t%.8g\t%.8g\t%d\t%.8g\t%.8g\t%.8g\n",rx[i],ry[i],rz[i],type[i]
				,px[i],py[i],pz[i]);
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
		fprintf(file,"%.8g\t%.8g\t%.8g\t%d\t%.8g\t%.8g\t%.8g\n",rx[i],ry[i],rz[i],type[i]
				,px[i],py[i],pz[i]);
	}
	fclose(file);

	return;
}

void readSnapShot(char *fileName){
	int i;
	double dummy;
	FILE *readFile;

	readFile = fopen(fileName, "rb");
	fscanf(readFile, "%lf",&(L));
	fscanf(readFile, "%lf",&(u));
	fscanf(readFile, "%lf",&(T));
	fscanf(readFile, "%lf",&(dummy));
	for (i=0;i<N;i++){
		fscanf(readFile, "%lf",&(rx[i]));
		fscanf(readFile, "%lf",&(ry[i]));
		fscanf(readFile, "%lf",&(rz[i]));
		fscanf(readFile, "%d",&(type[i]));
		px[i] = 0.0;	
		py[i] = 0.0;	
		pz[i] = 0.0;
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
	for (i=0;i<3;i++){
		fscanf(file,"%lf",&dummy);
	}
	for (i=0;i<N;i++){
		fscanf(file, "%lf",&(rx[i]));
		fscanf(file, "%lf",&(ry[i]));
		fscanf(file, "%lf",&(rz[i]));
		fscanf(file, "%d",&(type[i]));
		fscanf(file, "%lf",&(px[i]));
		fscanf(file, "%lf",&(py[i]));
		fscanf(file, "%lf",&(pz[i]));
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
		fscanf(file, "%lf",&(rz[i]));
		fscanf(file, "%d",&(type[i]));
		fscanf(file, "%lf",&(px[i]));
		fscanf(file, "%lf",&(py[i]));
		fscanf(file, "%lf",&(pz[i]));
		printf("%g\t%g\t%g\t%d\t%g\t%g\t%g\n",rx[i],ry[i],rz[i],type[i],px[i],py[i],pz[i]);
	}
	fclose(file);
	initializeSystem();

	return;
}

void fixDrift(){
	double Px,Py,Pz;
	int i;
	Px = 0.0; Py = 0.0; Pz = 0.0; 
	for (i=0;i<N;i++){
		Px += px[i]; Py += py[i]; Pz += pz[i];
	}
	Px = Px/DOUBLE_N; Py = Py/DOUBLE_N; Pz = Pz/DOUBLE_N;
	for (i=0;i<N;i++){
		px[i] -= Px; py[i] -= Py; pz[i] -= Pz;
	}

	return;
}

void applyPBC(){
	int i;

	for(i=0;i<N;i++){
		if (rx[i] >= LENGTH)
			rx[i] = rx[i] - LENGTH;
		else if (rx[i] < 0.0)
			rx[i] = rx[i] + LENGTH;

		if (ry[i] >= LENGTH)
			ry[i] = ry[i] - LENGTH;
		else if (ry[i] < 0.0)
			ry[i] = ry[i] + LENGTH;

		if (rz[i] >= LENGTH)
			rz[i] = rz[i] - LENGTH;
		else if (rz[i] < 0.0)
			rz[i] = rz[i] + LENGTH;
	}
	return;
}

void adjustStartingTemperature(){
	int i;
	double xi_T;

	kinetic = 0.0;
	for(i=0;i<N;i++){
		kinetic += px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i];
	} 
	kinetic = 0.50*kinetic;
	inst_T = (2.0*kinetic)/(3.0*DOUBLE_N - 4.0);
	xi_T = sqrt(T/inst_T);

	for(i=0;i<N;i++){
		px[i] *= xi_T;
		py[i] *= xi_T;
		pz[i] *= xi_T;
	}
	return;
}
//initiallize system - stuff for cell subdivision and nebz list.
void initializeSystem(){
	double rcutSq,sr2,sr6,sr12;
	int i,j;

	invL = 1.0/L;
	V = L*L*L;
	cellCutOff = 1.4*sqrt(CUTOFF_SQRD) + DR;
	listCutOffSqrd = cellCutOff*cellCutOff;
	strain = 0.0; // notice that the strain set to be zero as the minimizer is for general condition
	leDisplacement = fmod(strain,1.0);
	if (leDisplacement < -0.5)
		leDisplacement += 1.0;
	else if (leDisplacement > 0.5)
		leDisplacement -= 1.0;

	applyPBC();
	updateNebzLists();
	calculateForces();
	nebListCounter = 0;
	maxD = 0.0;
	fixDrift();

	for(i=0;i<N;i++){
		rxUnFolded[i] = rx[i];
		ryUnFolded[i] = ry[i];
		rzUnFolded[i] = rz[i];
	}	
	adjustStartingTemperature();

	return;
}
//initiallize system - stuff for cell subdivision and nebz list.
void initializeSystemRestart(){
	double rcutSq,sr2,sr6,sr12;
	int i,j;

	invL = 1.0/L;
	V = L*L*L;
	cellCutOff = 1.4*sqrt(CUTOFF_SQRD) + DR;
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
		rzUnFolded[i] = rz[i];
	}	
	//adjustStartingTemperature();

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
		fprintf(outCoordFile,"%g\t%g\t%g\t%d\n",rxUnFolded[k],ryUnFolded[k],rzUnFolded[k],type[k]);
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

void updateNebzLists(){
	int i,j,x,y,z,current,m,m2,m3;
	int a,b,c,k,l,numHere,numThere,target,w,xIndex,yIndex;
	double dx,dy,dz,r2,invCellSize;
	double rxi,ryi,rzi;

	nebListCounter = 0;
	maxD = 0.0;

	if ( N < 512 ){
		//set all neb lists to zero
		for (i=0; i<N; i++)
			numOfNebz[i] = 0;

		for (i=0;i<N-1;i++){
			rxi = rx[i];
			ryi = ry[i];
			rzi = rz[i];
			for (j=i+1; j<N; j++){
				dx = rx[j] - rxi;
				dy = ry[j] - ryi;
				dz = rz[j] - rzi;
				// mess due to LEES EDWARDS periodic boundary conditions 
				if ( dz >= HALF_LENGTH ){
					dz -= LENGTH;
				}
				else if ( dz < -HALF_LENGTH ){
					dz += LENGTH;
				}
				if ( dx >= HALF_LENGTH )
					dx -= LENGTH;
				else if ( dx < -HALF_LENGTH )
					dx += LENGTH;
				if ( dy >= HALF_LENGTH )
					dy -= LENGTH;
				else if ( dy < -HALF_LENGTH )
					dy += LENGTH;
				// end of mess 
				r2 = L*L*( dx*dx + dy*dy + dz*dz );
				if (r2 < listCutOffSqrd){
					nebz[i][numOfNebz[i]] = j;
					nebz[j][numOfNebz[j]] = i;
					numOfNebz[i]++;
					numOfNebz[j]++;
				}
			}
		}
	}
	else{
		m =(int)(L/cellCutOff);
		m3 = m*m*m; 		//the length of cells and numInCell...
		m2 = m*m;
		invCellSize = (double)m; 	// for reduced coordinates...
		//set number in each cell to zero
		for (i=0; i<m3; i++)
			numInCell[i] = 0;

		//********cells sorting********
		for (i=0;i<N;i++){
			current = m2*(int)(rz[i]*invCellSize)+m*(int)(ry[i]*invCellSize)+(int)(rx[i]*invCellSize);
			cells[current][numInCell[current]] = i;
			numInCell[current]++;
		}

		//set all neb lists to zero
		for (i=0; i<N; i++)
			numOfNebz[i] = 0;

		for (z=0;z<m;z++){
			for (y=0;y<m;y++){
				for (x=0;x<m;x++){
					current = m2*z + m*y + x; //this cell index.
					numHere = numInCell[current];
					//first check interactions within a cell
					for (k=0;k<numHere-1;k++){
						for (l=k+1;l<numHere;l++){
							i = cells[current][k];
							j = cells[current][l];
							dx = rx[j] - rx[i];
							dy = ry[j] - ry[i];
							dz = rz[j] - rz[i];
							r2 = L*L*( dx*dx + dy*dy + dz*dz );
							if (r2 < listCutOffSqrd){
								nebz[i][numOfNebz[i]] = j;
								nebz[j][numOfNebz[j]] = i;
								numOfNebz[i]++;
								numOfNebz[j]++;
							}
						}
					}

					//shift[13][3] = {{-1,-1,1},{-1,0,1},{-1,1,1},{0,-1,1},{0,0,1},
					//{0,1,1},{1,-1,1},{1,0,1},{1,1,1},{-1,-1,0},{-1,0,0},{-1,1,0},{0,1,0}}
					//now there are 13 nebz that need to be checked.
					for (w=0;w<13;w++){ 
						a = x + shift[w][0];
						if (a<0) a = m-1;
						else if (a==m) a = 0;
						b = y + shift[w][1];
						if (b<0) b = m-1;
						else if (b==m) b = 0;
						c = z + shift[w][2];
						if (c<0) c = m-1;
						else if (c==m) c = 0;
						target = a + m*b + m2*c;
						numThere = numInCell[target];
						for (k=0;k<numHere;k++){
							for (l=0;l<numThere;l++){
								i = cells[current][k];
								j = cells[target][l];
								dx = rx[j] - rx[i];
								dy = ry[j] - ry[i];
								dz = rz[j] - rz[i];
								// mess due to periodic boundary conditions 
								if ( dx >= HALF_LENGTH )
									dx -= LENGTH;
								else if ( dx < -HALF_LENGTH )
									dx += LENGTH;
								if ( dy >= HALF_LENGTH )
									dy -= LENGTH;
								else if ( dy < -HALF_LENGTH )
									dy += LENGTH;
								if ( dz >= HALF_LENGTH )
									dz -= LENGTH;
								else if ( dz < -HALF_LENGTH )
									dz += LENGTH;
								// end of mess 
								r2 = L*L*( dx*dx + dy*dy + dz*dz );
								if (r2 < listCutOffSqrd){
									nebz[i][numOfNebz[i]] = j;
									nebz[j][numOfNebz[j]] = i;
									numOfNebz[i]++;
									numOfNebz[j]++;
								}
							}
						}
					}
				}
			}
		}
	}

	return;
}

void calculateForces(){
	int i,j,m,k,typeI;
	double g,dx,dy,dz,rxi,ryi,rzi,r2OverSigma2,invSigmaSqrd;
	double sigmaOverR12,temp;

	//first set 'em all to zero
	u = 0.0;
	for (i=0;i<N;i++){
		fx[i] = 0.0;
		fy[i] = 0.0;
		fz[i] = 0.0;
	}

	for (i=0;i<N-1;i++){
		rxi = rx[i];
		ryi = ry[i];
		rzi = rz[i];
		typeI = type[i];
		m = numOfNebz[i];
		for (k=0;k<m;k++){
			j = nebz[i][k];
			if (j > i){
				dx = rx[j] - rxi;
				dy = ry[j] - ryi;
				dz = rz[j] - rzi;
				// mess due to LEES EDWARDS periodic boundary conditions 
				if ( dz >= HALF_LENGTH )
					dz -= LENGTH;
				else if ( dz < -HALF_LENGTH )
					dz += LENGTH;
				if ( dx >= HALF_LENGTH )
					dx -= LENGTH;
				else if ( dx < -HALF_LENGTH )
					dx += LENGTH;
				if ( dy >= HALF_LENGTH )
					dy -= LENGTH;
				else if ( dy < -HALF_LENGTH )
					dy += LENGTH;
				// end of mess 
				invSigmaSqrd = invSizeSqrd[typeI][type[j]];
				r2OverSigma2 = L*L*( dx*dx + dy*dy + dz*dz )*invSigmaSqrd;
				if (r2OverSigma2 < CUTOFF_SQRD){
					sigmaOverR12 = 1.0/(r2OverSigma2*r2OverSigma2*r2OverSigma2*r2OverSigma2*r2OverSigma2*r2OverSigma2);
					g = (N_0*sigmaOverR12 - r2OverSigma2*4.0*C_4  - 2.0*C_2)*invSigmaSqrd;
					temp = L*dx*g;
					fx[j] += temp;
					fx[i] -= temp;
					temp = L*dy*g;
					fy[j] += temp;
					fy[i] -= temp;
					temp = L*dz*g;
					fz[j] += temp;
					fz[i] -= temp;
					u += r2OverSigma2*(sigmaOverR12 + C_2 + r2OverSigma2*C_4) + C_0;
				}
			}
		}
	}
	return;
}
/*********************************** MINIMIZER FUNCTIONS *************************************************/
void grad(double *point, double *res){
	int i,k,j,m,typeI;
	double rxi,ryi,rzi,dx,dy,dz,invSigmaSqrd, halfL;
	double interactionFactor,sigma2OverR2,sigma6OverR6;
	double r2OverSigma2,g,temp,sr2,sr6,sr8,rijSq;

	halfL = 0.5*L;

	//set to zero
	for (i=0;i<3*N;i++)
		res[i] = 0.0;

	for (i=0;i<N-1;i++){
		rxi = point[3*i];
		ryi = point[3*i+1];
		rzi = point[3*i+2];
		typeI = type[i];
		m = numOfNebz[i];
		for (k=0;k<m;k++){
			j = nebz[i][k];
			if (j > i){
				dx = point[3*j] - rxi;
				dy = point[3*j+1] - ryi;
				dz = point[3*j+2] - rzi;
				// mess due to LEES EDWARDS periodic boundary conditions 
				if ( dz >= halfL){
					dz -= L;
					dx -= leDisplacement*L;
				}
				else if ( dz < -halfL ){
					dz += L;
					dx += leDisplacement*L;
				}
				if ( dx >= halfL )
					dx -= L;
				else if ( dx < -halfL )
					dx += L;
				if ( dy >= halfL )
					dy -= L;
				else if ( dy < -halfL )
					dy += L;
				// end of mess 
				invSigmaSqrd = invSizeSqrd[typeI][type[j]];
				rijSq = ( dx*dx + dy*dy + dz*dz );
				r2OverSigma2 = rijSq*invSigmaSqrd; //real space coordinates...
				if (r2OverSigma2 < CUTOFF_SQRD){
					sr2=1.0/r2OverSigma2;
					sr6=sr2*sr2*sr2;
					sr8=sr6*sr2;
					//g = fc[typeI][type[j]]*(sr8*(2.0*sr6-1.0)+f1[typeI][type[j]]
					//			  +f2[typeI][type[j]]);
					temp = dx*g;
					res[3*j] -= temp;
					res[3*i] += temp;
					temp = dy*g;
					res[3*j+1] -= temp;
					res[3*i+1] += temp;
					temp = dz*g;
					res[3*j+2] -= temp;
					res[3*i+2] += temp;
				}
			}
		}
	}
	return;
}

//take care of periodic boundary conditions - move to model specific file 
//also check if neb lists need updating...
void fix(double *point){
	int i,xIndex,yIndex,zIndex;
	double tmpd,dx,dy,dz,length,halfLength;

	maxD = 0.0; length = L; halfLength = 0.5*L;

	for (i=0;i<N;i++){
		xIndex = 3*i;
		yIndex = xIndex+1;
		zIndex = yIndex+1;

		if (point[zIndex] >= length){
			point[zIndex] -= length;
			point[xIndex] -= length*leDisplacement;
		}
		else if (point[zIndex] < 0.0){
			point[zIndex] += length;
			point[xIndex] += length*leDisplacement;
		}
		if (point[xIndex] >= length)
			point[xIndex] -= length;
		else if (point[xIndex] < 0.0)
			point[xIndex] += length;
		if (point[yIndex] >= length)
			point[yIndex] -= length;
		else if (point[yIndex] < 0.0)
			point[yIndex] += length;

		//after fixing boundary conditions, I can now look for the maximal distance
		//traveled by some particle
		dx = point[xIndex] - rx[i]*L;
		dy = point[yIndex] - ry[i]*L;
		dz = point[zIndex] - rz[i]*L;
		// mess due to LEES EDWARDS periodic boundary conditions 
		if ( dz >= halfLength ){
			dz -= length;
			dx -= length*leDisplacement;
		}
		else if ( dz < -halfLength ){
			dz += length;
			dx += length*leDisplacement;
		}
		if ( dx >= halfLength )
			dx -= length;
		else if ( dx < -halfLength )
			dx += length;
		if ( dy >= halfLength )
			dy -= length;
		else if ( dy < -halfLength )
			dy += length;
		// end of mess 
		tmpd = dx*dx + dy*dy + dz*dz;
		if ( tmpd > maxD )
			maxD = tmpd;
	}
	if ( 2.0*sqrt(maxD) > DR ){
		for (i=0;i<N;i++){
			rx[i] = point[3*i]*invL;
			ry[i] = point[3*i+1]*invL;
			rz[i] = point[3*i+2]*invL;
		}
		updateNebzLists();
	}
	//that's it.
	return;
}

void initializeMinimizer(int start){ 
	int i;
	if ( start==0 ) 
		lastx = LAST_X_DEFAULT;
	for (i=0;i<3*N;i++){
		if ( restart !=2 )
			g[i] = -q[i];
		q[i] = h[i] = g[i];
	}
	restart = 0;
	return;
}

double prod(double x, double *gradScratch){
	int i;
	double s;
	for (i=0;i<3*N;i++)
		scratch[i] = p[i] + x*q[i];
	fix(scratch);
	grad(scratch, gradScratch);
	for (s = 0.0, i=0; i<3*N; i++)
		s += gradScratch[i]*q[i];
	return s;
}

double linmin(){
	double x,y,s,t,m,tmpd,step;
	double *gx,*gy,*gUnused;
	int it,i;

	gx = mosesX; gy = mosesY;

	x = lastx/gtyp;
	s = prod(x, gx); //result goes into gx
	it = 0;
	if ( s<0.0 ){
		do{
			y = x*LINMIN_G1;
			t = prod(y, gy);
			if ( t >= 0.0 ) break ;
			x = y; s = t; gUnused = gx; gx = gy; gy = gUnused; 
			it++;
		}while (it < LINMIN_MAX_ITERATIONS);
	}
	else if ( s>0 ){		
		do{
			y = x*LINMIN_G3;
			t = prod(y, gy);
			if ( t <= 0.0 ) break ;
			x = y ; s = t ; gUnused = gx ; gx = gy ; gy = gUnused;
			it++;
		}while (it < LINMIN_MAX_ITERATIONS);
	}
	else{		/* hole in one s = 0.0 */
		t = 1.0; y = x;
	}

	if ( it >= LINMIN_MAX_ITERATIONS){
		printf("Warning! linmin overran, ");
		tmpd = prod(0.0, gy);
		if ( tmpd > 0.0 )
			restart = 1;
	}

	if ( s < 0.0 ) s = - s;
	if ( t < 0.0 ) t = - t;
	m = ( s + t );
	s /= m; t /= m;

	m =  s * y + t * x;
	/* evaluate the step length, not that it necessarily means anything */
	for ( step=0.0, i=0; i<3*N; i++ ){
		tmpd = m*q[i];
		p[i] += tmpd; 			/* this is the point where the parameter vector steps */
		step += fabs(tmpd);
		q[i] = s*gy[i] + t*gx[i];	/* send back the estimated gradient in q (NB not like linmin) */
	}

	fix(p);			/* take care of periodic boundary conditions */
	lastx = m*LINMIN_G2*gtyp;
	return ( step / (3.0*DOUBLE_N) ); 
}

int min(){
	int it,i,j,k;
	double gg, gam, dgg ,step, tmpd,maxGradSqr;

	lastx = LAST_X_DEFAULT;
	restart = 0;

	//initialize
	updateNebzLists();
	for (i=0;i<N;i++){
		p[3*i] = L*rx[i];
		p[3*i+1] = L*ry[i];
		p[3*i+2] = L*rz[i];
	}

	grad(p, q);
	initializeMinimizer(1);

	//main loop
	for (it = 0; it < MAX_ITERATIONS; it++){
		maxGradSqr = 0.0; gg=0.0;
		for (i=0;i<3*N;i++){
			tmpd = g[i]*g[i];
			gg += tmpd;
			if (tmpd > maxGradSqr)
				maxGradSqr = tmpd;
		}

		gtyp = sqrt( gg / (DOUBLE_N*3.0) );
		//fprintf(stderr, "\r");
		//fprintf(stderr, "iteration %d out of %d : gg = %6.3g tol = %6.3g: ", it , MAX_ITERATIONS , gg , TOL ) ;

		//printf("gtyp is %g\n",gtyp);
		if ( maxGradSqr < TOL ){
			for (i=0;i<N;i++){
				rx[i] = p[3*i]*invL;
				ry[i] = p[3*i+1]*invL;
				rz[i] = p[3*i+2]*invL;
			}
			updateNebzLists();
			calculateForces();
			return 1;
		}

		step = linmin();
		//printf("here?\n");
		if (restart)
			grad(p, q);
		if (restart){
			printf("\nRestarting optimizer\n"); 
			initializeMinimizer(0);
		}
		else{
			for ( dgg=0.0, i=0; i<3*N; i++)
				dgg += ( q[i]+g[i] )*q[i];
			gam = dgg/gg;
			if (gam > 3.0){
				//printf("gamma was too big, restarting optimizer\n");
				restart = 1;
				initializeMinimizer(0);
			}
			else{//this is my addition, since <gam> can't be too large, it fucks everything up.
				for ( tmpd=0.0, i=0; i<3*N; i++){
					g[i] = -q[i];                /* g stores (-) the most recent gradient */
					q[i] = h[i] = g[i] + gam*h[i];
					/* h stores q, the current line direction */
					/* check that the inner product of gradient and line search is < 0 */
					tmpd -= q[i]*g[i];
				}
				if ( tmpd > 0.0 ){
					restart = 2; 	/* signifies that g[j] = -q[j] is already done */
					printf("Restarting optimizer (2)\n");
					initializeMinimizer(0);
				}
			}
		}
	}
	return 0;
}
/*********************************** END OF MINIMIZER FUNCTIONS **********************************************/
/*this runs NVE and the update rules are 1. update coordinate first then update velocity by
  half time step and then calculate forces and then update velocity by another half time step*/
void advanceTimeNVE(){
	int i;
	double temp,dx,dy,dz,sum;

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

		dz = TIME_STEP*pz[i] + 0.5*TIME_STEP*TIME_STEP*fz[i];
		rzUnFolded[i] += dz*invL;
		temp = rz[i]  + dz*invL;
		if (temp >= LENGTH)
			rz[i] = temp - LENGTH;
		else if (temp < 0.0)
			rz[i] = temp + LENGTH;
		else
			rz[i] = temp;

		px[i] += HALF_TIME_STEP*fx[i];
		py[i] += HALF_TIME_STEP*fy[i];
		pz[i] += HALF_TIME_STEP*fz[i];

		temp = dx*dx + dy*dy + dz*dz;
		if (temp > maxD)
			maxD = temp;
	}
	nebListCounter++;
	if ( 2.0*((double)nebListCounter)*sqrt(maxD) > DR ){
		//fprintf(beFile,"%d\n",nebListCounter);
		updateNebzLists();
		nebListCounter = 0;
		maxD = 0.0;
	}

	calculateForces();
	kinetic = 0.0;
	for (i=0; i<N; i++){
		px[i] += HALF_TIME_STEP*fx[i];
		py[i] += HALF_TIME_STEP*fy[i];
		pz[i] += HALF_TIME_STEP*fz[i];
		kinetic += px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i];
	}
	kinetic *= 0.5;
	T = 2.0*kinetic/(3*DOUBLE_N-3.0);

	return;
}

void advanceTimeNVT(){
	int i;
	double temp,dx,dy,dz;
	double xi_T;

	for (i=0; i<N; i++){
		px[i] += HALF_TIME_STEP*fx[i];
		py[i] += HALF_TIME_STEP*fy[i];
		pz[i] += HALF_TIME_STEP*fz[i];

		dz = TIME_STEP*pz[i]; //real space displacement
		rzUnFolded[i] += dz*invL;
		temp = rz[i]  + dz*invL;
		if (temp >= LENGTH){
			rz[i] = temp - LENGTH;
		}
		else if (temp < 0.0){
			rz[i] = temp + LENGTH;
		}
		else
			rz[i] = temp;

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

		temp = dx*dx + dy*dy + dz*dz;
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
		pz[i] += HALF_TIME_STEP*fz[i];
		kinetic += px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i];
	}
	kinetic = 0.5*kinetic;

	inst_T = 2.0*kinetic/(3.0*DOUBLE_N-3.0);
	xi_T = sqrt(1.0 + TIME_STEP*(T/inst_T - 1.0)/TAU_T);
	//xi_T = sqrt(T/inst_T);
	for (i=0;i<N;i++){
		px[i] *= xi_T;
		py[i] *= xi_T;
		pz[i] *= xi_T;
	}
	return;
}

void advanceTimeNVTBrownClarke(){
	int i;
	double temp,dx,dy,dz;
	double pxTemp,pyTemp,pzTemp;
	double xi_T;

	for (i=0; i<N; i++){
		dz = TIME_STEP*pz[i]; //real space displacement
		rzUnFolded[i] += dz*invL;
		temp = rz[i]  + dz*invL;
		if (temp >= LENGTH){
			rz[i] = temp - LENGTH;
		}
		else if (temp < 0.0){
			rz[i] = temp + LENGTH;
		}
		else
			rz[i] = temp;

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

		temp = dx*dx + dy*dy + dz*dz;
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
		pzTemp = pz[i] + HALF_TIME_STEP*fz[i];
		kinetic += pxTemp*pxTemp + pyTemp*pyTemp + pzTemp*pzTemp;
	}
	kinetic = 0.5*kinetic;
	inst_T = 2.0*kinetic/(3.0*DOUBLE_N-4.0); //Notice the 4 : due to conservations
	xi_T = sqrt(T/inst_T);

	kinetic = 0.0;
	for (i=0;i<N;i++){
		px[i] = (2.0*xi_T-1.0)*px[i] + xi_T*TIME_STEP*fx[i];
		py[i] = (2.0*xi_T-1.0)*py[i] + xi_T*TIME_STEP*fy[i];
		pz[i] = (2.0*xi_T-1.0)*pz[i] + xi_T*TIME_STEP*fz[i];
		kinetic += px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i];
	}
	kinetic = 0.5*kinetic;
	inst_T = 2.0*kinetic/(3.0*DOUBLE_N-4.0); //Notice the 4 : due to conservations

	return;
}

void advanceTimeNVTWithImpurity(){
	int i;
	double temp,dx,dy,dz;
	double xi_T;

	for (i=0; i<N; i++){
		px[i] += HALF_TIME_STEP*fx[i];
		py[i] += HALF_TIME_STEP*fy[i];
		pz[i] += HALF_TIME_STEP*fz[i];

		dz = TIME_STEP*pz[i]; //real space displacement
		rzUnFolded[i] += dz*invL;
		temp = rz[i]  + dz*invL;
		if (temp >= LENGTH){
			rz[i] = temp - LENGTH;
		}
		else if (temp < 0.0){
			rz[i] = temp + LENGTH;
		}
		else
			rz[i] = temp;

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

		temp = dx*dx + dy*dy + dz*dz;
		if (temp > maxD)
			maxD = temp;
	}

	for(i=0; i<NUMBER_IMPURITY; i++){
		rx[impurityParticleIndex[i]] = rxImpurity[i];
		ry[impurityParticleIndex[i]] = ryImpurity[i];
		rz[impurityParticleIndex[i]] = rzImpurity[i];
		rxUnFolded[impurityParticleIndex[i]] = rxImpurity[i];
		ryUnFolded[impurityParticleIndex[i]] = ryImpurity[i];
		rzUnFolded[impurityParticleIndex[i]] = rzImpurity[i];

	} 	

	nebListCounter++;
	if ( 2.0*((double)nebListCounter)*sqrt(maxD) > DR )
		updateNebzLists();

	calculateForces();
	for (i=0; i<N; i++){
		px[i] += HALF_TIME_STEP*fx[i];
		py[i] += HALF_TIME_STEP*fy[i];
		pz[i] += HALF_TIME_STEP*fz[i];
	}

	for(i=0; i<NUMBER_IMPURITY; i++){
		px[impurityParticleIndex[i]] = 0.0;
		py[impurityParticleIndex[i]] = 0.0;
		pz[impurityParticleIndex[i]] = 0.0;
	}

	kinetic = 0.0;
	for(i=0; i<N; i++){
		kinetic += px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i];
	} 	
	kinetic = 0.5*kinetic;
	inst_T = 2.0*kinetic/(3.0*(DOUBLE_N - NUMBER_IMPURITY)-1.0);
	xi_T = sqrt(1.0 + TIME_STEP*(T/inst_T - 1.0)/TAU_T);
	//xi_T = sqrt(T/inst_T);
	for (i=0;i<N;i++){
		px[i] *= xi_T;
		py[i] *= xi_T;
		pz[i] *= xi_T;
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
		T -= TIME_STEP*rate;
		i++;
	}
	t1 = time(0);
	printf("\n\n Time in seconds: %d\n", t1-t0);

	return;
}

int prepareSavingArray(int runLength, int numOfOrigins, double factor){
	int i,j,k;
	int linearInterval;
	int maximalInterval;
	int offset,current,index;
	char nextStepIndexFileName[128];
	FILE *outFile;

	linearInterval = runLength/numOfOrigins;
	//maximalInterval = runLength/DIVIDE_TO;
	maximalInterval = runLength;

	printf("linearInterval = %d\n",linearInterval);
	printf("maximalInterval = %d\n",maximalInterval);

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

	sprintf(nextStepIndexFileName,"nextIndex_%d_%.2f_%.3d.dat",N,T,jobNo);
	outFile = fopen(nextStepIndexFileName,"wb");
	for (i=0;i<current;i++)
		fprintf(outFile,"%d\n",nextStepIndex[i]);
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
	sprintf(restartFileName,"%.3d/restartEq_%d_%.3f_%.3d.dat",jobNo,N,T,jobNo);

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
		rzUnFolded[i] = rz[i];
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

	sprintf(restartFileName,"%.3d/restartPostEq_%d_%.3f_%.3d.dat",jobNo,N,T,jobNo);

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
		rzUnFolded[i] = rz[i];
	}
	t1 = time(0);
	printf("Post-equilibration at T = %g - Time in seconds: %d\n\n\n",T, t1-t0);

	return;
}

void run(double duration){
	int k,i,steps,between,current,start,numOfFrames;
	int t0,t1;      /*for timing purposes...*/
	char stuffFileName[128], dataFileName[128];
	FILE *stuffFile, *dataFile;
	char restartFileName[128],tempFileName[128],moveDataFile[128];
	FILE *tempFile;
	double dummy;
	int intDummy;

	sprintf(tempFileName,"%.3d/tempData.dat",jobNo);
	sprintf(restartFileName,"%.3d/restartRun_%d_%.3f_%.3d.dat",jobNo,N,T,jobNo);
	sprintf(stuffFileName,"%.3d/stuff_%d_%.3f_%.3d.dat",jobNo,N,T,jobNo);
	sprintf(dataFileName,"%.3d/data_%d_%.3f_%.3d.dat",jobNo,N,T,jobNo);

	t0 = time(0);
	steps = (int)(duration/TIME_STEP);
	between = (int)(1.0/TIME_STEP);

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
				fscanf(dataFile,"%lf %lf %lf %d",&(rx[i]),&(ry[i]),&(rz[i]),&(type[i]));
				fprintf(tempFile,"%.10g\t%.10g\t%.10g\t%d\n",rx[i],ry[i],rz[i],type[i]);
			}
		}
		fclose(dataFile);
		fclose(tempFile);
		sprintf(moveDataFile,"mv %.3d/tempData.dat %.3d/data_%d_%.3f_%.3d.dat",jobNo,jobNo,N,T,jobNo);
		system(moveDataFile);
		dataFile = fopen(dataFileName,"ab"); // lets see whether the writing is done properly
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
			fprintf(stuffFile,"%d\t%.10g\t%.10g\t%.10g\n",k,u,kinetic,inst_T);
			for (i=0;i<N;i++)
				fprintf(dataFile,"%.10g\t%.10g\t%.10g\t%d\n",rxUnFolded[i],ryUnFolded[i],rzUnFolded[i],type[i]);
			current++;
		}
		if(!(k%((int)(steps/10))))
			saveCompleteRestartState(restartFileName,k);

		if ( !(k%65536) )
			fixDrift();
	}

	printf("Last index of this run =%d\n",current);

	t1 = time(0);
	printf("data collection at T = %g - Time in seconds: %d\n\n\n",T, t1-t0);
	fclose(stuffFile);
	fclose(dataFile);
	return;
}


void selectRandomImpurity(){
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

	sprintf(impurityFileName,"%.3d/impurityFileIndex_%d_%.3d.dat",jobNo,N,jobNo);
	impurityFile = fopen(impurityFileName,"wb"); 	
	printf("number of impurity = %d\n",NUMBER_IMPURITY);
	for(i=0; i<NUMBER_IMPURITY; i++){
		impurityParticleIndex[i] = indices[i];
		fprintf(impurityFile,"%d\n",impurityParticleIndex[i]);	
		rxImpurity[i] = rx[impurityParticleIndex[i]];
		ryImpurity[i] = ry[impurityParticleIndex[i]];
		rzImpurity[i] = rz[impurityParticleIndex[i]];
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
		rzImpurity[i] = rz[impurityParticleIndex[i]];
	}
	fclose(file);
	return;
}

int main(int argc,char *argv[]){
	int runLength,numOfOrigins,nextStepIndexSize,restartFlag;
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
	jobNo = serial*NUM_TASK + TASK_ID;

	for_ran1 = -(int)(time(0)+jobNo);
	printf("random no seed = %d\n",for_ran1);
	ran1(&for_ran1); //initialize the random number
	srand( jobNo + 100 - for_ran1);
	sprintf(dir,"mkdir -p %.3d",jobNo);
	system(dir);

	duration = 300000.0;
	durationPostPin = 0;
	prodDuration = 300000.0;

	/*********** initialization ************/
	T = 0.520;
	T_org = T;
	sprintf(impurityFileName,"%.3d/impurityFileIndex_%d_%.3d.dat",jobNo,N,jobNo);

	runLength = (int)(prodDuration/TIME_STEP);
	factor = 1.4;
	numOfOrigins = 200;

	restartEq=0;
	restartPostEq=0;
	restartRun=0;
	restartFlag = 3;
	/* 0 = 0 => fresh run */
	/* 0 = 1 => restart equilibration */
	/* 0 = 2 => restart post-equilibration */
	/* 0 = 3 => restart production */

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
		nextStepIndex = (int *)malloc(sizeof(int)*nextStepIndexSize);
		prepareSavingArray(runLength,numOfOrigins,factor);
		/**************run*****************/
		run(prodDuration);
		/**********************************/
	}

	free(nextStepIndex);
	ierr = MPI_Finalize ( );

	return 0;
}
