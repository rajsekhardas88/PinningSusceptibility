#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>


#define	N						    9842
#define	DOUBLE_N					9842.0
#define	RHO						    0.81

/* stuff that we don't touch usually */
#define	TIME_STEP					0.005
#define	DELTA(a,b)					((a)==(b)?1.0:0.0)
#define	MAX_NEBZ					64
#define	MAX_CELLS					(N/2)
#define	X_COMP						0
#define Y_COMP 						1
#define	Z_COMP						2
#define	CUTOFF_SQRD					2.5*2.5
#define	C_0						    -0.806140903539923
#define	C_2						    0.7
#define	C_4						    -0.156300219287607
#define	N_0						    10.0
#define	DIVIDE_TO					2
#define	PARAMETER					0.09

/* constant constants... :P*******************************************************/
#define	LENGTH						1.0
#define	HALF_LENGTH					0.5

#define	TWO_PI						6.28318530717958
#define	UNI						    ((double)rand()/((double)RAND_MAX + 1.0))


/*globals*/
double L; 			/*	length of <square> box 	*/
double invL;			/* inverse length of box	*/
double V;			/*	volume			*/
double T;			/*	temperature			*/
double t;			/*	time				*/
double u;			/*	potential energy 		*/
double kinetic;		/*	kinetic energy 		*/
double sxz;
double xi_T;

double rx[N];		/*	x component of position	*/
double ry[N];		/*	y component of position	*/
double fx[N];		/*	x component of position	*/
double fy[N];		/*	y component of position	*/
long int type[N]; 		/*	for binary system	*/

long int serial,origin,numOfFrames;
long int *frameIndexArray;
double *allDataX;
double *allDataY;
double *allDataZ;
long int *allDataType;
long int *pinned;
long int *pinned_new;
long int numOfImp;

const double invSizeSqrd[2][2] = {{1.00,1.00},{1.00,1.00}}; // [0][0] = small-small; [0][1] = [1][0] = large-small; [1][1] = large-large;

/*void readDataFiles(){
  int k,dummy,index,consecutive;
  char frameIndexFileName[128];
  char dataFileName[128];
  char impFileName[128];
  FILE *frameIndexFile, *dataFile, *impFile;

  sprintf(impFileName,"impurityFileIndex_%d_%.3d.dat",N,serial);
  impFile = fopen(impFileName,"rb");
  numOfImp = 0;
  while ( fscanf(impFile,"%d",&dummy) != EOF )
  numOfImp++;
//printf("number of impurity particles %d\n",numOfImp);	
fclose(impFile);

pinned = (int *)malloc(sizeof(int)*numOfImp);
impFile = fopen(impFileName,"rb");
numOfImp = 0;
while ( fscanf(impFile,"%d",&(pinned[numOfImp])) != EOF )
numOfImp++;
printf("number of impurity particles %d\n",numOfImp);	
fclose(impFile);



sprintf(frameIndexFileName,"../nextIndex_%d_%.2f_%.3d.dat",N,T,serial);
frameIndexFile = fopen(frameIndexFileName,"rb");
numOfFrames = 0;
while ( fscanf(frameIndexFile,"%d",&dummy) != EOF )
numOfFrames++;
rewind(frameIndexFile);
frameIndexArray = ( int *)malloc(sizeof( int)*numOfFrames);

for (k=0;k<numOfFrames;k++)
fscanf(frameIndexFile,"%d",&(frameIndexArray[k]));
fclose(frameIndexFile);

index = 0;
while ( frameIndexArray[index]+1 == frameIndexArray[index+1] )
index++;
dummy = index;
consecutive = 0;
while ( consecutive < dummy ){
if ( frameIndexArray[index]+1 == frameIndexArray[index+1] )
consecutive++;
else
consecutive = 0;
index++;
}
printf("origin? %d\n", frameIndexArray[index - dummy]);
origin = frameIndexArray[index - dummy];


sprintf(dataFileName,"data_%d_%.3f_%.3d.dat",N,T,serial);
dataFile = fopen(dataFileName,"rb");
allDataX = (double *)malloc(sizeof(double)*numOfFrames*N);
allDataY = (double *)malloc(sizeof(double)*numOfFrames*N);
allDataType = (int *)malloc(sizeof(int)*numOfFrames*N);
k = 0;
while ( fscanf(dataFile,"%lf %lf %d",&(allDataX[k]),&(allDataY[k]),&(allDataType[k])) != EOF )
k++;
fclose(dataFile);
return;
}

void loadFrame(int frameToLoad, double *rx, double *ry, int *type){
int i,j;

for (j=0,i=frameToLoad*N; i < (frameToLoad+1)*N; i++){
rx[j] = allDataX[i];
ry[j] = allDataY[i];
type[j] = allDataType[i];
j++;
}
return;
}
*/

void readDataFiles(){
	long int k,dummy,index,consecutive,j,i;
	char dummyC[10];
	double dummyD;
	char frameIndexFileName[128];
	char dataFileName[128];
	FILE *frameIndexFile, *dataFile;

	//sprintf(frameIndexFileName,"nextIndex_%d_%.2f.dat",N,T);
	sprintf(frameIndexFileName,"tmp.txt");
	frameIndexFile = fopen(frameIndexFileName,"rb");
	numOfFrames = 0;
	while ( fscanf(frameIndexFile,"%d",&dummy) != EOF )
		numOfFrames++;
	numOfFrames = numOfFrames - 1;
	rewind(frameIndexFile);
	frameIndexArray = (long int *)malloc(sizeof(long int)*numOfFrames);

	for (k=0;k<numOfFrames;k++)
		fscanf(frameIndexFile,"%d",&(frameIndexArray[k]));
	fclose(frameIndexFile);

	index = 0;
	while ( frameIndexArray[index]+1 == frameIndexArray[index+1] )
		index++;
	dummy = index;
	consecutive = 0;
	while ( consecutive < dummy ){
		if ( frameIndexArray[index]+1 == frameIndexArray[index+1] )
			consecutive++;
		else
			consecutive = 0;
		index++;
	}
	printf("origin? %d and No of Frames %d\n", frameIndexArray[index - dummy],numOfFrames);
	origin = frameIndexArray[index - dummy];


	//sprintf(dataFileName,"data_%d_%.2f_%.3d.dat",N,T,serial);
	sprintf(dataFileName,"dump.liq");
	dataFile = fopen(dataFileName,"rb");
	allDataX = (double *)malloc(sizeof(double)*numOfFrames*N);
	allDataY = (double *)malloc(sizeof(double)*numOfFrames*N);
	allDataZ = (double *)malloc(sizeof(double)*numOfFrames*N);
	allDataType = (long int *)malloc(sizeof(long int)*numOfFrames*N);
	k = 0;
	for(j=0;j<numOfFrames;j++){
		fscanf(dataFile,"%s %s",&dummyC, &dummyC);
		fscanf(dataFile,"%lf",&dummyD);	
		fscanf(dataFile,"%s %s %s %s",&dummyC, &dummyC, &dummyC, &dummyC);
		fscanf(dataFile,"%lf",&dummyD);	
		fscanf(dataFile,"%s %s %s %s %s %s",&dummyC, &dummyC, &dummyC, &dummyC,&dummyC,&dummyC);
		fscanf(dataFile,"%lf %lf",&dummyD, &dummyD);	
		fscanf(dataFile,"%lf %lf",&dummyD, &dummyD);	
		fscanf(dataFile,"%lf %lf",&dummyD, &dummyD);	
		fscanf(dataFile,"%s %s %s %s %s %s",&dummyC, &dummyC, &dummyC, &dummyC,&dummyC,&dummyC);

		for(i=0;i<N;i++){
			fscanf(dataFile,"%lf %lf %lf %d",&(allDataX[k]),&(allDataY[k]),&(allDataZ[k]),&(allDataType[k]));
			allDataX[k] = (allDataX[k] + 0.00*L)*invL;	
			allDataY[k] = (allDataY[k] + 0.00*L)*invL;	
			allDataZ[k] = (allDataZ[k] + 0.00*L)*invL;	
			k++;
		}
	}
	fclose(dataFile);
	return;
}

void loadFrame(long int frameToLoad, double *rx, double *ry, double *rz, long int *type){
	long long int i,j;

	for (j=0,i=frameToLoad*N; i < (frameToLoad+1)*N; i++){
		rx[j] = allDataX[i];
		ry[j] = allDataY[i];
		rz[j] = allDataZ[i];
		type[j] = allDataType[i];
		j++;
	}
	return;
}



void gridOverlap(double factor, double parameter ){
	FILE *cFile;
	char cFileName[128];
	int k,i,current,other,originNumber,offset,maximalInterval;
	long int currentFrame,otherFrame,timeIndex,maxTimeIndex;
	long int numOfOrigins,numOfIntervals;
	double *rx_0, *ry_0,*rz_0, *rx_t, *ry_t, *rz_t, *qcBox;
	double *rx0, *ry0, *rz0, *rxt,*ryt, *rzt;
	long int *type_0, *type_t, *counters;
	double *c, *delta_t;
	double sum,sum1,dx,dy,dz,r2,invCellSize,sL,sL2,r;
	double tau,shift;
	double pinnedX,pinnedY,pinnedZ;
	long int tauFrameNo,lastFrame;
	long int j,ll0,llt;
	tau = 4.0*10850.50;
	tauFrameNo = (long int)(tau/TIME_STEP);
	current = 0;
	currentFrame = frameIndexArray[current];
	otherFrame = currentFrame + tauFrameNo; 
	other = current;
	while (other < numOfFrames && frameIndexArray[other] <= otherFrame)
		other++;
	lastFrame = other;


	sL2 = 1.60*1.60*1.60;
	sL = pow(sL2,1.0/3.0);
	long int m = (long int)((L)/sL) ;
	sL = L/m;
	sL = sL/L;
	invCellSize = 1.0/sL;
	long int m3 = m*m*m;
	long int m2 = m*m;
	double n01[m3], nt1[m3];
	long int k1,aa,bb,cc,boxN,bin;
	qcBox = (double *)malloc(sizeof(double)*m3);
	rx0 = (double *)malloc(sizeof(double)*N);
	rxt = (double *)malloc(sizeof(double)*N);
	ryt = (double *)malloc(sizeof(double)*N);
	ry0 = (double *)malloc(sizeof(double)*N);
	rz0 = (double *)malloc(sizeof(double)*N);
	rzt = (double *)malloc(sizeof(double)*N);

	rx_0 = (double *)malloc(sizeof(double)*N);
	ry_0 = (double *)malloc(sizeof(double)*N);
	rx_t = (double *)malloc(sizeof(double)*N);
	ry_t = (double *)malloc(sizeof(double)*N);
	rz_t = (double *)malloc(sizeof(double)*N);
	rz_0 = (double *)malloc(sizeof(double)*N);
	type_0 = (long int *)malloc(sizeof(long int)*N);
	type_t = (long int *)malloc(sizeof(long int)*N);
	current = 0; originNumber = 0;
	do{
		while (current<numOfFrames && frameIndexArray[current] != originNumber*origin)
			current++;
		if (current<numOfFrames)
			originNumber++;
	}while(current<numOfFrames);
	maximalInterval = originNumber*origin/DIVIDE_TO;
	numOfIntervals = 10+(int)( log((double)maximalInterval)/log(factor) ); //just to make sure
	c = (double *)malloc(sizeof(double)*numOfFrames);
	delta_t = (double *)malloc(sizeof(double)*numOfFrames);
	for (k=0;k<numOfFrames;k++){
		c[k] = 0.0;
		delta_t[k] = 0.0;
	}


	long int l,sft;
	bin = 120;
	double dr = 1.0/bin;
	// double pinnedX,pinnedY,pinnedZ;

	/* find out the ternary particles*/
	/*_______________________________________________________________________*/

	shift = 0.1;
	current = 0;
	loadFrame(current,rx_0,ry_0,rz_0,type_0);
	numOfImp = 0;
	for(i=0;i<N;i++){
		if(type_0[i]==3){
			numOfImp++;
		}
	}
	printf("Number of impFile = %d\n",numOfImp);
	pinned = (long int *)malloc(sizeof(long int)*numOfImp);
	numOfImp = 0;
	for(i=0;i<N;i++){
		if(type_0[i]==3){
			pinned[numOfImp] =i;
			numOfImp++;
		}
	}

	/*_ Check the Ternary particles which act as pinned particles *********************/
	/*__________________________________________________________________________________*/

	for(i = 0;i<m3;i++){
		qcBox[i] =0.0;

	}
	for(l=0;l<numOfImp;l++){
		for (sft =0;sft<7;sft++){
			current = 0;
			loadFrame(current,rx_0,ry_0,rz_0,type_0);
			if(sft ==0){
				for(i=0;i<N;i++){
					rx_0[i] = rx_0[i]+shift/L;
				}
			}
			if(sft ==1){
				for(i=0;i<N;i++){
					rx_0[i] = rx_0[i]-shift/L;
				}
			}
			if(sft ==2){
				for(i=0;i<N;i++){
					ry_0[i] = ry_0[i]+shift/L;
				}
			}
			if(sft ==3){
				for(i=0;i<N;i++){
					ry_0[i] = ry_0[i]-shift/L;
				}
			}
			if(sft ==4){
				for(i=0;i<N;i++){
					rx_0[i] = rx_0[i];
					ry_0[i]= ry_0[i];
					rz_0[i]= rz_0[i];
				}
			}
			if(sft ==5){
				for(i=0;i<N;i++){
					rz_0[i] = rz_0[i]-shift/L;
				}
			}
			if(sft ==6){
				for(i=0;i<N;i++){
					rz_0[i] = rz_0[i]+shift/L;
				}
			}

			for(i=0;i<N;i++){
				if (rx_0[i] >= LENGTH){
					while(rx_0[i] >= LENGTH)
						rx_0[i] = rx_0[i] - LENGTH;
				}
				else if (rx_0[i] < 0.0){
					while(rx_0[i] < 0.0)
						rx_0[i] = rx_0[i] + LENGTH;
				}

				if (ry_0[i] >= LENGTH){
					while(ry_0[i] >= LENGTH)
						ry_0[i] = ry_0[i] - LENGTH;
				}
				else if (ry_0[i] < 0.0){
					while(ry_0[i] < 0.0)
						ry_0[i] = ry_0[i] + LENGTH;
				}
				if (rz_0[i] >= LENGTH){
					while(rz_0[i] >= LENGTH)
						rz_0[i] = rz_0[i] - LENGTH;
				}
				else if (rz_0[i] < 0.0){
					while(rz_0[i] < 0.0)
						rz_0[i] = rz_0[i] + LENGTH;
				}
			}
			rx0[l] = rx_0[pinned[l]];
			ry0[l] = ry_0[pinned[l]];
			rz0[l] = rz_0[pinned[l]];
			for(i=0;i<m3;i++){
				n01[i] = 0.0;
			}
			aa = (long int)(rx0[l]*invCellSize);
			bb = (long int)(ry0[l]*invCellSize);
			cc = (long int)(rz0[l]*invCellSize);
			boxN = m2*(cc) + m*(bb) + aa;
			ll0 = boxN; 
			n01[boxN] +=1.0; 
			//sum1 = 0.0;
			//for(i=0;i<m3;i++){
			// sum1 = sum1 + n01[i];
			//}
			loadFrame(lastFrame,rx_t,ry_t,rz_t,type_t);
			if(sft ==0){
				for(i=0;i<N;i++){
					rx_t[i] = rx_t[i]+shift/L;
				}
			}
			if(sft ==1){
				for(i=0;i<N;i++){
					rx_t[i] = rx_t[i]-shift/L;
				}
			}
			if(sft ==2){
				for(i=0;i<N;i++){
					ry_t[i] = ry_t[i]+shift/L;
				}
			}
			if(sft ==3){
				for(i=0;i<N;i++){
					ry_t[i] = ry_t[i]-shift/L;
				}
			}
			if(sft ==4){
				for(i=0;i<N;i++){
					rx_t[i] = rx_t[i];
					ry_t[i]= ry_t[i];
					rz_t[i]= rz_t[i];
				}
			}
			if(sft ==5){
				for(i=0;i<N;i++){
					rz_t[i] = rz_t[i]-shift/L;
				}
			}
			if(sft ==6){
				for(i=0;i<N;i++){
					rz_t[i] = rz_t[i]+shift/L;
				}
			}
			for(i=0;i<N;i++){
				if (rx_t[i] >= LENGTH){
					while(rx_t[i] >= LENGTH)
						rx_t[i] = rx_t[i] - LENGTH;
				}
				else if (rx_t[i] < 0.0){
					while(rx_t[i] < 0.0)
						rx_t[i] = rx_t[i] + LENGTH;
				}

				if (ry_t[i] >= LENGTH){
					while(ry_t[i] >= LENGTH)
						ry_t[i] = ry_t[i] - LENGTH;
				}
				else if (ry_t[i] < 0.0){
					while(ry_t[i] < 0.0)
						ry_t[i] = ry_t[i] + LENGTH;
				}
				if (rz_t[i] >= LENGTH){
					while(rz_t[i] >= LENGTH)
						rz_t[i] = rz_t[i] - LENGTH;
				}
				else if (rz_t[i] < 0.0){
					while(rz_t[i] < 0.0)
						rz_t[i] = rz_t[i] + LENGTH;
				}
			}
			rxt[l] = rx_t[pinned[l]];
			ryt[l] = ry_t[pinned[l]];
			rzt[l] = rz_t[pinned[l]];
			for(i=0;i<m3;i++){
				nt1[i] = 0.0;
			}
			aa = (long int)(rxt[l]*invCellSize);
			bb = (long int)(ryt[l]*invCellSize);
			cc = (long int)(rzt[l]*invCellSize);
			boxN = m2*(cc) + m*(bb) + aa;
			llt = boxN;
			nt1[boxN] +=1.0; 
			qcBox[l] = qcBox[l]+ (n01[ll0]*nt1[ll0]);
		}
	}
	l = 0;
	for(i=0;i<numOfImp;i++){
		if(qcBox[i]/7.0>=0.850){
			l +=1;
		}
	}
	long int numOfImpReal = l;	
	printf("Number of real impFile = %d\n",numOfImpReal);
	pinned_new = (long int *)malloc(sizeof(long int)*numOfImpReal);
	l = 0;
	for(i=0;i<numOfImp;i++){
		if(qcBox[i]/7.0>=0.850){
			pinned_new[l] = pinned[i];
			l +=1;
		}
	}



	sL2 = 0.8*0.8*0.8;
	sL = pow(sL2,1.0/3.0);
	m = (long int)((L)/sL) ;
	sL = L/m;
	sL = sL/L;
	invCellSize = 1.0/sL;
	m3 = m*m*m;
	m2 = m*m;
	double n0[m3], nt[m3];
	shift = 0.1;
	sprintf(cFileName,"spatialOverlap_%d_%.3f_%.3d.dat",N,T,serial);
	cFile = fopen(cFileName,"wb");
	// for (sft =0;sft<5;sft++){
	for (l =0;l<numOfImpReal;l++){
		int w = 4;
		current = 0;
		loadFrame(current,rx_0,ry_0,rz_0,type_0);
		for(i=0;i<N;i++){
			if (rx_0[i] >= LENGTH){
				while(rx_0[i] >= LENGTH)
					rx_0[i] = rx_0[i] - LENGTH;
			}
			else if (rx_0[i] < 0.0){
				while(rx_0[i] < 0.0)
					rx_0[i] = rx_0[i] + LENGTH;
			}

			if (ry_0[i] >= LENGTH){
				while(ry_0[i] >= LENGTH)
					ry_0[i] = ry_0[i] - LENGTH;
			}
			else if (ry_0[i] < 0.0){
				while(ry_0[i] < 0.0)
					ry_0[i] = ry_0[i] + LENGTH;
			}
			if (rz_0[i] >= LENGTH){
				while(rz_0[i] >= LENGTH)
					rz_0[i] = rz_0[i] - LENGTH;
			}
			else if (rz_0[i] < 0.0){
				while(rz_0[i] < 0.0)
					rz_0[i] = rz_0[i] + LENGTH;
			}
		}

		for(i=0;i<m3;i++)
			n0[i] = 0.0;
		pinnedX = rx_0[pinned_new[l]];
		pinnedY = ry_0[pinned_new[l]];
		pinnedZ = rz_0[pinned_new[l]];
		k1 = 0;
		for(i=0;i<N;i++){
			if(type_t[i]!=3){
				dx = rx_0[i] - pinnedX;
				dy = ry_0[i] - pinnedY;
				dz = rz_0[i] - pinnedZ;
				while ( dx >= HALF_LENGTH ){
					dx -= LENGTH;
				}
				while ( dx < -HALF_LENGTH ){
					dx += LENGTH;
				}
				while ( dy >= HALF_LENGTH ){
					dy -= LENGTH;
				}
				while ( dy < -HALF_LENGTH ){
					dy += LENGTH;  
				}
				while ( dz >= HALF_LENGTH ){
					dz -= LENGTH;
				}
				while ( dz < -HALF_LENGTH ){
					dz += LENGTH;  
				}     
				r2 =( dx*dx + dy*dy + dz*dz);
				r = sqrt(r2);
				if(r > w*dr && r<= w*dr+sL){
					rx0[k1] = rx_0[i] ;
					ry0[k1] = ry_0[i] ;
					rz0[k1] = rz_0[i] ;
					k1 = k1 + 1;
				}
			}
		}	
		for(i = 0;i<k1;i++){
			aa = (long int)(rx0[i]*invCellSize);
			bb = (long int)(ry0[i]*invCellSize);
			cc = (long int)(rz0[i]*invCellSize);
			boxN = m2*cc + m*(bb) + aa;
			n0[boxN] +=1; 
		}



		sum1 = 0.0;
		for(i=0;i<m3;i++)
			sum1 = sum1+n0[i];
		current = 0; originNumber = 0;
		maxTimeIndex = 0;	
		currentFrame = frameIndexArray[current];
		other = current;
		for(k=1;k<(numOfFrames);k++){
			loadFrame(k,rx_t,ry_t,rz_t,type_t);
			for(i=0;i<N;i++){
				if (rx_t[i] >= LENGTH){
					while(rx_t[i] >= LENGTH)
						rx_t[i] = rx_t[i] - LENGTH;
				}
				else if (rx_t[i] < 0.0){
					while(rx_t[i] < 0.0)
						rx_t[i] = rx_t[i] + LENGTH;
				}

				if (ry_t[i] >= LENGTH){
					while(ry_t[i] >= LENGTH)
						ry_t[i] = ry_t[i] - LENGTH;
				}
				else if (ry_t[i] < 0.0){
					while(ry_t[i] < 0.0)
						ry_t[i] = ry_t[i] + LENGTH;
				}
				if (rz_t[i] >= LENGTH){
					while(rz_t[i] >= LENGTH)
						rz_t[i] = rz_t[i] - LENGTH;
				}
				else if (rz_t[i] < 0.0){
					while(rz_t[i] < 0.0)
						rz_t[i] = rz_t[i] + LENGTH;
				}
			}

			for(i=0;i<m3;i++){
				nt[i] = 0.0;
			}
			k1 = 0;
			for(i=0;i<N;i++){
				if(type_t[i]!=3){
					dx = rx_t[i] - pinnedX;
					dy = ry_t[i] - pinnedY;
					dz = rz_t[i] - pinnedZ;
					while ( dx >= HALF_LENGTH ){
						dx -= LENGTH;
					}
					while ( dx < -HALF_LENGTH ){
						dx += LENGTH;
					}
					while ( dy >= HALF_LENGTH ){
						dy -= LENGTH;
					}
					while ( dy < -HALF_LENGTH ){
						dy += LENGTH;  
					}
					while ( dz >= HALF_LENGTH ){
						dz -= LENGTH;
					}
					while ( dz < -HALF_LENGTH ){
						dz += LENGTH;  
					}                
					r2 =( dx*dx + dy*dy +dz*dz);
					r = sqrt(r2);
					if(r > w*dr && r<= w*dr+sL){
						rxt[k1] = rx_t[i] ;
						ryt[k1] = ry_t[i] ;
						rzt[k1] = rz_t[i] ;
						k1 = k1 + 1;
					}
				}
			}
			for(i=0;i<k1;i++){
				aa = (long int)(rxt[i]*invCellSize);
				bb = (long int)(ryt[i]*invCellSize);
				cc = (long int)(rzt[i]*invCellSize);
				boxN = m2*cc + m*(bb) + aa;
				if (boxN >m3){
					printf("out of bound box found %d\n",boxN);
				}
				nt[boxN] +=1; 
			}


			sum = 0.0;
			for (i=0;i<m3;i++){
				sum = sum + n0[i]*nt[i];
			}
			if(sum1>0.0)
				sum = sum/(sum1);
			c[k] = c[k]+sum;
			otherFrame = frameIndexArray[k];

			delta_t[k] = TIME_STEP*(double)(otherFrame);
		}

	}
	//}
	for(k=1;k<numOfFrames;k++){
		fprintf(cFile,"%.10g\t%.10g\n",delta_t[k],c[k]/(numOfImpReal));

	}

	fclose(cFile);free(delta_t);free(c);
	free(rx_0); free(rx_t); free(rx0); free(rxt); 
	free(ry_0); free(ry_t); free(ry0); free(ryt);
	free(rz_0); free(rz_t); free(rz0); free(rzt);
	free(type_0); free(type_t);free(qcBox);
	return;
}

int main(int n,char **inputStrings){
	double factor = 1.40000000;
	sscanf(inputStrings[1],"%d",&serial);
	L = pow(DOUBLE_N/RHO,1.0/3.0);
	invL = 1.0/L;
	T = 0.520;


	readDataFiles();

	gridOverlap(factor,PARAMETER);
	printf("are we here?\n");

	free(pinned);
	free(allDataX);
	free(allDataY);
	free(allDataType);
	free(frameIndexArray);
	return 0;
}


