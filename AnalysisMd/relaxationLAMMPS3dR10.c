#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#define	N					    9842
#define	DOUBLE_N				    9842.0
#define	RHO					    0.81

/* stuff that we don't touch usually */
#define	TIME_STEP				0.005
#define	DELTA(a,b)				((a)==(b)?1.0:0.0)
#define	MAX_NEBZ				64
#define	MAX_CELLS				(N/2)
#define	X_COMP					0
#define Y_COMP	 				1
#define	Z_COMP					2
#define	DIVIDE_TO				2
#define	PARAMETER				0.09

/* constant constants... :P*******************************************************/
#define	LENGTH					1.0
#define	HALF_LENGTH				0.5
#define	TWO_PI					6.28318530717958
#define	UNI					    ((double)rand()/((double)RAND_MAX + 1.0))



/*globals*/
double L; 			/*	length of <square> box 	*/
double invL;			/* inverse length of box	*/
double V;			/*	volume			*/
double T;			/*	temperature			*/

double rx[N];		/*	x component of position	*/
double ry[N];		/*	y component of position	*/
double rz[N];
int type[N]; 		/*	for binary system	*/

int serial,origin,numOfFrames;
int *frameIndexArray;
double *allDataX;
double *allDataY;
double *allDataZ;
int *allDataType;

void readDataFiles(){
    int k,dummy,index,consecutive,j,i;
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
    frameIndexArray = (int *)malloc(sizeof(int)*numOfFrames);

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
    allDataType = (int *)malloc(sizeof(int)*numOfFrames*N);
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
    printf("Data file reading is done");
    fclose(dataFile);
    return;
}

void loadFrame(int frameToLoad, double *rx, double *ry, double *rz, int *type){
    int i,j;

    for (j=0,i=frameToLoad*N; i < (frameToLoad+1)*N; i++){
        //printf("%g\n",allDataX[i]);
        rx[j] = allDataX[i];
        ry[j] = allDataY[i];
        rz[j] = allDataZ[i];
        type[j] = allDataType[i];
        j++;
    }
    return;
}

void meanSquareDisplacement(double factor){
    FILE *msdFile;
    char msdFileName[128];
    int k,i,current,other,originNumber,offset,maximalInterval;
    int currentFrame,otherFrame,timeIndex;
    int numOfOrigins,numOfIntervals;
    double *rx_0, *ry_0, *rz_0, *rx_t, *ry_t, *rz_t;
    int *type_0, *type_t, *counters;
    double *msd, *delta_t;
    double sum,dx,dy,dz,r2;

    //count number of origins

    rx_0 = (double *)malloc(sizeof(double)*N);
    ry_0 = (double *)malloc(sizeof(double)*N);
    rz_0 = (double *)malloc(sizeof(double)*N);
    rx_t = (double *)malloc(sizeof(double)*N);
    ry_t = (double *)malloc(sizeof(double)*N);
    rz_t = (double *)malloc(sizeof(double)*N);
    type_0 = (int *)malloc(sizeof(int)*N);
    type_t = (int *)malloc(sizeof(int)*N);

    current = 0; originNumber = 0;
    do{
        while (current<numOfFrames && frameIndexArray[current] != originNumber*origin+1)
            current++;
        //printf("originNumber is %d\n",originNumber);
        if (current<numOfFrames)
            originNumber++;
    }while(current<numOfFrames);

    maximalInterval = originNumber*origin/DIVIDE_TO;
    printf("numOfIntervals = %d\n",maximalInterval);
    numOfIntervals = 10+(int)( log((double)maximalInterval)/log(factor) ); //just to make sure
    msd = (double *)malloc(sizeof(double)*numOfIntervals);
    counters = (int *)malloc(sizeof(int)*numOfIntervals);
    delta_t = (double *)malloc(sizeof(double)*numOfIntervals);
    for (k=0;k<numOfIntervals;k++){
        msd[k] = 0.0;
        counters[k] = 0;
        delta_t[k] = 0.0;
    }

    current = 0; originNumber = 0;
    do{
        currentFrame = frameIndexArray[current];
        loadFrame(current,rx_0,ry_0, rz_0,type_0);
        //printf("loaded frame number %d which is %d\n",current,frameIndexArray[current]);

        other = current;
        offset = 1; //the smallest interval
        timeIndex = 1;
        while (other < numOfFrames && offset < maximalInterval){
            otherFrame = currentFrame + offset; //this is the frame to look for

            //look for that frame
            while (other < numOfFrames && frameIndexArray[other] != otherFrame)
                other++;
            if (other < numOfFrames){
                //at this point other and current have the indices of the frames.
                //printf("comparing %d with %d  (indices %d with %d)\n",frameIndexArray[current],frameIndexArray[other],current,other);
                loadFrame(other,rx_t,ry_t,rz_t,type_t);
                for (sum = 0.0,i=0;i<N;i++){
                    dx = rx_t[i] - rx_0[i];
                    dy = ry_t[i] - ry_0[i];
                    dz = rz_t[i] - rz_0[i];
                    // mess due to periodic boundary conditions
                    //		if ( dy >= HALF_LENGTH )
                    //			dy -= LENGTH;
                    //		else if ( dy < -HALF_LENGTH )
                    //			dy += LENGTH;
                    //		if ( dx >= HALF_LENGTH )
                    //			dx -= LENGTH;
                    //		else if ( dx < -HALF_LENGTH )
                    //			dx += LENGTH;
                    //		if ( dz >= HALF_LENGTH )
                    //			dz -= LENGTH;
                    //		else if ( dz < -HALF_LENGTH )
                    //			dz += LENGTH;
                    // end of mess 
                    r2 = L*L*( dx*dx + dy*dy + dz*dz );
                    sum += r2;
                }
                sum = sum/DOUBLE_N;
                msd[timeIndex] += sum;
                counters[timeIndex]++;
                delta_t[timeIndex] = TIME_STEP*(double)(otherFrame - currentFrame);
                //advance offset and timeIndex
                if ( (int)(offset*factor) == offset )
                    offset++;
                else
                    offset = (int)(((double)offset)*factor);
                timeIndex++;
            }
        }
        //go to next origin
        originNumber++;
        while (current<numOfFrames && frameIndexArray[current] != originNumber*origin)
            current++;
        //printf("originNumber is %d\n",originNumber);
    } while (current<numOfFrames);

    sprintf(msdFileName,"msd_%d_%.4f_%.3d.dat",N,T,serial);
    msdFile = fopen(msdFileName,"wb");
    fprintf(msdFile,"0.0\t0.0\n");
    printf("numOfIntervals = %d\n",numOfIntervals);
    for (k=1;k<numOfIntervals;k++){
        if ( counters[k] ){
            fprintf(msdFile,"%.10g\t%.10g\n",delta_t[k],msd[k]/(double)counters[k]);
            printf("%.10g\t%.10g\n",delta_t[k],msd[k]/(double)counters[k]);
        }
    }
    fclose(msdFile);


    free(msd); free(delta_t); free(counters);
    free(rx_0); free(rx_t);
    free(ry_0); free(ry_t);
    free(rz_0); free(rz_t);
    free(type_0); free(type_t);

    return;
}

void overlapCorrelationFunction(double factor, double parameter ){
    FILE *cFile;
    char cFileName[128];
    int k,i,current,other,originNumber,offset,maximalInterval;
    int currentFrame,otherFrame,timeIndex,maxTimeIndex;
    int numOfOrigins,numOfIntervals;
    double *rx_0, *ry_0, *rz_0, *rx_t, *ry_t, *rz_t;
    int *type_0, *type_t, *counters;
    double *c, *delta_t, *X4;
    double sum,dx,dy,dz,r2;

    //count number of origins

    rx_0 = (double *)malloc(sizeof(double)*N);
    ry_0 = (double *)malloc(sizeof(double)*N);
    rz_0 = (double *)malloc(sizeof(double)*N);
    rx_t = (double *)malloc(sizeof(double)*N);
    ry_t = (double *)malloc(sizeof(double)*N);
    rz_t = (double *)malloc(sizeof(double)*N);
    type_0 = (int *)malloc(sizeof(int)*N);
    type_t = (int *)malloc(sizeof(int)*N);

    current = 0; originNumber = 0;
    do{
        while (current<numOfFrames && frameIndexArray[current] != originNumber*origin+1)
            current++;
        //printf("originNumber is %d\n",originNumber);
        if (current<numOfFrames)
            originNumber++;
    }while(current<numOfFrames);

    maximalInterval = originNumber*origin/DIVIDE_TO;
    numOfIntervals = 10+(int)( log((double)maximalInterval)/log(factor) ); //just to make sure
    c = (double *)malloc(sizeof(double)*numOfIntervals);
    X4 = (double *)malloc(sizeof(double)*numOfIntervals);
    counters = (int *)malloc(sizeof(int)*numOfIntervals);
    delta_t = (double *)malloc(sizeof(double)*numOfIntervals);
    for (k=0;k<numOfIntervals;k++){
        c[k] = 0.0;
        X4[k] = 0.0;
        counters[k] = 0;
        delta_t[k] = 0.0;
    }

    current = 0; originNumber = 0;
    maxTimeIndex = 0;
    do{
        currentFrame = frameIndexArray[current];
        loadFrame(current,rx_0,ry_0,rz_0,type_0);
        //printf("loaded frame number %d which is %d\n",current,frameIndexArray[current]);

        other = current;
        offset = 1; //the smallest interval
        timeIndex = 1;
        while (other < numOfFrames && offset < maximalInterval){
            otherFrame = currentFrame + offset; //this is the frame to look for

            //look for that frame
            while (other < numOfFrames && frameIndexArray[other] != otherFrame)
                other++;
            if (other < numOfFrames){
                //at this point other and current have the indices of the frames.
                //printf("comparing %d with %d  (indices %d with %d)\n",frameIndexArray[current],frameIndexArray[other],current,other);
                loadFrame(other,rx_t,ry_t,rz_t,type_t);
                for (sum = 0.0,i=0;i<N;i++){
                    dx = rx_t[i] - rx_0[i];
                    dy = ry_t[i] - ry_0[i];
                    dz = rz_t[i] - rz_0[i];
                    // mess due to periodic boundary conditions 
                    if ( dy >= HALF_LENGTH )
                        dy -= LENGTH;
                    else if ( dy < -HALF_LENGTH )
                        dy += LENGTH;
                    if ( dx >= HALF_LENGTH )
                        dx -= LENGTH;
                    else if ( dx < -HALF_LENGTH )
                        dx += LENGTH;
                    if ( dz >= HALF_LENGTH )
                        dz -= LENGTH;
                    else if ( dz < -HALF_LENGTH )
                        dz += LENGTH;
                    // end of mess 
                    r2 = L*L*( dx*dx + dy*dy + dz*dz);
                    if (r2 < parameter)
                        sum += 1.0;
                }
                sum = sum/DOUBLE_N;
                c[timeIndex] += sum;
                X4[timeIndex] += sum*sum;
                counters[timeIndex]++;
                delta_t[timeIndex] = TIME_STEP*(double)(otherFrame - currentFrame);
                //advance offset and timeIndex
                if ( (int)(offset*factor) == offset )
                    offset++;
                else
                    offset = (int)(((double)offset)*factor);
                timeIndex++;
                if(maxTimeIndex < timeIndex)
                    maxTimeIndex = timeIndex;
            }
        }
        //go to next origin
        originNumber++;
        while (current<numOfFrames && frameIndexArray[current] != originNumber*origin)
            current++;
        //printf("originNumber is %d\n",originNumber);
    } while (current<numOfFrames);


    printf("maxTimeIndex = %d\n",maxTimeIndex);

    maxTimeIndex++;
    // calculating the longest time point 
    current = 0;
    currentFrame = frameIndexArray[current];
    loadFrame(current,rx_0,ry_0, rz_0,type_0);
    other = numOfFrames-10;
    otherFrame = frameIndexArray[other];
    printf("last frame = %d\n",otherFrame);
    loadFrame(other,rx_t,ry_t,rz_t,type_t);
    for (sum = 0.0,i=0;i<N;i++){
        dx = rx_t[i] - rx_0[i];
        dy = ry_t[i] - ry_0[i];
        dz = rz_t[i] - rz_0[i];
        // mess due to periodic boundary conditions 
        if ( dy >= HALF_LENGTH )
            dy -= LENGTH;
        else if ( dy < -HALF_LENGTH )
            dy += LENGTH;
        if ( dx >= HALF_LENGTH )
            dx -= LENGTH;
        else if ( dx < -HALF_LENGTH )
            dx += LENGTH;
        if ( dz >= HALF_LENGTH )
            dz -= LENGTH;
        else if ( dz < -HALF_LENGTH )
            dz += LENGTH;

        // end of mess 
        r2 = L*L*( dx*dx + dy*dy + dz*dz);
        if (r2 < parameter)
            sum += 1.0;
    }
    sum = sum/DOUBLE_N;
    c[maxTimeIndex] += sum;
    X4[maxTimeIndex] += sum*sum;
    counters[maxTimeIndex]++;
    delta_t[maxTimeIndex] = TIME_STEP*(double)(otherFrame - currentFrame);
    printf("numOfIntervals = %d, time index = %d, time = %.10g\n",numOfIntervals,timeIndex,delta_t[maxTimeIndex]);
    printf("counter = %d\n",counters[maxTimeIndex]);
    // -----------------------------------------------------------------

    sprintf(cFileName,"relax_%d_%.4f_%.3d.dat",N,T,serial);
    cFile = fopen(cFileName,"wb");
    fprintf(cFile,"0.0\t1.0\t0.0\n");
    for (k=1;k<numOfIntervals;k++){
        if ( counters[k] )
            fprintf(cFile,"%.10g\t%.10g\t%.10g\n",delta_t[k],
                    c[k]/(double)counters[k],
                    X4[k]/(double)counters[k] - c[k]*c[k]/((double)(counters[k]*counters[k])));
    }
    fclose(cFile);


    free(c); free(X4); free(delta_t); free(counters);
    free(rx_0); free(rx_t);
    free(ry_0); free(ry_t);
    free(rz_0); free(rz_t);
    free(type_0); free(type_t);
    return;
}
void calculateGr(){
    int i,j,maxBin,index,current,originNumber;
    char snapFileName[128];
    char grFileName[128];
    FILE *grFile;
    int *hist;
    double bin,r,dx,dy,dz,r2,norm;
    double *rx_0, *ry_0, *rz_0;
    int *type_0;
    int NCC;
    double RHOCC;
    rx_0 = (double *)malloc(sizeof(double)*N);
    ry_0 = (double *)malloc(sizeof(double)*N);
    rz_0 = (double *)malloc(sizeof(double)*N);
    type_0 = (int *)malloc(sizeof(int)*N);

    bin = 0.01;
    maxBin = (int)(0.50*L/bin);
    hist = (int *)malloc(sizeof(int)*maxBin);

    for(i=0;i<maxBin;i++)
        hist[i] = 0;

    current = 0; originNumber = 0;
    loadFrame(current,rx_0,ry_0,rz_0,type_0);
    NCC=0;
    for(i=0;i<N;i++){
      if(type_0[i]==3)
          NCC++;
    }
    printf("NCC = %d\n",NCC);
    RHOCC = (double)(NCC)/(L*L*L);
    do{
        if(current<numOfFrames && frameIndexArray[current] == originNumber*origin +1){
            loadFrame(current,rx_0,ry_0,rz_0,type_0);

            for(i=0;i<N-1;i++){
                if(type_0[i]==3){
                    for(j=i+1;j<N;j++){
                        if(type_0[j]==1 || type_0[j]==2){
                            dx = rx_0[j] - rx_0[i];
                            dy = ry_0[j] - ry_0[i];
                            dz = rz_0[j] - rz_0[i];
                            if ( dx >= HALF_LENGTH )
                                dx -= LENGTH;
                            else if ( dx < -HALF_LENGTH )
                                dx += LENGTH;
                            if ( dz >= HALF_LENGTH )
                                dz -= LENGTH;
                            else if ( dz < -HALF_LENGTH )
                                dz += LENGTH;
                            if ( dy >= HALF_LENGTH )
                                dy -= LENGTH;
                            else if ( dy < -HALF_LENGTH )
                                dy += LENGTH;
                            r2 = L*L*( dx*dx + dy*dy + dz*dz );
                            if(sqrt(r2)<0.50)printf("problem in frame = %d\n",current);
                            index = (int)(sqrt(r2)/bin);
                            if( index < maxBin )
                                hist[index] += 2;
                        }
                    }
                }
            }
            originNumber++;
        }
        current++;
    }while(current<numOfFrames);

    sprintf(grFileName,"gr_%d_%.2f_%.3d.dat",N,T,serial);
    grFile = fopen(grFileName,"wb");
    //norm = 2.0*TWO_PI*RHO*bin*DOUBLE_N*(double)(originNumber);
    norm = 2.0*TWO_PI*RHOCC*bin*(double)(N-NCC)*(double)(originNumber);
    for(i=0;i<maxBin;i++){
        r = (double)((i+1)*bin);
        fprintf(grFile,"%.10g\t%.10g\n",r,((double)hist[i])/(norm*r*r) );
    }
    fclose(grFile);
    free(hist);

    free(rx_0);free(ry_0);free(rz_0);free(type_0);

    return;
}





int main(int n,char **inputStrings){
    double factor = 1.40000000;
    sscanf(inputStrings[1],"%d",&serial);
    L = pow(DOUBLE_N/RHO,1.0/3.0);
    T = 0.520; 
    invL = 1.0/L;
    printf("box size is %.3f and inverse is %.10f\n",L,invL);

    readDataFiles();

    meanSquareDisplacement(factor);
    overlapCorrelationFunction(factor,PARAMETER);
    calculateGr();
    printf("are we here?\n");


    free(allDataX);
    free(allDataY);
    free(allDataType);
    free(frameIndexArray);
    return 0;
}


