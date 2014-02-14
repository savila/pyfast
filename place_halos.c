//#include <Python.h>
#include <stdio.h>     
#include <stdlib.h>     
#include <time.h>



						             			   		      
void place_halos(long NHalosTot, double *HaloMass, int NCells, int *FirstHaloInCell,  int NTotPart, double *PartX, double *PartY, double *PartZ, double L, double *HaloX, double *HaloY, double *HaloZ){


	int i,j,k,lin_ijk,rnd, starthalo, endhalo;
	int *NPartPerCell, *count, *excluded;
	long ihalo,ilong;
	long **ListOfPart;
	double invL = 1./L;	


	excluded = (int *) calloc(NTotPart, sizeof(int));

	NPartPerCell = (int *) calloc(NCells*NCells*NCells,sizeof(int));
	count = (int *) calloc(NCells*NCells*NCells,sizeof(int));

	srand (time(NULL));
	fprintf(stdout,"RAND_MAX=%d\n",RAND_MAX);


	//Assign particles to grid  ------------------------------------
	
	//count partclies per node
	for (ilong=0;ilong<NTotPart;ilong++) {
		i = (int) (invL * PartX[ilong]*NCells);
		j = (int) (invL * PartY[ilong]*NCells);
		k = (int) (invL * PartZ[ilong]*NCells);
				
		lin_ijk = i+j*NCells+k*NCells*NCells;
		NPartPerCell[lin_ijk]++;
	}
	//Alloc Enough Memory	
	ListOfPart = (long **) calloc(NCells*NCells*NCells,sizeof(long *));
	for (i=0;i<NCells;i++){
		for (j=0;j<NCells;j++){
			for (k=0;k<NCells;k++){
				lin_ijk = i+j*NCells+k*NCells*NCells;
				ListOfPart[lin_ijk] = (long *) calloc(NPartPerCell[lin_ijk],sizeof(long *));
			}		
		}
	}

	for (ilong=0;ilong<NTotPart;ilong++) {
		i = (int) (invL * PartX[ilong]*NCells);
		j = (int) (invL * PartY[ilong]*NCells);
		k = (int) (invL * PartZ[ilong]*NCells);
		
		lin_ijk = i+j*NCells+k*NCells*NCells;
		ListOfPart[lin_ijk][count[lin_ijk]] = ilong;
		count[lin_ijk]++;
	}
	//-----------------------------------  Particles assigned to grid




	//Walk the grid
	//OPENMP
	for (i=0;i<NCells;i++){
		for (j=0;j<NCells;j++){
			for (k=0;k<NCells;k++){
				lin_ijk = i+j*NCells+k*NCells*NCells;
				starthalo = FirstHaloInCell[lin_ijk];

				if (lin_ijk=NCells*NCells*NCells-1)
					endhalo=NHalosTot;
				else
					endhalo = FirstHaloInCell[lin_ijk+1];
				for(ihalo = starthalo; ihalo<endhalo; ihalo++){
					//Now we actually place them!	
					do {				
						rnd = (int)  NPartPerCell[lin_ijk] * ((double)rand()/(RAND_MAX+1.0));		
					} while(excluded[ListOfPart[lin_ijk][rnd]]==1);
					HaloX[ihalo] = PartX[ListOfPart[lin_ijk][rnd]];
					HaloY[ihalo] = PartY[ListOfPart[lin_ijk][rnd]];
					HaloZ[ihalo] = PartZ[ListOfPart[lin_ijk][rnd]];
					excluded[ListOfPart[lin_ijk][rnd]]=1;
							
				}
			}		
		}
	}

	free(count); free(NPartPerCell); free(ListOfPart);
}
