/*
 * Main.cpp
 *
 *  Created on: Jun 9, 2016
 *      Author: Song Chen
 */
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <sys/stat.h>
#include <sys/types.h>

#include "Wpofd.h"



int main(void){

	int Nbins=350;//38;//350;
	double DataArray[3][Nbins];

	FILE *fpr;
	//fpr=fopen("./out/my-hist-stacking-Kmag.txt","r");
	//fpr=fopen("./out/Jon-hist-stacking.txt","r");
	//fpr=fopen("./out/my-hist-stacking-85and008_1.txt","r");
//        fpr=fopen("./out/my-hist-PD-cluster_85and002.txt","r");
	fpr=fopen("./out/my-hist-stacking-85and002_1.txt","r");
	     double temp1,temp2;
	     double temp3;
	     int Ncount=0;
	     double totalpixel=0;

		  while ((Ncount<Nbins)&&(fscanf(fpr,"%lf  %lf  %lf\n", &temp1, &temp2, &temp3)!=EOF)){

			  if(Ncount<10)printf(" Read: %e %e %e\n",temp1, temp2, temp3);


			  DataArray[0][Ncount]=temp1;
			  DataArray[1][Ncount]=temp2;
			  DataArray[2][Ncount]=temp3;
			  totalpixel+=temp3;
		      Ncount++;


		  }
		  fclose(fpr);


		  printf("Ncount =%d\n", Ncount);
		  printf("Totalpixels %e\n",totalpixel);


		 struct PD_params ParamsArray;
                 
		 ParamsArray.inbetween_size=pow(2,19); //Size of FFTW space

		 ParamsArray.d_max=0.001; //in unit  Jy
		 ParamsArray.d_min=pow(10,-8); //pixel min flux without noise

		 //source_max and source_min
		 //is only used in the case that
		 //p or stk model only constrain part of the dnds function
		 //(for instance the faint tail)


		 ParamsArray.PSFresultionFWHM=6.0;//6.0;
		 //ParamsArray.pixelsize=3600.0/10000;
		 ParamsArray.pixelsize=1; //is only relevant for the beam resultion 
		 ParamsArray.sigma_noise=16.2e-6; 

		 ParamsArray.totalpixel=(int) 12960000;//(int)totalpixel;

		 Set_beam(&ParamsArray);

		 double result[Ncount];
		 double LogLH;	


		 ParamsArray.p_model.interplot_length=10;
		 // interplot the dN/dS with {{log10x} {log10y}} per str
		 double arrary[2][ParamsArray.p_model.interplot_length]={{-8.074,-7.33,-6.73,-6.01,-4.87,-4.09,-3.55,-3.13,-2.53,-1.87},{16.7,15.84,15.1,14.12,12.24,10.68,9.462,8.602,7.424,6.267}};
		 ParamsArray.p_model.interplot_pointer=(double *) malloc(sizeof(double)*ParamsArray.p_model.interplot_length*2);
		 memcpy(ParamsArray.p_model.interplot_pointer, arrary, sizeof(double)*ParamsArray.p_model.interplot_length*2);


		 ParamsArray.stk_model.interplot_length=10;
		 ParamsArray.stk_model.interplot_pointer=(double *) malloc(sizeof(double)*ParamsArray.stk_model.interplot_length*2);
		 memcpy(ParamsArray.stk_model.interplot_pointer, arrary, sizeof(double)*ParamsArray.stk_model.interplot_length*2);


		 enum LH_type mLH_type=Stk; // choose either 'Stk' or 'PD'

		 //----------if run PD-----------------------
		 if(mLH_type==PD){
		 printf("running PofD !!\n");
		 ParamsArray.p_model.dnds_type=powerlaw;//powerlaw;
		 ParamsArray.p_model.Smax=85e-6; //model Smax NOT equal to the inject Smax
		 ParamsArray.p_model.Smin=0.02e-6; //model Smin

		 ParamsArray.p_model.number_of_break=1;
		 
		 ParamsArray.p_model.C=4000;
		 ParamsArray.p_model.S0=15e-6;
		 ParamsArray.p_model.alpha=-1.75;
		 ParamsArray.p_model.beta=-2.0;
		 ParamsArray.p_model.S1=19e-6;
		 ParamsArray.p_model.gamma=-2.5;

		 LogLH=ComputePD_LH(Ncount, (double *) DataArray,result, &ParamsArray );
		 //---------end run PD-----------------------
		 }else if(mLH_type==Stk){
		 printf("running Stacking !!\n");
		 //----------if run Stk---------------------
		 ParamsArray.p_model.dnds_type=interpolation;
		 ParamsArray.p_model.Smax=85e-6; //model Smax NOT equal to the inject Smax
		 ParamsArray.p_model.Smin=0.02e-6; //model Smin

		 ParamsArray.p_model.number_of_break=1;

		 ParamsArray.p_model.C=4000;
		 ParamsArray.p_model.S0=15e-6;
		 ParamsArray.p_model.alpha=-1.75;
		 ParamsArray.p_model.beta=-2.0;
		 ParamsArray.p_model.S1=19e-6;
		 ParamsArray.p_model.gamma=-2.5;

		 ParamsArray.stk_model.dnds_type=powerlaw;//powerlaw;		 

		 ParamsArray.stk_model.Smax=85e-6;
		 ParamsArray.stk_model.Smin=0.2e-6;

		 ParamsArray.stk_model.number_of_break=1;

		 ParamsArray.stk_model.C=pow(10,3.00);
		 ParamsArray.stk_model.S0=25e-6;
		 ParamsArray.stk_model.alpha=-1.85;
		 ParamsArray.stk_model.beta=-2.5;	 
		 //ParamsArray.stk_model.S1=19e-6;
		 //ParamsArray.stk_model.gamma=-2.5;
		


		 Compute_pofd_elements(&ParamsArray); 
		// this is function should be run just once at initialization(if confusion noise is fixed) 
		// you need to look at the old stacking code to see where this function is called


		LogLH=ComputeStk_LH(Ncount, (double *) DataArray,result, &ParamsArray );
		

 		free_pofd_elements(&ParamsArray);
		 }
 
    	 

     	printf("LogLH is %e\n",LogLH);


    	 
    	 free_beam(&ParamsArray);

    	FILE *fpw;
    	fpw=fopen("./result_low_density.txt","w");
    	for(int i=0;i<Ncount;i++){
    		fprintf(fpw,"%e   %e   %e \n",(DataArray[0][i]+DataArray[1][i])*0.5,DataArray[2][i],result[i]);
    	}
    	fclose(fpw);


	return 0;
}


