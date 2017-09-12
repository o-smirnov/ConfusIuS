//============================================================================
// Name        : Cpofd.cpp
// Author      : Song Chen
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "Wpofd.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/stat.h>
#include <sys/types.h>

#include <math.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_randist.h>



#define PI 3.1415926
#define float_error 1e-6

#define Error_I -111.1

#define amp_I 1.0
#define amp_II 1.0




int findbin(double flux, double Delta, int shift){
	return int(flux/Delta)+shift;
}

double interp( double xi, int Ndata, double* data ){

	double x[2][Ndata];
	memcpy(x, data, sizeof(double)*2*Ndata);



	if(xi>=x[0][0]&&xi<=x[0][Ndata-1]){

	gsl_interp_accel *acc
	      = gsl_interp_accel_alloc ();
	    gsl_spline *spline
	      = gsl_spline_alloc (gsl_interp_cspline, Ndata);

	    gsl_spline_init (spline, x[0], x[1], Ndata);


	       double yi = gsl_spline_eval (spline, xi, acc);
	        //printf ("%g %g\n", xi, yi);

	    gsl_spline_free (spline);
	    gsl_interp_accel_free (acc);
	    return yi;
	}else{
		printf("error: The x=%e should not beyond data[%e~%e] ! \n",xi,x[0][0],x[0][Ndata-1]);
		exit(0);
	}

}

double interp_adv( double xi, int Ndata, double* data, int kmin, int kmax ){

	double x[2][Ndata];
    // when Ndata too large this line will have segment error!

	memcpy(x, data, sizeof(double)*2*Ndata);

	if(xi>=x[0][kmin]&&xi<=x[0][kmax]){

	gsl_interp_accel *acc
	      = gsl_interp_accel_alloc ();
	    gsl_spline *spline
	      = gsl_spline_alloc (gsl_interp_cspline, kmax-kmin+1);

	    gsl_spline_init (spline, &x[0][kmin], &x[1][kmin], kmax-kmin+1);


	       double yi = gsl_spline_eval (spline, xi, acc);
	        //printf ("%g %g\n", xi, yi);

	    gsl_spline_free (spline);
	    gsl_interp_accel_free (acc);
	    return yi;
	}else{
		printf("Find error in interp_adv, input: N=%d, kmin=%d, kmax=%d\n",Ndata,kmin,kmax);
		printf("error: The x=%e should not beyond data[%e~%e] ! \n",xi,x[0][kmin],x[0][kmax]);
		return Error_I;
	}

}
double np_roll(double * in, double *out, int Nmax, int shift){

for(int i=0; i<Nmax; i++){
	if(i<Nmax-1-shift+1){
		out[i]=in[i+shift];
	}else{
		out[i]=in[i+shift-(Nmax-1)-1];
	}

}

	return 0.0;
}




double input_fw(int Mod,int Nmax, double dmin, double Delta, double * in)
{
	if (Mod==FFTW_FORWARD){
	   for(int j=0;j<Nmax;j++){

	            	   if(j==0){
	            		   in[j]=0.0+dmin;
	            	   }else {
	            		   in[j]=dmin+j*Delta;
	            	   }


	  }
	}else if(Mod==FFTW_BACKWARD){
		for(int j=0;j<Nmax;j++){


				    	   if(j==0){
				    		   in[j]=0.0+dmin;

				    	   }else if(j<0.5*Nmax+1) {
				    		   in[j]=dmin+j*1.0/(Nmax*Delta);

				    	   }else{
				    		   in[j]=-(Nmax-j)*1.0/(Nmax*Delta)-dmin;
				    	   }


				    }

	}else{
		printf("error: wrong input !\n");
		return -1.0;

	}


	return 0.0;
}




double output_res( int Nbins, double * xbin_min, double * xbin_max, double * data, double Delta, int N, double * out, double sigma_noise, double * result,double totalpixel){



	double * flux=new double[N];

	double * outshift;

	outshift =(double *) malloc(sizeof(double)*N);

	double sigma=10*sigma_noise;


	if(input_fw(FFTW_FORWARD,N,-sigma, Delta,flux)==-1.0) exit(0);

	int xshift=int((sigma/Delta));

	int shift=N-1-xshift;

	//printf("shift %d\n",shift);
	np_roll(out,outshift,N,shift);
	//printf("save mode: FORWARD !!\n");


	int Ncut= findbin(xbin_max[Nbins-1]*1e-6,Delta,xshift)+1;
//    printf("Ncut %d\n",Ncut);

	double marray[2][Ncut+2];// why +2
	//memcpy(array[0],flux, sizeof(double)*N);
	//memcpy(array[1], out, sizeof(double)*N);

//	double check_first=0.0;



	for(int i=0;i<Ncut+2;i++){
		marray[0][i]=flux[i];
		marray[1][i]=outshift[i]/(N*Delta)/amp_II;
//		check_first+=marray[1][i]*Delta;

	}


//	printf("check_first =%e  delta %e\n",check_first, check_first-1);

    int bin_split=10;
    double bin_integration_element=0.0;

    double bin_delta;
    //double chi2=0.0;
    double loglike=0.0;
//    double normalization_check=0.0;



//    printf("total pixels %d \n",totalpixel);

    int kmin,kmax;

    for(int i=0;i<Nbins;i++){


    	bin_integration_element=0.0;
    	bin_delta=(xbin_max[i]-xbin_min[i])/(bin_split);

    	if(xbin_min[i]>0){
    	  kmin=findbin(xbin_min[i]*1e-6,Delta,xshift);
    	}else{
    	  kmin=findbin(xbin_min[i]*1e-6,Delta,xshift)-1;
    	}
    	if(xbin_max[i]>0){
    	  kmax=findbin(xbin_max[i]*1e-6,Delta,xshift)+1;
    	}else{
    	  kmax=findbin(xbin_max[i]*1e-6,Delta,xshift);
    	}

    	// to make interpolt smooth, we include more data points
    	if((kmin-10)>0){
    		kmin=kmin-10;
    	}else{
    		kmin=0;
    	}

    	if((kmax+10)<Ncut+2){
    		kmax=kmax+10;
    	}else{
    		kmax=Ncut+1;
    	}


    	for(int j=0;j<bin_split;j++){
    		double inputflux=(xbin_min[i]+(0.5+j)*bin_delta)*1e-6;
    	//	printf("come in 3  inputflux %lf !\n",inputflux);
    		double temp=interp_adv( inputflux,Ncut+2,(double *)marray, kmin, kmax);

    		if((temp<Error_I+1.0)&&(temp>Error_I-1.0)){
    		    printf("error: sigma %lf xshift %d \n",sigma,xshift);
    			exit(0);
    		}else{
//    			double P_0=1.0/sqrt(2*PI)*exp(-0.5*pow(inputflux/sigma_noise,2))*exp(-nbar);
    			bin_integration_element+=temp*bin_delta*1e-6;
    		}
    	}

    	// covariance mertix ??
    	result[i]=bin_integration_element*totalpixel;


//    	normalization_check+=bin_integration_element;


    	   	 double lambda=bin_integration_element*totalpixel;

    	   	if(data[i]==0.0){

    	   		loglike+=-lambda;
    	   	      //printf("loglike= %f\n",loglike);
    	   	}else{
    	    		loglike+=data[i]*log(lambda)-lambda-data[i]*log(data[i])+data[i];
    	   	      //printf("loglike= %f\n",loglike);
    	   	}

    }



/*
    int cs=gsl_stats_max_index(data,1,Nbins);
    int aa=gsl_stats_max_index(outshift,1,N);
    printf("cs=%d   aa=%d\n ",cs,aa);

    printf("peak Hist at %e   peak theory at %e \n",(xbin_max[cs]+xbin_min[cs])*0.5*1.0E-6,flux[aa]);
    double CMB=(xbin_max[cs]+xbin_min[cs])*0.5*1e-6-flux[aa];
    printf("background intensity= %e\n ",CMB);
    printf("index %lf \n",CMB/Delta);
*/

    free(flux);
	free(outshift);
	//free(marray);//forget to free

//	printf("normalization check = %e\n",normalization_check);


	return loglike;
}


double Set_beam(void * ParamsArray){

	struct PD_params *p
				    = (struct PD_params *) ParamsArray;

	double FHWM=p->PSFresultionFWHM;
	double pixelsize_asec=p->pixelsize;

	int  b_fwhmpx = (int)FHWM/pixelsize_asec;   //BEAM FWHM IN PIXELS

	int size=b_fwhmpx*4+1;
	double * bx = new double[size]; //CREATE X ARRAY TO CALC GAUSSIAN, PROBABLY 2 OR 3 TIMES THE SIZE
	for(int i=0;i<size;i++){
		bx[i] = i-(size-1)*1.0/2;

	}
	double * beam = new double[size*size];
	p->m_beam_size=size;
	p->m_beam=(double *) malloc(sizeof(double)*p->m_beam_size*p->m_beam_size);

	printf("Beam FWHM in Pixels is %d\n",b_fwhmpx);
	printf("Beam size is: %d x %d\n",p->m_beam_size,p->m_beam_size);

	for(int i=0;i<size;i++){
		for(int j=0;j<size;j++){
	    //beam[j+i*size] = exp(-(pow(bx[i],2)+pow(bx[j],2))*1.0/(2.*pow(b_fwhmpx/2.3548,2)));
	    beam[j+i*size] = exp(-(pow(bx[i],2)+pow(bx[j],2))*1.0/(2.*pow(b_fwhmpx/2.3548,2)));
	    //it is not a normalized Gaussian 1.0/(2.*PI*pow(b_fwhmpx/2.3548,2))*
		}
	}

	memcpy(p->m_beam, beam, sizeof(double)*size*size);

	delete[] beam;
	delete[] bx;
	return 0.0;
}


double Total_number_of_stk_Intgrand(double flux, void * model_pointer){
	struct dNdS_model_params * model
			    = (struct dNdS_model_params *) model_pointer;
	double dnds=0.0;


if(model->dnds_type==interpolation){
	if((flux>model->Smin)&&(flux<=model->Smax)){
		dnds=pow(10,interp(log10(flux),model->interplot_length,(double *)model->interplot_pointer));
	}else{
		dnds=0.0;
	}
}else if(model->dnds_type==powerlaw){

	if(model->number_of_break==1){

		if((flux>model->Smin)&&(flux<=model->S0)){
			dnds=model->C*pow(flux,model->alpha);
		}else if((flux>model->S0)&&(flux<model->Smax)){
			dnds=model->C*pow(model->S0,model->alpha-model->beta)*pow(flux,model->beta);
	//	}else if ((flux>model->S1)&&(flux<model->Smax)){
	//		dnds=model->C*pow(model->S0,model->alpha-model->beta)*pow(model->S1,model->beta-model->gamma)*pow(flux,model->gamma);
		}else {
			dnds=0.0;
		}
	}else if(model->number_of_break==2){
		if((flux>model->Smin)&&(flux<=model->S0)){
			dnds=model->C*pow(flux,model->alpha);
		}else if((flux>model->S0)&&(flux<=model->S1)){
			dnds=model->C*pow(model->S0,model->alpha-model->beta)*pow(flux,model->beta);
		}else if ((flux>model->S1)&&(flux<=model->Smax)){
			dnds=model->C*pow(model->S0,model->alpha-model->beta)*pow(model->S1,model->beta-model->gamma)*pow(flux,model->gamma);
		}else {
			dnds=0.0;
		}

	}else{
		printf("error: powerlaw dnds number of break is wrong !\n");
		exit(0);
	}
}else{
	printf("error: unknown dnds model !\n");
	exit(0);
}

	
	return dnds/4.25452e10;

}

double Get_total_number_of_stk(dNdS_model_params * model){

	 double result,error;
			  gsl_integration_workspace * w
			      = gsl_integration_workspace_alloc (1000);
			  gsl_function F;
			  F.function = &Total_number_of_stk_Intgrand;
			  F.params = model;

			  gsl_integration_qags (&F, model->Smin, model->Smax, 0, 1e-4, 1000, w, &result, &error);

			  gsl_integration_workspace_free (w);

	//printf("The total number of stacking sources is  %e with error %e  \n", result*pow(3600,2), error*pow(3600,2));
	return result;

}


/*
double Get_dNdS_per_px(double flux, double px_size, dNdS_model_params * model){

	double dnds=0.0;


//	double arrary[2][8]={{log10x,-6.7,-6.3,-5.533,-4.766,-4.0,-3.25,-1.9},{log10y,15.0620,14.4342,13.4563,12.1586,10.2764,8.55,6.314}};
if(model->dnds_type==interpolation){
	if(flux>=0.0)dnds=pow(10,interp(log10(flux),model->interplot_length,(double *)model->interplot_pointer));
}else if(model->dnds_type==powerlaw){
	if((flux>model->Smin)&&(flux<=model->S0)){
		dnds=model->C*pow(flux,model->alpha);
	}else if((flux>model->S0)&&(flux<=model->S1)){
		dnds=model->C*pow(model->S0,model->alpha-model->beta)*pow(flux,model->beta);
	}else if ((flux>model->S1)&&(flux<=model->Smax)){
		dnds=model->C*pow(model->S0,model->alpha-model->beta)*pow(model->S1,model->beta-model->gamma)*pow(flux,model->gamma);
	}else {
		dnds=0.0;
	}
}else{
	printf("error: unknown dnds model !\n");
	exit(0);
}

	double pxsz_sr = pow(px_size,2)/4.25452e10;
	return dnds*pxsz_sr;

}
*/
double Get_dNdS_per_px(double flux, double px_size, dNdS_model_params * model){

	double dnds=0.0;

//	double arrary[2][8]={{log10x,-6.7,-6.3,-5.533,-4.766,-4.0,-3.25,-1.9},{log10y,15.0620,14.4342,13.4563,12.1586,10.2764,8.55,6.314}};
if(model->dnds_type==interpolation){
	
	if((flux>model->Smin)&&(flux<=model->Smax)){
		dnds=pow(10,interp(log10(flux),model->interplot_length,(double *)model->interplot_pointer));
	}else{
	  	dnds=0.0;
	}
}else if((model->dnds_type==powerlaw)&&(model->number_of_break==1)){

	
	  if((flux>model->Smin)&&(flux<=model->S0)){
		dnds=model->C*pow(flux,model->alpha);
	  }else if((flux>model->S0)&&(flux<model->Smax)){
		dnds=model->C*pow(model->S0,model->alpha-model->beta)*pow(flux,model->beta);
	  }else {
		dnds=0.0;
	  }
	
}else if((model->dnds_type==powerlaw)&&(model->number_of_break==2)){
	  if((flux>model->Smin)&&(flux<=model->S0)){
		dnds=model->C*pow(flux,model->alpha);
	  }else if((flux>model->S0)&&(flux<model->S1)){
		dnds=model->C*pow(model->S0,model->alpha-model->beta)*pow(flux,model->beta);
	  }else if ((flux>model->S1)&&(flux<model->Smax)){
		dnds=model->C*pow(model->S0,model->alpha-model->beta)*pow(model->S1,model->beta-model->gamma)*pow(flux,model->gamma);
	  }else {
		dnds=0.0;
	  
	  }
	
}else{
	printf("error: unknown dnds model !\n");
	exit(0);
}

	double pxsz_sr = pow(px_size,2)/4.25452e10;
	return dnds*pxsz_sr;

}


double Get_Rx(double x, void * params){

	struct PD_params *p
			    = (struct PD_params *) params;



	int Ntsize=p->m_beam_size*p->m_beam_size;
	double * flux=new double[Ntsize];
    double * newy=new double[Ntsize];


    int check_count=0;
	for(int i=0;i<Ntsize;i++){

		  flux[i]=x/p->m_beam[i];

		  if((flux[i]>=p->p_model.Smin)&&(flux[i]<=p->p_model.Smax)){
			   newy[i]=Get_dNdS_per_px(flux[i],p->pixelsize,&(p->p_model));
		  		}else{
		  			newy[i]=-100.0;
		  			check_count++;
		  		}
	}



	if(check_count==Ntsize){
	//	printf("error: Get_Rx( %f )  flux > source_max !  %f \n",log10(x),Get_dNdS_per_px(x,1.5667));
    //   exit(0);
		delete[] flux;
		delete[] newy;
		return 0.0;
	}


	double rx=0.0;
	for(int i=0;i<Ntsize;i++){

	if((newy[i]>0.0)&&(p->m_beam[i]>1.0e-7)) rx=rx+newy[i]/p->m_beam[i];

	}


delete[] flux;
delete[] newy;



	return rx;
}


double ComputePD_LH(int Nbins, double * DataArray,double *result,  void * ParamsArray ) {

	struct PD_params *p
			    = (struct PD_params *) ParamsArray;


	fftw_complex *inbetween;
	int N=p->inbetween_size;
	

		double *f, *in;

		fftw_plan p1;

		f =(double *) malloc(sizeof(double)*N);


		inbetween =(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*N);

		//p->inbetween =(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*N);

		in =(double *) malloc(sizeof(double)*N);


		//make a plan "estimate(guess) or measure(find an opptimal way)"
	    p1 = fftw_plan_dft_r2c_1d(N, in, inbetween, FFTW_ESTIMATE);




		//double Delta=(p->d_max)/(N-1);
		double Delta=(p->d_max)/(N-1);

	    //Set_beam(ParamsArray);
	    //it changes because of 1.0 show or not show in beam

	    if(input_fw(FFTW_FORWARD,N,0.0,Delta,f)==-1.0) exit(0);


	    double dx=Delta;


	    for (int i = 0; i < N; i++)
	    		{
	    		   //pixel min without noise
	    	       if(f[i]< p->d_min ){
	    				in[i]=0.0;
	    		   }else{

	    		       in[i]=Get_Rx(f[i],ParamsArray)*dx*amp_I;
	    		   }

			
	    		}




		fftw_execute(p1); // execute the plan


     double x[3][Nbins];
     	memcpy(x, DataArray, sizeof(double)*3*Nbins);



	double  *w, *out;

	fftw_plan p2;

	w =(double *) malloc(sizeof(double)*N);

	out =(double *) malloc(sizeof(double)*N);

	p2 = fftw_plan_dft_c2r_1d(N, inbetween, out, FFTW_ESTIMATE);
	//FFTW_MEASURE



	double nbar=inbetween[0][0];


	double noise=p->sigma_noise; //unit Jy


    if(input_fw(FFTW_BACKWARD,N,0.0,Delta,w)==-1.0) exit(0);

    double temp_real;
		  for (int i = 0; i < N; i++){
			  temp_real=inbetween[i][0];
		inbetween[i][0]=cos(inbetween[i][1]/amp_I)*exp(temp_real/amp_I-nbar/amp_I-0.5*pow(noise*w[i]*2*PI,2))*amp_II;
		inbetween[i][1]=sin(inbetween[i][1]/amp_I)*exp(temp_real/amp_I-nbar/amp_I-0.5*pow(noise*w[i]*2*PI,2))*amp_II;
				  }


	fftw_execute(p2); // execute the plan


//	printf("fftw plans finished ! \n");



	double LogLH=0.0;
	LogLH=output_res(Nbins, x[0],x[1],x[2],Delta,N,out,p->sigma_noise, result,p->totalpixel);
//	LogLH=output_res(Nbins, x[0],x[1],x[2],Delta,N,out,p->sigma_noise, result,nbar_stk);
//	LogLH=output_res(Nbins, x[0],x[1],x[2],Delta,N,out,p->sigma_noise, result,1.0);



	fftw_destroy_plan(p1);
	fftw_destroy_plan(p2);




	free(in);
	free(f);
	free(out);
	free(w);
	fftw_free(inbetween);



	return LogLH;

}



double ComputeStk_LH(int Nbins, double * DataArray,double *result,  void * ParamsArray ) {

	struct PD_params *p
			    = (struct PD_params *) ParamsArray;


     //int N=pow(2,18);
	int N=p->inbetween_size;



     double x[3][Nbins];
     	memcpy(x, DataArray, sizeof(double)*3*Nbins);


	fftw_complex *inbetween,*inbetweenStacking;

	double *f, *w, *out,*inStacking;

	fftw_plan p2,p3;

	w =(double *) malloc(sizeof(double)*N);
	f =(double *) malloc(sizeof(double)*N);


	inbetween =(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*N);
	inbetweenStacking =(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*N);

//	in =(double *) malloc(sizeof(double)*N);
	inStacking =(double *) malloc(sizeof(double)*N);

	out =(double *) malloc(sizeof(double)*N);

	//make a plan "estimate(guess) or measure(find an opptimal way)"
//    p1 = fftw_plan_dft_r2c_1d(N, in, inbetween, FFTW_ESTIMATE);
	p2 = fftw_plan_dft_c2r_1d(N, inbetween, out, FFTW_ESTIMATE);

	p3 = fftw_plan_dft_r2c_1d(N, inStacking, inbetweenStacking, FFTW_ESTIMATE);
	//FFTW_MEASURE




	//double Delta=(p->d_max)/(N-1);
	double Delta=(p->d_max)/(N-1);

    //Set_beam(ParamsArray);
    //it changes because of 1.0 show or not show in beam

    if(input_fw(FFTW_FORWARD,N,0.0,Delta,f)==-1.0) exit(0);


    double dx=Delta;


    for (int i = 0; i < N; i++)
    		{

    		  if((f[i]>p->stk_model.Smin)&&(f[i]< p->stk_model.Smax)){
    			  inStacking[i]=dx*Get_dNdS_per_px(f[i],p->pixelsize,&(p->stk_model))*pow(3600.0/p->pixelsize,2);
    		  }else{
    			  inStacking[i]=0.0;
    		  }
/*
    	       if(f[i]< p->d_min ){
    				in[i]=0.0;
    		   }else{
    		       in[i]=Get_Rx(f[i],ParamsArray)*dx*amp_I;
    		   }
*/

    		}

    printf("f[0]=%e inStacking[0]=%lf\n",f[0],inStacking[0]);
    printf("f[1]=%e inStacking[1]=%lf\n",f[1],inStacking[1]);



//	fftw_execute(p1); // execute the plan

	fftw_execute(p3); // execute the plan


	double nbar_stk=inbetweenStacking[0][0];
//	double nbar_stk=1.0;

	double nbar=p->inbetween[0][0];


//	printf("nbar=%lf nbar_stk=%lf\n",nbar, nbar_stk);

	double normalization=exp(nbar);

	double noise=p->sigma_noise; //unit Jy


    if(input_fw(FFTW_BACKWARD,N,0.0,Delta,w)==-1.0) exit(0);

		  for (int i = 0; i < N; i++){

			double  R_pd=cos(p->inbetween[i][1])*exp(p->inbetween[i][0]);
			double  I_pd=sin(p->inbetween[i][1])*exp(p->inbetween[i][0]);

			double R_stk=inbetweenStacking[i][0]/nbar_stk*exp(-0.5*pow(noise*w[i]*2*PI,2));
			double I_stk=inbetweenStacking[i][1]/nbar_stk*exp(-0.5*pow(noise*w[i]*2*PI,2));


			
			//Stacking
			inbetween[i][0]=(R_pd*R_stk-I_pd*I_stk)/normalization;
			inbetween[i][1]=(R_pd*I_stk+R_stk*I_pd)/normalization;

			//P(D)
			//inbetween[i][0]=cos(p->inbetween[i][1])*exp(p->inbetween[i][0]-0.5*pow(noise*w[i]*2*PI,2));
			//inbetween[i][1]=sin(p->inbetween[i][1])*exp(p->inbetween[i][0]-0.5*pow(noise*w[i]*2*PI,2));
		  }


	fftw_execute(p2); // execute the plan


//	printf("fftw plans finished ! \n");

/*	double unity_check=0.0;
	for(int i=0;i<N;i++){
		unity_check+=out[i]/N;///(N*Delta)
	}


	//printf("unity check =%e\n",unity_check);
*/

        double my_stk_total_number=Get_total_number_of_stk(&(p->stk_model));
	printf("nbar_stk  %lf     integral %lf  \n",nbar_stk, my_stk_total_number*pow(3600,2));

	double LogLH=0.0;
//	LogLH=output_res(Nbins, x[0],x[1],x[2],Delta,N,out,p->sigma_noise, result,p->totalpixel);
//	LogLH=output_res(Nbins, x[0],x[1],x[2],Delta,N,out,p->sigma_noise, result,nbar_stk);
	LogLH=output_res(Nbins, x[0],x[1],x[2],Delta,N,out,p->sigma_noise, result,my_stk_total_number*pow(3600,2));
//	LogLH=output_res(Nbins, x[0],x[1],x[2],Delta,N,out,p->sigma_noise, result,1.0);



//	fftw_destroy_plan(p1);
	fftw_destroy_plan(p2);
	fftw_destroy_plan(p3);




//	free(in);
	free(inStacking);
	free(f);
	free(out);
	free(w);
	fftw_free(inbetween);
	fftw_free(inbetweenStacking);



	return LogLH;

}
void free_beam(void * ParamsArray){
	struct PD_params *p
		    = (struct PD_params *) ParamsArray;

    free(p->m_beam);
    printf(" Beam freed.\n");

}
double Compute_pofd_elements(void * ParamsArray){
   printf("Computing the P(D) elements for stacking ...");

	struct PD_params *p
			    = (struct PD_params *) ParamsArray;

	fftw_complex *inbetween;
	int N=p->inbetween_size;
	

		double *f, *in;

		fftw_plan p1;

		f =(double *) malloc(sizeof(double)*N);


		inbetween =(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*N);

		p->inbetween =(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*N);

		in =(double *) malloc(sizeof(double)*N);


		//make a plan "estimate(guess) or measure(find an opptimal way)"
	    p1 = fftw_plan_dft_r2c_1d(N, in, inbetween, FFTW_ESTIMATE);




		//double Delta=(p->d_max)/(N-1);
		double Delta=(p->d_max)/(N-1);

	    //Set_beam(ParamsArray);
	    //it changes because of 1.0 show or not show in beam

	    if(input_fw(FFTW_FORWARD,N,0.0,Delta,f)==-1.0) exit(0);


	    double dx=Delta;


	    for (int i = 0; i < N; i++)
	    		{
	    		   //pixel min without noise
	    	       if(f[i]< p->d_min ){
	    				in[i]=0.0;
	    		   }else{
	    		       in[i]=Get_Rx(f[i],ParamsArray)*dx*amp_I;
	    		   }


	    		}




		fftw_execute(p1); // execute the plan

		//copy inbetween to the structure
		memcpy( p->inbetween,inbetween, sizeof(fftw_complex)*N);



		fftw_destroy_plan(p1);
		free(in);
		free(f);
		fftw_free(inbetween);

		printf("Done!\n");
		return 0.0;

}
void free_pofd_elements(void * ParamsArray){

	struct PD_params *p
		    = (struct PD_params *) ParamsArray;

	fftw_free(p->inbetween);
	printf(" P(D) elements freed.\n");
}
