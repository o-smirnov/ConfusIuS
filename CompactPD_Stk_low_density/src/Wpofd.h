/*
 * Wpofd.h
 *
 *  Created on: Jun 9, 2016
 *      Author: Song Chen
 */

#ifndef WPOFD_H_
#define WPOFD_H_

#include <fftw3.h>

enum dnds_type {interpolation,powerlaw};
enum LH_type {PD,Stk};

struct dNdS_model_params{
	enum dnds_type dnds_type;

	int number_of_break;

	double C;
	double alpha;
	double beta;
	double gamma;

	double S0;
	double S1;
	double Smin;
	double Smax;

	int interplot_length;
	double * interplot_pointer;


};


struct PD_params{
  
  double totalpixel;
  double d_max;
  double d_min;
  //double source_max;
  //double source_min;
  double PSFresultionFWHM;
  double pixelsize;
  double sigma_noise;
  //double flux_mean_shift;

  double * m_beam;
  int m_beam_size;

  //for pofd elements
  fftw_complex * inbetween;
  int inbetween_size;

  struct dNdS_model_params p_model;
  struct dNdS_model_params stk_model;



  };


double ComputePD_LH(int Nbins, double * DataArray,double *result, void * ParamsArray );
double ComputeStk_LH(int Nbins, double * DataArray,double *result, void * ParamsArray );

double Set_beam(void * ParamsArray);
void free_beam(void * ParamsArray);
double Compute_pofd_elements(void * ParamsArray);
void free_pofd_elements(void * ParamsArray);
#endif /* WPOFD_H_ */
