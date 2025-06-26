#pragma once
#include<armadillo>
#include<math.h>
#include<iostream>
//polyharmonic spline k = |x-y|^(2c+1)

extern int Multi_c;

double Spline_Kernel(const double x);
double Spline_Kernel_2p(const double* p1, const double* p2);		//(1,1)
void Spline_Gradient_Kernel_2p(const double* p1, const double* p2, double* G);	// dk/dx (2,1)
void Spline_Hessian_Kernel_2p(const double* p1, const double* p2, double* H);	// dk/dxdx  -(2,2)
//mean curvature
double Spine_MeanCurvature_Kernal_2p(const double* p1, const double* p2, const double* g1);	//(3,1)
void Spline_Gradient_MeanCurvature_Kernal_2p(const double* p1, const double* p2, const double* g1, double* re);		//(3,2)
double Spline_MeanCurvature_MeanCurvature_Kernal_2p(const double* p1, const double* p2, const double* g1, const double* g2);//(3,3)