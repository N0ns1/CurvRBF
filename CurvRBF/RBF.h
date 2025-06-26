#pragma once
#include<armadillo>
#include"surfacer/ImplicitedSurfacing.h"
#include<vector>
#include<chrono>
#include<iostream>
#include<cmath>
#include<algorithm>
#include<string>

using namespace std;

enum RBF_Kernal     //core function
{
    Polyharmonic_Spline,
    Gaussian,
    C4Compact
};

enum RBF_Method
{
    BasicRBF,
    HRBF,
    HRBF_LowOrderTerm,
    Curvature_HRBF
};

class RBF
{
public:
    RBF(vector<double>& v);
    RBF(vector<double>& v, vector<double>& vn);
    RBF(vector<double>& v, vector<double>& vn, vector<double>& h);

    void Set_Kernal_Func(RBF_Kernal f);
    void Set_RBF_Method(RBF_Method m);
    void Clean_Points();

    double (*Kernal_Function)(const double x);
    double (*Kernal_Function_2p)(const double* p1, const double* p2);
    void (*Kernal_Gradient_Function_2p)(const double* p1, const double* p2, double* G);
    void (*Kernal_Hessian_Function_2p)(const double* p1, const double* p2, double* H);
    //mean curvature
    double (*Kernal_MeanCurvature_Function_2p)(const double* p1, const double* p2, const double* g1);	//(3,1)
    void (*Kernal_Gradient_MeanCurvature_Function_2p)(const double* p1, const double* p2, const double* g1, double* re);		//(3,2)
    double (*Kernal_MeanCurvature_MeanCurvature_Function_2p)(const double* p1, const double* p2, const double* g1, const double* g2);

public:
    RBF_Method method;
    vector<double> Points;
    vector<double> Gradients;
    vector<double> MeanCurvture;

public:
    arma::mat A;        //interpolation matrix

    arma::mat K00;      //base matrix
    arma::mat K01;      //gradient matrix
    arma::mat K10;      //gradient matrix ^T
    arma::mat K11;      //Hessian matrix

    arma::mat K02;      //mean curvature matrix
    arma::mat K20;      //mean curvature matrix ^T
    arma::mat K12;      //mean curvature gradient matrix
    arma::mat K21;      //mean curvature gradient matrix ^T
    arma::mat K22;      //mean curvature mean curvature matrix

    arma::mat P01;        //low order term  matrix
    arma::mat P10;        //low order term  matrix ^T

    arma::mat C;          //constraint
    arma::mat Coef;       //solution coefficient

    arma::mat Lambda_s;     //position approximation    diagonal matrix  
    arma::mat Lambda_g;     //gradient approximation    diagonal matrix 

public:
    //method
    void Basic_Solution();
    int Gradient_Iteration_Solution(double factor1, double factor2, vector<int> ids);

    //solution process
    void Unitize_Gradient();
    void Build_Linear_System();
    void Solve_Linear_System();
    static double Dist_Function(const R3Pt& in_pt);
   
    //configuration approximation matrix
    void Build_Lam_s(vector<double>& ls);
    void Build_Lam_g(vector<double>& lg);

    /*method: HRBF*/
    //base
    void Build_RBF_Matrix_K00();
    void Build_HRBF_Matrix_K01_K10();
    void Build_HRBF_Matrix_K11();

    //Basic RBF without LowOrderTerm
    void Build_RBF_Matrix_A();
    void Build_RBF_Constraint_C();
    double Get_RBF_FuncValue(const double* p);

    double Get_Duchon_Energy(); //for HRBF
    //HRBF without LowOrderTerm
    void Build_HRBF_Matrix_A();
    void Build_HRBF_Constraint_C();
    double Get_HRBF_FuncValue(const double* p);
    vector<double> Get_HRBF_Gradient(const double* v);
    double Get_HRBF_MeanCurvature(const double* p, const double* g);
    void Get_HRBF_MeanCurvature_Allpts();
    //HRBF with LowOrderTerm
    void Build_HRBF_Matrix_A_with_LowOrderTerm();
    void Build_HRBF_Low_Order_Term_Degree1_P01_P10();
    void Build_HRBF_Constraint_C_with_LowOrderTerm();
    double Get_HRBF_FuncValue_LowOrderTerm(const double* p);
    vector<double> Get_HRBF_Gradient_LowOrderTerm(const double* v);
    // mean curvature HRBF without LowOrderTerm
    void Build_MeanCurvature_HRBF_Matrix_A();
    void Build_MeanCurvature_HRBF_Constraint_C();
    void Build_MeanCurvature_HRBF_Matrix_K20_K02();
    void Build_MeanCurvature_HRBF_Matrix_K21_K12();
    void Build_MeanCurvature_HRBF_Matrix_K22();
    double Get_MeanCurvature_HRBF_FuncValue(const double* p);
    vector<double> Get_MeanCurvature_HRBF_Gradient(const double* v);
    double Get_CurvRBF_MeanCurvature(const double* p, const double* g);
    void Get_CurvRBF_MeanCurvature_Allpts();

    //off-contact
    int isoffContact;
    vector<int> offContact_list;
    void SetOffContact(vector<int> ids);

    //multi-values
    vector<int> StartIds;
    vector<double> IoValues;
    int SetMultiValuesSolution;

public:
    int num_iter;
    vector<double>finalMesh_v;
    vector<uint>finalMesh_fv;
    void SetThis();
    void Surfacing(int n_voxels);

    void ParallelCalDistFunc(vector<float>& buffer,
        int idim, int jdim, int kdim,
        double dx, double dy, double dz, double* offset);
    void Surfacing_DMC_GPU_DataDirect(const string& outpath,
        int idim, int jdim, int kdim);  // Grid spacing

    void ParallelCalDistFunc_MultiValues(vector<float>& buffer,
        int idim, int jdim, int kdim,
        double dx, double dy, double dz, double* offset, double IoValue);
    void Surfacing_DMC_GPU_DataDirect_MultiValues(const string& outpath,
        int idim, int jdim, int kdim, vector<double>& IoValues);  // Grid spacing

public:
    void Write_Surface(string fname);
    void Write_XYZN(string fname);
    void Write_XYZNC(string fname);
    void Write_Multi_XYZNC(string fname);
    void Write_info(string fname);
    void Save_Interpolation_Matrix(string& path);
    void Save_Constraint_Matrix(string& path);

public:
    double solutionTime;
    double GenerateScalarFieldTime;
    double ComputeIsosurfaceTime;

};

static RBF* s_RBF;
