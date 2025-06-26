#include"RBF.h"
#include"kernelFunction.h"
#include"utility.h"
#include"reader.h"

// DMC_GPU Project files
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <vector>
#include <array>
#include <map>
#include ".\surfacer_Dmc_GPU\UniformGrid.h"
#include ".\surfacer_Dmc_GPU\DualMarchingCubes.h"

typedef std::chrono::high_resolution_clock Clock;

RBF::RBF(vector<double>& v):Points(v)
{
	Clean_Points();
	SetMultiValuesSolution = 0;
}

RBF::RBF(vector<double>& v, vector<double>& vn) :Points(v), Gradients(vn)
{
	Clean_Points();
	SetMultiValuesSolution = 0;
}

RBF::RBF(vector<double>& v, vector<double>& vn, vector<double>& h) :Points(v), Gradients(vn), MeanCurvture(h)
{
	Clean_Points();
	SetMultiValuesSolution = 0;
}

void RBF::Unitize_Gradient()
{
	int n = Gradients.size() / 3;

	vector<double> vec(3);

	for (int i = 0; i < n; ++i)
	{
		vec[0] = Gradients[i * 3];
		vec[1] = Gradients[i * 3 + 1];
		vec[2] = Gradients[i * 3 + 2];

		normalize<double>(vec, 3);

		Gradients[i * 3] = vec[0];
		Gradients[i * 3 + 1] = vec[1];
		Gradients[i * 3 + 2] = vec[2];
	}
}

void RBF::Set_Kernal_Func(RBF_Kernal f)
{
	switch (f)
	{
	case Polyharmonic_Spline:
		cout << "Kernal function: Polyharmonic Spline" << endl;
		Kernal_Function = Spline_Kernel;
		Kernal_Function_2p = Spline_Kernel_2p;
		Kernal_Gradient_Function_2p = Spline_Gradient_Kernel_2p;
		Kernal_Hessian_Function_2p = Spline_Hessian_Kernel_2p;
		Kernal_MeanCurvature_Function_2p = Spine_MeanCurvature_Kernal_2p;
		Kernal_Gradient_MeanCurvature_Function_2p = Spline_Gradient_MeanCurvature_Kernal_2p;
		Kernal_MeanCurvature_MeanCurvature_Function_2p = Spline_MeanCurvature_MeanCurvature_Kernal_2p;
		break;
	default:
		break;
	}
}

void RBF::Set_RBF_Method(RBF_Method m)
{
	this->method = m;
}

void RBF::Clean_Points()
{
	int n = Points.size() / 3;
	int k = 0;
	for (int i = 0; i < n - 1; ++i)
	{
		for (int j = i + 1; j < n; ++j)
		{
			if (Points[i * 3] == Points[j * 3] &&
				Points[i * 3 + 1] == Points[j * 3 + 1] &&
				Points[i * 3 + 2] == Points[j * 3 + 2])
			{
				cout << "Clean duplicate points : " << i << ' ' << j << endl;
				Points.erase(Points.begin() + j * 3, Points.begin() + j * 3 + 3);
				if (!Gradients.empty())
				{
					Gradients.erase(Gradients.begin() + j * 3, Gradients.begin() + j * 3 + 3);
				}
				--n;
				--j;
				++k;
			}
		}
	}
	cout << "Number of duplicate points : " << k << endl;
}

void RBF::Basic_Solution()
{

	Unitize_Gradient();
	Build_Linear_System();		//Get matrix A and C
	if (isoffContact == 1)	SetOffContact(offContact_list);
	auto t1 = Clock::now();
	Solve_Linear_System();
	auto t2 = Clock::now();
	solutionTime = std::chrono::nanoseconds(t2 - t1).count() / 1e9;
	cout << "Solution completed. Time: " << solutionTime << 's' << endl;
	SetThis();
}


void RBF::Build_Linear_System()
{
	cout << "Build linear system...";
	switch (method)
	{
	case BasicRBF:
		cout << "BasicRBF" << endl;
		Multi_c = 2;
		Build_RBF_Matrix_A();
		Build_RBF_Constraint_C(); 
		break;
	case HRBF_LowOrderTerm:
		cout << "Hermite RBF with Low Order Term" << endl;
		Multi_c = 1;
		Build_HRBF_Matrix_A_with_LowOrderTerm();
		Build_HRBF_Constraint_C_with_LowOrderTerm();
		break;
	case HRBF:
		cout << "Hermite RBF" << endl;
		Multi_c = 2;
		Build_HRBF_Matrix_A();
		Build_HRBF_Constraint_C();
		break;
	case Curvature_HRBF:
		cout << "Mean Curvature HRBF" << endl;
		Multi_c = 2;
		Build_MeanCurvature_HRBF_Matrix_A();
		Build_MeanCurvature_HRBF_Constraint_C();
		break;
	default:
		break;
	}
}

void RBF::Solve_Linear_System()
{
	/*cout << "Coefficient matrix is symmetric?   ";
	bool is_sym = A.is_symmetric();
	cout << is_sym << endl;*/

	cout << "Solving coefficient..." << endl;
	Coef.clear();

	Coef = arma::inv(A) * C;
	//Coef = arma::solve(A, C);
}

int RBF::Gradient_Iteration_Solution(double factor1, double factor2, vector<int> ids)
{
	
	Unitize_Gradient();
	int num = Points.size() / 3;

	//set lambda g
	vector<double> lg(3 * num, 0.0);
	for (int i : ids)
	{
		lg[i * 3] = factor1;
		lg[i * 3 + 1] = factor1;
		lg[i * 3 + 2] = factor1;
	}
	Build_Lam_g(lg);
	SetThis();
	double* pts = Points.data();
	double* gts = Gradients.data();
	int flag = 1;
	int j = 1;
	do
	{
		//Unitize_Gradient();
		flag = 1;
		cout << "Iteration: " << j << "-----------------------" << endl;
		Build_Linear_System();		//Get matrix A and C
		A.submat(num, num, 4 * num - 1, 4 * num - 1) += Lambda_g;

		if (isoffContact == 1) SetOffContact(ids);

		auto t1 = Clock::now();
		Solve_Linear_System();
		auto t2 = Clock::now();
		solutionTime = std::chrono::nanoseconds(t2 - t1).count() / 1e9;
		

		for (int i : ids)
		{
			vector<double> g = Get_MeanCurvature_HRBF_Gradient(pts + i * 3);
			double* new_gradient = g.data();
			double* old_gradient = gts + i * 3;
			cout << "intial gradient:"
				<< old_gradient[0] << ' ' << old_gradient[1] << ' ' << old_gradient[2] << endl;
			cout << "iter gradient:"
				<< new_gradient[0] << ' ' << new_gradient[1] << ' ' << new_gradient[2] << endl;

			double dotProduct = dot<double>(old_gradient, new_gradient, 3);
			double old_norm = length<double>(old_gradient, 3);
			double new_norm = length<double>(new_gradient, 3);

			double dLength = -fabs(new_norm - old_norm) / max(new_norm, old_norm);		// [-1, 0]
			double similarity = (dotProduct / (new_norm * old_norm) - 1.0) * 0.5;	// [-1, 0]

			if (fabs(similarity) > 1e-6 || fabs(dLength) > 1e-2)
			{
				gts[i * 3] = new_gradient[0];
				gts[i * 3 + 1] = new_gradient[1];
				gts[i * 3 + 2] = new_gradient[2];

				Lambda_g(i * 3, i * 3) = factor1 * (1 - exp(factor2 * similarity)) +
					factor1 * (1 - exp(factor2 * dLength));
				Lambda_g(i * 3 + 1, i * 3 + 1) = factor1 * (1 - exp(factor2 * similarity)) +
					factor1 * (1 - exp(factor2 * dLength));
				Lambda_g(i * 3 + 2, i * 3 + 2) = factor1 * (1 - exp(factor2 * similarity)) +
					factor1 * (1 - exp(factor2 * dLength));

				flag = 0;
				cout << j << ":" << similarity << "   " << dLength
					<< ' ' << Lambda_g(i * 3, i * 3) << endl;
			}
			else
			{
				cout << j << ":" << similarity << "   " << dLength
					<< ' ' << Lambda_g(i * 3, i * 3) << endl;
				cout << "Final gradient:"
					<< new_gradient[0] << ' ' << new_gradient[1] << ' ' << new_gradient[2] << endl << endl;
				Lambda_g(i * 3, i * 3) = 0;
				Lambda_g(i * 3 + 1, i * 3 + 1) = 0;
				Lambda_g(i * 3 + 2, i * 3 + 2) = 0;
			}
			
		}
		j++;

	} while (flag == 0);
	

	
	cout << "Solution completed. Time: " << solutionTime << 's' << endl;
	num_iter = j;
	return j;
}

void RBF::Build_HRBF_Matrix_A_with_LowOrderTerm()
{
	int n = Points.size() / 3;
	A.clear();
	A.set_size(4 * n + 4, 4 * n + 4);
	A.zeros();
	Build_RBF_Matrix_K00();
	Build_HRBF_Matrix_K01_K10();
	Build_HRBF_Matrix_K11();
	Build_HRBF_Low_Order_Term_Degree1_P01_P10();

	A.submat(0, 0, n - 1, n - 1) = K00;
	A.submat(0, n, n - 1, 4 * n - 1) = K01;
	A.submat(n, 0, 4 * n - 1, n - 1) = K10;
	A.submat(n, n, 4 * n - 1, 4 * n - 1) = K11;
	A.submat(0, 4 * n, 4 * n - 1, 4 * n + 3) = P01;
	A.submat(4 * n, 0, 4 * n + 3, 4 * n - 1) = P10;
}

void RBF::Build_HRBF_Matrix_A()
{
	int n = Points.size() / 3;
	A.clear();
	A.set_size(4 * n, 4 * n);
	A.zeros();
	Build_RBF_Matrix_K00();
	Build_HRBF_Matrix_K01_K10();
	Build_HRBF_Matrix_K11();


	A.submat(0, 0, n - 1, n - 1) = K00;
	A.submat(0, n, n - 1, 4 * n - 1) = K01;
	A.submat(n, 0, 4 * n - 1, n - 1) = K10;
	A.submat(n, n, 4 * n - 1, 4 * n - 1) = K11;
}


void RBF::Build_RBF_Matrix_K00()
{
	int n = Points.size() / 3;
	K00.clear();
	K00.set_size(n, n);
	double* pts = Points.data();

	for (int i = 0; i < n; ++i)
	{     
		for (int j = i; j < n; ++j)
		{
			K00(i, j) = K00(j, i) = Kernal_Function_2p(pts + i * 3, pts + j * 3);
		}
	}
}

void RBF::Build_HRBF_Matrix_K01_K10()
{
	int n = Points.size() / 3;
	K01.clear();
	K10.clear();
	K01.set_size(n, 3 * n);
	K10.set_size(3 * n, n);
	double* pts = Points.data();

	double* G = new double[3];
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			Kernal_Gradient_Function_2p(pts + i * 3, pts + j * 3, G);

			for (int k = 0; k < 3; ++k)
			{
				K01(i, j * 3 + k) = -G[k];       
				K10(i * 3 + k, j) = G[k];
			}

		}
	}
	delete[] G;
}

void RBF::Build_HRBF_Matrix_K11()
{
	int n = Points.size() / 3;
	K11.clear();
	K11.set_size(3 * n, 3 * n);
	double* pts = Points.data();

	double* H = new double[9];
	for (int i = 0; i < n; ++i)				
	{
		for (int j = i; j < n; ++j)
		{
			Kernal_Hessian_Function_2p(pts + i * 3, pts + j * 3, H);
			for (int k = 0; k < 3; ++k)
				for (int l = 0; l < 3; ++l)
					K11(k + i * 3, l + j * 3) = K11(l + j * 3, k + i * 3) = -H[k * 3 + l];
		}
	}
	delete[] H;
}

void RBF::Build_HRBF_Low_Order_Term_Degree1_P01_P10()
{
	int n = Points.size() / 3;
	P01.clear();
	P10.clear();
	P01.set_size(4 * n, 4);
	P10.set_size(4, 4 * n);
	P01.zeros();
	P10.zeros();

	for (int i = 0; i < n; ++i)
	{
		P01(i, 0) = Points[i * 3];
		P01(i, 1) = Points[i * 3 + 1];
		P01(i, 2) = Points[i * 3 + 2];
		P01(i, 3) = 1;
		P01.submat(n + i * 3, 0, n + i * 3 + 2, 2).eye();
	}
	P10 = P01.t();
}

void RBF::Build_HRBF_Constraint_C_with_LowOrderTerm()
{
	int n = Points.size() / 3;
	C.clear();
	C.set_size(4 * n + 4, 1);
	C.zeros();
	for (int i = 0; i < n; ++i)
	{
		C(n + i * 3, 0) = Gradients[i * 3];
		C(n + i * 3 + 1, 0) = Gradients[i * 3 + 1];
		C(n + i * 3 + 2, 0) = Gradients[i * 3 + 2];
	}
	if (SetMultiValuesSolution == 1)
	{
		for (int i = 0; i < StartIds.size(); ++i)
		{
			int begin = StartIds[i];
			int end;
			if (i != StartIds.size() - 1)	end = StartIds[i + 1];
			else  end = n;
			for (int j = begin; j < end; ++j)
			{
				C(j, 0) = IoValues[i];
			}
		}
	}
}

void RBF::Build_HRBF_Constraint_C()
{
	int n = Points.size() / 3;
	C.clear();
	C.set_size(4 * n, 1);
	C.zeros();
	for (int i = 0; i < n; ++i)
	{
		C(n + i * 3, 0) = Gradients[i * 3];
		C(n + i * 3 + 1, 0) = Gradients[i * 3 + 1];
		C(n + i * 3 + 2, 0) = Gradients[i * 3 + 2];
	}
	if (SetMultiValuesSolution == 1)
	{
		for (int i = 0; i < StartIds.size(); ++i)
		{
			int begin = StartIds[i];
			int end;
			if (i != StartIds.size() - 1)	end = StartIds[i + 1];
			else  end = n;
			for (int j = begin; j < end; ++j)
			{
				C(j, 0) = IoValues[i];
			}
		}
	}
}

double RBF::Get_HRBF_FuncValue_LowOrderTerm(const double* v)
{
	int n = Points.size() / 3;
	arma::mat a(1, 4 * n + 4);
	arma::mat re;
	a.zeros();
	double* pts = Points.data();

	for (int i = 0; i < n; ++i)
	{
		a(0, i) = Kernal_Function_2p(v, pts + i * 3);
		double G[3];
		Kernal_Gradient_Function_2p(v, pts + i * 3, G);
		a(0, n + i * 3) = -G[0];
		a(0, n + i * 3 + 1) = -G[1];
		a(0, n + i * 3 + 2) = -G[2];
	}
	a(0, 4 * n) = *v;
	a(0, 4 * n + 1) = *(v + 1);
	a(0, 4 * n + 2) = *(v + 2);
	a(0, 4 * n + 3) = 1;

	//re = a.submat(0, 0, 0, 4 * n - 1) * Coef;
	re = a * Coef;

	return re(0, 0);
}

double RBF::Get_HRBF_FuncValue(const double* v)
{
	int n = Points.size() / 3;
	arma::mat a(1, 4 * n);
	arma::mat re;
	a.zeros();
	double* pts = Points.data();

	for (int i = 0; i < n; ++i)
	{
		a(0, i) = Kernal_Function_2p(v, pts + i * 3);
		double G[3];
		Kernal_Gradient_Function_2p(v, pts + i * 3, G);
		a(0, n + i * 3) = -G[0];
		a(0, n + i * 3 + 1) = -G[1];
		a(0, n + i * 3 + 2) = -G[2];
	}

	re = a * Coef;

	return re(0, 0);
}

vector<double> RBF::Get_HRBF_Gradient(const double* v)
{
	int n = Points.size() / 3;
	arma::mat a(3, 4 * n);
	arma::mat re;
	a.zeros();
	double* pts = Points.data();
	for (int i = 0; i < n; ++i)
	{
		double G[3];
		Kernal_Gradient_Function_2p(v, pts + i * 3, G);
		a(0, i) = G[0];
		a(1, i) = G[1];
		a(2, i) = G[2];

		double H[9];
		Kernal_Hessian_Function_2p(v, pts + i * 3, H);
		for (int j = 0; j < 3; ++j)
		{
			for (int k = 0; k < 3; ++k)
			{
				a(j, n + i * 3 + k) = -H[j * 3 + k];
			}
		}
	}

	re = a * Coef;
	vector<double> g = { re(0,0), re(1,0), re(2,0) };

	return g;
}

vector<double> RBF::Get_HRBF_Gradient_LowOrderTerm(const double* v)
{
	int n = Points.size() / 3;
	arma::mat a(3, 4 * n + 4);
	arma::mat re;
	a.zeros();
	double* pts = Points.data();
	for (int i = 0; i < n; ++i)
	{
		double G[3];
		Kernal_Gradient_Function_2p(v, pts + i * 3, G);
		a(0, i) = G[0];
		a(1, i) = G[1];
		a(2, i) = G[2];

		double H[9];
		Kernal_Hessian_Function_2p(v, pts + i * 3, H);
		for (int j = 0; j < 3; ++j)
		{
			for (int k = 0; k < 3; ++k)
			{
				a(j, n + i * 3 + k) = -H[j * 3 + k];
			}
		}
	}

	a(0, 4 * n) = 1;
	a(1, 4 * n + 1) = 1;
	a(2, 4 * n + 2) = 1;

	re = a * Coef;
	vector<double> g = { re(0,0), re(1,0), re(2,0) };

	return g;
}

double RBF::Get_HRBF_MeanCurvature(const double* p, const double* g)
{
	double re1, re2;
	re1 = re2 = 0.0;
	double* pts = Points.data();

	int all_num = this->Points.size() / 3;
	
	for (int i = 0; i < all_num; ++i)
	{
		re1 += Kernal_MeanCurvature_Function_2p(p, pts + i * 3, g) * Coef[i];

		double b[3] = { 0, 0, 0 };
		Kernal_Gradient_MeanCurvature_Function_2p(p, pts + i * 3, g, b);
		for (int j = 0; j < 3; ++j)
		{
			re2 += Coef[all_num + i * 3 + j] * b[j];
		}
	}
	return (re1 + re2);
}

void RBF::Get_HRBF_MeanCurvature_Allpts()
{
	int num_p = Points.size() / 3;
	double* pts = Points.data();
	double* gts = Gradients.data();
	MeanCurvture.resize(num_p);
	for (int i = 0; i < num_p; ++i)
	{
		if (gts[i * 3] == 0 && gts[i * 3 + 1] == 0 && gts[i * 3 + 2] == 0)
		{
			vector<double> g = Get_HRBF_Gradient(pts + i * 3);
			gts[i * 3] = g[0];
			gts[i * 3 + 1] = g[1];
			gts[i * 3 + 2] = g[2];
		}
		MeanCurvture[i] = Get_HRBF_MeanCurvature(pts + i * 3, gts + i * 3);
	}
}

double RBF::Get_Duchon_Energy()
{
	int n = Points.size() / 3;
	arma::mat L = Coef.submat(0, 0, 4 * n - 1, 0);
	arma::mat M = A.submat(0, 0, 4 * n - 1, 4 * n - 1);
	arma::mat re = L.t() * M * L;
	return re(0, 0);
}

void RBF::SetThis()
{
	s_RBF = this;
}

double RBF::Dist_Function(const R3Pt& in_pt)
{
	switch (s_RBF->method)
	{
	case BasicRBF:
		return s_RBF->Get_RBF_FuncValue(&(in_pt[0]));
		break;
	case HRBF_LowOrderTerm:
		return s_RBF->Get_HRBF_FuncValue_LowOrderTerm(&(in_pt[0]));
		break;
	case HRBF:
		return s_RBF->Get_HRBF_FuncValue(&(in_pt[0]));
		break;
	case Curvature_HRBF:
		return s_RBF->Get_MeanCurvature_HRBF_FuncValue(&(in_pt[0]));
		break;
	default:
		break;
	}
}

void RBF::Surfacing(int n_voxels)
{
	cout << "surfacing..." << endl;
	Surfacer sf;
	ComputeIsosurfaceTime = sf.Surfacing_Implicit(Points, n_voxels, false, Dist_Function);
	sf.WriteSurface(finalMesh_v, finalMesh_fv);
}

void RBF::Write_Surface(string fname)
{
	//writeObjFile(fname,finalMesh_v,finalMesh_fv);

	writePLYFile_VF(fname, finalMesh_v, finalMesh_fv);
}

void RBF::Save_Interpolation_Matrix(string& path)
{
	ofstream outer(path, ios::out);
	for (int i = 0; i < A.n_rows; ++i)
	{
		for (int j = 0; j < A.n_cols; ++j)
		{
			outer << setprecision(16) << A(i, j);
			if (j != A.n_cols - 1)
				outer << " ";
		}
		outer << endl;
	}
}

void RBF::Save_Constraint_Matrix(string& path)
{
	ofstream outer(path, ios::out);
	for (int i = 0; i < C.n_rows; ++i)
	{
		outer << setprecision(16) << C(i, 0) << endl;
	}
}

void RBF::Build_Lam_s(vector<double>& ls)
{
	Lambda_s.clear();
	int n = Points.size() / 3;
	Lambda_s.resize(n, n);
	Lambda_s.zeros();
	for (int i = 0; i < n; ++i)
	{
		Lambda_s(i, i) = ls[i];
	}
}

void RBF::Build_Lam_g(vector<double>& lg)
{
	Lambda_g.clear();
	int n = Gradients.size();
	Lambda_g.resize(n, n);
	Lambda_g.zeros();
	for (int i = 0; i < n; ++i)
	//for (int i = 0; i < 3; ++i)
	{
		Lambda_g(i, i) = lg[i];
	}
}

void RBF::Write_XYZN(string fname)
{
	writeXYZnormal(fname, this->Points, this->Gradients, 3);
}

void RBF::Write_XYZNC(string fname)
{
	writeXYZnormalCurvature(fname, this->Points, this->Gradients, this->MeanCurvture);
}

void RBF::Write_Multi_XYZNC(string fname)
{
	writeXYZnormalCurvature_Multi(fname, this->Points, this->Gradients, this->MeanCurvture,
		this->StartIds, this->IoValues);
}

void RBF::Build_MeanCurvature_HRBF_Matrix_A()
{
	int n = Points.size() / 3;
	A.clear();
	A.set_size(5 * n, 5 * n);
	A.zeros();
	Build_RBF_Matrix_K00();
	Build_HRBF_Matrix_K01_K10();
	Build_HRBF_Matrix_K11();
	Build_MeanCurvature_HRBF_Matrix_K20_K02();
    Build_MeanCurvature_HRBF_Matrix_K21_K12();
	Build_MeanCurvature_HRBF_Matrix_K22();


	A.submat(0, 0, n - 1, n - 1) = K00;
	A.submat(0, n, n - 1, 4 * n - 1) = K01;
	A.submat(n, 0, 4 * n - 1, n - 1) = K10;
	A.submat(n, n, 4 * n - 1, 4 * n - 1) = K11;
	A.submat(0, 4 * n, n - 1, 5 * n - 1) = K02;
	A.submat(4 * n, 0, 5 * n - 1, n - 1) = K20;
	A.submat(n, 4 * n, 4 * n - 1, 5 * n - 1) = K12;
	A.submat(4 * n, n, 5 * n - 1, 4 * n - 1) = K21;
	A.submat(4 * n, 4 * n, 5 * n - 1, 5 * n - 1) = K22;

	/*string outfile2 = "D:\\Files\\Research_Group\\Data\\GeoData\\jinchuanmodel\\cyx\\1kq\\Matrix.txt";
	ofstream outer2(outfile2, ios::out);
	for (int i = 0; i < A.n_rows; i++)
	{
		for (int j = 0; j < A.n_cols; j++)
		{
			outer2 << setprecision(16) << A(i, j);
			if (j != A.n_cols - 1)
				outer2 << " ";
		}
		outer2 << endl;
	}
	outer2.close();*/

	/*cout << "\nCoefficient matrix is symmetric?   ";
	bool is_sym = A.is_symmetric();
	cout << is_sym << endl;*/
}

void RBF::Build_MeanCurvature_HRBF_Constraint_C()
{
	int n = Points.size() / 3;
	C.clear();
	C.set_size(5 * n, 1);
	C.zeros();
	for (int i = 0; i < n; ++i)
	{
		C(n + i * 3, 0) = Gradients[i * 3];
		C(n + i * 3 + 1, 0) = Gradients[i * 3 + 1];
		C(n + i * 3 + 2, 0) = Gradients[i * 3 + 2];
		C(4 * n + i, 0) = MeanCurvture[i];
	}
	if (SetMultiValuesSolution == 1)
	{
		for (int i = 0; i < StartIds.size(); ++i)
		{
			int begin = StartIds[i];
			int end;
			if (i != StartIds.size() - 1)	end = StartIds[i + 1];
			else  end = n;
			for (int j = begin; j < end; ++j)
			{
				C(j, 0) = IoValues[i];
			}
		}
	}
}

void RBF::Build_MeanCurvature_HRBF_Matrix_K20_K02()
{
	int n = Points.size() / 3;
	K02.clear();
	K20.clear();
	K02.set_size(n, n);
	K20.set_size(n, n);
	double* pts = Points.data();
	double* gts = Gradients.data();

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			//(2,0)  act on the first x
			double re = Kernal_MeanCurvature_Function_2p(pts + i * 3, pts + j * 3, gts + i * 3); 
			K20(i, j) = K02(j, i) = re;
		}
	}
}

void RBF::Build_MeanCurvature_HRBF_Matrix_K21_K12()
{
	int n = Points.size() / 3;
	K12.clear();
	K21.clear();
	K12.set_size(3 * n, n);
	K21.set_size(n, 3 * n);
	double* pts = Points.data();
	double* gts = Gradients.data();

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			double g[3];
			//(2,1)  act on the first x
			Kernal_Gradient_MeanCurvature_Function_2p(pts + i * 3, pts + j * 3, gts + i * 3, g);
			for (int k = 0; k < 3; ++k)
			{
				K12(j * 3 + k, i) = K21(i, j * 3 + k) = g[k];
			}

		}
	}
}
void RBF::Build_MeanCurvature_HRBF_Matrix_K22()
{
	int n = Points.size() / 3;
	K22.clear();
	K22.set_size(n, n);
	double* pts = Points.data();
	double* gts = Gradients.data();

	for (int i = 0; i < n; ++i)
	{
		for (int j = i; j < n; ++j)
		{
			double re = Kernal_MeanCurvature_MeanCurvature_Function_2p(pts + i * 3, pts + j * 3, 
																	   gts + i * 3, gts + j * 3);
			K22(i, j) = K22(j, i) = re;
		}
	}
}

double RBF::Get_MeanCurvature_HRBF_FuncValue(const double* v)
{
	int n = Points.size() / 3;
	arma::mat a(1, 5 * n);
	arma::mat re;
	a.zeros();
	double* pts = Points.data();
	double* gts = Gradients.data();

	for (int i = 0; i < n; ++i)
	{
		a(0, i) = Kernal_Function_2p(v, pts + i * 3);
		double G[3];
		Kernal_Gradient_Function_2p(v, pts + i * 3, G);
		a(0, n + i * 3) = -G[0];
		a(0, n + i * 3 + 1) = -G[1];
		a(0, n + i * 3 + 2) = -G[2];
		a(0, 4 * n + i) = Kernal_MeanCurvature_Function_2p(pts + i * 3, v, gts + i * 3);
	}
	if (isoffContact == 1)
	{
		arma::uvec rows_to_remove = arma::conv_to<arma::uvec>::from(offContact_list);
		rows_to_remove = arma::sort(rows_to_remove, "descend");
		a.shed_cols(rows_to_remove);
	}
	re = a * Coef;

	return re(0, 0);
}

vector<double> RBF::Get_MeanCurvature_HRBF_Gradient(const double* v)
{
	int n = Points.size() / 3;
	arma::mat a(3, 5 * n);
	arma::mat re;
	a.zeros();
	double* pts = Points.data();
	double* gts = Gradients.data();
	for (int i = 0; i < n; ++i)
	{
		double G[3];
		Kernal_Gradient_Function_2p(v, pts + i * 3, G);
		a(0, i) = G[0];
		a(1, i) = G[1];
		a(2, i) = G[2];

		double H[9];
		Kernal_Hessian_Function_2p(v, pts + i * 3, H);
		for (int j = 0; j < 3; ++j)
		{
			for (int k = 0; k < 3; ++k)
			{
				a(j, n + i * 3 + k) = -H[j * 3 + k];
			}
		}

		double re[3];
		Kernal_Gradient_MeanCurvature_Function_2p(pts + i * 3, v, gts + i * 3, re);
		a(0, 4 * n + i) = re[0];
		a(1, 4 * n + i) = re[1];
		a(2, 4 * n + i) = re[2];
	}
	if (isoffContact == 1)
	{
		arma::uvec rows_to_remove = arma::conv_to<arma::uvec>::from(offContact_list);
		rows_to_remove = arma::sort(rows_to_remove, "descend");
		a.shed_cols(rows_to_remove);
	}

	re = a * Coef;
	vector<double> g = { re(0,0), re(1,0), re(2,0) };

	return g;
}

void RBF::Surfacing_DMC_GPU_DataDirect(const string& outpath,
	int idim, int jdim, int kdim)  
{
	auto t1 = Clock::now();
	//generate scalar field
	cout << "...Generate scalar field" << endl;
	double scale = 1.5;
	int vn = Points.size() / 3;
	double maxXYZ[3] = { Points[0], Points[1],Points[2] };
	double minXYZ[3] = { Points[0], Points[1],Points[2] };
	for (int i = 1; i < vn; ++i)
	{
		double x = Points[i * 3];
		double y = Points[i * 3 + 1];
		double z = Points[i * 3 + 2];
		if (maxXYZ[0] < x)	maxXYZ[0] = x;
		if (maxXYZ[1] < y)	maxXYZ[1] = y;
		if (maxXYZ[2] < z)	maxXYZ[2] = z;
		if (minXYZ[0] > x)	minXYZ[0] = x;
		if (minXYZ[1] > y)	minXYZ[1] = y;
		if (minXYZ[2] > z)	minXYZ[2] = z;
	}
	float X_compensation = (maxXYZ[0] - minXYZ[0]) * (scale - 1) / 2;
	float Y_compensation = (maxXYZ[1] - minXYZ[1]) * (scale - 1) / 2;
	float Z_compensation = (maxXYZ[2] - minXYZ[2]) * (scale - 1) / 2;
	minXYZ[0] = minXYZ[0] - X_compensation;	maxXYZ[0] = maxXYZ[0] + X_compensation;
	minXYZ[1] = minXYZ[1] - Y_compensation; maxXYZ[1] = maxXYZ[1] + Y_compensation;
	minXYZ[2] = minXYZ[2] - Z_compensation; maxXYZ[2] = maxXYZ[2] + Z_compensation;

	float dx = (maxXYZ[0] - minXYZ[0]) / idim;
	float dy = (maxXYZ[1] - minXYZ[1]) / jdim;
	float dz = (maxXYZ[2] - minXYZ[2]) / kdim;


	vector<float> buffer(idim * jdim * kdim);

	ParallelCalDistFunc(buffer, idim, jdim, kdim, dx, dy, dz, minXYZ);

	cout << "...Generate done" << endl;
	auto t2 = Clock::now();
	GenerateScalarFieldTime = (std::chrono::nanoseconds(t2 - t1).count() / 1e9);

	//GPU DMC
	t1 = Clock::now();
	using SurfaceCase = p_mc::UniformGrid::SurfaceCase;
	std::map<std::string, int> config;
	config["valence"] = 0; // computes vertex valences
	config["element-quality"] = 0; // computes element quality
	config["p3X3YColor"] = 1; // color based mesh simplification
	config["p3X3YOld"] = 1;  // simplify isolated elements, i.e. no neighbor with same valence pattern
	config["p3333"] = 1; // simplify vertex valence pattern 3333
	config["halfedge-datastructure"] = 0; // computes halfedge data structure for quad only mesh
	config["non-manifold"] = 0; // compute number of non-manifold edges
	//std::array<int, 3> dim{300, 300, 300};
	p_mc::DualMarchingCubes dmc;

	// Init: generate a synthetic data set
	std::cout << " ... init" << std::endl;
	//dmc.init(dim, SurfaceCase::TwoHoledTorus);
	dmc.init(idim, jdim, kdim, dx, dy, dz, buffer, minXYZ);
	const float i0 = 0.0;

	// mesh
	using Vertex = p_mc::DualMarchingCubes::Vertex;
	using Normal = p_mc::DualMarchingCubes::Normal;
	using Halfedge = p_mc::DualMarchingCubes::Halfedge;
	using HalfedgeFace = p_mc::DualMarchingCubes::HalfedgeFace;
	using HalfedgeVertex = p_mc::DualMarchingCubes::HalfedgeVertex;
	using Triangle = p_mc::DualMarchingCubes::Triangle;
	using Quadrilateral = p_mc::DualMarchingCubes::Quadrilateral;

	std::vector<Vertex> v;
	std::vector<Normal> n;
	std::vector<Triangle> t;
	std::vector<Quadrilateral> q;
	std::vector<Halfedge> h;
	std::vector<HalfedgeFace> hf;
	std::vector<HalfedgeVertex> hv;

	// compute iso-surface
	std::cout << " ... compute iso-surface" << std::endl;
	dmc.dualMC(i0, v, n, t, q, h, hf, hv, config);

	t2 = Clock::now();
	ComputeIsosurfaceTime = (std::chrono::nanoseconds(t2 - t1).count() / 1e9);

	cout << "Generate Scalar FieldTime: " << GenerateScalarFieldTime << endl;
	cout << "Compute Isosurface Time: " << ComputeIsosurfaceTime << endl;

	// write mesh in obj file format
	std::cout << " ... writing obj\n";
	std::ofstream objF;
	// write triangle mesh
	objF.open(outpath + "_tris.obj");
	if (!objF.is_open()) {
		std::cout << "ERROR: can't open output file " << std::endl;
	}
	else {
		objF << "#Dual Marching Cubes, triangle mesh\n";
		const int nr_v = v.size();
		objF << setprecision(18);
		for (int i = 0; i < nr_v; i++)
		{
			objF << "v " << v[i][0] << " " << v[i][1] << " " << v[i][2] << std::endl;
		}
		for (int i = 0; i < nr_v; i++)
		{
			objF << "vn " << n[i][0] << " " << n[i][1] << " " << n[i][2] << std::endl;
		}
		const int nr_t = t.size();
		for (int i = 0; i < nr_t; i++)
		{
			objF << "f " << (t[i][0] + 1) << "//" << (t[i][0] + 1)
				<< " " << (t[i][1] + 1) << "//" << (t[i][1] + 1)
				<< " " << (t[i][2] + 1) << "//" << (t[i][2] + 1) << std::endl;
		}
		objF.close();
	}
	// write quad mesh
	/*objF.open(outpath + "_quad.obj");
	if (!objF.is_open()) {
		std::cout << "ERROR: can't open output file " << std::endl;
	}
	else {
		objF << "#Dual Marching Cubes, quad only mesh\n";
		const int nr_v = v.size();
		objF << setprecision(18);
		for (int i = 0; i < nr_v; i++)
		{
			objF << "v " << v[i][0] << " " << v[i][1] << " " << v[i][2] << std::endl;
		}
		for (int i = 0; i < nr_v; i++)
		{
			objF << "vn " << n[i][0] << " " << n[i][1] << " " << n[i][2] << std::endl;
		}
		const int nr_q = q.size();
		for (int i = 0; i < nr_q; i++)
		{
			objF << "f " << (q[i][0] + 1) << "//" << (q[i][0] + 1)
				<< " " << (q[i][1] + 1) << "//" << (q[i][1] + 1)
				<< " " << (q[i][2] + 1) << "//" << (q[i][2] + 1)
				<< " " << (q[i][3] + 1) << "//" << (q[i][3] + 1) << std::endl;
		}
		objF.close();
	}*/
	return;
}

void RBF::ParallelCalDistFunc(vector<float>& buffer,
	int idim, int jdim, int kdim,
	double dx, double dy, double dz, double* offset)
{
	const size_t total_size = idim * jdim * kdim;
	buffer.resize(total_size);

	const unsigned num_threads = std::max(1u, std::thread::hardware_concurrency());
	std::vector<std::thread> workers;
	workers.reserve(num_threads);

	const int base_chunk = kdim / num_threads;
	int extra_chunk = kdim % num_threads;

	auto worker_func = [&](int z_start, int z_end) {
		const size_t layer_stride = idim * jdim;
		for (int z = z_start; z < z_end; ++z) {
			const double z_coord = z * dz + offset[2];
			const size_t z_offset = z * layer_stride;

			for (int y = 0; y < jdim; ++y) {
				const double y_coord = y * dy + offset[1];
				const size_t y_offset = y * idim;

				for (int x = 0; x < idim; ++x) {
					const R3Pt pt(x * dx + offset[0], y_coord, z_coord);
					const size_t idx = z_offset + y_offset + x;
					buffer[idx] = static_cast<float>(Dist_Function(pt));
				}
			}
		}
	};

	int z_start = 0;
	for (unsigned t = 0; t < num_threads; ++t) {
		int z_end = z_start + base_chunk;
		if (t < extra_chunk) z_end += 1;

		workers.emplace_back(worker_func, z_start, z_end);
		z_start = z_end;
	}

	for (auto& th : workers) {
		if (th.joinable()) th.join();
	}
}

void RBF::ParallelCalDistFunc_MultiValues(vector<float>& buffer,
	int idim, int jdim, int kdim,
	double dx, double dy, double dz, double* offset,
	double IoValue)
{
	const size_t total_size = idim * jdim * kdim;
	buffer.resize(total_size);

	const unsigned num_threads = std::max(1u, std::thread::hardware_concurrency());
	std::vector<std::thread> workers;
	workers.reserve(num_threads);

	const int base_chunk = kdim / num_threads;
	int extra_chunk = kdim % num_threads;

	auto worker_func = [&](int z_start, int z_end) {
		const size_t layer_stride = idim * jdim;
		for (int z = z_start; z < z_end; ++z) {
			const double z_coord = z * dz + offset[2];
			const size_t z_offset = z * layer_stride;

			for (int y = 0; y < jdim; ++y) {
				const double y_coord = y * dy + offset[1];
				const size_t y_offset = y * idim;

				for (int x = 0; x < idim; ++x) {
					const R3Pt pt(x * dx + offset[0], y_coord, z_coord);
					const size_t idx = z_offset + y_offset + x;
					buffer[idx] = static_cast<float>(Dist_Function(pt) - IoValue);
				}
			}
		}
	};

	int z_start = 0;
	for (unsigned t = 0; t < num_threads; ++t) {
		int z_end = z_start + base_chunk;
		if (t < extra_chunk) z_end += 1;

		workers.emplace_back(worker_func, z_start, z_end);
		z_start = z_end;
	}

	for (auto& th : workers) {
		if (th.joinable()) th.join();
	}
}

void RBF::Surfacing_DMC_GPU_DataDirect_MultiValues(const string& outpath,
	int idim, int jdim, int kdim, vector<double>& IoValues)  // Grid spacing
{
	auto t1 = Clock::now();
	//generate scalar field
	cout << "...Generate scalar field" << endl;
	double scale = 1.00;
	int vn = Points.size() / 3;
	double maxXYZ[3] = { Points[0], Points[1],Points[2] };
	double minXYZ[3] = { Points[0], Points[1],Points[2] };
	for (int i = 1; i < vn; ++i)
	{
		double x = Points[i * 3];
		double y = Points[i * 3 + 1];
		double z = Points[i * 3 + 2];
		if (maxXYZ[0] < x)	maxXYZ[0] = x;
		if (maxXYZ[1] < y)	maxXYZ[1] = y;
		if (maxXYZ[2] < z)	maxXYZ[2] = z;
		if (minXYZ[0] > x)	minXYZ[0] = x;
		if (minXYZ[1] > y)	minXYZ[1] = y;
		if (minXYZ[2] > z)	minXYZ[2] = z;
	}
	float X_compensation = (maxXYZ[0] - minXYZ[0]) * (scale * 0.95 - 1) / 2;
	float Y_compensation = (maxXYZ[1] - minXYZ[1]) * (scale - 1) / 2;
	float Z_compensation = (maxXYZ[2] - minXYZ[2]) * (scale * 1.1 - 1) / 2;
	minXYZ[0] = minXYZ[0] - X_compensation;	maxXYZ[0] = maxXYZ[0] + X_compensation;
	minXYZ[1] = minXYZ[1] - Y_compensation; maxXYZ[1] = maxXYZ[1] + Y_compensation;
	minXYZ[2] = minXYZ[2] - Z_compensation; maxXYZ[2] = maxXYZ[2] + Z_compensation;
	float dx = (maxXYZ[0] - minXYZ[0]) / idim;
	float dy = (maxXYZ[1] - minXYZ[1]) / jdim;
	float dz = (maxXYZ[2] - minXYZ[2]) / kdim;

	vector<float> buffer(idim * jdim * kdim);

	ParallelCalDistFunc(buffer, idim, jdim, kdim, dx, dy, dz, minXYZ);

	/*string path = R"(D:\desktop\PaperExperiment_r2\diceng\scalarField\scalar.txt)";
	fstream writer;
	writer.open(path, ios::out);
	const size_t layer_stride = idim * jdim;
	for (int z = 0; z < kdim; ++z) {
		const double z_coord = z * dz + minXYZ[2];
		const size_t z_offset = z * layer_stride;

		for (int y = 0; y < jdim; ++y) {
			const double y_coord = y * dy + minXYZ[1];
			const size_t y_offset = y * idim;

			for (int x = 0; x < idim; ++x) {
				const R3Pt pt(x * dx + minXYZ[0], y_coord, z_coord);
				const size_t idx = z_offset + y_offset + x;
				writer << x * dx + minXYZ[0] << ' ' << y_coord << ' ' << z_coord << ' ' << buffer[idx] << endl;
			}
		}
	}*/
	auto t2 = Clock::now();
	GenerateScalarFieldTime = (std::chrono::nanoseconds(t2 - t1).count() / 1e9);
	cout << "...Generate done" << endl;

	//IoValues.push_back(IoValues[2] - 100.0);
	IoValues.push_back(IoValues[5] - 100.0);
	IoValues.push_back(IoValues[6] - 150.0);
	for (double Value : IoValues)
	{
		//GPU DMC
		t1 = Clock::now();
		using SurfaceCase = p_mc::UniformGrid::SurfaceCase;
		std::map<std::string, int> config;
		config["valence"] = 0; // computes vertex valences
		config["element-quality"] = 0; // computes element quality
		config["p3X3YColor"] = 1; // color based mesh simplification
		config["p3X3YOld"] = 1;  // simplify isolated elements, i.e. no neighbor with same valence pattern
		config["p3333"] = 1; // simplify vertex valence pattern 3333
		config["halfedge-datastructure"] = 0; // computes halfedge data structure for quad only mesh
		config["non-manifold"] = 0; // compute number of non-manifold edges
		//std::array<int, 3> dim{300, 300, 300};
		p_mc::DualMarchingCubes dmc;

		// Init: generate a synthetic data set
		std::cout << " ... init" << std::endl;
		//dmc.init(dim, SurfaceCase::TwoHoledTorus);
		dmc.init(idim, jdim, kdim, dx, dy, dz, buffer, minXYZ);
		const float i0 = Value;

		// mesh
		using Vertex = p_mc::DualMarchingCubes::Vertex;
		using Normal = p_mc::DualMarchingCubes::Normal;
		using Halfedge = p_mc::DualMarchingCubes::Halfedge;
		using HalfedgeFace = p_mc::DualMarchingCubes::HalfedgeFace;
		using HalfedgeVertex = p_mc::DualMarchingCubes::HalfedgeVertex;
		using Triangle = p_mc::DualMarchingCubes::Triangle;
		using Quadrilateral = p_mc::DualMarchingCubes::Quadrilateral;

		std::vector<Vertex> v;
		std::vector<Normal> n;
		std::vector<Triangle> t;
		std::vector<Quadrilateral> q;
		std::vector<Halfedge> h;
		std::vector<HalfedgeFace> hf;
		std::vector<HalfedgeVertex> hv;

		// compute iso-surface
		std::cout << " ... compute iso-surface" << std::endl;
		dmc.dualMC(i0, v, n, t, q, h, hf, hv, config);

		t2 = Clock::now();
		ComputeIsosurfaceTime = (std::chrono::nanoseconds(t2 - t1).count() / 1e9);

		cout << "Generate Scalar FieldTime: " << GenerateScalarFieldTime << endl;
		cout << "Compute Isosurface Time: " << ComputeIsosurfaceTime << endl;

		// write mesh in obj file format
		std::cout << " ... writing obj\n";
		std::ofstream objF;
		// write triangle mesh
		objF.open(outpath  + to_string(Value) + "_tris.obj");
		if (!objF.is_open()) {
			std::cout << "ERROR: can't open output file " << std::endl;
		}
		else {
			objF << "#Dual Marching Cubes, triangle mesh\n";
			const int nr_v = v.size();
			objF << setprecision(18);
			for (int i = 0; i < nr_v; i++)
			{
				objF << "v " << v[i][0] << " " << v[i][1] << " " << v[i][2] << std::endl;
			}
			for (int i = 0; i < nr_v; i++)
			{
				objF << "vn " << n[i][0] << " " << n[i][1] << " " << n[i][2] << std::endl;
			}
			const int nr_t = t.size();
			for (int i = 0; i < nr_t; i++)
			{
				objF << "f " << (t[i][0] + 1) << "//" << (t[i][0] + 1)
					<< " " << (t[i][1] + 1) << "//" << (t[i][1] + 1)
					<< " " << (t[i][2] + 1) << "//" << (t[i][2] + 1) << std::endl;
			}
			objF.close();
		}
		// write quad mesh
		/*objF.open(outpath + "_quad.obj");
		if (!objF.is_open()) {
			std::cout << "ERROR: can't open output file " << std::endl;
		}
		else {
			objF << "#Dual Marching Cubes, quad only mesh\n";
			const int nr_v = v.size();
			objF << setprecision(18);
			for (int i = 0; i < nr_v; i++)
			{
				objF << "v " << v[i][0] << " " << v[i][1] << " " << v[i][2] << std::endl;
			}
			for (int i = 0; i < nr_v; i++)
			{
				objF << "vn " << n[i][0] << " " << n[i][1] << " " << n[i][2] << std::endl;
			}
			const int nr_q = q.size();
			for (int i = 0; i < nr_q; i++)
			{
				objF << "f " << (q[i][0] + 1) << "//" << (q[i][0] + 1)
					<< " " << (q[i][1] + 1) << "//" << (q[i][1] + 1)
					<< " " << (q[i][2] + 1) << "//" << (q[i][2] + 1)
					<< " " << (q[i][3] + 1) << "//" << (q[i][3] + 1) << std::endl;
			}
			objF.close();
		}*/
	}
	
	return;
}

double RBF::Get_CurvRBF_MeanCurvature(const double* p, const double* g)
{
	double re1, re2, re3;
	re1 = re2 = re3 = 0.0;
	double* pts = Points.data();
	double* gts = Gradients.data();

	int all_num = this->Points.size() / 3;

	for (int i = 0; i < all_num; ++i)
	{
		re1 += Kernal_MeanCurvature_Function_2p(p, pts + i * 3, g) * Coef[i];

		double b[3] = { 0, 0, 0 };
		Kernal_Gradient_MeanCurvature_Function_2p(p, pts + i * 3, g, b);
		for (int j = 0; j < 3; ++j)
		{
			re2 += Coef[all_num + i * 3 + j] * b[j];
		}

		re3 += Kernal_MeanCurvature_MeanCurvature_Function_2p(p, pts + i * 3, g, gts + i * 3) * Coef[all_num * 4 + i];
	}
	return (re1 + re2 + re3);
}

void RBF::Get_CurvRBF_MeanCurvature_Allpts()
{
	int num_p = Points.size() / 3;
	double* pts = Points.data();
	double* gts = Gradients.data();
	MeanCurvture.resize(num_p);
	for (int i = 0; i < num_p; ++i)
	{
		MeanCurvture[i] = Get_CurvRBF_MeanCurvature(pts + i * 3, gts + i * 3);
	}
}

void RBF::Write_info(string filename)
{

	ofstream outer(filename.data(), ofstream::out);
	if (!outer.good()) {
		cout << "Can not create output file " << filename << endl;
		return ;
	}
	outer << setprecision(6);
	outer << "Solution Time: " << solutionTime << endl;
	outer << "Generate Scalar Field Time: " << GenerateScalarFieldTime << endl;
	outer << "Compute Isosurface Time: " << ComputeIsosurfaceTime << endl;
	if(num_iter)	outer << "Number of gradient iteration: " << num_iter << endl;
	outer.close();
	return;

}

void RBF::SetOffContact(vector<int> ids)
{
	arma::uvec rows_to_remove = arma::conv_to<arma::uvec>::from(ids);
	rows_to_remove = arma::sort(rows_to_remove, "descend");
	A.shed_rows(rows_to_remove);
	A.shed_cols(rows_to_remove);
	C.shed_rows(rows_to_remove);
}

void RBF::Build_RBF_Matrix_A()
{
	int n = Points.size() / 3;
	A.clear();
	A.set_size(n , n);
	A.zeros();
	Build_RBF_Matrix_K00();

	A.submat(0, 0, n - 1, n - 1) = K00;
}

void RBF::Build_RBF_Constraint_C()
{
	int n = Points.size() / 3;
	C.clear();
	C.set_size(n, 1);
	C.zeros();
	if (SetMultiValuesSolution == 1)
	{
		for (int i = 0; i < StartIds.size(); ++i)
		{
			int begin = StartIds[i];
			int end;
			if (i != StartIds.size() - 1)	end = StartIds[i + 1];
			else  end = n;
			for (int j = begin; j < end; ++j)
			{
				C(j, 0) = IoValues[i];
			}
		}
	}
}

double RBF::Get_RBF_FuncValue(const double* v)
{
	int n = Points.size() / 3;
	arma::mat a(1, n);
	arma::mat re;
	a.zeros();
	double* pts = Points.data();

	for (int i = 0; i < n; ++i)
	{
		a(0, i) = Kernal_Function_2p(v, pts + i * 3);
	}

	re = a * Coef;

	return re(0, 0);
}