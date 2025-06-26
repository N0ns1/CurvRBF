#include"kernelFunction.h"
#include"utility.h"


using namespace std;
int Multi_c = 1;
inline double Spline_Kernel(const double x)
{
	return pow(x, 2 * Multi_c + 1);
}

double Spline_Kernel_2p(const double* p1, const double* p2)
{
	double diff[3] = { *p1 - *p2, *(p1+1) - *(p2+1), *(p1+2) - *(p2+2) };
	double len = length<double>(diff, 3);
	return Spline_Kernel(len);
}

void Spline_Gradient_Kernel_2p(const double* p1, const double* p2, double* G)
{
	double diff[3] = { *p1 - *p2, *(p1 + 1) - *(p2 + 1), *(p1 + 2) - *(p2 + 2) };
	double len = length<double>(diff, 3);
	for (int i = 0; i < 3; ++i)
	{
		G[i] = (2 * Multi_c + 1) * pow(len, 2 * Multi_c - 1) * diff[i];
	}
}

void Spline_Hessian_Kernel_2p(const double* p1, const double* p2, double* H)
{
	double diff[3] = { *p1 - *p2, *(p1 + 1) - *(p2 + 1), *(p1 + 2) - *(p2 + 2) };
	double len = length<double>(diff, 3);

	if (len < 1e-8)
	{
		for (int i = 0; i < 9; ++i)
			H[i] = 0;
	}
	else
	{
		for (int i = 0; i < 3; ++i)
			for (int j = i; j < 3; ++j)
			{
				if (i == j)
				{
					H[i * 3 + j] = ((2 * Multi_c + 1) * pow(len, 2 * Multi_c - 1) +
						(4 * pow(Multi_c, 2) - 1) * pow(diff[i], 2) * pow(len, 2 * Multi_c - 3));
				}
				else
				{
					H[i * 3 + j] = H[j * 3 + i] = ((4 * pow(Multi_c, 2) - 1) * diff[i] * diff[j] * pow(len, 2 * Multi_c - 3));
				}
			}
	}
}


double Spine_MeanCurvature_Kernal_2p(const double* p1, const double* p2, const double* g1)
{
	double diff[3] = { *p1 - *p2, *(p1 + 1) - *(p2 + 1), *(p1 + 2) - *(p2 + 2) };
	double len_dist = length<double>(diff, 3);

	if (len_dist < 1e-8)
	{
		return 0;
	}
	else
	{
		double temp[3];
		double temp_re = (2 * Multi_c + 1) * pow(len_dist, 2 * Multi_c - 1);
		for (int i = 0; i < 3; ++i)
		{
			temp[i] = (temp_re +
				(4 * pow(Multi_c, 2) - 1) * pow(diff[i], 2) * pow(len_dist, 2 * Multi_c - 3));
		}

		double norm[3];
		norm[0] = *g1;
		norm[1] = *(g1 + 1);
		norm[2] = *(g1 + 2);
		double Len = length<double>(norm, 3);
		normalize<double>(norm, 3);

		arma::mat H(3, 3);
		double h[9];
		Spline_Hessian_Kernel_2p(p1, p2, h);
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				H(i, j) = h[i * 3 + j];
			}
		}
		double term1 = temp[0] + temp[1] + temp[2];

		double term2 = pow(norm[0], 2) * H(0, 0) +
			pow(norm[1], 2) * H(1, 1) +
			pow(norm[2], 2) * H(2, 2) +
			2 * norm[0] * norm[1] * H(0, 1) +
			2 * norm[0] * norm[2] * H(0, 2) +
			2 * norm[1] * norm[2] * H(1, 2);

		return -(term1 - term2) / (2 * Len);
	}
}

void Spline_Gradient_MeanCurvature_Kernal_2p(const double* p1, const double* p2, const double* g1, double* re)
{
	double diff[3] = { *p1 - *p2, *(p1 + 1) - *(p2 + 1), *(p1 + 2) - *(p2 + 2) };
	double len_dist = length<double>(diff, 3);
	double fxxx, fyyy, fzzz, fxxy, fxxz, fyyx, fyyz, fzzx, fzzy, fxyz;
	if (len_dist < 1e-8)
	{
		fxxx = 0;
		fyyy = 0;
		fzzz = 0;

		fxxy = 0;
		fxxz = 0;
		fyyx = 0;
		fyyz = 0;
		fzzx = 0;
		fzzy = 0;

		fxyz = 0;
	}
	else
	{
		int temp_c = Multi_c;
		//temp_c = 3;
		fxxx = -(4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * pow(diff[0], 3) * pow(len_dist, 2 * temp_c - 5) -
			3 * (4 * pow(temp_c, 2) - 1) * diff[0] * pow(len_dist, 2 * temp_c - 3);
		fyyy = -(4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * pow(diff[1], 3) * pow(len_dist, 2 * temp_c - 5) -
			3 * (4 * pow(temp_c, 2) - 1) * diff[1] * pow(len_dist, 2 * temp_c - 3);
		fzzz = -(4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * pow(diff[2], 3) * pow(len_dist, 2 * temp_c - 5) -
			3 * (4 * pow(temp_c, 2) - 1) * diff[2] * pow(len_dist, 2 * temp_c - 3);

		fxxy = -(4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * pow(diff[0], 2) * diff[1] * pow(len_dist, 2 * temp_c - 5) -
			(4 * pow(temp_c, 2) - 1) * diff[1] * pow(len_dist, 2 * temp_c - 3);
		fxxz = -(4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * pow(diff[0], 2) * diff[2] * pow(len_dist, 2 * temp_c - 5) -
			(4 * pow(temp_c, 2) - 1) * diff[2] * pow(len_dist, 2 * temp_c - 3);
		fyyx = -(4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * pow(diff[1], 2) * diff[0] * pow(len_dist, 2 * temp_c - 5) -
			(4 * pow(temp_c, 2) - 1) * diff[0] * pow(len_dist, 2 * temp_c - 3);
		fyyz = -(4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * pow(diff[1], 2) * diff[2] * pow(len_dist, 2 * temp_c - 5) -
			(4 * pow(temp_c, 2) - 1) * diff[2] * pow(len_dist, 2 * temp_c - 3);
		fzzx = -(4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * pow(diff[2], 2) * diff[0] * pow(len_dist, 2 * temp_c - 5) -
			(4 * pow(temp_c, 2) - 1) * diff[0] * pow(len_dist, 2 * temp_c - 3);
		fzzy = -(4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * pow(diff[2], 2) * diff[1] * pow(len_dist, 2 * temp_c - 5) -
			(4 * pow(temp_c, 2) - 1) * diff[1] * pow(len_dist, 2 * temp_c - 3);

		fxyz = -(4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * diff[0] * diff[1] * diff[2] * pow(len_dist, 2 * temp_c - 5);

	}

	double term1[3];
	term1[0] = fxxx + fyyx + fzzx;
	term1[1] = fxxy + fyyy + fzzy;
	term1[2] = fxxz + fyyz + fzzz;

	double norm[3];
	norm[0] = *g1;
	norm[1] = *(g1 + 1);
	norm[2] = *(g1 + 2);
	double Len = length<double>(norm, 3);
	normalize<double>(norm, 3);

	double term2[3];
	term2[0] = pow(norm[0], 2) * fxxx + pow(norm[1], 2) * fyyx + pow(norm[2], 2) * fzzx + 2 * norm[0] * norm[1] * fxxy + 2 * norm[0] * norm[2] * fxxz + 2 * norm[1] * norm[2] * fxyz;
	term2[1] = pow(norm[0], 2) * fxxy + pow(norm[1], 2) * fyyy + pow(norm[2], 2) * fzzy + 2 * norm[0] * norm[1] * fyyx + 2 * norm[0] * norm[2] * fxyz + 2 * norm[1] * norm[2] * fyyz;
	term2[2] = pow(norm[0], 2) * fxxz + pow(norm[1], 2) * fyyz + pow(norm[2], 2) * fzzz + 2 * norm[0] * norm[1] * fxyz + 2 * norm[0] * norm[2] * fzzx + 2 * norm[1] * norm[2] * fzzy;
	for (int i = 0; i < 3; i++)
	{
		re[i] = -(term1[i] - term2[i]) / (2 * Len);
	}
}

double Spline_MeanCurvature_MeanCurvature_Kernal_2p(const double* p1, const double* p2, const double* g1, const double* g2)
{
	double diff[3] = { *p1 - *p2, *(p1 + 1) - *(p2 + 1), *(p1 + 2) - *(p2 + 2) };
	double len_dist = length<double>(diff, 3);

	double fxxxx, fyyyy, fzzzz,
		fxxyy, fxxzz, fyyzz,
		fxxyz, fyyxz, fzzxy,
		fxxxy, fxxxz, fyyyx, fyyyz, fzzzx, fzzzy;

	if (len_dist < 1e-8)
	{
		fxxxx = 0;
		fyyyy = 0;
		fzzzz = 0;

		fxxyy = 0;
		fxxzz = 0;
		fyyzz = 0;

		fxxyz = 0;
		fyyxz = 0;
		fzzxy = 0;

		fxxxy = 0;
		fxxxz = 0;
		fyyyx = 0;
		fyyyz = 0;
		fzzzx = 0;
		fzzzy = 0;
	}
	else
	{
		int temp_c = Multi_c;
		fxxxx = ((4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * (2 * temp_c - 5) * pow(diff[0], 4) * pow(len_dist, 2 * temp_c - 7) +
			6 * (4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * pow(diff[0], 2) * pow(len_dist, 2 * temp_c - 5) +
			3 * (4 * pow(temp_c, 2) - 1) * pow(len_dist, 2 * temp_c - 3));
		fyyyy = ((4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * (2 * temp_c - 5) * pow(diff[1], 4) * pow(len_dist, 2 * temp_c - 7) +
			6 * (4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * pow(diff[1], 2) * pow(len_dist, 2 * temp_c - 5) +
			3 * (4 * pow(temp_c, 2) - 1) * pow(len_dist, 2 * temp_c - 3));
		fzzzz = ((4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * (2 * temp_c - 5) * pow(diff[2], 4) * pow(len_dist, 2 * temp_c - 7) +
			6 * (4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * pow(diff[2], 2) * pow(len_dist, 2 * temp_c - 5) +
			3 * (4 * pow(temp_c, 2) - 1) * pow(len_dist, 2 * temp_c - 3));

		fxxyy = ((4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * (2 * temp_c - 5) * pow(diff[0], 2) * pow(diff[1], 2) * pow(len_dist, 2 * temp_c - 7) +
			(4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * (pow(diff[0], 2) + pow(diff[1], 2)) * pow(len_dist, 2 * temp_c - 5) +
			(4 * pow(temp_c, 2) - 1) * pow(len_dist, 2 * temp_c - 3));
		fxxzz = ((4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * (2 * temp_c - 5) * pow(diff[0], 2) * pow(diff[2], 2) * pow(len_dist, 2 * temp_c - 7) +
			(4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * (pow(diff[0], 2) + pow(diff[2], 2)) * pow(len_dist, 2 * temp_c - 5) +
			(4 * pow(temp_c, 2) - 1) * pow(len_dist, 2 * temp_c - 3));
		fyyzz = ((4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * (2 * temp_c - 5) * pow(diff[1], 2) * pow(diff[2], 2) * pow(len_dist, 2 * temp_c - 7) +
			(4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * (pow(diff[1], 2) + pow(diff[2], 2)) * pow(len_dist, 2 * temp_c - 5) +
			(4 * pow(temp_c, 2) - 1) * pow(len_dist, 2 * temp_c - 3));

		fxxyz = ((4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * (2 * temp_c - 5) * pow(diff[0], 2) * diff[1] * diff[2] * pow(len_dist, 2 * temp_c - 7) +
			(4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * diff[1] * diff[2] * pow(len_dist, 2 * temp_c - 5));
		fyyxz = ((4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * (2 * temp_c - 5) * pow(diff[1], 2) * diff[0] * diff[2] * pow(len_dist, 2 * temp_c - 7) +
			(4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * diff[0] * diff[2] * pow(len_dist, 2 * temp_c - 5));
		fzzxy = ((4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * (2 * temp_c - 5) * pow(diff[2], 2) * diff[0] * diff[1] * pow(len_dist, 2 * temp_c - 7) +
			(4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * diff[0] * diff[1] * pow(len_dist, 2 * temp_c - 5));

		fxxxy = ((4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * (2 * temp_c - 5) * pow(diff[0], 3) * diff[1] * pow(len_dist, 2 * temp_c - 7) +
			3 * (4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * diff[0] * diff[1] * pow(len_dist, 2 * temp_c - 5));
		fxxxz = ((4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * (2 * temp_c - 5) * pow(diff[0], 3) * diff[2] * pow(len_dist, 2 * temp_c - 7) +
			3 * (4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * diff[0] * diff[2] * pow(len_dist, 2 * temp_c - 5));
		fyyyx = ((4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * (2 * temp_c - 5) * pow(diff[1], 3) * diff[0] * pow(len_dist, 2 * temp_c - 7) +
			3 * (4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * diff[1] * diff[0] * pow(len_dist, 2 * temp_c - 5));
		fyyyz = ((4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * (2 * temp_c - 5) * pow(diff[1], 3) * diff[2] * pow(len_dist, 2 * temp_c - 7) +
			3 * (4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * diff[1] * diff[2] * pow(len_dist, 2 * temp_c - 5));
		fzzzx = ((4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * (2 * temp_c - 5) * pow(diff[2], 3) * diff[0] * pow(len_dist, 2 * temp_c - 7) +
			3 * (4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * diff[2] * diff[0] * pow(len_dist, 2 * temp_c - 5));
		fzzzy = ((4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * (2 * temp_c - 5) * pow(diff[2], 3) * diff[1] * pow(len_dist, 2 * temp_c - 7) +
			3 * (4 * pow(temp_c, 2) - 1) * (2 * temp_c - 3) * diff[2] * diff[1] * pow(len_dist, 2 * temp_c - 5));
	}


	arma::vec norm(3);
	arma::vec norm_x(3);

	norm[0] = g2[0];
	norm[1] = g2[1];
	norm[2] = g2[2];

	norm_x[0] = g1[0];
	norm_x[1] = g1[1];
	norm_x[2] = g1[2];

	double Len = length<double>(g2, 3);
	double Len_x = length<double>(g1, 3);
	norm = norm / Len;
	norm_x = norm_x / Len_x;
	//cout << "normal length: " << Len << ' ' << Len_x << endl;

	double term1 = fxxxx + fyyyy + fzzzz + 2 * fxxyy + 2 * fxxzz + 2 * fyyzz;

	arma::mat H(3, 3);
	H(0, 0) = fxxxx + fxxyy + fxxzz;
	H(1, 1) = fxxyy + fyyyy + fyyzz;
	H(2, 2) = fxxzz + fyyzz + fzzzz;
	H(0, 1) = H(1, 0) = fxxxy + fyyyx + fzzxy;
	H(0, 2) = H(2, 0) = fxxxz + fyyxz + fzzzx;
	H(1, 2) = H(2, 1) = fxxyz + fyyyz + fzzzy;
	arma::mat temp1, temp2;
	temp1 = norm.t() * H * norm;
	temp2 = norm_x.t() * H * norm_x;
	double term2 = temp1(0, 0) + temp2(0, 0);

	//arma::mat H(3, 3);

	H(0, 0) = pow(norm[0], 2) * fxxxx + pow(norm[1], 2) * fxxyy + pow(norm[2], 2) * fxxzz +
		2 * norm[0] * norm[1] * fxxxy + 2 * norm[0] * norm[2] * fxxxz + 2 * norm[1] * norm[2] * fxxyz;
	H(1, 1) = pow(norm[0], 2) * fxxyy + pow(norm[1], 2) * fyyyy + pow(norm[2], 2) * fyyzz +
		2 * norm[0] * norm[1] * fyyyx + 2 * norm[0] * norm[2] * fyyxz + 2 * norm[1] * norm[2] * fyyyz;
	H(2, 2) = pow(norm[0], 2) * fxxzz + pow(norm[1], 2) * fyyzz + pow(norm[2], 2) * fzzzz +
		2 * norm[0] * norm[1] * fzzxy + 2 * norm[0] * norm[2] * fzzzx + 2 * norm[1] * norm[2] * fzzzy;
	H(0, 1) = H(1, 0) = pow(norm[0], 2) * fxxxy + pow(norm[1], 2) * fyyyx + pow(norm[2], 2) * fzzxy +
		2 * norm[0] * norm[1] * fxxyy + 2 * norm[0] * norm[2] * fxxyz + 2 * norm[1] * norm[2] * fyyxz;
	H(0, 2) = H(2, 0) = pow(norm[0], 2) * fxxxz + pow(norm[1], 2) * fyyxz + pow(norm[2], 2) * fzzzx +
		2 * norm[0] * norm[1] * fxxyz + 2 * norm[0] * norm[2] * fxxzz + 2 * norm[1] * norm[2] * fzzxy;
	H(1, 2) = H(2, 1) = pow(norm[0], 2) * fxxyz + pow(norm[1], 2) * fyyyz + pow(norm[2], 2) * fzzzy +
		2 * norm[0] * norm[1] * fyyxz + 2 * norm[0] * norm[2] * fzzxy + 2 * norm[1] * norm[2] * fyyzz;

	arma::mat term3 = norm_x.t() * H * norm_x;
	return (term1 - term2 + term3(0, 0)) / (4 * Len * Len_x);
}