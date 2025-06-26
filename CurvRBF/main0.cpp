#include"RBF.h"
#include"reader.h"
#include"utility.h"
#include<iostream>
#include<format>

#include <random>
#include <limits>

using namespace std;

double torusMeanCurvature(double* p);
void PerturbVector(
    std::vector<double>& H,
    double strength,
    unsigned int seed,
    bool use_gaussian,
    double clamp_min,
    double clamp_max);

int main()
{
    int turn = 3;

    //basic HRBF
    if (turn == 0) 
    {
        string outpath, inpath, pcname, ext;
        string inpathName = R"(C:\Users\71708\Desktop\CurvRBF\PaperExperiment_r2\ErrorReconstruction\4kq_NEW\sam\p400.xyz)";
        outpath = R"(C:\Users\71708\Desktop\CurvRBF\PaperExperiment_r2\ErrorReconstruction\4kq_NEW\sam\)";
      
        if (outpath.empty())
            SplitFileName(inpathName, outpath, pcname, ext);
        else
            SplitFileName(inpathName, inpath, pcname, ext);
        cout << "input file: " << inpathName << endl;
        cout << "output path: " << outpath << endl;

        vector<double> Vs;
        vector<double> Gs;
        readXYZnormal(inpathName, Vs, Gs);
        cout << "number of points: " << Vs.size() / 3 << endl;

        RBF HRBF_Core(Vs, Gs);

        //HRBF_Core.Set_RBF_Method(HRBF_LowOrderTerm);
        HRBF_Core.Set_RBF_Method(HRBF);
        HRBF_Core.Set_Kernal_Func(Polyharmonic_Spline);
        HRBF_Core.Basic_Solution();

        //DMC GPU
        int n_voxels = 200;
        string objPath = outpath + pcname + "_HRBF";
        HRBF_Core.Surfacing_DMC_GPU_DataDirect(objPath, n_voxels, n_voxels, n_voxels);

        //Tracking-based MC
        //int n_voxel_line = 100;
        //HRBF_Core.Surfacing(n_voxel_line);
        //HRBF_Core.Write_Surface(outpath + pcname + "_surface");
        return 1;
    }

    //mean curvature HRBF
    if (turn == 1)
    {
        string outpath, inpath, pcname, ext;
        string inpathName = R"(C:\Users\71708\Desktop\CurvRBF\PaperExperiment_r2\ErrorReconstruction\torus_test\torus_cur.xyz)";
        outpath = R"(C:\Users\71708\Desktop\CurvRBF\PaperExperiment_r2\ErrorReconstruction\torus_test\)";

        if (outpath.empty())
            SplitFileName(inpathName, outpath, pcname, ext);
        else
            SplitFileName(inpathName, inpath, pcname, ext);
        cout << "input file: " << inpathName << endl;
        cout << "output path: " << outpath << endl;

        vector<double> Vs;
        vector<double> Gs;
        vector<double> h;
        readXYZnormalCurvature(inpathName, Vs, Gs, h);
        cout << "number of points: " << Vs.size() / 3 << endl;

        //h[1851] *= 2;
        //h[764] *= 2;
        //for (int i = 192; i < 216; ++i)
        //{
        //    h[i] = -0.01;
        //}
        //h[777] = h[777] * -1;
        //h[778] = h[778] * -1;
        //h[98] = h[98] * -1;
        //h[2507] = h[2507] * 3;
        //PerturbVector(h, 0.005, 12345, true, -0.03, 0.03);

        RBF HRBF_Core(Vs, Gs, h);
;
        HRBF_Core.Set_RBF_Method(Curvature_HRBF);
        HRBF_Core.Set_Kernal_Func(Polyharmonic_Spline);
        HRBF_Core.Basic_Solution();

        int n_voxels = 200;
        string objPath = outpath + pcname + "_CurvRBF";
        HRBF_Core.Surfacing_DMC_GPU_DataDirect(objPath, n_voxels, n_voxels, n_voxels);
        
        //int n_voxel_line = 100;
        //HRBF_Core.Surfacing(n_voxel_line);
        //HRBF_Core.Write_Surface(outpath + pcname + "_surface");

        HRBF_Core.Write_info(outpath + pcname + "_info.txt");
        return 1;
    }

    //two pass
    if (turn == 2)
    {
        /*-----------------The first pass--------------*/
        string outpath, inpath, pcname, ext;
        string inpathName = R"(C:\Users\71708\Desktop\CurvRBF\PaperExperiment_r2\support\1kq\1kq1.xyz)";
        outpath = R"()";

        if (outpath.empty())
            SplitFileName(inpathName, outpath, pcname, ext);
        else
            SplitFileName(inpathName, inpath, pcname, ext);
        cout << "input file: " << inpathName << endl;
        cout << "output path: " << outpath << endl;

        vector<double> Vs;
        vector<double> Gs;
        readXYZnormal(inpathName, Vs, Gs);
        cout << "number of points: " << Vs.size() / 3 << endl;
        RBF HRBF_0(Vs, Gs);
        
        HRBF_0.Set_RBF_Method(HRBF);
        HRBF_0.Set_Kernal_Func(Polyharmonic_Spline);
        HRBF_0.Basic_Solution();

        HRBF_0.Get_HRBF_MeanCurvature_Allpts();
        string writePath = outpath + pcname + "_xyznc.txt";
        HRBF_0.Write_XYZNC(writePath);

        int n_voxels = 150;
        string objPath = outpath + pcname + "_HRBF";
        //HRBF_0.Surfacing_DMC_GPU_DataDirect(objPath, n_voxels, n_voxels, n_voxels);
        //HRBF_0.Write_info(outpath + pcname + "_info.txt");
        //Tracking-based MC
        //int n_voxel_line = 100;
        //HRBF_0.Surfacing(n_voxel_line);
        //HRBF_0.Write_Surface(outpath + pcname + "_surface");

        /*-------------The second pass--------------*/
        int contin = 0;
        cout << "Press 1 to continue if xyznc is ready.\n";
        cin >> contin;
        if (contin)
        {
            vector<double> Vs;
            vector<double> Gs;
            vector<double> h;
            readXYZnormalCurvature(writePath, Vs, Gs, h);
            cout << "number of points: " << Vs.size() / 3 << endl;

            RBF CurvRBF(Vs, Gs, h);
            
            CurvRBF.Set_RBF_Method(Curvature_HRBF);
            CurvRBF.Set_Kernal_Func(Polyharmonic_Spline);
            CurvRBF.Basic_Solution();

            /*CurvRBF.Get_CurvRBF_MeanCurvature_Allpts();
            string writePath = outpath + pcname + "_xyznc2.txt";
            CurvRBF.Write_XYZNC(writePath);*/

            int n_voxels = 150;
            string objPath = outpath + pcname + "_CurvRBF";
            CurvRBF.Surfacing_DMC_GPU_DataDirect(objPath, n_voxels, n_voxels, n_voxels);
            CurvRBF.Write_info(outpath + pcname + "_info.txt");
        }

        return 1;
    }

    //two pass with gradient iteration
    if (turn == 3)
    {
        /*-----------------The first pass--------------*/
        string outpath, inpath, pcname, ext;
        string inpathName = R"(C:\Users\71708\Desktop\CurvRBF\PaperExperiment_r2\support\1kq\1kq1_150deg.xyz)";
        outpath = R"()";

        if (outpath.empty())
            SplitFileName(inpathName, outpath, pcname, ext);
        else
            SplitFileName(inpathName, inpath, pcname, ext);
        cout << "input file: " << inpathName << endl;
        cout << "output path: " << outpath << endl;

        vector<double> Vs;
        vector<double> Gs;
        readXYZnormal(inpathName, Vs, Gs);
        cout << "number of points: " << Vs.size() / 3 << endl;
        RBF HRBF_0(Vs, Gs);

        HRBF_0.Set_RBF_Method(HRBF);
        HRBF_0.Set_Kernal_Func(Polyharmonic_Spline);
        HRBF_0.Basic_Solution();

        HRBF_0.Get_HRBF_MeanCurvature_Allpts();
        string writePath = outpath + pcname + "_xyznc.txt";
        HRBF_0.Write_XYZNC(writePath);
        int n_voxels = 200;
        string objPath = outpath + pcname + "_HRBF";
        //HRBF_0.Surfacing_DMC_GPU_DataDirect(objPath, n_voxels, n_voxels, n_voxels);

        /*-------------The second pass--------------*/
        int contin = 0;
        cout << "Press 1 to continue if xyznc is ready.\n";
        cin >> contin;
        if (contin)
        {
            vector<double> Vs;
            vector<double> Gs;
            vector<double> h;
            readXYZnormalCurvature(writePath, Vs, Gs, h);
            int n = Vs.size() / 3;
            cout << "number of points: " << n << endl;

            //Calculate information of new points
            vector<int> ids;
            for (int i = n - 1; i < n; ++i)
            {
                ids.push_back(i);
            }
            for (int i : ids)
            {
                vector<double> gradient = HRBF_0.Get_HRBF_Gradient(Vs.data() + i * 3);
                Gs[i * 3] = gradient[0]; Gs[i * 3 + 1] = gradient[1]; Gs[i * 3 + 2] = gradient[2];
                h[i] = HRBF_0.Get_HRBF_MeanCurvature(Vs.data() + i * 3, Gs.data() + i * 3);
                cout << "Target position:" << Vs[i * 3] << ' ' << Vs[i * 3 + 1] << ' ' << Vs[i * 3 + 2] << endl;
                cout << "Initial gradient:" << Gs[i * 3] << ' ' << Gs[i * 3 + 1] << ' ' << Gs[i * 3 + 2] << endl;
                cout << "Initial mean curvature:" << h[i] << endl;
                cout << "Enter new mean curvature:";
                double f;
                cin >> f;
                h[i] *= f;
            }
            //h[1037] *= -5000;

            //solve
            RBF CurvRBF(Vs, Gs, h);
            CurvRBF.Set_RBF_Method(Curvature_HRBF);
            CurvRBF.Set_Kernal_Func(Polyharmonic_Spline);
            //double factor1 = 1e7;
            double factor1 = 0;
            double factor2 = 5;

            CurvRBF.isoffContact = 0;
            CurvRBF.offContact_list = ids;
            //CurvRBF.Basic_Solution();
            CurvRBF.Gradient_Iteration_Solution(factor1, factor2, ids);


            CurvRBF.Write_XYZNC(writePath);

            cout << "iosurface: " << CurvRBF.Get_MeanCurvature_HRBF_FuncValue(Vs.data() + (n - 1) * 3) << endl;

            int n_voxels = 200;
            string objPath = outpath + pcname + "_CurvRBF";
            CurvRBF.Surfacing_DMC_GPU_DataDirect(objPath, n_voxels, n_voxels, n_voxels);
            CurvRBF.Write_info(outpath + pcname + "_info.txt");
        }
        return 1;
    }

    //compare difference 
    if (turn == 4)
    {
        string referenceModelPointsPath = R"(C:\Users\71708\Desktop\CurvRBF\PaperExperiment_r2\support\3kq\initial_10x.txt)"; //原始模型控制点
        //string testModelPointsPath = R"(D:\desktop\data\1\machineLearning\samplePoints.xyz)";  //待比较模型控制点
        string testModel = R"(C:\Users\71708\Desktop\CurvRBF\PaperExperiment_r2\support\3kq\opt_10x.ply)";	 //待比较模型
        string outPath = R"(C:\Users\71708\Desktop\CurvRBF\PaperExperiment_r2\support\3kq\opt_10x_differ2.ply)";
        double mean_diff = 0.0;

        //Get Reference model implicit func
        vector<double> Vs, Gs;
        readXYZnormal(referenceModelPointsPath, Vs, Gs);
        RBF Rerfer_HRBF(Vs, Gs);
        Rerfer_HRBF.Set_RBF_Method(HRBF);
        Rerfer_HRBF.Set_Kernal_Func(Polyharmonic_Spline);
        Rerfer_HRBF.Basic_Solution();

        //Rerfer_HRBF.Surfacing(100);
        //Rerfer_HRBF.Write_Surface(outpath + pcname + "_surface");

        //toy implicit function
        auto Torus = [](R3Pt& p)
        {
            double r = 1.0;
            double R = 4.0;
            double x = p[0];
            double y = p[1];
            double z = p[2];
            return (sqrt(x * x + y * y) - R) * (sqrt(x * x + y * y) - R) + z * z - r * r;
        };
        auto Grad_Torus = [](R3Pt& p)
        {
            double r = 1.0;
            double R = 4.0;
            double x = p[0];
            double y = p[1];
            double z = p[2];
            double dx = 2 * x - 2 * R * x / sqrt(x * x + y * y);
            double dy = 2 * y - 2 * R * y / sqrt(x * x + y * y);
            double dz = 2 * z;
            vector<double> g = { dx,dy,dz };
            return g;
        };

        //get points from test model surface
        vector<double> pos;
        readPLYPos(testModel, pos);

        //calculate difference
        vector<double> diff;
        diff.reserve(pos.size() / 3);
        for (int i = 0; i < pos.size() / 3; ++i)
        {
            R3Pt in_pt(pos[i * 3], pos[i * 3 + 1], pos[i * 3 + 2]);

            //double d = abs(Rerfer_HRBF.Dist_Function(in_pt));
            double d = Rerfer_HRBF.Dist_Function(in_pt);
            //double d = abs(Torus(in_pt));

            vector<double> Gradient = Rerfer_HRBF.Get_HRBF_Gradient(pos.data() + i * 3);
            //vector<double> Gradient = Grad_Torus(in_pt);

            double L = length<double>(Gradient, 3);
            
            diff.push_back(d / 1);
        }
        WriteXYZQ(testModel, outPath, diff);
        double sum = 0.0;
        for (double e : diff)
        {
            sum += e;
        }
        cout << "mean difference: " << sum / diff.size() << endl;
    }

    //计算初始法向量
    if (turn == 5)
    {
        string outpath, inpath, pcname, ext;
        string inpathName = R"(D:\desktop\毕业论文\实验\missSamplePoints.xyz)";
        outpath = R"(D:\desktop\毕业论文\实验\)";

        if (outpath.empty())
            SplitFileName(inpathName, outpath, pcname, ext);
        else
            SplitFileName(inpathName, inpath, pcname, ext);
        cout << "input file: " << inpathName << endl;
        cout << "output path: " << outpath << endl;

        vector<double> Vs;
        vector<double> Gs;
        readXYZnormal(inpathName, Vs, Gs);
        cout << "number of points: " << Vs.size() / 3 << endl;

        RBF HRBF_Core(Vs, Gs);
        HRBF_Core.Set_RBF_Method(HRBF);
        HRBF_Core.Set_Kernal_Func(Polyharmonic_Spline);
        HRBF_Core.Basic_Solution();

        string TargetInpathName = R"(D:\desktop\毕业论文\实验\targetPoints.xyz)";
        string TargetOutpathName = R"(D:\desktop\毕业论文\实验\targetPoints_initial.xyz)";
        readXYZnormal(TargetInpathName, Vs, Gs);
        int n = Vs.size() / 3;
        vector<double> re_Gs;
        for (int i = 0; i < n; ++i)
        {
            double* p = Vs.data() + i * 3;
            vector<double> re = HRBF_Core.Get_HRBF_Gradient(p);
            re_Gs.insert(re_Gs.end(), re.begin(), re.end());
        }
        writeXYZnormal(TargetOutpathName, Vs, re_Gs);


    }

    //计算圆环面上曲率
    if (turn == 6)
    {
        vector<double> Vs;
        vector<double> Gs;
        string TargetInpathName = R"(D:\desktop\毕业论文\实验\targetPoints_initial.xyz)";
        string TargetOutpathName = R"(D:\desktop\毕业论文\实验\targetPoints_initial_curvature.txt)";
        readXYZnormal(TargetInpathName, Vs, Gs);
        int n = Vs.size() / 3;
        vector<double> re_c;
        for (int i = 0; i < n; ++i)
        {
            double* p = Vs.data() + i * 3;
            double re = torusMeanCurvature(p);
            re_c.push_back(re);
        }
        writeXYZnormalCurvature(TargetOutpathName, Vs, Gs, re_c);
    }

    if (turn == 7)
    {
        vector<double> Vs;
        vector<double> Gs;
        vector<double> h;
        string inpathName = R"(D:\desktop\毕业论文\实验\miss_target_points.txt)";
        readXYZnormalCurvature(inpathName, Vs, Gs, h);
        cout << "number of points: " << Vs.size() / 3 << endl;

        //Calculate information of new points
        vector<int> ids = { 287, 288, 289, 290,291,292,293,294,295 };

        //solve
        RBF CurvRBF(Vs, Gs, h);
        CurvRBF.Set_RBF_Method(Curvature_HRBF);
        CurvRBF.Set_Kernal_Func(Polyharmonic_Spline);
        //double factor1 = 10;
        //double factor2 = 5;
        double factor1 = 0.01;
        double factor2 = 5;
        CurvRBF.isoffContact = 0;
        CurvRBF.offContact_list = ids;
        CurvRBF.Gradient_Iteration_Solution(factor1, factor2, ids);

    }

    //MutiXYZ_first_pass
    if (turn == 8)
    {
        string inpathName = R"(D:\desktop\PaperExperiment_r2\diceng\scalarField\)";
        string outpath = R"(D:\desktop\PaperExperiment_r2\diceng\scalarField\)";

        cout << "input file: " << inpathName << endl;
        cout << "output path: " << outpath << endl;

        vector<double> Vs;
        vector<double> Gs;
        vector<int> StartIds;
        vector<double> IoValues;
        readMultiXYZnormal(inpathName, Vs, Gs, StartIds, IoValues);
        cout << "number of points: " << Vs.size() / 3 << endl;

        //IoValues = { 0, -798.252, -1954.86, -2996.39, -4094.3, -5485.67, -5948.01, -6899.25 };
        RBF HRBF_Core(Vs, Gs);

        HRBF_Core.Set_RBF_Method(HRBF_LowOrderTerm);
        //HRBF_Core.Set_RBF_Method(HRBF);
        HRBF_Core.Set_Kernal_Func(Polyharmonic_Spline);
        HRBF_Core.SetMultiValuesSolution = 1;
        HRBF_Core.IoValues = IoValues;
        HRBF_Core.StartIds = StartIds;
        HRBF_Core.Basic_Solution();

        //DMC GPU
        int n_voxels = 200;
        string objPath = outpath + "CurvRBF";
        HRBF_Core.Surfacing_DMC_GPU_DataDirect_MultiValues(objPath, n_voxels, n_voxels, n_voxels, IoValues);

        HRBF_Core.Get_HRBF_MeanCurvature_Allpts();
        HRBF_Core.Write_Multi_XYZNC(outpath);

        return 1;
    }
    
    //MutiXYZ_second_pass
    if (turn == 9)
    {
        string inpathName = R"(D:\desktop\PaperExperiment_r2\diceng\off-contact\onePoint\FirstPass\)";
        string outpath = R"(D:\desktop\PaperExperiment_r2\diceng\off-contact\onePoint\SecondPass\0.01\)";

        cout << "input file: " << inpathName << endl;
        cout << "output path: " << outpath << endl;

        vector<double> Vs;
        vector<double> Gs;
        vector<double> H;
        vector<int> StartIds;
        vector<double> IoValues;
        readMultiXYZnormalCurvature(inpathName, Vs, Gs, H, StartIds, IoValues);
        int n = Vs.size() / 3;
        cout << "number of points: " << n << endl;

        //PerturbVector(H, 0.02, 12345, true, -0.1, 0.1);
        //vector<int> ids = { StartIds[5] + 74, StartIds[5] + 75, StartIds[5] + 76, StartIds[5] + 77, StartIds[5] + 78,
        //StartIds[5] + 79, StartIds[5] + 80, StartIds[5] + 81, StartIds[5] + 82, StartIds[5] + 83 };
        //for (int i : ids)
        //{
        //    cout << "Target Points: " << Vs[i * 3] << ' ' << Vs[i * 3 + 1] << ' ' << Vs[i * 3 + 2] << endl;
        //    cout << "Target Points Curvature: " << H[i] << endl;
        //    double s;
        //    //cin >> s;
        //    H[i] = -0.005;
        //}
        vector<int> ids;
        for (int i = 1; i < 3; ++i)
        {
            ids.push_back(n - i);
        }
        for (int i : ids)
        {
            H[i] = 0.005;
        }


        RBF HRBF_Core(Vs, Gs, H);
        
        //HRBF_Core.Set_RBF_Method(HRBF_LowOrderTerm);
        HRBF_Core.Set_RBF_Method(Curvature_HRBF);
        HRBF_Core.Set_Kernal_Func(Polyharmonic_Spline);
        HRBF_Core.SetMultiValuesSolution = 1;
        HRBF_Core.IoValues = IoValues;
        HRBF_Core.StartIds = StartIds;
        HRBF_Core.isoffContact = 1;
        HRBF_Core.offContact_list = ids;
        HRBF_Core.Basic_Solution();

        //DMC GPU
        int n_voxels = 200;
        string objPath = outpath + "CurvRBF";
        HRBF_Core.Surfacing_DMC_GPU_DataDirect_MultiValues(objPath, n_voxels, n_voxels, n_voxels, IoValues);

        return 1;
    }

    //Multi-stra gradient iteration
    if (turn == 10)
    {
        /*-----------------The first pass--------------*/
        string inpathName = R"(D:\desktop\PaperExperiment_r2\diceng\off-contact\onePoint\sour\)";
        string outpath = R"(D:\desktop\PaperExperiment_r2\diceng\off-contact\onePoint\FirstPass\)";
        string outpath_second = R"(D:\desktop\PaperExperiment_r2\diceng\off-contact\onePoint\SecondPass\sour\)";

        cout << "input file: " << inpathName << endl;
        cout << "output path: " << outpath << endl;

        vector<double> Vs;
        vector<double> Gs;
        vector<int> StartIds;
        vector<double> IoValues;
        readMultiXYZnormal(inpathName, Vs, Gs, StartIds, IoValues);
        cout << "number of points: " << Vs.size() / 3 << endl;

        RBF HRBF_0(Vs, Gs);
        HRBF_0.Set_RBF_Method(HRBF);
        HRBF_0.Set_Kernal_Func(Polyharmonic_Spline);
        HRBF_0.SetMultiValuesSolution = 1;
        HRBF_0.IoValues = IoValues;
        HRBF_0.StartIds = StartIds;
        HRBF_0.Basic_Solution();

        //DMC GPU
        //int n_voxels = 200;
        //string objPath = outpath + "CurvRBF";
        //HRBF_Core.Surfacing_DMC_GPU_DataDirect_MultiValues(objPath, n_voxels, n_voxels, n_voxels, IoValues);

        HRBF_0.Get_HRBF_MeanCurvature_Allpts();
        HRBF_0.Write_Multi_XYZNC(outpath);

        /*-------------The second pass--------------*/
        int contin = 0;
        cout << "Press 1 to continue if xyznc is ready.\n";
        cin >> contin;
        if (contin)
        {
            vector<double> Vs;
            vector<double> Gs;
            vector<double> H;
            vector<int> StartIds;
            vector<double> IoValues;
            readMultiXYZnormalCurvature(outpath, Vs, Gs, H, StartIds, IoValues);
            int n = Vs.size() / 3;
            cout << "number of points: " << n << endl;

            //Calculate information of new points
            vector<int> ids;
            for (int i = 1; i < 3; ++i)
            {
                ids.push_back(n - i);
            }
            for (int i : ids)
            {
                vector<double> gradient = HRBF_0.Get_HRBF_Gradient(Vs.data() + i * 3);
                Gs[i * 3] = gradient[0]; Gs[i * 3 + 1] = gradient[1]; Gs[i * 3 + 2] = gradient[2];
                H[i] = HRBF_0.Get_HRBF_MeanCurvature(Vs.data() + i * 3, Gs.data() + i * 3);
                cout << "Target position:" << Vs[i * 3] << ' ' << Vs[i * 3 + 1] << ' ' << Vs[i * 3 + 2] << endl;
                cout << "Initial gradient:" << Gs[i * 3] << ' ' << Gs[i * 3 + 1] << ' ' << Gs[i * 3 + 2] << endl;
                cout << "Initial mean curvature:" << H[i] << endl;
                cout << "Enter new mean curvature:";
                double f;
                //cin >> f;
                //H[i] = 0.005;
            }

            RBF HRBF_Core(Vs, Gs, H);
            HRBF_Core.Set_RBF_Method(Curvature_HRBF);
            HRBF_Core.Set_Kernal_Func(Polyharmonic_Spline);
            double factor1 = 1e4;
            double factor2 = 5;
            HRBF_Core.SetMultiValuesSolution = 1;
            HRBF_Core.IoValues = IoValues;
            HRBF_Core.StartIds = StartIds;
            HRBF_Core.isoffContact = 1;
            HRBF_Core.offContact_list = ids;
            HRBF_Core.Gradient_Iteration_Solution(factor1, factor2, ids);
            HRBF_Core.Write_Multi_XYZNC(outpath_second);

            //DMC GPU
            int n_voxels = 200;
            string objPath = outpath_second + "CurvRBF";
            HRBF_Core.Surfacing_DMC_GPU_DataDirect_MultiValues(objPath, n_voxels, n_voxels, n_voxels, IoValues);
            HRBF_Core.Write_info(outpath_second + "_info.txt");
        }

        return 1;
    }

    //计算点的平均曲率
    if (turn == 11)
    {
        string outpath, inpath, pcname, ext;
        string inpathName = R"(C:\Users\71708\Desktop\CurvRBF\PaperExperiment_r2\ErrorReconstruction\4kq_NEW\sam\4kq.xyz)";   //模型的控制点

        if (outpath.empty())
            SplitFileName(inpathName, outpath, pcname, ext);
        else
            SplitFileName(inpathName, inpath, pcname, ext);
        cout << "input file: " << inpathName << endl;
        cout << "output path: " << outpath << endl;

        vector<double> Vs;
        vector<double> Gs;
        readXYZnormal(inpathName, Vs, Gs);
        cout << "number of points: " << Vs.size() / 3 << endl;

        RBF HRBF_Core(Vs, Gs);
        HRBF_Core.Set_RBF_Method(HRBF);
        HRBF_Core.Set_Kernal_Func(Polyharmonic_Spline);
        HRBF_Core.Basic_Solution();

        string PointsInpathName = R"(C:\Users\71708\Desktop\CurvRBF\PaperExperiment_r2\ErrorReconstruction\4kq_NEW\sam\p400.xyz)";     //待计算曲率的控制点
        string PointsOutpathName = R"(C:\Users\71708\Desktop\CurvRBF\PaperExperiment_r2\ErrorReconstruction\4kq_NEW\sam\p400_cur.txt)";
        readXYZnormal(PointsInpathName, Vs, Gs);

        int n = Vs.size() / 3;
        vector<double> h_list;
        h_list.reserve(n);

        for (int i = 0; i < n; ++i)
        {
            double* p = Vs.data() + i * 3;
            double* g = Gs.data() + i * 3;
            double h = HRBF_Core.Get_HRBF_MeanCurvature(p, g);
            h_list.push_back(h);
        }
        writeXYZnormalCurvature(PointsOutpathName, Vs, Gs, h_list);

    }

}

double torusMeanCurvature(double* p) 
{
    double x = p[0];
    double y = p[1];
    double z = p[2];
    double R = 2;
    double r = 1;
    // 计算环面隐函数的平均曲率
    double s = std::sqrt(x * x + y * y);
    if (s == 0.0) return 0.0; // 避免除以零（环面上s >= R-r > 0）

    // 一阶偏导数
    double Fx = 2.0 * (s - R) * x / s;
    double Fy = 2.0 * (s - R) * y / s;
    double Fz = 2.0 * z;

    // 二阶偏导数
    double Fxx = 2.0 - (2 * R * y * y) / pow(s, 3);
    double Fyy = 2.0 - (2 * R * x * x) / pow(s, 3);
    double Fzz = 2.0;
    double Fxy = 2.0 * R* x * y / pow(s, 3);
    double Fxz = 0;
    double Fyz = 0;

    arma::vec G;
    G.set_size(3);
    G(0) = Fx; G(1) = Fy; G(2) = Fz;
    arma::mat H;
    H.set_size(3, 3);
    H(0, 0) = Fxx; H(1, 1) = Fyy; H(2, 2) = Fzz;
    H(0, 1) = H(1, 0) = Fxy;
    H(0, 2) = H(2, 0) = Fxz;
    H(1, 2) = H(2, 1) = Fyz;
    double Lap = Fxx + Fyy + Fzz;
    double Len_g = arma::norm(G, 2);
    // 计算分子和分母
    arma::mat re = G.t() * H * G;
    double numerator = re(0,0) - Len_g * Len_g * Lap;
    double denominator = 2.0 * Len_g * Len_g * Len_g;

    return numerator / denominator;
}


//Perturb the values in the array
void PerturbVector(
    std::vector<double>& H,
    double strength,
    unsigned int seed,
    bool use_gaussian,
    double clamp_min,
    double clamp_max)
{
    if (H.empty() || strength <= 0.0) return;

    // Initialize stochastic engine (using hardware entropy source as seed)
    std::random_device rd;
    std::mt19937_64 engine(seed != 0 ? seed : rd());

    // Select Distribution Type
    if (use_gaussian)
    {
        std::normal_distribution<double> dist(0.0, strength);
        for (double& h : H)
        {
            h += dist(engine);
            h = std::clamp(h, clamp_min, clamp_max);
        }
    }
    else
    {
        std::uniform_real_distribution<double> dist(-strength, strength);
        for (double& h : H)
        {
            h += dist(engine);
            h = std::clamp(h, clamp_min, clamp_max);
        }
    }
}




