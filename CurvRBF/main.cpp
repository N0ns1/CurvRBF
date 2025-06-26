#include"RBF.h"
#include"reader.h"
#include"utility.h"
#include<iostream>
#include<format>

const double DEFAULT_FACTOR_A = 1.0;
const double DEFAULT_FACTOR_B = 5.0;
const int DEFAULT_N_VOXELS = 200;

void print_usage(const char* prog_name)
{
    std::cerr << "Usage: " << prog_name << " <input_xyz_file> [output_path]" << std::endl;
    std::cerr << "Example: " << prog_name << " C:\\data\\points.xyz C:\\data\\output\\" << std::endl;
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        print_usage(argv[0]);
        return 1; 
    }

    std::string inpathName = argv[1];
    std::string outpath, temp;
    std::string pcname, ext;

    if (argc > 2) 
    {
        SplitFileName(inpathName, temp, pcname, ext);
        outpath = argv[2];
        outpath = outpath + "\\";
    }
    else 
    {
        SplitFileName(inpathName, outpath, pcname, ext);
    }

    std::cout << "Input file: " << inpathName << std::endl;
    std::cout << "Output path: " << outpath << std::endl;

    /*----------------First pass--------------*/
    std::vector<double> Vs;
    std::vector<double> Gs;
    if (!readXYZnormal(inpathName, Vs, Gs)) 
    { 
        std::cerr << "Error: Failed to read input file: " << inpathName << std::endl;
        return 1;
    }
    std::cout << "Number of points: " << Vs.size() / 3 << std::endl;

    RBF HRBF_0(Vs, Gs);

    HRBF_0.Set_RBF_Method(HRBF);
    HRBF_0.Set_Kernal_Func(Polyharmonic_Spline);
    HRBF_0.Basic_Solution();

    HRBF_0.Get_HRBF_MeanCurvature_Allpts();
    std::string intermediatePath = outpath + pcname + "_xyznc.txt";
    HRBF_0.Write_XYZNC(intermediatePath);
    std::cout << "First pass complete. Intermediate file written to: " << intermediatePath << std::endl;

    /*-------------Second pass--------------*/
    std::vector<double> Vs_2;
    std::vector<double> Gs_2;
    std::vector<double> h;
    readXYZnormalCurvature(intermediatePath, Vs_2, Gs_2, h);

    vector<int> OffContact_ids;
    OffContact_ids.clear();
    vector<int> Target_ids;

    cout << "Start adding mean curvature control points...\n";
    int contin = 1;
    while (contin)
    {
        cout << "Add new points? Enter 1 or 0\n";
        int isNewPoint;
        cin >> isNewPoint;
        if (isNewPoint)
        {
            double x, y, z;
            cout << "Enter the XYZ coordinates of the points, separated by spaces: \n";
            cin >> x >> y >> z;
            Vs_2.insert(Vs_2.end(), { x,y,z });
            int np = Vs_2.size() / 3;
            vector<double> gradient = HRBF_0.Get_HRBF_Gradient(Vs_2.data() + (np - 1) * 3);
            Gs_2.insert(Gs_2.end(), gradient.begin(), gradient.end());
            h.push_back(HRBF_0.Get_HRBF_MeanCurvature(Vs_2.data() + (np - 1) * 3, Gs_2.data() + (np - 1) * 3));
            cout << "The initial mean curvature at this point is " << -h[h.size() - 1] << endl;
            cout << "Please enter the modification multiplier: ";
            double f;
            cin >> f;
            h[h.size() - 1] *= f;
            cout << "Enter the constraint type (0 for on-contact and 1 for off-contact):";
            int isOff;
            cin >> isOff;
            if (isOff)   OffContact_ids.push_back(h.size() - 1);
            Target_ids.push_back(h.size() - 1);
        }
        else
        {
            cout << "Modify existing point mean curvature? Enter 1 or 0\n";
            int isModify;
            cin >> isModify;
            if (isModify)
            {
                cout << "Enter the id of an existing point (count from 0): \n";
                int id;
                cin >> id;
                cout << "The initial mean curvature at this point is " << -h[id] << endl;
                cout << "Please enter the modification multiplier: ";
                double f;
                cin >> f;
                h[id] *= f;
                cout << "Enter the constraint type (0 for on-contact and 1 for off-contact):";
                int isOff;
                cin >> isOff;
                if (isOff)   OffContact_ids.push_back(id);
                Target_ids.push_back(id);
            }
            else
            {
                break;
            }
        }
    }


    //solve
    RBF CurvRBF(Vs_2, Gs_2, h);
    CurvRBF.Set_RBF_Method(Curvature_HRBF);
    CurvRBF.Set_Kernal_Func(Polyharmonic_Spline);
    CurvRBF.isoffContact = 1;
    CurvRBF.offContact_list = OffContact_ids;
    //CurvRBF.Basic_Solution();
    CurvRBF.Gradient_Iteration_Solution(DEFAULT_FACTOR_A, DEFAULT_FACTOR_B, OffContact_ids);
    CurvRBF.Write_XYZNC(intermediatePath);

    ////DMC GPU
    string objPath = outpath + pcname + "_CurvRBF";
    CurvRBF.Surfacing_DMC_GPU_DataDirect(objPath, DEFAULT_N_VOXELS, DEFAULT_N_VOXELS, DEFAULT_N_VOXELS);
    CurvRBF.Write_info(outpath + pcname + "_info.txt");

    //Tracking-based MC
    //CurvRBF.Surfacing(DEFAULT_N_VOXELS);
    //CurvRBF.Write_Surface(outpath + pcname + "_surface");
    
    return 1;

}
