#include"reader.h"
#include <filesystem>

bool readXYZ(string filename, vector<double>& v) 
{
    ifstream file(filename, ios::in);
    char buf[1024];
    if (!file.is_open())
    {
        cout << "fail to read file!" << endl;
        exit(0);
    }
    while (file.getline(buf, 1024))
    {
        char delims[] = " ";
        double result;
        char* re = NULL;
        char* next_token = NULL;

        re = strtok_s(buf, delims, &next_token);
        sscanf_s(re, "%lf", &result);
        v.push_back(result);
        for (int i = 0; i < 2; ++i)
        {
            re = strtok_s(NULL, delims, &next_token);
            sscanf_s(re, "%lf", &result);
            v.push_back(result);
        }
    }
    return true;
}

bool readXYZnormal(string filename, vector<double>& v, vector<double>& vn, int dim) 
{
    ifstream file(filename, ios::in);
    char buf[1024];
    if (!file.is_open()) 
    { 
        return false;
    }
    v.clear();
    vn.clear();
    while (file.getline(buf, 1024))
    {
        char delims[] = " ";
        double result;
        char* re = NULL;
        char* next_token = NULL;

        re = strtok_s(buf, delims, &next_token);
        sscanf_s(re, "%lf", &result);
        v.push_back(result);
        for (int i = 0; i < 2; ++i)
        {
            re = strtok_s(NULL, delims, &next_token);
            sscanf_s(re, "%lf", &result);
            v.push_back(result);
        }

        for (int i = 0; i < 3; i++)
        {
            re = strtok_s(NULL, delims, &next_token);
            if (re != NULL)
            {
                sscanf_s(re, "%lf", &result);
                vn.push_back(result);
            }
            else
            {
                vn.push_back(0);
            }

        }
    }
    return true;
}


bool readMultiXYZnormal(string filename, vector<double>& v, vector<double>& vn,
    vector<int>& fileStartIndices, vector<double>& IoValues)
{
    vector<filesystem::path> xyzFiles;
    for (const auto& entry : filesystem::directory_iterator(filename)) 
    {
        if (entry.is_regular_file() && entry.path().extension() == ".xyz") 
        {
            xyzFiles.push_back(entry.path());
        }
    }
    sort(xyzFiles.begin(), xyzFiles.end());

    for (const auto& filePath : xyzFiles)
    {
        ifstream inFile(filePath);
        if (!inFile)
        {
            cerr << "Can not open the file: " << filePath << std::endl;
            continue;
        }

        fileStartIndices.push_back(v.size() / 3);
        string line;
        while (getline(inFile, line))
        {
            if (line.empty()) continue;

            std::istringstream iss(line);
            double x, y, z, nx, ny, nz;
            if (!(iss >> x >> y >> z >> nx >> ny >> nz))
            {
                std::cerr << "The format is error: " << filePath << std::endl;
                break;
            }
            v.push_back(x);
            v.push_back(y);
            v.push_back(z);

            vn.push_back(nx);
            vn.push_back(ny);
            vn.push_back(nz);
        }
        inFile.close();
    }
    vector<double> centralXYZ(fileStartIndices.size() * 3, 0);
    int n = v.size() / 3;
    for (int i = 0; i < fileStartIndices.size(); ++i)
    {
        int begin = fileStartIndices[i];
        int end;
        if (i != fileStartIndices.size() - 1)	end = fileStartIndices[i + 1];
        else  end = n;
        for (int j = begin; j < end; ++j)
        {
            centralXYZ[i * 3] += v[j * 3];
            centralXYZ[i * 3 + 1] += v[j * 3 + 1];
            centralXYZ[i * 3 + 2] += v[j * 3 + 2];
        }
        centralXYZ[i * 3] /= (end - begin);
        centralXYZ[i * 3 + 1] /= (end - begin);
        centralXYZ[i * 3 + 2] /= (end - begin);
    }
    //for (double x : centralXYZ)
    //{
    //    cout << x << ' ';
    //}
    IoValues.push_back(0);
    for (int i = 1; i < fileStartIndices.size(); ++i)
    {
        double dx = centralXYZ[i * 3] - centralXYZ[0];
        double dy = centralXYZ[i * 3 + 1] - centralXYZ[1];
        double dz = centralXYZ[i * 3 + 2] - centralXYZ[2];
        double d = sqrt(dx * dx + dy * dy + dz * dz);
        IoValues.push_back(-d);
    }

    

    return true;
}

bool readXYZnormalCurvature(string filename, vector<double>& v, vector<double>& vn, vector<double>& h, int dim)
{
    ifstream file(filename, ios::in);
    char buf[1024];
    if (!file.is_open())
    {
        cout << "fail to read file!" << endl;
        exit(0);
    }
    while (file.getline(buf, 1024))
    {
        char delims[] = " ";
        double result;
        char* re = NULL;
        char* next_token = NULL;

        re = strtok_s(buf, delims, &next_token);
        sscanf_s(re, "%lf", &result);
        v.push_back(result);
        for (int i = 0; i < 2; ++i)
        {
            re = strtok_s(NULL, delims, &next_token);
            sscanf_s(re, "%lf", &result);
            v.push_back(result);
        }

        for (int i = 0; i < 3; i++)
        {
            re = strtok_s(NULL, delims, &next_token);
            if (re != NULL)
            {
                sscanf_s(re, "%lf", &result);
                vn.push_back(result);
            }
            else
            {
                vn.push_back(0);
            }
        }

        re = strtok_s(NULL, delims, &next_token);
        if (re != NULL)
        {
            sscanf_s(re, "%lf", &result);
            h.push_back(result);
        }
        else
        {
            h.push_back(-9999);
        }
    }
    return true;
}

bool readMultiXYZnormalCurvature(string filename, vector<double>& v, vector<double>& vn, vector<double>& h,
    vector<int>& fileStartIndices, vector<double>& IoValues)
{
    vector<filesystem::path> xyzFiles;
    for (const auto& entry : filesystem::directory_iterator(filename))
    {
        if (entry.is_regular_file() && entry.path().extension() == ".txt")
        {
            xyzFiles.push_back(entry.path());
        }
    }
    sort(xyzFiles.begin(), xyzFiles.end());

    for (const auto& filePath : xyzFiles)
    {
        ifstream inFile(filePath);
        if (!inFile)
        {
            cerr << "Can not open the file: " << filePath << std::endl;
            continue;
        }

        fileStartIndices.push_back(v.size() / 3);
        string line;
        while (getline(inFile, line))
        {
            if (line.empty()) continue;

            std::istringstream iss(line);
            double x, y, z, nx, ny, nz, cur;
            if (!(iss >> x >> y >> z >> nx >> ny >> nz >> cur))
            {
                std::cerr << "The format is error: " << filePath << std::endl;
                break;
            }
            v.push_back(x);
            v.push_back(y);
            v.push_back(z);

            vn.push_back(nx);
            vn.push_back(ny);
            vn.push_back(nz);

            h.push_back(cur);
        }
        inFile.close();
    }
    vector<double> centralXYZ(fileStartIndices.size() * 3, 0);
    int n = v.size() / 3;
    for (int i = 0; i < fileStartIndices.size(); ++i)
    {
        int begin = fileStartIndices[i];
        int end;
        if (i != fileStartIndices.size() - 1)	end = fileStartIndices[i + 1];
        else  end = n;
        for (int j = begin; j < end; ++j)
        {
            centralXYZ[i * 3] += v[j * 3];
            centralXYZ[i * 3 + 1] += v[j * 3 + 1];
            centralXYZ[i * 3 + 2] += v[j * 3 + 2];
        }
        centralXYZ[i * 3] /= (end - begin);
        centralXYZ[i * 3 + 1] /= (end - begin);
        centralXYZ[i * 3 + 2] /= (end - begin);
    }

    IoValues.push_back(0);
    for (int i = 1; i < fileStartIndices.size(); ++i)
    {
        double dx = centralXYZ[i * 3] - centralXYZ[0];
        double dy = centralXYZ[i * 3 + 1] - centralXYZ[1];
        double dz = centralXYZ[i * 3 + 2] - centralXYZ[2];
        double d = sqrt(dx * dx + dy * dy + dz * dz);
        IoValues.push_back(-d);
    }

    return true;
}

bool writeXYZ(string filename, vector<double>& v, int dim)
{

    int npt = v.size() / dim;
    ofstream outer(filename.data(), ofstream::out);
    if (!outer.good()) 
    {
        cout << "Can not create output file " << filename << endl;
        return false;
    }
    outer << setprecision(18);
    if (dim == 3)
    {
        for (int i = 0; i < npt; ++i)
        {
            double* p_v = v.data() + i * 3;
            outer << p_v[0] << ' ' << p_v[1] << ' ' << p_v[2] << endl;
        }
    }
    else if (dim == 2)
    {
        for (int i = 0; i < npt; ++i)
        {
            double* p_v = v.data() + i * 2;
            outer << p_v[0] << ' ' << p_v[1] << endl;
        }
    }
   
    outer.close();
    return true;

}

bool writeXYZnormal(string filename, vector<double>& v, vector<double>& vn, int dim) 
{

    int npt = v.size() / 3;
    ofstream outer(filename.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not create output file " << filename << endl;
        return false;
    }
    outer << setprecision(18);
    if (dim == 3)
    {
        for (int i = 0; i < npt; ++i) 
        {
            double* p_v = v.data() + i * 3;
            double* p_vn = vn.data() + i * 3;
            outer  << p_v[0] << ' ' << p_v[1] << ' ' << p_v[2] << ' ';
            outer  << p_vn[0] << ' ' << p_vn[1] << ' ' << p_vn[2] << endl;
        }
    }
    else if (dim == 2)
    {
        for (int i = 0; i < npt; ++i)
        {
            double* p_v = v.data() + i * 2;
            double* p_vn = vn.data() + i * 2;
            outer  << p_v[0] << ' ' << p_v[1] << ' ';
            outer << p_vn[0] << ' ' << p_vn[1] << endl;
        }
    }
    
    outer.close();
    return true;

}

bool writeXYZnormalCurvature(string filename, vector<double>& v, vector<double>& vn, vector<double>& h, int dim)
{
    int npt = v.size() / 3;
    ofstream outer(filename.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not create output file " << filename << endl;
        return false;
    }
    outer << setprecision(18);
    if (dim == 3)
    {
        for (int i = 0; i < npt; ++i)
        {
            double* p_v = v.data() + i * 3;
            double* p_vn = vn.data() + i * 3;
            outer << p_v[0] << ' ' << p_v[1] << ' ' << p_v[2] << ' ';
            outer << p_vn[0] << ' ' << p_vn[1] << ' ' << p_vn[2] << ' ' << h[i] << endl;
        }
    }
    else if (dim == 2)
    {
        for (int i = 0; i < npt; ++i)
        {
            double* p_v = v.data() + i * 2;
            double* p_vn = vn.data() + i * 2;
            outer << p_v[0] << ' ' << p_v[1] << ' ';
            outer << p_vn[0] << ' ' << p_vn[1] << endl;
        }
    }

    outer.close();
    return true;
}

bool writeXYZnormalCurvature_Multi(string filename, vector<double>& v, vector<double>& vn, vector<double>& h,
    vector<int>& StartIds, vector<double>& IoValues)
{
    int npt = v.size() / 3;
    for (int i = 0; i < StartIds.size(); ++i)
    {
        int begin = StartIds[i];
        int end;
        if (i != StartIds.size() - 1)	end = StartIds[i + 1];
        else  end = npt;
        string outpath = filename + to_string(i) + "CurvRBF_" + to_string(IoValues[i]) + "_xyznc.txt";

        ofstream outer(outpath.data(), ofstream::out);
        if (!outer.good()) {
            cout << "Can not create output file " << outpath << endl;
            return false;
        }
        outer << setprecision(18);
        for (int j = begin; j < end; ++j)
        {  
            outer << v[j * 3] << ' ' << v[j * 3 + 1] << ' ' << v[j * 3 + 2] << ' ';
            outer << vn[j * 3] << ' ' << vn[j * 3 + 1] << ' ' << vn[j * 3 + 2] << ' ' << h[j] << endl; 
        }
        outer.close();
    }
    return true;
}

bool writeObjFile(string filename, const vector<double>& vertices, const vector<unsigned int>& faces2vertices) 
{
    filename = filename + ".obj";
    ofstream outer(filename.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not create output Obj file " << filename << endl;
        return false;
    }

    outer << setprecision(18);
    int n_vertices = vertices.size() / 3;
    int n_faces = faces2vertices.size() / 3;
    for (int i = 0; i < n_vertices; ++i) {
        auto p_v = vertices.data() + i * 3;
        outer << "v " << p_v[0] << " " << p_v[1] << " " << p_v[2] << endl;
    }

    for (int i = 0; i < n_faces; ++i) {
        auto p_fv = faces2vertices.data() + i * 3;
        outer << "f " << p_fv[0] + 1 << " " << p_fv[1] + 1 << " " << p_fv[2] + 1 << endl;
    }


    outer.close();
    cout << "saving finish: " << filename << endl;
    return true;
}

bool writePLYFile_VF(string filename, const vector<double>& vertices, const vector<unsigned int>& faces2vertices)
{
    filename = filename + ".ply";
    ofstream outer(filename.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not create output PLY file " << filename << endl;
        return false;
    }

    outer << setprecision(18);
    int n_vertices = vertices.size() / 3;
    int n_faces = faces2vertices.size() / 3;
    outer << "ply" << endl;
    outer << "format ascii 1.0" << endl;
    outer << "element vertex " << n_vertices << endl;
    outer << "property float x" << endl;
    outer << "property float y" << endl;
    outer << "property float z" << endl;
    outer << "element face " << n_faces << endl;
    outer << "property list uchar int vertex_indices" << endl;
    outer << "end_header" << endl;

    for (int i = 0; i < n_vertices; ++i) {
        auto p_v = vertices.data() + i * 3;
        for (int j = 0; j < 3; ++j)outer << setprecision(16) << p_v[j] << " ";
        outer << endl;
    }

    for (int i = 0; i < n_faces; ++i) {
        auto p_fv = faces2vertices.data() + i * 3;
        outer << "3 ";
        for (int j = 0; j < 3; ++j)outer << setprecision(16) << p_fv[j] << " ";
        outer << endl;
    }
    outer.close();
    cout << "saving finish: " << filename << endl;
    return true;
}

void SplitFileName(const string& fullfilename, string& filepath, string& filename, string& extname) 
{
    int pos;

    // Find the last dot for extension
    pos = fullfilename.find_last_of('.');
    if (pos != string::npos) 
    {
        extname = fullfilename.substr(pos);  // Extract the extension including the dot
    }
    else 
    {
        extname.clear();  
    }

    pos = fullfilename.find_last_of("\\/");
    if (pos != string::npos) 
    {
        filepath = fullfilename.substr(0, pos + 1);
        filename = fullfilename.substr(pos + 1, (fullfilename.length() - pos - 1) - extname.length());
    }
    else 
    {
        filepath.clear();  // No path
        filename = fullfilename.substr(0, fullfilename.length() - extname.length());  
    }
}

void writeMat(string filename, vector<vector<double>>& m)
{
    ofstream writer(filename.data(), ofstream::out);
    int n_col = m.size();
    int n_row = m[0].size();
    for (int i = 0;i<n_row;++i)
    {
        for (int j = 0; j < n_col; ++j)
        {
            writer << m[j][i] << ",";
        }
        writer << endl;
    }
    writer.close(); 
}

void writeVec(string filename, vector<double>& v)
{
    ofstream writer(filename.data(), ofstream::out);
    int n = v.size();
    for (int i = 0; i < n; ++i)
    {
        writer << v[i] << endl;
    }
    writer.close();
}

bool readPLYPos(string path, vector<double>& pos)
{
    ifstream inputFile(path);
    if (!inputFile)
    {
        cerr << "can't open the file." << endl;
        return 1;
    }

    string line;

    bool isDataSection = false;

    while (getline(inputFile, line))
    {
        if (isDataSection)
        {
            float x, y, z;
            if (sscanf_s(line.c_str(), "%f %f %f", &x, &y, &z) == 3)
            {
                if (x != 3)
                {
                    pos.push_back(x);
                    pos.push_back(y);
                    pos.push_back(z);
                    //cout << x << " " << y << ' ' << z << endl;
                }
                else break;
            }
        }
        else if (line == "end_header")
        {
            isDataSection = true;
        }
    }

    inputFile.close();
    return 0;
}

bool WriteXYZQ(string ReadPath, string OutPath, vector<double>& diff)
{
    ifstream inputFile(ReadPath);
    ofstream outputFile(OutPath);

    string line;
    int i = 0;

    bool isDataSection = false;
    while (getline(inputFile, line))
    {
        if (isDataSection)
        {
            float x, y, z;
            if (sscanf_s(line.c_str(), "%f %f %f", &x, &y, &z) == 3)
            {
                if (x != 3)
                {
                    line += to_string(diff[i++]);
                    outputFile << line << endl;
                    //cout <<line<< endl;
                }
                else
                {
                    outputFile << line << endl;
                }
            }
        }
        else if (line == "end_header")
        {
            isDataSection = true;
            outputFile << line << endl;
        }
        else if (line == "property float z")
        {
            outputFile << line << endl;
            outputFile << "property float quality" << endl;
        }
        else
        {
            outputFile << line << endl;
        }
    }

    inputFile.close();
    outputFile.close();

    return 0;
}