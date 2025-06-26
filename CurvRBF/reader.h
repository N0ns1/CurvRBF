#pragma once
#include<vector>
#include<string>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
using namespace std;

bool readXYZ(string filename, vector<double>& v);
bool readXYZnormal(string filename, vector<double>& v, vector<double>& vn, int dim = 3);
bool readXYZnormalCurvature(string filename, vector<double>& v, vector<double>& vn, vector<double>& h, int dim = 3);
bool readMultiXYZnormal(string filename, vector<double>& v, vector<double>& vn, vector<int>& fileStartIndices, vector<double>& IoValues);
bool readMultiXYZnormalCurvature(string filename, vector<double>& v, vector<double>& vn, vector<double>& h, 
	vector<int>& fileStartIndices, vector<double>& IoValues);
bool readPLYPos(string path, vector<double>& pos);

bool writeXYZ(string filename, vector<double>& v, int dim = 3);
bool writeXYZnormal(string filename, vector<double>& v, vector<double>& vn, int dim = 3);
bool writeXYZnormalCurvature(string filename, vector<double>& v, vector<double>& vn, vector<double>& h, int dim = 3);
bool writeXYZnormalCurvature_Multi(string filename, vector<double>& v, vector<double>& vn, vector<double>& h, 
	vector<int>& StartIds, vector<double>& IoValues);

bool WriteXYZQ(string ReadPath, string OutPath, vector<double>& diff);

bool writeObjFile(string filename, const vector<double>& vertices, const vector<unsigned int>& faces2vertices);
bool writePLYFile_VF(string filename, const vector<double>& vertices, const vector<unsigned int>& faces2vertices);

void SplitFileName(const string& fullfilename, string& filepath, string& filename, string& extname);

void writeMat(string filename, vector<vector<double>>& m);

void writeVec(string filename, vector<double>& v);