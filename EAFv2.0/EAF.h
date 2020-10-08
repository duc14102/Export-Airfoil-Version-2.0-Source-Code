#pragma once
#include "Vector_matrix.h"
vector<string> print_database();
string input_airfoil(vector<string>data);
void export_airfoil_data(string NACA);
void export_4_digit(string NACA, int N, string namefile);
void export_5_digit(string NACA, int N, string namefile);
void convert_airfoil_data();
void convert_parsec(vector <double> x0, vector <double> y0, int N, string name, string filename_out);
vectorMT Parsec(vector <double>x0, vector <double> y0);
double nrsovler(vector <double> ca, vector <double> cx, double res, int maxloop);
double evalf(vector <double>ca, vector <double>cx, double x);