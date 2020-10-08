#include "Header.h"
#include "EAF.h"
void export_airfoil_data(string NACA)
{
	int N = 200;
	do
	{
		if (N < 10)
			cout << "\nMinimum number of points for running this function is 10\n";
		cout << "Enter the number of points: ";
		cin >> N;
	} while (N<10);
	string namefile;
	cout << "Enter the name of file: ";
	cin >> namefile;
	if (NACA.size() == 4)
	{
		export_4_digit(NACA, N, namefile);
	}
	else if(NACA.size() == 5)
	{
		export_5_digit(NACA, N, namefile);
	}
	cout << "\nExport to file \""<<namefile<<"\"...Done!\n" << endl;
}
void export_4_digit(string NACA, int N, string namefile)
{
	// http://www.airfoiltools.com/airfoil/naca4digit
	double M, P, XX;
	int temp;
	double a0 = 0.2969, a1 = -0.126, a2 = -0.3516, a3 = 0.2843, a4 = -0.1015;
	// or a4 = -0.1036 for a closed trailing edge
	stringstream geek(NACA);
	geek >> temp;
	M = (temp / 1000) / 100.0;
	P = ((temp - (temp / 1000) * 1000) / 100) / 10.0;
	XX = (temp % 100) / 100.0;
	// cout<< fixed << setprecision(2) << "M = " << M << "\tP = " << P << "\tXX = " << XX << endl;
	double beta = 0, theta, x, yc, yt;
	vector <double> xu, xl, yu, yl;
	for (int i = 0; i < N/2 + 1; i++)
	{
		x = (1 - cos(beta)) / 2;
		yt = XX / 0.2 * (a0 * pow(x, 0.5) + a1 * x + a2 * pow(x, 2) + a3 * pow(x, 3) + a4 * pow(x, 4));
		if (x < P)
		{
			theta = atan(2 * M / (P * P) * (P - x));
			yc = M / (P * P) * (2 * P * x - x * x);
			xu.push_back(x - yt * sin(theta));
			yu.push_back(yc + yt * cos(theta));
			xl.push_back(x + yt * sin(theta));
			yl.push_back(yc - yt * cos(theta));
		}
		else
		{
			theta = atan(2 * M / ((1 - P) * (1 - P)) * (P - x));
			yc = M / ((1 - P) * (1 - P)) * (1 - 2 * P + 2 * P * x - x * x);
			xu.push_back(x - yt * sin(theta));
			yu.push_back(yc + yt * cos(theta));
			xl.push_back(x + yt * sin(theta));
			yl.push_back(yc - yt * cos(theta));
		}
		beta += pi / (N / 2);
	}

	ofstream outfile(namefile);
	outfile << "NACA " << NACA;
	for (int i = xu.size() - 1; i > 0; i--)
	{
		outfile << "\n" << fixed << setprecision(6) << xu[i] << "     " << yu[i];
	}
	for (int i = 0; i < xl.size(); i++)
	{
		outfile << "\n" << fixed << setprecision(6) << xl[i] << "     " << yl[i];
	}
	outfile.close();
}
void export_5_digit(string NACA, int N, string namefile)
{
	// http://www.airfoiltools.com/airfoil/naca5digit
	double L, P, Q, XX;
	int temp, Digits;
	double a0 = 0.2969, a1 = -0.126, a2 = -0.3516, a3 = 0.2843, a4 = -0.1015;
	// or a4 = -0.1036 for a closed trailing edge
	stringstream geek(NACA);
	geek >> temp;
	L = (temp / 10000) * 3.0 /20;
	P = ((temp / 1000) % 10) / 20.0;
	Q = (temp / 100) % 10;
	XX = (temp % 100) / 100.0;
	Digits = (temp / 100) % 100;
	// cout << fixed << setprecision(2) << "L = " << L << "\tP = " << P << "\tQ = " << Q << "\tXX = " << XX << "\tDigits = " << Digits << endl;
	double beta = 0, theta, x, yc, yt;
	double r, k1, k2_k1=-1;
	vector <double> xu, xl, yu, yl;
	switch (Digits)
	{
	case 10:
		r = 0.0580;
		k1 = 361.4;
		break;
	case 20:
		r = 0.1260;
		k1 = 51.64;
		break;
	case 30:
		r = 0.2025;
		k1 = 15.957;
		break;
	case 40:
		r = 0.29;
		k1 = 6.643;
		break;
	case 50:
		r = 0.391;
		k1 = 3.23;
		break;
	case 21:
		r = 0.13;
		k1 = 51.99;
		k2_k1 = 0.000764;
		break;
	case 31:
		r = 0.217;
		k1 = 15.793;
		k2_k1 = 0.00677;
		break;
	case 41:
		r = 0.318;
		k1 = 6.52;
		k2_k1 = 0.0303;
		break;
	case 51:
		r = 0.441;
		k1 = 3.191;
		k2_k1 = 0.1355;
		break;
	default:
		cout << "\nError digits" << endl;
		break;
	}
	// cout << "r = " << r << "\tk1 = " << k1 << "\tk2_1 = " << k2_k1 << endl;
	for (int i = 0; i < N / 2 + 1; i++)
	{
		x = (1 - cos(beta)) / 2;
		yt = XX / 0.2 * (a0 * pow(x, 0.5) + a1 * x + a2 * pow(x, 2) + a3 * pow(x, 3) + a4 * pow(x, 4));
		if (x < r)
		{
			if (Q == 0) {
				theta = atan(k1 / 6 * (3 * x * x - 6 * r * x + r * r * (3 - r)));
				yc = k1 / 6 * (pow(x, 3) - 3 * r * x * x + r * r * (3 - r) * x);
			}
			else if(Q == 1){
				theta = atan(k1 / 6 * (3 * (x - r) * (x - r) - k2_k1 * pow(1 - r, 3) - pow(r, 3)));
				yc = k1 / 6 * (pow(x - r, 3) - k2_k1 * pow(1 - r, 3) * x - pow(r, 3) * x + pow(r, 3));
			}
			xu.push_back(x - yt * sin(theta));
			yu.push_back(yc + yt * cos(theta));
			xl.push_back(x + yt * sin(theta));
			yl.push_back(yc - yt * cos(theta));
		}
		else
		{
			if (Q == 0) {
				theta = atan(-k1 / 6 * pow(r, 3));
				yc = k1 * pow(r, 3) / 6 * (1 - x);
			}
			else if (Q == 1) {
				theta = atan(k1 / 6 * (3*k2_k1 *pow(x-r,2) - k2_k1 * pow(1 - r, 3) - pow(r, 3)));
				yc = k1 / 6 * (k2_k1 * pow(x - r, 3) - k2_k1 * pow(1 - r, 3) * x - pow(r, 3) * x + pow(r, 3));
			}
			xu.push_back(x - yt * sin(theta));
			yu.push_back(yc + yt * cos(theta));
			xl.push_back(x + yt * sin(theta));
			yl.push_back(yc - yt * cos(theta));
		}
		beta += pi / (N / 2);
	}

	ofstream outfile(namefile);
	outfile << "NACA " << NACA;
	for (int i = xu.size() - 1; i > 0; i--)
	{
		outfile << "\n" << fixed << setprecision(6) << xu[i] << "     " << yu[i];
	}
	for (int i = 0; i < xl.size(); i++)
	{
		outfile << "\n" << fixed << setprecision(6) << xl[i] << "     " << yl[i];
	}
	outfile.close();
}