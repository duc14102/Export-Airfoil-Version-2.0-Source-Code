#include "Header.h"
#include "EAF.h"
#include "Matrix.h"
#include "Vector_matrix.h"

void convert_airfoil_data() {
	string filename_in, filename_out;
	string name_airfoil;
	bool check = false;
	int N;
	vector <double> x0, y0;

	do
	{
		cout << "\nEnter the filename of the input file: ";
		cin >> filename_in;
		ifstream infile(filename_in);
		if (!infile) cout << "\tFile name error";
		else
		{
			check = true;

			double temp;
			getline(infile, name_airfoil);
			while (!infile.eof())
			{
				infile >> temp;
				x0.push_back(fabs(temp));
				infile >> temp;
				y0.push_back(temp);
			}
			infile.close();

			cout << "\tInput successful!\n";
		}
	} while (!check);

	cout << "\nEnter the filename of the output file : ";
	cin >> filename_out;
	cout << "Enter the number of points: ";
	cin >> N;

	convert_parsec(x0,y0,N,name_airfoil,filename_out);

    cout << "\nConvert to file \"" << filename_out << "\"...Done!\n" << endl;
}

void convert_parsec(vector <double> x0, vector <double> y0, int N, string name, string filename_out) {
    vectorMT a = Parsec(x0, y0);
    //cout << "p:\n" << p;
    vectorMT x(N / 2 + 1), yu(N / 2 + 1), yl(N / 2 + 1);
    double beta = 0;
    for (int i = 1; i <= N / 2 + 1; i++)
    {
        x(i) = (1 - cos(beta)) / 2;
        //coordinate array
        yu(i) = a(1) * pow(x(i), 0.5) + a(2) * pow(x(i), 1.5) + a(3) * pow(x(i), 2.5) + a(4) * pow(x(i), 3.5) + a(5) * pow(x(i), 4.5) + a(6) * pow(x(i), 5.5);
        yl(i) = a(7) * pow(x(i), 0.5) + a(8) * pow(x(i), 1.5) + a(9) * pow(x(i), 2.5) + a(10) * pow(x(i), 3.5) + a(11) * pow(x(i), 4.5) + a(12) * pow(x(i), 5.5);
        beta += pi / (N / 2);
    }
    ofstream outfile(filename_out);
    outfile << name;
    for (int i = rows(x); i >= 2; i--)
    {
        outfile << "\n" << fixed << setprecision(6) << x(i) << "     " << yu(i);
    }
    for (int i = 1; i <= rows(x); i++)
    {
        outfile << "\n" << fixed << setprecision(6) << x(i) << "     " << yl(i);
    }
    outfile.close();
}

vectorMT Parsec(vector <double>x0, vector <double> y0)
{
    int row = x0.size() / 2 + 1;
    int colum = 6;
    matrix Aup(row, colum), bup(row, 1);
    matrix Alo(row, colum), blo(row, 1);

    for (int i = 1; i <= row; i++)
    {
        for (int j = 1; j <= colum; j++)
        {
            Aup(i, j) = pow(fabs(x0[i-1]), j - 0.5);
            Alo(i, j) = pow(fabs(x0[row + i - 1 - 1]), j - 0.5);
        }
    }

    for (int i = 1; i <= row; i++)
    {
        bup(i, 1) = y0[i - 1];
        blo(i, 1) = y0[row + i - 1 - 1];
    }

    vectorMT WXup = LeastSquares(Aup, bup);
    //cout << "up\n" << WXup;
    vectorMT WXlo = LeastSquares(Alo, blo);
    //cout << "lo\n" << WXlo;
    vectorMT a(12);
    for (int i = 1; i <= rows(WXup); i++)
    {
        a(i) = WXup(i);
        a(i+6) = WXlo(i);
    }
    return a;
}