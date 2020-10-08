#include "Header.h"
#include "EAF.h"
string input_airfoil(vector<string>data)
{
	string NACA;
	bool out = false;
	do
	{
		cout << "\nEnter the name of the airfoil: ";
		cin >> NACA;
		for (int i = 0; i < data.size(); i++)
		{
			if (NACA == data[i])
			{
				out = true;
				break;
			}
		}
	} while (!out);
	return NACA;
}