
#include "Header.h"
#include "EAF.h"

int main()
{
    int mode;
    string check;
    do
    {
        system("cls");
        cout << "\tCopyright (TM) 2020 by Nguyen Van Duc - All Rights Reserved" << endl;
        cout << "\t\tEAF v.2.0  |  Ha Noi - 12/09/2020\n" << endl;

        cout << "Mode:\n";
        cout << "Export airfoil: \t1\n";
        cout << "Convert airfoil:\t2\n\n";

        do
        {
            cout << "Enter the mode: ";
            cin >> mode;
        } while (mode!=0 && mode!=1 && mode!=2);

        if (mode == 1) {
            export_airfoil_data(input_airfoil(print_database()));
            system("pause");
        }
        else if (mode == 2) {
            convert_airfoil_data();
            system("pause");
        }
        
        cout << "\nEnter 0 to exit or any key to continue: ";
        cin >> check;
    } while (check!="0");

    return 0;
}

