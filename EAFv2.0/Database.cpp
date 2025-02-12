#include "Header.h"
vector<string> print_database()
{
    vector<string>data4digit = { "0006","0008","0009","0010","0012","0015","0018","0021","0024","1408","1410",
                                 "1412","2408","2410","2411","2412","2414","2415","2418","2421","2424","4412",
                                 "4415","4418","4421","4424","6409","6412" };
    vector<string>data5digit = { "21012","22112","22012","23012","23015","23018","23021",
                                 "23024","23112","24012","24112","25012","25112" };
    vector<string>data = { "0006","0008","0009","0010","0012","0015","0018","0021","0024","1408","1410",
                           "1412","2408","2410","2411","2412","2414","2415","2418","2421","2424","4412",
                           "4415","4418","4421","4424","6409","6412","21012","22112","22012","23012",
                           "23015","23018","23021","23024","23112","24012","24112","25012","25112" };
	cout << "\n\t\t\tAirfoil database" << endl;
    int t=0;
	cout << "\nNACA 4 digits:"<<endl;
    for (int i = 0; i < data4digit.size(); i++)
    {
        cout << "NACA " << data4digit[i] << "\t";
        t++;
        if (t == 5) {
            t = 0;
            cout << endl;
        }
    }
    cout << endl;
    t = 0;
    cout << "\nNACA 5 digits:" << endl;
    for (int i = 0; i < data5digit.size(); i++)
    {
        cout << "NACA " << data5digit[i] << "\t";
        t++;
        if (t == 5) {
            t = 0;
            cout << endl;
        }
    }
    cout << endl;
    return data;
}