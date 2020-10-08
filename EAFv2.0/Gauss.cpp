#include "Header.h"
#include "EAF.h"
double nrsovler(vector <double> ca, vector <double> cx, double res, int maxloop)
{
    vector <double> dca;
    vector <double> dcx;
    for (int i = 0; i < 6; i++)
    {
        dca.push_back(ca[i] * cx[i]);
        dcx.push_back(cx[i] - 1);
    }
    double x0 = 1;
    double x1 = x0 - evalf(ca, cx, x0) / evalf(dca, dcx, x0);
    x0 = x1;
    // return          x1;
    int iCnt = 0;
    while (fabs(evalf(ca, cx, x1)) > res && iCnt < maxloop)
    {
        x1 = x0 - evalf(ca, cx, x0) / evalf(dca, dcx, x0);
        // cout<<"x1 = "<<x1<<endl;
        x0 = x1;
        iCnt++;
    }
    return x1;
}
double evalf(vector <double>ca, vector <double>cx, double x)
{
    double val = 0;
    for (int i = 0; i < ca.size(); i++)
    {
        if (ca[i] != 0)
        {
            val += ca[i] * pow(x, cx[i]);
        }
    }
    return val;
}