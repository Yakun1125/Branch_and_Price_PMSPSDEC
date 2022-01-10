#pragma once
#include <string>
using namespace std;

class Reader
{
private:
    int m;
    int n;
    double e0;
    double* aTime; //basic processing time
    double* hTime; //deteriorating threshold
    double* bTime; //extra processing time
    double* dTime; //due time
    double* qTime; //unit energy consumption
    double Alpha; //weight factor
    double Beta;

public:
    Reader(string add);
    ~Reader();
    int getm();
    int getn();
    double gete();
    double a(int num);
    double h(int num);
    double b(int num);
    double d(int num);
    double q(int num);
    double alpha();
    double beta();
};