#pragma once
#include <vector>
#include "Reader.h"
using namespace std;

class Evaluate
{
private:
    Reader* Data;

public:
    Evaluate(Reader* Data);
    void Construction(int* schedule);
    double Calculation(int* schedule);
    void Output(vector<double>& cost, vector<vector<int>>& schedule, int* skd);
};