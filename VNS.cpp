#include "Reader.h"
#include "VNS.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <random>

using namespace std;

VNS::VNS(int Itermax, int Shakeiter, int Tabutenure, int Neighbormax, Reader* data, Evaluate* eval)
{
    IterMax = Itermax;
    ShakeIter = Shakeiter;
    TabuTenure = Tabutenure;
    NeighborMax = Neighbormax;
    Data = data;
    Sol = eval;
}

bool VNS::NS1(int* skd)
{
    bool opt = false;
    int NS1Iter = 0;
    int* sub_schedule = new int[Data->getn() + 1];
    int* sub_skd = new int[Data->getn() + 1];
    copy(skd, skd + Data->getn(), sub_skd);
    copy(skd, skd + Data->getn(), sub_schedule);
    default_random_engine random(time(NULL));
    uniform_int_distribution<int> dis1(0, Data->getn() - 1);
    double Best = Sol->Calculation(skd);
    while (NS1Iter < 10)
    {
        double LocalBest = 10000000.0;
        int* temp = new int[Data->getn()];
        copy(sub_schedule, sub_schedule + Data->getn(), temp);
        int TabuNS1[100][100] = { 0 };
        int Iter = 0;
        while (Iter < NeighborMax)
        {
            copy(temp, temp + Data->getn(), sub_skd);
            int i = dis1(random);
            int j = 0;
            while (true)
            {
                j = dis1(random);
                if (i != j)
                {
                    break;
                }
            }
            int mid = sub_skd[i];
            sub_skd[i] = sub_skd[j];
            sub_skd[j] = mid;
            double TempVal = Sol->Calculation(sub_skd);
            if (TempVal < Best)
            {
                Best = TempVal;
                copy(sub_skd, sub_skd + Data->getn(), sub_schedule);
                copy(sub_skd, sub_skd + Data->getn(), skd);
                TabuNS1[i][j] += Iter + TabuTenure;
                TabuNS1[j][i] += Iter + TabuTenure;
                opt = true;
            }
            else if (TempVal < LocalBest && TabuNS1[i][j] < Iter)
            {
                LocalBest = TempVal;
                if (TempVal < Sol->Calculation(sub_schedule))
                {
                    copy(sub_skd, sub_skd + Data->getn(), sub_schedule);
                    TabuNS1[i][j] += Iter + TabuTenure;
                    TabuNS1[j][i] += Iter + TabuTenure;
                }
            }
            Iter += 1;
        }
        delete[] temp;
        temp = nullptr;
        NS1Iter += 1;
    }
    delete[] sub_schedule;
    sub_schedule = nullptr;
    delete[] sub_skd;
    sub_skd = nullptr;
    return opt;
}

bool VNS::NS2(int* skd)
{
    int NS2Iter = 0;
    bool opt = false;
    int* sub_schedule = new int[Data->getn() + 1];
    int* sub_skd = new int[Data->getn() + 1];
    copy(skd, skd + Data->getn(), sub_schedule);
    copy(skd, skd + Data->getn(), sub_skd);
    double Best = Sol->Calculation(skd);
    default_random_engine random(time(NULL));
    uniform_int_distribution<int> dis1(0, Data->getn() - 1);
    while (NS2Iter < 10)
    {
        int Iter = 0;
        int TabuNS2[100][100] = { 0 };
        double LocalBest = 10000000.0;
        int* temp = new int[Data->getn()];
        copy(sub_schedule, sub_schedule + Data->getn(), temp);
        while (Iter < NeighborMax)
        {
            copy(temp, temp + Data->getn(), sub_skd);
            int i = dis1(random);
            int j = 0;
            while (true)
            {
                j = dis1(random);
                if (i != j)
                {
                    break;
                }
            }
            int mid = sub_skd[i];
            while (i < Data->getn() - 1)
            {
                sub_skd[i] = sub_skd[i + 1];
                i++;
            }
            int m = Data->getn() - 1;
            while (m >= j)
            {
                sub_skd[m + 1] = sub_skd[m];
                m--;
            }
            sub_skd[j] = mid;
            double TempVal = Sol->Calculation(sub_skd);
            if (TempVal < Best)
            {
                Best = TempVal;
                copy(sub_skd, sub_skd + Data->getn(), sub_schedule);
                copy(sub_skd, sub_skd + Data->getn(), skd);
                TabuNS2[i][j] += Iter + TabuTenure;
                TabuNS2[j][i] += Iter + TabuTenure;
                opt = true;
            }
            else if (TempVal < LocalBest && TabuNS2[i][j] < Iter)
            {
                LocalBest = TempVal;
                if (TempVal < Sol->Calculation(sub_schedule))
                {
                    copy(sub_skd, sub_skd + Data->getn(), sub_schedule);
                    TabuNS2[i][j] += Iter + TabuTenure;
                    TabuNS2[j][i] += Iter + TabuTenure;
                }
            }
            Iter += 1;
        }
        delete[] temp;
        temp = nullptr;
        NS2Iter += 1;
    }
    delete[] sub_schedule;
    sub_schedule = nullptr;
    delete[] sub_skd;
    sub_skd = nullptr;
    return opt;
}

bool VNS::NS3(int* skd)
{ //2-Exchange
    bool opt = false;
    int* sub_schedule = new int[Data->getn() + 1];
    copy(skd, skd + Data->getn(), sub_schedule);
    double Best = Sol->Calculation(skd);
    int i = 0;
    default_random_engine random(time(NULL));
    uniform_int_distribution<int> dis1(0, Data->getn() - 1);
    while (i < Data->getn())
    {
        int j = i + 1;
        while (j < Data->getn())
        {
            int p = 101;
            int q = 101;
            while (true)
            {
                if (p != i && q != j && p != j && q != i && p != q)
                {
                    break;
                }
                p = dis1(random);
                q = dis1(random);
            }
            int temp1 = sub_schedule[p];
            int temp2 = sub_schedule[q];
            sub_schedule[p] = sub_schedule[i];
            sub_schedule[i] = temp1;
            sub_schedule[q] = sub_schedule[j];
            sub_schedule[j] = temp2;
            double tempV = Sol->Calculation(sub_schedule);
            if (tempV < Best)
            {
                copy(sub_schedule, sub_schedule + Data->getn(), skd);
                i = Data->getn();
                j = Data->getn();
                opt = true;
            }
            j += 1;
        }
        i += 1;
    }
    delete[] sub_schedule;
    sub_schedule = nullptr;
    return opt;
}

bool VNS::NS4(int* skd)
{ //two-opt-swap
    bool opt = false;
    int Iter = 0;
    int* temp = new int[Data->getn()];
    copy(skd, skd + Data->getn(), temp);
    double Best = Sol->Calculation(skd);
    default_random_engine random(time(NULL));
    uniform_int_distribution<int> dis1(0, Data->getn() - 1);
    double LocalBest = 1000000.0;
    while (Iter < NeighborMax)
    {
        int i = dis1(random);
        int j = 0;
        while (true)
        {
            j = dis1(random);
            if (i != j)
            {
                break;
            }
        }
        if (j < i)
        {
            int mid = i;
            i = j;
            j = mid;
        }
        int mid = 0;
        for (int m = i; m < (i + j) / 2 + 1; m++)
        {
            mid = temp[i + j - m];
            temp[i + j - m] = temp[m];
            temp[m] = mid;
        }
        double tempV = Sol->Calculation(temp);
        if (tempV < LocalBest)
        {
            LocalBest = tempV;
        }
        if (tempV < Best)
        {
            Best = tempV;
            copy(temp, temp + Data->getn(), skd);
            opt = true;
        }
        Iter += 1;
    }
    if ((LocalBest - Best) / Best < 0.05)
    {
        copy(temp, temp + Data->getn(), skd);
    }
    delete[] temp;
    temp = nullptr;
    return opt;
}

void VNS::VNSearch(int* skd)
{
    int Iteration = 0;
    bool opt = false;
    while (Iteration < IterMax)
    {
        while (opt)
        {
            opt = NS4(skd);
        }
        opt = true;
        while (opt)
        {
            opt = NS1(skd);
        }
        opt = true;
        while (opt)
        {
            opt = NS2(skd);
        }
        opt = true;
        while (opt)
        {
            opt = NS3(skd);
        }
        opt = true;
        Iteration += 1;
    }
}