#include "Evaluation.h"
#include "Reader.h"
#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

Evaluate::Evaluate(Reader* data)
{
    Data = data;
}

void Evaluate::Construction(int* skd)
{
    double val = 10000.0;
    int start = 0;
    for (int i = 0; i < Data->getn(); i++)
    {
        if (Data->d(i) < val)
        {
            val = Data->d(i);
            start = i;
        }
    }
    int cur_point = 0;
    *(skd + cur_point) = start;
    cur_point += 1;
    double m1 = 0.6;
    double m2 = 0.2;
    double m3 = 0.15;
    double m4 = 0.05;
    vector<double> seq(Data->getn());
    for (int i = 0; i < Data->getn(); i++)
    {
        seq[i] = m1 * Data->d(i) - m2 * Data->q(i) + m3 * Data->h(i) + m4 * (Data->a(i) + Data->b(i));
    }
    for (int i = 0; i < Data->getn(); i++)
    {
        int pos = min_element(seq.begin(), seq.end()) - seq.begin();
        if (pos != start)
        {
            *(skd + cur_point) = pos;
            cur_point += 1;
        }
        seq[pos] = 10000;
    }
    for (int i = 0; i < Data->getn(); i++)
    {
        cout << *(skd + i) << ", ";
    }
    cout << endl;
}

double Evaluate::Calculation(int* skd)
{
    int point = 0;
    vector<double> start(Data->getm());
    vector<double> cost(Data->getm());
    double totalcost = 0;
    for (int i = 0; i < Data->getm(); i++)
    {
        start[i] = Data->a(skd[i]);
        cost[i] = Data->beta() * Data->a(skd[i]) * Data->q(skd[i]);
    }
    int i = Data->getm();
    double tardiness = 0.0;
    double energy = 0.0;
    do
    {
        point = min_element(start.begin(), start.end()) - start.begin();
        if (start[point] <= Data->h(skd[i]))
        {
            start[point] = start[point] + Data->a(skd[i]);
            tardiness = Data->alpha() * max(start[point] - Data->d(skd[i]), 0.0);
            energy = Data->beta() * Data->a(skd[i]) * Data->q(skd[i]);
            cost[point] = cost[point] + energy + tardiness;
            i += 1;
        }
        else
        {
            if (Data->alpha() * max(start[point] + Data->b(skd[i]) + Data->a(skd[i]) - Data->d(skd[i]), 0.0) + Data->beta() * Data->b(skd[i]) * Data->q(skd[i]) <= Data->beta() * Data->gete() * (start[point] - Data->h(skd[i])) + Data->alpha() * max(start[point] + Data->a(skd[i]) - Data->d(skd[i]), 0.0))
            {
                start[point] = start[point] + Data->a(skd[i]) + Data->b(skd[i]);
                tardiness = Data->alpha() * max(start[point] - Data->d(skd[i]), 0.0);
                energy = Data->beta() * Data->q(skd[i]) * (Data->a(skd[i]) + Data->b(skd[i]));
                cost[point] = cost[point] + energy + tardiness;
                i += 1;
            }
            else
            {
                energy = Data->beta() * (Data->gete() * (start[point] - Data->h(skd[i])) + Data->q(skd[i]) * Data->a(skd[i]));
                start[point] = start[point] + Data->a(skd[i]);
                tardiness = Data->alpha() * max(start[point] - Data->d(skd[i]), 0.0);
                cost[point] = cost[point] + energy + tardiness;
                i += 1;
            }
        }
    } while (i != Data->getn());
    for (int j = 0; j < Data->getm(); j++)
    {
        totalcost += cost[j];
    }
    return totalcost;
}

void Evaluate::Output(vector<double>& cost, vector<vector<int>>& schedule, int* skd) {
    int point = 0;
    vector<double> start(Data->getm());
    for (int i = 0; i < Data->getm(); i++)
    {
        schedule[i].push_back(skd[i]);
        start[i] = Data->a(skd[i]);
        cost[i] = Data->beta() * Data->a(skd[i]) * Data->q(skd[i]);
    }
    int i = Data->getm();
    double tardiness = 0.0;
    double energy = 0.0;
    do
    {
        point = min_element(start.begin(), start.end()) - start.begin();
        schedule[point].push_back(skd[i]);
        if (start[point] <= Data->h(skd[i]))
        {
            start[point] = start[point] + Data->a(skd[i]);
            tardiness = Data->alpha() * max(start[point] - Data->d(skd[i]), 0.0);
            energy = Data->beta() * Data->a(skd[i]) * Data->q(skd[i]);
            cost[point] = cost[point] + energy + tardiness;
            i += 1;
        }
        else
        {
            if (Data->alpha() * max(start[point] + Data->b(skd[i]) + Data->a(skd[i]) - Data->d(skd[i]), 0.0) + Data->beta() * Data->b(skd[i]) * Data->q(skd[i]) <= Data->beta() * Data->gete() * (start[point] - Data->h(skd[i])) + Data->alpha() * max(start[point] + Data->a(skd[i]) - Data->d(skd[i]), 0.0))
            {
                start[point] = start[point] + Data->a(skd[i]) + Data->b(skd[i]);
                tardiness = Data->alpha() * max(start[point] - Data->d(skd[i]), 0.0);
                energy = Data->beta() * Data->q(skd[i]) * (Data->a(skd[i]) + Data->b(skd[i]));
                cost[point] = cost[point] + energy + tardiness;
                i += 1;
            }
            else
            {
                energy = Data->beta() * (Data->gete() * (start[point] - Data->h(skd[i])) + Data->q(skd[i]) * Data->a(skd[i]));
                start[point] = start[point] + Data->a(skd[i]);
                tardiness = Data->alpha() * max(start[point] - Data->d(skd[i]), 0.0);
                cost[point] = cost[point] + energy + tardiness;
                i += 1;
            }
        }
    } while (i != Data->getn());

}