#pragma once
#include "Reader.h"
#include "LinkTree.h"
#include "ColGen.h"
#include <ilcplex/ilocplex.h>
#include <vector>

using namespace std;

class BandB {
private:
	Reader* Data;
	ColumnGen* Col;
	double lowerbound;
	double upperbound;
public:
	BandB(Reader* Data, ColumnGen* Col, double lowerbound, double upperbound);
	bool runBB(vector<vector<int>>& schedule, vector<vector<int>>& best_schedule, vector<double>& skd_value, int& sol_size, LinkTree* branch, int** branch_info, int depth);
	void Branching(int** branch_info, bool recur, LinkTree* branch);
	double Single_Eval(vector<int>& new_skd);
};