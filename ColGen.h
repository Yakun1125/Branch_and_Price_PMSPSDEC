#pragma once
#include "Reader.h"
#include "Evaluation.h"
#include "VNS.h"
#include <ilcplex/ilocplex.h>

using namespace std;

class ColumnGen{
private:
	Reader* Data;
public:
	ColumnGen(Reader* Data);
	double RunColGen(vector<vector<int>>& schedule, vector<double>& skd_value, vector<double>& getY, int& sol_size, int** branch_info);
	void Initialize(Evaluate* sol, VNS* search, vector<vector<int>>& skd, vector<double>& cost);
	void Add_Column(IloEnv& env, IloNumVarArray& var, IloRangeArray& cons, IloObjective& obj, vector<vector<int>>& schedule, vector<double>& skd_value, int& sol_size, int sol_num, IloCplex& SubP, IloIntVarArray& a, IloNumVarArray& T, IloNumVarArray& E, IloNumVarArray& S);
};