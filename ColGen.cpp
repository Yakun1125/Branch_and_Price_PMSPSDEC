#include "ColGen.h"
#include "Reader.h"
#include "Evaluation.h"
#include "VNS.h"
#include <ilcplex/ilocplex.h>
#include <algorithm>

ColumnGen::ColumnGen(Reader* data) {
	Data = data;
}

double ColumnGen::RunColGen(vector<vector<int>>& schedule, vector<double>& skd_value, vector<double>& getY, int& sol_size, int** branch_info) {
	bool cycle = true;

	IloEnv env;
	try {
		IloNumVarArray var(env);
		IloModel Mp(env);
		IloObjective obj = IloAdd(Mp, IloMinimize(env));
		IloRangeArray cons = IloAdd(Mp, IloRangeArray(env, Data->getn() + 1, 1, 1));
		cons[Data->getn()].setBounds(0, Data->getm());
		for (int m = 0; m < sol_size; m++) {
			IloNumColumn col = obj(skd_value[m]);
			for (int n = 0; n < schedule[m].size(); n++) {
				col += cons[schedule[m][n]](1);
			}
			col += cons[Data->getn()](1);
			var.add(IloNumVar(col, 0, 1));
			col.end();
		}
		IloCplex LMP(Mp);
		LMP.setOut(env.getNullStream());
		//LMP.exportModel("MasterProblem.lp");

		/// Define Subproblem///
		IloModel Sub(env);
		IloObjective ReducedCost = IloAdd(Sub, IloMinimize(env));

		IloArray<IloIntVarArray> x(env, Data->getn());
		IloIntVarArray a(env, Data->getn());
		IloNumVarArray LiT(env, Data->getn());
		IloNumVarArray LiE(env, Data->getn());
		IloNumVarArray S(env, Data->getn());
		IloIntVarArray Z(env, Data->getn());
		IloIntVarArray U(env, Data->getn());
		IloNumVarArray P(env, Data->getn());
		IloNumVarArray T(env, Data->getn());
		IloNumVarArray C(env, Data->getn());
		IloNumVarArray e(env, Data->getn());
		IloNumVarArray E(env, Data->getn());

		for (int k = 0; k < Data->getn(); k++) {
			x[k] = IloIntVarArray(env, Data->getn(), 0, 1);
			a[k] = IloIntVar(env, 0, 1);
			LiT[k] = IloNumVar(env, 0, 50000);
			LiE[k] = IloNumVar(env, 0, 50000);
			S[k] = IloNumVar(env, 0, 50000);
			Z[k] = IloIntVar(env, 0, 1);
			U[k] = IloIntVar(env, 0, 1);
			P[k] = IloNumVar(env, Data->a(k), Data->a(k) + Data->b(k));
			T[k] = IloNumVar(env, 0, 50000);
			C[k] = IloNumVar(env, 0, 50000);
			e[k] = IloNumVar(env, 0, 50000);
			E[k] = IloNumVar(env, 0, 50000);
		}
		double M = 50000;

		for (int i = 0; i < Data->getn(); i++) {
			for (int j = 0; j < Data->getn(); j++) {
				if (i != j) {
					Sub.add(a[i] + a[j] - x[i][j] - x[j][i] <= 1);
				}
				Sub.add(C[i] - S[j] + M * x[i][j] <= M);
			}

			Sub.add(S[i] + P[i] <= C[i]);
			Sub.add(T[i] + Data->d(i) >= C[i]);
			Sub.add(S[i] + M * Z[i] <= Data->h(i) + M);
			Sub.add(U[i] + Z[i] <= 1);
			Sub.add(Data->b(i) * U[i] + P[i] >= Data->a(i) + Data->b(i));
			Sub.add(Data->gete() * S[i] + M * U[i] - e[i] <= Data->gete() * Data->h(i) + M);
			Sub.add(Data->q(i) * P[i] + e[i] <= E[i]);
			Sub.add(LiT[i] <= M * a[i]);
			Sub.add(T[i] + M * a[i] - LiT[i] <= M);
			Sub.add(LiE[i] <= M * a[i]);
			Sub.add(E[i] + M * a[i] - LiE[i] <= M);
			Sub.add(LiT[i] <= T[i]);
			Sub.add(LiE[i] <= E[i]);

		}
		for (int i = 0; i < Data->getn(); i++) {
			for (int j = 0; j < Data->getn(); j++) {
				if (branch_info[i][j] == 1) {//forbidden job (j+i+1) in the front of job i
					Sub.add(x[i][j] == 0);
				}
				else if (branch_info[i][j] == 2) {
					Sub.add(x[i][j] == 1);
					IloIntVar bi(env, 0, 1);
					Sub.add(a[i] + a[j] == 2 * bi);
				}
			}
		}
		IloCplex SubP(Sub);
		//SubP.setParam(IloCplex::Param::Emphasis::Numerical, true);
		SubP.setOut(env.getNullStream());

		IloNumArray pi(env, Data->getn() + 1);
		double CGvalue = -1e10;
		while (cycle) {
			cycle = false;
			if (!LMP.solve()) {
				cout << "Master Problem Infeasible..........................." << endl;
				return 1e10;
			}
			CGvalue = LMP.getObjValue();
			IloExpr SubObj(env);
			for (int n = 0; n < Data->getn(); n++) {
				pi[n] = LMP.getDual(cons[n]);
				SubObj += Data->alpha() * LiT[n] + Data->beta() * LiE[n] - pi[n] * a[n];
			}
			pi[Data->getn()] = LMP.getDual(cons[Data->getn()]);
			SubObj += -pi[Data->getn()];
			ReducedCost.setExpr(SubObj); SubObj.end();
			SubP.solve();
			if (SubP.getValue(ReducedCost) < -1e-05) {
				cycle = true;
				if ((sol_size + SubP.getSolnPoolNsolns()) > skd_value.size()) {
					schedule.resize(sol_size + SubP.getSolnPoolNsolns());
					for (int m = sol_size; m < sol_size + SubP.getSolnPoolNsolns(); m++) {
						schedule[m].reserve(Data->getn());
					}
					skd_value.resize(sol_size + SubP.getSolnPoolNsolns());
				}
				Add_Column(env, var, cons, obj, schedule, skd_value, sol_size, SubP.getSolnPoolNsolns(), SubP, a, T, E, S);
			}
			cout << "Pricing Problem: " << SubP.getValue(ReducedCost) << "   LMP: " << CGvalue << "   Skd Size:" << sol_size << endl;
		}
		getY.resize(sol_size);
		for (int n = 0; n < sol_size; n++) {
			getY[n] = LMP.getValue(var[n]);
		}
		return CGvalue;
	}
	catch (IloException& e) {
		cerr << "Concert exception caught: " << e << endl;
	}
	env.end();
	return 1e10;
}

void ColumnGen::Add_Column(IloEnv& env, IloNumVarArray& var, IloRangeArray& cons, IloObjective& obj, vector<vector<int>>& schedule, vector<double>& skd_value, int& sol_size, int sol_num, IloCplex& SubP, IloIntVarArray& a, IloNumVarArray& T, IloNumVarArray& E, IloNumVarArray& S) {
	for (int i = 0; i < sol_num; i++) {
		if (SubP.getObjValue(i) < -1e-05) {
			IloNumColumn col(env);
			vector<double> Start(Data->getn(),1000000);
			int count = 0;
			for (int j = 0; j < Data->getn(); j++) {
				if (SubP.getValue(a[j], i) > 0.5) {
					skd_value[sol_size] += Data->alpha() * SubP.getValue(T[j], i) + Data->beta() * SubP.getValue(E[j], i);
					Start[j] = SubP.getValue(S[j], i);
					col += cons[j](1);
					count += 1;
				}
			}
			for (int j = 0; j < count; j++) {
				int job = min_element(Start.begin(), Start.end()) - Start.begin();
				schedule[sol_size].push_back(job);
				Start[job] = 1000000;
			}
			vector<double>().swap(Start);
			col += cons[Data->getn()](1);
			col += obj(skd_value[sol_size]);
			var.add(IloNumVar(col, 0, 1)); col.end();
			sol_size += 1;
		}
	}
}

void ColumnGen::Initialize(Evaluate* sol, VNS* search, vector<vector<int>>& skd, vector<double>& cost) {
	int* Schedule = new int[Data->getn()];

	sol->Construction(Schedule);
	//search->VNSearch(Schedule);
    sol->Output(cost, skd, Schedule);

	delete[] Schedule; Schedule = nullptr;

}