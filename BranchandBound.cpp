#include"BranchBound.h"
#include"Reader.h"
#include"LinkTree.h"
#include<cmath>
#include<algorithm>

BandB::BandB(Reader* data, ColumnGen* col, double Lowerbound, double Upperbound) {
	Data = data;
	Col = col;
	lowerbound = Lowerbound;
	upperbound = Upperbound;
}

bool BandB::runBB(vector<vector<int>>& schedule, vector<vector<int>>& best_schedule, vector<double>& skd_value, int& sol_size, LinkTree* branch, int** branch_info, int depth) {
	if (upperbound - lowerbound < 1e-12) {
		return true;
	}
	if (branch == NULL) {
		LinkTree* newNode = new LinkTree;
		newNode->father = NULL;
		newNode->son = NULL;
		newNode->toplevel = true;
		newNode->Node1 = -1;
		newNode->Node2 = -1;
		newNode->is_branch = 0;
		branch = newNode;
	}
	vector<double> getY(sol_size, -1);
	double CGvalue = Col->RunColGen(schedule, skd_value, getY, sol_size, branch_info);
	if (CGvalue > 1e10 - 1) {
		cout << "master problem infeasible, revise the initialization" << endl;//this may happen when branch forbidden arc
		return true;
	}
	branch->Node_lowerbound = CGvalue;

	if ((branch->father != NULL) && (branch->father->son != NULL) && branch->father->toplevel) {
		lowerbound = min(branch->Node_lowerbound, branch->father->son->Node_lowerbound);
		branch->toplevel = true;
	}
	else if (branch->father == NULL) {//root node
		lowerbound = CGvalue;
	}

	if (branch->Node_lowerbound > upperbound) {
		cout << "branch pruning | lower bound " << lowerbound << " | upper bound " << upperbound << " | gap "
			<< ((upperbound - lowerbound) / upperbound) << " | skd size " << sol_size << " | depth " << depth << endl;
		return true;
	}
	else {
		bool int_feasible = true;
		double** arc_val = new double* [Data->getn()]();
		for (int i = 0; i < Data->getn() ; i++) {
			arc_val[i] = new double[Data->getn()]();
		}
		for (int i = 0; i < sol_size; i++) {
			if (getY[i] > 1e-6) {
				for (int j = 0; j < schedule[i].size() - 1; j++) {
					for (int k = j + 1; k < schedule[i].size(); k++) {
						arc_val[schedule[i][j]][schedule[i][k]] += getY[i];
					}
				}
			}
		}
		double min_pair = 1e6;
		int BestNode1 = -1;
		int BestNode2 = -1;
		int branch_dire = 2;
		for (int i = 0; i < Data->getn(); i++) {
			for (int j = 0; j < Data->getn(); j++) {
				if ((arc_val[i][j] > 1e-6) && ((arc_val[i][j] < 0.999999) || (arc_val[i][j] > 1.0000001))) {
					int_feasible = false;
					double temp = abs(arc_val[i][j] - 0.5);
					if (temp < min_pair) {
						BestNode1 = i;
						BestNode2 = j;
						min_pair = temp;
					}
				}
			}
		}
		for (int i = 0; i < Data->getn() ; i++) {
			delete[] arc_val[i]; arc_val[i] = nullptr;
		}
		delete[] arc_val; arc_val = nullptr;

		if (int_feasible) {
			if (branch->Node_lowerbound <= upperbound) {
				upperbound = branch->Node_lowerbound;
				int best_num = 0;
				for (int m = 0; m < sol_size; m++) {
					if (getY[m] > 1e-6) {
						best_schedule[best_num].clear();
						best_schedule[best_num].assign(schedule[m].begin(), schedule[m].end());
						best_num += 1;
					}
				}
				cout << "incumbant feasible soluton | lower bound " << lowerbound << " | upper bound " << upperbound << " | gap "
					<< ((upperbound - lowerbound) / upperbound) << " | skd size " << sol_size << " | depth " << depth << endl;
			}
			else {
				cout << "integer feasible solution | lower bound " << lowerbound << " | upper bound " << upperbound << " | gap "
					<< ((upperbound - lowerbound) / upperbound) << " | skd size " << sol_size << " | depth " << depth << endl;
			}
			return true;
		}
		else {
			cout << "integer infeasible solution | lower bound " << lowerbound << " | upper bound " << upperbound << " | gap "
				<< ((upperbound - lowerbound) / upperbound) << " | skd size " << sol_size << " | depth " << depth << endl;

			// forbide Node1->Node2
			int node_sol_size = 0;
			vector<vector<int>> node_skd;
			node_skd.reserve(sol_size);
			vector<double> node_skd_value;
			node_skd_value.reserve(sol_size);
			for (int i = 0; i < sol_size; i++) {
				bool accept = true;
				for (int j = 0; j < schedule[i].size() - 1; j++) {
					if (schedule[i][j] == BestNode1) {
						for (int k = j + 1; k < schedule[i].size(); k++) {
							if (schedule[i][k] == BestNode2) {
								vector<int> new_skd;
								new_skd.reserve(Data->getn());
								new_skd.assign(schedule[i].begin(), schedule[i].end());
								swap(new_skd[j], new_skd[k]);
								node_skd.push_back(new_skd);
								node_skd_value.push_back(Single_Eval(new_skd));
								node_sol_size += 1;
								vector<int>().swap(new_skd);
								accept = false;
								break;
							}
						}
					}
				}
				if (accept) {
					node_skd.push_back(schedule[i]);
					node_skd_value.push_back(skd_value[i]);
					node_sol_size += 1;
				}
			}
			LinkTree* newNode1 = new  LinkTree;
			newNode1->father = branch;
			newNode1->Node1 = BestNode1;
			newNode1->Node2 = BestNode2;
			newNode1->toplevel = false;
			newNode1->is_branch = 1;// branch_dire;
			newNode1->Node_lowerbound = -1e10;
			newNode1->son = NULL;
			Branching(branch_info, false, newNode1);

			bool ok = true;
			ok = runBB(node_skd, best_schedule, node_skd_value, node_sol_size, newNode1, branch_info, depth + 1);
			vector<vector<int>>().swap(node_skd);
			vector<double>().swap(node_skd_value);
			if (!ok) {
				return false;
			}

			branch->son = newNode1;

			LinkTree* newNode2 = new LinkTree;
			newNode2->father = branch;
			newNode2->Node1 = BestNode1;
			newNode2->Node2 = BestNode2;
			newNode2->toplevel = false;
			newNode2->is_branch = 2;
			newNode2->Node_lowerbound = -1e10;
			newNode2->son = NULL;

			for (int i = 0; i < Data->getn(); i++) {
				for (int j = 0; j < Data->getn(); j++) {
					branch_info[i][j] = 0;
				}
			}
			Branching(branch_info, true, newNode2);
			// set Node1->Node2
			int node_sol_size2 = 0;
			vector<vector<int>> node_skd2;
			node_skd2.reserve(sol_size);
			vector<double> node_skd_value2;
			node_skd_value2.reserve(sol_size);

			for (int i = 0; i < sol_size; i++) {
				auto is_Node1 = find(schedule[i].begin(), schedule[i].end(), BestNode1);
				int pos1 = is_Node1 - schedule[i].begin();
				auto is_Node2 = find(schedule[i].begin(), schedule[i].end(), BestNode2);
				int pos2 = is_Node2 - schedule[i].begin();
				if (is_Node1 != schedule[i].end() && is_Node2 != schedule[i].end()) {
					if (pos1<pos2) {
						node_skd2.push_back(schedule[i]);
						node_skd_value2.push_back(skd_value[i]);
						node_sol_size2 += 1;
					}
					else {
						vector<int> new_skd;
						new_skd.reserve(Data->getn());
						new_skd.assign(schedule[i].begin(), schedule[i].end());
						swap(new_skd[pos1], new_skd[pos2]);
						node_skd2.push_back(new_skd);
						node_skd_value2.push_back(Single_Eval(new_skd));
						node_sol_size2 += 1;
						vector<int>().swap(new_skd);
					}
				}
				else if (is_Node1 != schedule[i].end() && is_Node2 == schedule[i].end()) {
					vector<int> new_skd;
					new_skd.reserve(Data->getn());
					new_skd.assign(schedule[i].begin(), schedule[i].end());
					new_skd.insert(new_skd.begin() + pos1 + 1, BestNode2);
					node_skd2.push_back(new_skd);
					node_skd_value2.push_back(Single_Eval(new_skd));
					node_sol_size2 += 1;
					vector<int>().swap(new_skd);
				}
				else if (is_Node2 != schedule[i].end() && is_Node1 == schedule[i].end()) {
					vector<int> new_skd;
					new_skd.reserve(Data->getn());
					new_skd.assign(schedule[i].begin(), schedule[i].end());
					new_skd.insert(new_skd.begin() + pos2, BestNode1);
					node_skd2.push_back(new_skd);
					node_skd_value2.push_back(Single_Eval(new_skd));
					node_sol_size2 += 1;
					vector<int>().swap(new_skd);
				}
				else {
					node_skd2.push_back(schedule[i]);
					node_skd_value2.push_back(skd_value[i]);
					node_sol_size2 += 1;
				}
			}
			ok = runBB(node_skd2, best_schedule, node_skd_value2, node_sol_size2, newNode2, branch_info, depth + 1);
			vector<vector<int>>().swap(node_skd2);
			vector<double>().swap(node_skd_value2);

			branch->Node_lowerbound = min(newNode1->Node_lowerbound, newNode2->Node_lowerbound);
			return ok;
		}
	}
	return false;
}

void BandB::Branching(int** branch_info, bool recur, LinkTree* branch) {
	if (branch->father != NULL) {
		branch_info[branch->Node1][branch->Node2] = branch->is_branch;
		if (recur) {
			Branching(branch_info, true, branch->father);
		}
	}
}

double BandB::Single_Eval(vector<int>& new_skd) {
	double cost = 0;
	double start = 0;
	cost += Data->beta() * Data->q(new_skd[0]) * Data->a(new_skd[0]);
	start += Data->a(new_skd[0]);
	for (int i = 1; i < new_skd.size(); i++) {
		if (start > Data->h(new_skd[i])) {
			if (Data->alpha() * max(start + Data->b(new_skd[i]) + Data->a(new_skd[i]) - Data->d(new_skd[i]), 0.0) + Data->beta() * Data->b(new_skd[i]) * Data->q(new_skd[i]) <= Data->beta() * Data->gete() * (start - Data->h(new_skd[i])) + Data->alpha() * max(start + Data->a(new_skd[i]) - Data->d(new_skd[i]), 0.0)) {
				start += Data->a(new_skd[i]) + Data->b(new_skd[i]);
				cost += Data->alpha() * max(start - Data->d(new_skd[i]), 0.0) + Data->beta() * Data->q(new_skd[i]) * (Data->a(new_skd[i]) + Data->b(new_skd[i]));
			}
			else {
				cost += Data->beta() * (Data->gete() * (start - Data->h(new_skd[i])) + Data->q(new_skd[i]) * Data->a(new_skd[i]));
				start += Data->a(new_skd[i]);
				cost += Data->alpha() * max(start - Data->d(new_skd[i]), 0.0);
			}
		}
		else {
			start += Data->a(new_skd[i]);
			cost += Data->alpha() * max(start - Data->d(new_skd[i]), 0.0) + Data->beta() * Data->a(new_skd[i]) * Data->q(new_skd[i]);
		}
	}
	return cost;
}