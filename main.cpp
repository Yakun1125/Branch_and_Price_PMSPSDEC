
#include <ilcplex/ilocplex.h>
#include "Evaluation.h"
#include "Reader.h"
#include "VNS.h"
#include "ColGen.h"
#include "LinkTree.h"
#include "BranchBound.h"
#include <vector>


ILOSTLBEGIN


int
main(void) {

    Reader* Data = new Reader("C:\\Users\\WYKge\\Desktop\\BranchandPrice_PMSPSDEC\\2_8.txt");
    Evaluate* sol = new Evaluate(Data);
    VNS* Search = new VNS(10, 10, 2, 200, Data, sol);
    ColumnGen* Col = new ColumnGen(Data);
    int sol_size = Data->getm();
    vector<vector<int>> schedule;
    vector<vector<int>> best_skd(Data->getm(), vector<int>(Data->getn(),-1));
    schedule.resize(sol_size);
    for (int m = 0; m < sol_size; m++) {
        schedule[m].reserve(Data->getn());
    }
    vector<double> skd_value;
    skd_value.resize(sol_size);

    int** branch_info = new int* [Data->getn()];
    for (int n = 0; n < Data->getn(); n++) {
        branch_info[n] = new int[Data->getn()]();
    }

    Col->Initialize(sol, Search, schedule, skd_value);
    double BestBound = 0;
    for (int m = 0; m < Data->getm(); m++) {
        BestBound += skd_value[m];
        best_skd[m].clear();
        best_skd[m].assign(schedule[m].begin(), schedule[m].end());
    }

    for (int m = 0; m < Data->getm(); m++) {
        for (int n = 0; n < best_skd[m].size(); n++) {
            cout << best_skd[m][n] << ", ";
        }
        cout << endl;
    }
    cout << "Upperbound " << BestBound << endl;

    //vector<double> storeY(sol_size, -1);
    //Col->RunColGen(schedule, skd_value, storeY, sol_size, branch_info);

    int depth = 0;
    BandB* BranchandPrice = new BandB(Data, Col, -1e10, BestBound);
    BranchandPrice->runBB(schedule, best_skd, skd_value, sol_size, NULL, branch_info, depth);

    for (int m = 0; m < Data->getm(); m++) {
        for (int n = 0; n < best_skd[m].size(); n++) {
            cout << best_skd[m][n] << ", ";
        }
        cout << endl;
    }

    delete Col; Col = nullptr;
    delete Search; Search = nullptr;
    delete sol; sol = nullptr;
    for (int n = 0; n < Data->getn(); n++) {
        delete[] branch_info[n]; branch_info[n] = nullptr;
    }
    delete[] branch_info;
    delete BranchandPrice; BranchandPrice = nullptr;
    delete Data; Data = nullptr;
    return 0;

}