#include "Reader.h"
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

Reader::Reader(string add) {
    ifstream fin(add);
    double temp = 0;
    fin >> m;
    fin >> n;
    fin >> e0;
    aTime = new double[n];
    hTime = new double[n];
    bTime = new double[n];
    dTime = new double[n];
    qTime = new double[n];
    for (int i = 0; i < n; i++) {
        fin >> aTime[i];
    }
    for (int i = 0; i < n; i++) {
        fin >> hTime[i];
    }
    for (int i = 0; i < n; i++) {
        fin >> bTime[i];
    }
    for (int i = 0; i < n; i++) {
        fin >> qTime[i];
    }
    for (int i = 0; i < n; i++) {
        fin >> dTime[i];
    }
    fin >> Alpha;
    fin >> Beta;
    fin.close();
}
Reader::~Reader() {
    delete[] aTime; aTime = nullptr;
    delete[] hTime; hTime = nullptr;
    delete[] bTime; bTime = nullptr;
    delete[] dTime; dTime = nullptr;
    delete[] qTime; qTime = nullptr;
}

int Reader::getm() {
    return m;
}

int Reader::getn() {
    return n;
}

double Reader::gete() {
    return e0;
}

double Reader::a(int num) {
    return aTime[num];
}

double Reader::h(int num) {
    return hTime[num];
}

double Reader::b(int num) {
    return bTime[num];
}

double Reader::d(int num) {
    return dTime[num];
}

double Reader::q(int num) {
    return qTime[num];
}

double Reader::alpha() {
    return Alpha;
}

double Reader::beta() {
    return Beta;
}