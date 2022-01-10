#pragma once
struct LinkTree{
	LinkTree* father;
	LinkTree* son;
	int Node1;
	int Node2;
	int is_branch;
	double Node_lowerbound;
	bool toplevel;
};