#pragma once
#include "components.h"
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <utility>
#include <map>
#include <algorithm>

class circuit
{
private:
	std::map<int, component*> component_map;
	std::vector<matrix> component_array;
	int number_of_qbits;
	matrix qbit_state_vector;
	std::vector<std::string> draw_array;
public:
	circuit() = default;
	circuit(int n) {number_of_qbits = n;}
	~circuit(){}
	void initialise();
	void add(component *c, int qbit_number);
	void combine();
	void run();
	void print_state() {std::cout<<qbit_state_vector<<std::endl;}
	void draw();
};