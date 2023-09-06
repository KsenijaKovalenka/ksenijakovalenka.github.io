#include "circuit.h"
#include "components.h"
#include "matrix.h"
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>

void circuit::initialise()
{
	qbit_state_vector = matrix(pow(2, number_of_qbits),1);
	qbit_state_vector(1,1) = 1;
	draw_array.resize(number_of_qbits, "|0>  ");
}

void circuit::add(component *c, int qbit_number)
{
	component_map.insert(std::make_pair(qbit_number, c));
	std::cout<<"added components"<<std::endl;
}

void circuit::combine()
{
	// define identity for empty gates
	Identity I;
	size_t start;
	matrix combined_component;
	// check the first qubit and add a corresponding component
	if (component_map.count(1)) {
		combined_component = component_map[1]->get_matrix();
		// add component to the diagram
		draw_array[0] += component_map[1]->draw();
	} else {
		// assume identity if nothing specified
		combined_component = I.get_matrix();
		// add component to the diagram
		draw_array[0] += I.draw();
	}
	// check if the first component is applied to two gates
	if (combined_component.get_cols()==2) {start = 2;}
	else {start = 3; draw_array[1] += "----------------";}
	// perform a reccursive tensor product on the rest of the components
	for (size_t i{start}; i<=number_of_qbits; i++) {
		if (component_map.count(i)) {
			matrix single_component{component_map[i]->get_matrix()};
			combined_component = combined_component.tensor_product(single_component);
			// add component to the diagram
			draw_array[i - 1] += component_map[i]->draw();
			// check if the component is applied to two gates
			if (single_component.get_cols()==4) {i++; draw_array[i - 1] += "----------------";}
		} else {
			// assume identity if nothing specified
			combined_component = combined_component.tensor_product(I.get_matrix());
			// add component to the diagram
			draw_array[i - 1] += I.draw();
		}
	}
	// add the newly specified matrix to the circuit array
	component_array.push_back(combined_component);
	// clean the component map
	component_map.clear();
}

void circuit::run()
{
	// apply gates in the component array one by one
	for (auto element : component_array) {
		matrix new_vector{element * qbit_state_vector};
		qbit_state_vector = new_vector;
	}
}

void circuit::draw()
{
	for (const std::string& comp : draw_array) {
        std::cout << comp << '\n';
    }
}