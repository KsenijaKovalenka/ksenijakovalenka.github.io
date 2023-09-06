#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>

#include "components.h"
#include "matrix.h"
#include "circuit.h"


bool add_more(std::string object)
{
	// ask for any extra components/layers
	std::string continue_response;
	std::cout<<"Would you like to add another "<<object<<"? [y/n]"<<std::endl;
	std::cin>>continue_response;
	while(1) {
		if (continue_response=="y") {return true;}
		else if (continue_response=="n") {return false;}
		else {
			std::cout<<"Please enter y or n"<<std::endl;
			std::cin>>continue_response;
			}
	}
}

template<typename T>
T get_value(T min, T max) 
{
    T value;
    while (true) {
        std::cout << "Enter a value between " << min << " and " << max << ": ";
        if (!(std::cin >> value)) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "Invalid input. Please try again." << std::endl;
            continue;
        }
        if (value >= min && value <= max) {
            break; // Valid input, exit the loop
        } else {
            std::cout << "Value is out of range. Please try again." << std::endl;
        }
    }
    return value;
}

int main()
{
	// options to have only an example or only user input
	bool example{true};
	bool user{true};

	// instantiate component objects
	RotationX Rx(M_PI);
	RotationY Ry(M_PI);
	RotationZ Rz(M_PI);
	PhaseShift P(M_PI);
	ControlledNot CNOT;
	PauliX X;
	PauliY Y;
	PauliZ Z;
	Identity I;
	Hadamard H;
	SWAP S;

	if (example) {
		std::cout<<"..................."<<std::endl;
		std::cout<<"..EXAMPLE CIRCUIT.."<<std::endl;
		std::cout<<"..................."<<std::endl;

		circuit my_circuit_example(3);
		my_circuit_example.initialise();
		std::cout<<"initial state"<<std::endl;
		my_circuit_example.print_state();
		my_circuit_example.add(&X, 1);
		my_circuit_example.combine();
		std::cout<<"here"<<std::endl;

		my_circuit_example.add(&CNOT, 1);
		my_circuit_example.add(&X, 3);
		my_circuit_example.combine();
		std::cout<<"here"<<std::endl;

		my_circuit_example.add(&Ry, 2);
		my_circuit_example.combine();

		my_circuit_example.run();
		std::cout<<"final state"<<std::endl;
		my_circuit_example.print_state();
		my_circuit_example.draw();
		std::cout<<".........................."<<std::endl;
		std::cout<<"..END OF EXAMPLE CIRCUIT.."<<std::endl;
		std::cout<<".........................."<<std::endl;
	}

	if (user) {

		// flags to complete the layer and the whole circuit
		bool layer_progress{true};
		bool circuit_progress{true};

		// create a component library
		std::map<std::string, component*> component_library;
    	component_library.insert(std::make_pair("rx", &Rx));
		component_library.insert(std::make_pair("ry", &Ry));
		component_library.insert(std::make_pair("rz", &Rz));
		component_library.insert(std::make_pair("p", &P));
		component_library.insert(std::make_pair("cnot", &CNOT));
		component_library.insert(std::make_pair("x", &X));
		component_library.insert(std::make_pair("y", &Y));
		component_library.insert(std::make_pair("z", &Z));
		component_library.insert(std::make_pair("i", &I));
		component_library.insert(std::make_pair("h", &H));
		component_library.insert(std::make_pair("s", &S));

		// ask to initialise the circuit
		std::cout<<std::endl;
		std::cout<<"..................."<<std::endl;
		std::cout<<"..  USER  INPUT  .."<<std::endl;
		std::cout<<"..................."<<std::endl;
		std::cout<<"let's construct the circuit!"<<std::endl;
		std::cout<<"Please enter the number of qubits in the circuit:"<<std::endl;
		std::cout<<"(1-10)"<<std::endl;

		int qbit_number;
		qbit_number = get_value<int>(1, 6);
		circuit my_circuit(qbit_number);
		my_circuit.initialise();
		std::cout<<"Initialsised a circuit with "<<qbit_number<<" qubits."<<std::endl;
		std::cout<<std::endl;

		std::cout<<"Here're the available gates:"<<std::endl;
		std::cout<<"---------------------------------------------------"<<std::endl;
		std::cout<<"        Gate name        |           key           "<<std::endl;
		std::cout<<"---------------------------------------------------"<<std::endl;
		std::cout<<"       X Rotation        |            rx           "<<std::endl;
		std::cout<<"       Y Rotation        |            ry           "<<std::endl;
		std::cout<<"       Z Rotation        |            rz           "<<std::endl;
		std::cout<<"       Phase Shift       |            p            "<<std::endl;
		std::cout<<"       CNOT              |            cnot         "<<std::endl;
		std::cout<<"       Pauli X           |            x            "<<std::endl;
		std::cout<<"       Pauli Y           |            y            "<<std::endl;
		std::cout<<"       Pauli Z           |            z            "<<std::endl;
		std::cout<<"       Identity          |            i            "<<std::endl;
		std::cout<<"       Hadamard          |            h            "<<std::endl;
		std::cout<<"       SWAP              |            s            "<<std::endl;
		std::cout<<"---------------------------------------------------"<<std::endl;
		
	
		// outside loop for completing the circuit
		while(circuit_progress){
			layer_progress = circuit_progress;
			std::vector<int> occupied_qubits;
			// inside loop for comleting a layer
			std::string key;
			while(layer_progress){
				// check if there are available qubits
				if (occupied_qubits.size() >= qbit_number) {
					std::cout<<"no qubits available in this layer"<<std::endl;
					break;
				}
				std::cout<<"Please pick the gate by entering the key:"<<std::endl;
				std::cin>>key;
				while(1){
					// check if component with the entered key exists
					if (!component_library.count(key)) {
						std::cin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
						std::cout<<"wrong input, try again: "<<std::endl;
						std::cout<<"your value should be one of the keys"<<std::endl;
						std::cin>>key;
					} else if ((key == "cnot" || key == "s") && qbit_number <= 1) {  // check if there are two gates in the circuit for 2-qubit gates
						std::cin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
						std::cout<<"two qubit gates are unavailable because there is only one qubit in your circuit"<<std::endl;
						std::cin>>key;
					} else break;
				}

				// some gates require a parameter
				std::string custom;
				double gate_parameter = 0.0;
				if (key == "rx" || key == "ry" || key == "rz" || key == "p") {
					std::cout<<"The gate you've picked requires a parameter"<<std::endl;
					std::cout<<"Please type c for custom or p for default (pi) parameter"<<std::endl;
					// input parameter
					while(1){
						std::cin>>custom;
						if (custom == "c") {gate_parameter = get_value<double>(0, 2*M_PI); break;}
						else if (custom == "p") {break;}
						else {
							std::cin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
							std::cout<<"Wrong input, please type c or p"<<std::endl;
						}
					}
				}

				// dynamically create the component based on the selected key and parameter
				if (custom == "c") {
					component* new_component = nullptr;
					if (key == "rx") {
						component_library.insert(std::make_pair("rxc", new RotationX(gate_parameter)));
						key = "rxc";
					} else if (key == "ry") {
						component_library.insert(std::make_pair("ryc", new RotationX(gate_parameter)));
						key = "ryc";
					} else if (key == "rz") {
						component_library.insert(std::make_pair("rzc", new RotationX(gate_parameter)));
						key = "rzc";
					} else if (key == "p") {
						component_library.insert(std::make_pair("pc", new RotationX(gate_parameter)));
						key = "pc";
					}
				}
				// reset the custom key parameter
				custom = "p";
				bool break_layer = false;

				// gate picked, now the target qbit
				int which_qubit;
				while (1) {
					std::cout << "Which qubit do you want to apply the gate to?" << std::endl;
					// applying 2 qubit gates on the last qubit doesn't make sence
					if (key == "cnot" || key == "s") {
						which_qubit = get_value<int>(1, qbit_number-1);
					}
					else {
						which_qubit = get_value<int>(1, qbit_number);
					}

					// check if the qubit is available
					bool unavailable = std::any_of(occupied_qubits.begin(), occupied_qubits.end(), [which_qubit](int value) {
						return value == which_qubit;
					});

					// check if the 2nd qubit is available for 2 qubit gates
					if (key == "cnot" || key == "s") {
						unavailable = std::any_of(occupied_qubits.begin(), occupied_qubits.end(), [which_qubit](int value) {
							return value == (which_qubit + 1);
						});
						// check if there are still available qubits
						if (occupied_qubits.size() + 1 >= qbit_number) {
							std::cout<<"no TWO qubits available in this layer"<<std::endl;
							break_layer = true;
							break;
						}
					}

					if (unavailable) {
						std::cout << "The qubit you have picked already has a gate. Please pick another one." << std::endl;
						std::cin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
					} else {
						// mark qubit as unavailable
						if (key == "cnot" || key == "s") {
							occupied_qubits.push_back(which_qubit);
							occupied_qubits.push_back(which_qubit + 1);
						}
						else {occupied_qubits.push_back(which_qubit);}
						break;
					}
				}
				if (break_layer) {break;}

				// add the component to the layer and repeat
				my_circuit.add(component_library[key], which_qubit);

				// ask for any extra components in the layer
				layer_progress = add_more("component");
			}

			std::cout<<"layer completed"<<std::endl;
			// combine the layers to one matrix
			my_circuit.combine();

			if (key == "rx" || key == "ry" || key == "rz" || key == "p") {
				// Free dynamically allocated memory by deleting the components
				for (const auto& key_value : component_library) {
					if (key_value.first == "rxc" || key_value.first == "ryc" || key_value.first == "rzc" || key_value.first == "pc") {
						delete key_value.second;
					}
				}
			}
			// ask for any extra layers in the circuit
			circuit_progress = add_more("layer");
		}

		// run the created circuit and output the result along with circuit architecture
		my_circuit.run();
		std::cout<<"final state:"<<std::endl;
		my_circuit.print_state();

		std::cout<<"circuit diagram:"<<std::endl;
		my_circuit.draw();

		std::cout<<".........................."<<std::endl;
		std::cout<<"..  END OF USER INPUT   .."<<std::endl;
		std::cout<<".........................."<<std::endl;
	}
	
	return 0;
}