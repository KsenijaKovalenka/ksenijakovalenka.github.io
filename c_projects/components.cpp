#include "components.h"
#include "matrix.h"
#include <iomanip> 

void component::display(){std::cout<<gate<<std::endl;}

qbit::qbit(double theta, double phi) : angle{theta},phase{phi}
{
	qbit_vector(1,1) = cos(angle/2);
	qbit_vector(2,1) = sin(angle/2)*exp(phase*1i);
}
qbit::qbit(std::complex<double> z1, std::complex<double> z2)
{
	qbit_vector(1,1) = z1; 
	qbit_vector(2,1) = z2;
}

std::complex<double> qbit::get_zero_amplitude() {return qbit_vector(1,1);}
std::complex<double> qbit::get_one_amplitude() {return qbit_vector(2,1);}

void qbit::apply(component gate)
{
	matrix new_vector{gate.get_matrix()*qbit_vector};
	qbit_vector=new_vector;
}

void qbit::show()
{
	std::cout<<qbit_vector<<std::endl;
}

std::complex<double> qbit::dot(qbit q)
{
	return std::conj(qbit_vector(1,1)) * q.get_zero_amplitude() + std::conj(qbit_vector(2,1)) * q.get_one_amplitude();
}

// gate initialisatinons
RotationX::RotationX(double theta)
{
	angle = theta;
	std::ostringstream oss;
	oss << std::fixed << std::setprecision(0) << theta; // Set precision to 0 to remove decimals
    std::string trimmedtheta = oss.str();
	diagram = "--|  R_x(" + trimmedtheta + ")  |--" ;
	gate(1,1) = cos(angle/2);       gate(1,2) = -sin(angle/2) * 1i;
	gate(2,1) = -sin(angle/2) * 1i; gate(2,2) = cos(angle/2);
}

RotationY::RotationY(double theta)
{
	angle = theta;
	std::ostringstream oss;
    oss << std::fixed << std::setprecision(0) << theta; // Set precision to 0 to remove decimals
    std::string trimmedtheta = oss.str();
	diagram = "--|  R_y(" + trimmedtheta + ")  |--" ;
	gate(1,1) = cos(angle/2); gate(1,2) = -sin(angle/2);
	gate(2,1) = sin(angle/2); gate(2,2) = cos(angle/2);
}

RotationZ::RotationZ(double theta)
{
	angle = theta;
	std::ostringstream oss;
	oss << std::fixed << std::setprecision(0) << theta; // Set precision to 0 to remove decimals
    std::string trimmedtheta = oss.str();
	diagram = "--|  R_z(" + trimmedtheta + ")  |--" ;
	gate(1,1) = exp(-angle*1i); gate(1,2) = 0;
	gate(2,1) = 0             ; gate(2,2) = exp(angle*1i);
}

PhaseShift::PhaseShift(double phi)
{
	phase = phi;
	std::ostringstream oss;
	oss << std::fixed << std::setprecision(1) << phi; // Set precision to 0 to remove decimals
    std::string trimmedtheta = oss.str();
	diagram = "--|  P(" + trimmedtheta + ")  |--" ;
	gate(1,1) = 1; gate(1,2) = 0;
	gate(2,1) = 0; gate(2,2) = exp(phase*1i);
}

ControlledNot::ControlledNot()
{
	gate = matrix(4,4);
	gate(1,1) = 1; gate(1,2) = 0; gate(1,3) = 0; gate(1,4) = 0;
	gate(2,1) = 0; gate(2,2) = 1; gate(2,3) = 0; gate(2,4) = 0;
	gate(3,1) = 0; gate(3,2) = 0; gate(3,3) = 0; gate(3,4) = 1;
	gate(4,1) = 0; gate(4,2) = 0; gate(4,3) = 1; gate(4,4) = 0;
	diagram = "--|   CNOT   |--" ;
}

PauliX::PauliX()
{
	gate(1,1) = 0; gate(1,2) = 1;
	gate(2,1) = 1; gate(2,2) = 0;
	diagram = "--|    X     |--" ;
}

PauliY::PauliY()
{
	gate(1,1) = 0; gate(1,2) = -1i;
	gate(2,1) = 1i; gate(2,2) = 0;
	diagram = "--|    Y     |--" ;
}

PauliZ::PauliZ()
{
	gate(1,1) = 1; gate(1,2) = 0;
	gate(2,1) = 0; gate(2,2) = -1;
	diagram = "--|    Z     |--" ;
}

Identity::Identity()
{
	gate(1,1) = 1; gate(1,2) = 0;
	gate(2,1) = 0; gate(2,2) = 1;
	diagram = "----------------" ;
}

Hadamard::Hadamard()
{
	gate(1,1) = 1/sqrt(2); gate(1,2) = 1/sqrt(2);
	gate(2,1) = 1/sqrt(2); gate(2,2) = -1/sqrt(2);
	diagram = "--|    H     |--" ;
}

SWAP::SWAP()
{
	gate = matrix(4,4);
	gate(1,1) = 1; gate(1,2) = 0; gate(1,3) = 0; gate(1,4) = 0;
	gate(2,1) = 0; gate(2,2) = 0; gate(2,3) = 1; gate(2,4) = 0;
	gate(3,1) = 0; gate(3,2) = 1; gate(3,3) = 0; gate(3,4) = 0;
	gate(4,1) = 0; gate(4,2) = 0; gate(4,3) = 0; gate(4,4) = 1;
	diagram = "--|   SWAP   |--" ;
}