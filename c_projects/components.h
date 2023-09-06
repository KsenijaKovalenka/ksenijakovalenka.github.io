#pragma once
#include "matrix.h"

using namespace std::complex_literals;

class component
{
protected:
	matrix gate{matrix(2,2)};
public:
	component() = default;
	~component(){}
	matrix get_matrix(){return gate;}
	void display();
	virtual std::string draw(){return "--|  c  |--";}
};

class qbit
{
private:
	double angle{0};
	double phase{0};
	matrix qbit_vector{matrix(2,1)};
public:
	qbit(){}
	qbit(double theta, double phi);
	qbit(std::complex<double> z1, std::complex<double> z2);
	~qbit(){}
	std::complex<double> get_zero_amplitude();
	std::complex<double> get_one_amplitude();
	void apply(component gate);
	void show();
	std::complex<double> dot(qbit q);
	matrix get_vector() {return qbit_vector;}
};

// Universal set of gates
class RotationX : public component
{
private:
	double angle;
	std::string diagram;
public:
	RotationX(double theta);
	~RotationX(){}
	std::string draw(){return diagram;}
};

class RotationY : public component
{
private:
	double angle;
	std::string diagram;
public:
	RotationY(double theta);
	~RotationY(){}
	std::string draw(){return diagram;}
};

class RotationZ : public component
{
private:
	double angle;
	std::string diagram;
public:
	RotationZ(double theta);
	~RotationZ(){}
	std::string draw(){return diagram;}
};

class PhaseShift : public component
{
private:
	double phase;
	std::string diagram;
public:
	PhaseShift(double phi);
	~PhaseShift(){}
	std::string draw(){return diagram;}
};

// 2 qbit gate for entanglement
class ControlledNot : public component
{
private:
	std::string diagram;
public:
	ControlledNot();
	~ControlledNot(){}
	std::string draw(){return diagram;}
};

// Usefull shortcut gates 
class PauliX : public component
{
private:
	std::string diagram;
public:
	PauliX();
	PauliX(double angle);
	~PauliX(){}
	std::string draw(){return diagram;}
};

class PauliY : public component
{
private:
	std::string diagram;
public:
	PauliY();
	~PauliY(){}
	std::string draw(){return diagram;}
};

class PauliZ : public component
{
private:
	std::string diagram;
public:
	PauliZ();
	~PauliZ(){}
	std::string draw(){return diagram;}
};

class Identity : public component
{
private:
	std::string diagram;
public:
	Identity();
	~Identity(){}
	std::string draw(){return diagram;}
};

class Hadamard : public component
{
private:
	std::string diagram;
public:
	Hadamard();
	~Hadamard(){}
	std::string draw(){return diagram;}
};

class SWAP : public component
{
private:
	std::string diagram;
public:
	SWAP();
	~SWAP(){}
	std::string draw(){return diagram;}
};