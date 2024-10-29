---
layout: default
usemathjax: true
---

## Quantum Lattice Modelling with Machine Learning and Quantum Computing

> I’m going to have fun telling you about this absurdity, because I find it delightful. 
> 
> (R. P. Feynman, QED)

Solid state physics aims to understand and design specific material properties, including electron transport and magnetism. This goal requires a profound understanding of the electrons and their dynamics in the crystalline structures. The behaviour of a single electron can be relatively well described by the rules of quantum mechanics. However, a great challenge comes in considering many interacting electrons which behave collectively. The interesting phenomena usually come as a result of the interplay of many different physical properties. The search domain is therefore vast and extremely diverse; the development of new theoretical and computational techniques is vital to deepen our understanding of such systems.

Neural networks are an attractive tool to navigate the large domain of properties, exhibited by systems with many degrees of freedom. It is also more desirable to simulate quantum phenomena using structures which can replicate quantum mechanical behaviour. Recognising that the system of interest is quantum in its nature requires further embedding of the quantum algorithms to make the calculations faster and more effective.

| Code List    |                   | 
|:-------------|:------------------|
| Energy Band Gap Calculation | [_view on GitHub_](https://github.com/KsenijaKovalenka/surface_states)|
| Surface State Calculation | [_view on GitHub_](https://github.com/KsenijaKovalenka/Gap-Calculation-for-Topological-Phase-Transition) |
| Neural Network Architectures | [_view on GitHub_](https://github.com/KsenijaKovalenka/QNN) | 

### Topological Lattices: The Hole Truth

Topology unites spaces into sets that preserve some qualities – called topological invariants. In Euclidian space, this means that two spaces have the same topology if we can continuously deform one to get the other. Stretching and squishing preserve the topology, but not twisting, tearing or poking holes!

![p1-1](https://github.com/KsenijaKovalenka/ksenijakovalenka.github.io/assets/83072550/f4727759-1b8b-4231-87c1-68b7bbe10644)

### Wavefunction Topology

Non-trivial topology of momentum space of a Hamiltonian (in general, any parameter $R(t)$-space) leads to acquiring a geometric phase by the eigenstates when the system undergoes adiabatic transformations in a closed loop. The geometric phase is a quantised topological invariant. In 3D, the geometric phase is calculated using a closed surface integral of the Berry curvature, which is analogous to a magnetic field and reflects the ‘shape’ of the Hamiltonian (or its eigenvectors).

![p2-1](https://github.com/KsenijaKovalenka/ksenijakovalenka.github.io/assets/83072550/ef8ad456-466e-479b-9958-f178eb26b2fd)

### Phase Transition in BiTeI

BiTeI is a crystal undergoing a phase transition from a trivial to a topological insulator when subjected to hydrostatic pressure. Strong spin-orbit coupling in BiTeI leads to translational symmetry breaking, allowing non-trivial topology.


![p3-1](https://github.com/KsenijaKovalenka/ksenijakovalenka.github.io/assets/83072550/069f3e0b-2d96-4ca3-8b00-dc13f4633b0a)

When bands are crossed, the band inversion occurs, meaning that the states are mixing, and the topological invariant of the system changes. We see the jump of electron states corresponding to Bismuth’s outer shell states from conduction to valence band, effectively making a ‘hole’ in the k-space of Tellurium and Iodine electrons. This phenomenon inspired one of my artworks, which can be accessed [_here_](./art.html).

### Surface States

If the vacuum has a trivial topology and a crystal has a non-trivial topology, a transition between the two requires crossing the energy bands at the interface – edge of the crystal. This leads to the formation of conducting edge states in the material. This emerging property is protected from perturbations as small changes do not affect the topology.

![p4-1](https://github.com/KsenijaKovalenka/ksenijakovalenka.github.io/assets/83072550/473e7249-f11c-4171-9c1e-d7f9453e93c2)


### Learning the Transition

A couple of different neural network architectures were successful in predicting the quantum state of BiTeI when given a set of hoppings from the corresponding Hamiltonian. You can enjoy the process of binary classification of topological states performed by the layer of just 2 qubits [_here_](https://youtu.be/c4cicro8FFQ).

[back](./)
