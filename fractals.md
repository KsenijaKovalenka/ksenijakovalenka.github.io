---
layout: default
usemathjax: true
---

## Disorder-Driven Topological Insulator Phases in Two-Dimensional Materials

> Beautiful, damn hard, increasingly useful. That’s fractals. 
> 
> (B. Mandelbrot)

| Code List    |                   | 
|:-------------|:------------------|
| Large Matrix Diagonalisation | [_view on GitHub_](https://github.com/KsenijaKovalenka/ksenijakovalenka.github.io/tree/main/ksenija)|
| LCM Calculation | [_view on GitHub_](https://github.com/KsenijaKovalenka/ksenijakovalenka.github.io/tree/main/chern%20marker) |
| Multifractality Analysis | [_view on GitHub_](https://github.com/KsenijaKovalenka/ksenijakovalenka.github.io/tree/main/multifractals) | 

Topology has shown to be an important tool in classifying solid-state systems with peculiar behaviour. 
It is a global property of the system, hence exploring the effects of local perturbations on topological materials is of great interest. 
A way of calculating a topological invariant locally is presented here, shedding light on the relationship between topology and disorder.

### Haldane Model

Haldane Model is a graphene-like tight-binding model with the energy shift at the different sites of the unit cell, which is referred to as the massive 
Dirac Hamiltonian. In addition, there are complex next-nearest neighbour hoppings, resulting in electrons’ topological phase. Depending on the 
Dirac mass term M, the system can be a [trivial or topological insulator](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.61.2015). We added random onsite energies to the model, representing the 
disorder. 

![fractals1_page-0001](https://github.com/user-attachments/assets/63913806-3556-42f2-ba64-8d4659cc889d)

For a range of mass parameters, increasing disorder strength can lead to the series of phase transitions from trivial insulator to 
topological insulator to Anderson localisation (insulating) state. This can be seen from the LCM calculation for a finite slab of the material below. [ (package used) ](https://github.com/roberta-favata/spinv?tab=%20readme-ov-file,%20accessed:%202023-04-10)

![fractals3_page-0001](https://github.com/user-attachments/assets/1e000159-05c2-4b0a-a4d7-c6db8dec72f3)

A sum of LCM over the edges for a range of masses and disorder strength was computed to 
produce a phase diagram below.

![fractals2_page-0001](https://github.com/user-attachments/assets/9fbfe3c0-0758-45a4-aacc-7f9153677b43)

### Multifractality

There is a critical point in the Andeson localisation phase transition as the disorder strength is increased in the system. At this point, the wavefunction of the electrons at the Fermi level exhibits [multifractal scaling properties](https://en.wikipedia.org/wiki/Multifractal_system), which can be shown by calculating its dimensionality spectrum as shown below. This feature of the eigenstates might be the key to the emergence of topology in the presence of weak disorder due to the enhancement of correlation.

![fractals4_page-0001](https://github.com/user-attachments/assets/4c05193f-2e77-4521-af06-a0e461d97947)

[back](./)
