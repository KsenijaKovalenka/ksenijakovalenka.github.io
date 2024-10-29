---
layout: default
usemathjax: true
--- 

## About

I am **Ksenija Kovalenka**, a PhD student at **the University of Manchester**. Throughout my studies, I have developed a strong interest in quantum theory as something completely unfamiliar to my intuition and therefore extremely compelling. My current work focuses on the relation between [_topology and disorder_](./fractals.html) in two-dimensional lattices. I am exploring many-body properties of the multi-fractal states of electronic systems at the critical point of the phase transition. I also have a passion for quantum computation and machine learning, which by the nature of their novelty and intrinsic distinctiveness from other common computation methods are opening unexpected approaches to modern problems. My Master's work on [_'Modelling Quantum Lattices Using a Combination of Deep Learning and Quantum Computing Methods'_](./qnn-page.html) provides me with great experience and skillset allowing me to progress effectively in my fields of interest. This page is a portfolio showcasing my programming projects and other recent work. You can also have a look at the latest [_hackathon_](https://devpost.com/software/space-trash) project I completed with the computer science society at the University of Manchester, or access my CV [_here_](files/CV2.pdf).


## Quantum Lattices

The ability to rationally design and control materials with desirable properties has been a long-standing quest in condensed matter physics. Achieving this goal requires a deep understanding of the intricacies of quantum lattices and high-throughput methods to identify and characterise collective patterns emerging from such systems.

We investigate whether a combination of deep learning and quantum computing approaches can be effective predictors of these properties. This involves training deep learning models on existing datasets and exploring whether they remain good approximations when interpolating and extrapolating from such data. Various models representing quantum lattices with specific orderings will be developed and compared with conventional quantum methods to examine their performance and accuracy.

[_visit project page_](./qnn-page.html)

## C++ Projects

### Quantum Circuit Simulator

The created program allows one to construct circuits with up to 6 quantum bits (qubits) and apply the selection of single-qubit and two-qubit gates. An example circuit is included, which demonstrates all the basic features of the program. The available gates construct the complete basis and hence allow in principle any quantum computation within the qubit limit. Amongst other attributes, the source code is featuring a custom matrix class for gate representation, utilising smart pointers for clear ownership and safe memory management.

[_view on GitHub_](https://github.com/KsenijaKovalenka/ksenijakovalenka.github.io/tree/main/c_projects)

[_Report_](c_projects/report/quantum_circuit_simulator_report.pdf)

## Python Projects

### Neutron Transport

This project aims to investigate particle transport, specifically, penetration of neutrons through the materials of different thickness and composition. Transmission, absorption and reflection rated are calculated for three types of materials: water, lead, and graphite. The procedure includes simulating a path of a neutron through the material, assuming collisions with material particles with two possible outcomes: scattering or absorption, both having constant probabilities throughout the neutron's journey. Monte Carlo technique and randomness explored during this project. Special features include user input, neutron path plots and animation.

[_view on GitHub_](https://github.com/KsenijaKovalenka/ksenijakovalenka.github.io/tree/main/python_projects/programming_courses/neutrons_MC)

### Forced Oscillations

The aim of this project was to investigate how different numerical integration methods apply to a physical systems of forced harmonic oscillator. The equation of motion is $mx’’(t)+bx’(t)+kx(t) = F(t)$, where x is position and it's corresponding derivatives, F is the external force on the system, m is mass and k, b are spring and damping constants respectively. This differential equation can be solved analytically for specific functional forms $F(t)$, but not for all. Therefore, four different integration techniques are tested:
1. Euler 
2. Improved Euler 
3. Verlet 
4. Euler-Cromer 

[_view on GitHub_](https://github.com/KsenijaKovalenka/ksenijakovalenka.github.io/tree/main/python_projects/programming_courses/oscillations)

### Doppler Spectroscopy

Fitting routine on a 2-dimensional landscape to determine the mass of the exoplanet using transit method. The code takes in the data of the doppler shifts obtained from the star during a period of several years, converts it to the velocity of the star using a formula for Doppler's shift and fits a predicted curve to it by minimising the chi squared. A calculation is also performed to find the orbital distance and the mass of the exoplanet, causing the observed variation in the velocity of the star, and it's uncertainty. Plots representing the data and the fit are also created to aid the analysis.

[_view on GitHub_](https://github.com/KsenijaKovalenka/ksenijakovalenka.github.io/tree/main/python_projects/programming_courses/doppler)

### Measuring Drop Spreading Law

The task of the project was to analyse experimental data from the spreading of picolitre droplets on a flat substrate to determine the corresponding spreading law by performing a fit and finding the best parameters for the model. A detailed analysis and comparison of the different models was performed. A spreading law is a relationship between the speed of the contact line and the contact angle.

[_view on GitHub_](https://github.com/KsenijaKovalenka/ksenijakovalenka.github.io/tree/main/python_projects/programming_courses/spread_law)

## Selected Laboratory Projects

### Layer-by-Layer Assembly of van der Waals Heterostructures

Made a combined heterostructure of $WS_{2}$ and $InSe$ encapsulated in $hBN$, using the cantilever-assisted transfer technique. Performed characterisation with Raman spectroscopy and Atomic Force Microscopy (AFM).

![romanlab_page-0001](https://github.com/user-attachments/assets/71171d3f-6238-4d16-a991-88529f040c03)

### Nuclear Magnetic Resonance

Nuclear Magnetic Resonance phenomenon refers to the absorption of the radio waves by nuclei in the magnetic field. This technique was used to determine the gyromagnetic ratios of proton and fluorine, later to high precision. Spin-lattice relaxation times were investigated in aqueous glycerol and ferric nitrate solutions.

[_Report_](python_projects/lab/nmr/NMR_laboratory_report_Ksenija_Kovalenka.pdf)

### Z boson detection in the ATLAS experiment

The Standard Model of particle physics is aiming to describe fundamental building blocks of matter and their interactions with a relatively small number of particles. As with any other physics theory, its predictions have to be tested experimentally. The Large Hadron Collider (LHC) is the largest particle physics experiment built for this purpose. This experiment was aiming to measure the cross-section of the Z boson decay to 2 leptons using data from the [ATLAS experiment](https://arxiv.org/abs/1603.09222) at the LHC collected in 2016.

[_Report_](python_projects/lab/ATLAS/Ksenija_Kovalenka_ATLAS_laboratory_report.pdf)

[_Example Code_](https://github.com/KsenijaKovalenka/ksenijakovalenka.github.io/tree/main/python_projects/lab/ATLAS)

### Galactic Hydrogen

Hydrogen 1s ground state energy undergoes hyperfine splitting due to the interaction of the electron cloud spin and the nuclear spin. Transition in between the two energy levels produces long wavelength radiation penetrating through the dust clouds and the Earth’s atmosphere making it a good astronomical information source. Hydrogen (21cm) spectral scans were taken at a range of longitudes along the galactic plane of the Milky Way Galaxy using the '7m' radio telescope at Jodrell Bank Observatory, Manchester, UK. From these spectra, velocities of hydrogen clouds were obtained to plot the rotation curve, which is important in providing evidence for the existence of the dark matter.

[_Report_](python_projects/lab/galactic_hydrogen/Ksenija_and_Matthew_galactic_hydrogen_report.pdf)

### Cepheid Variables

A distance to the NGC 4258 galaxy was determined by looking at Hubble Space Telescope observations of twelve Cepheid variable stars in the outer parts of the galaxy. Phase folding computational technique was used to construct the light curves and determine the luminosity of stars. Comparing absolute and observed magnitudes allowed to determine the distance to each Cepheid and averaged out gave a value of the distance to the galaxy.

[_Report_](python_projects/lab/cepheid_variables/Ksenija_Kovalenka_Cepheid_Variables.pdf)
