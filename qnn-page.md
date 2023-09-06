---
layout: default
---

## Quantum Neural Networks

> I’m going to have fun telling you about this absurdity, because I find it delightful. 
> 
> (R. P. Feynman, QED)

Noisy Intermediate Scale Quantum (NISQ) hardware and algorithms are showing the promising development of quantum computation. Quantum machine learning in turn is a fast-accelerating field aspiring to create a powerful data-analysis tool for problems which are unreachable even for modern supercomputers. A lot of these potential problems lay in the field of solid-state physics as systems studied in this field are large collections of quantum object. This project aims to investigate the applications of machine learning and quantum computing for modelling the behaviour of different materials.

A classical neural network structure and behaviour was modified to exhibit quantum behaviour by mimicking a quantum system of spins in the magnetic field acting as a qubit circuit. The main manipulations were done on the activation function, which has the structure of the hydrogen radial orbit superposition. With certain amplitudes of different orbitals (determined using the Monte Carlo technique) good learning rate is achieved alongside with dispersion in the intermediate qubit position. Latest implementation includes the CNOT gate acting on 2 qubit system.

[back](./)

<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 3.2 Final//EN">
<html>
<head>
<title>Mathedemo</title>
<script type="text/x-mathjax-config">
  MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}});
</script>
<script type="text/javascript"
  src="http://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script>
</head>

<body>
<h2>Math in TeX notation</h2>

When $a \ne 0$, there are two solutions to \(ax^2 + bx + c = 0\) and they are
$$x = {-b \pm \sqrt{b^2-4ac} \over 2a}.$$

</body>
</html>
