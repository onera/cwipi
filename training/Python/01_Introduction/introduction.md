---
jupytext:
  text_representation:
    extension: '.md'
    format_name: myst
    format_version: '0.7'
    jupytext_version: 1.4.0+dev
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Inntroduction

Welcome to the CWIPI library introductory day !

The aim of the day is to give you : 
  - An overview of CWIPI 1.0 functionalities
  - An understanding on how CWIPI can be called to couple solvers written in C/C++, Fortran or Python.

This is an interactive course, so don't hesitate to interrupt us to ask your questions.

The training will take place in three stages:
  - General presentation of CWIPI
  - Exercise 1 : a basic coupling
  - Exercise 2 : an extension on time step manadgment with a deformed mesh

## What is CWIPI?

#### LGPL code coupling library in a massively parallel distributed setting:
- Provides tools to define code coupling algorithms to control each exchange between the codes
- Transfer of the interpoleted field through a geometric coupling interface

#### Version 0.x:
- since 2009
- dependency to FVM (LGPL, EDF R&D)

#### Version 1.0:
- released in 2023
- dependency to ParaDiGM (LPGL, ONERA)

#### Contributors:
- ONERA: equivalent of one full time developper (4 physical persons)
- CERFACS: post-doctoral researchers

#### Acess
- https://w3.onera.fr/cwipi/bibliotheque-couplage-cwipi (in french)
- on GitHub (in progress)

#### Documentation
- https://numerics.gitlab-pages.onera.net/coupling/cwipi/dev/index.html

# Core concepts

Before starting the hands-on exercise, we will spend some time on the philosophy to build a coupling using the CWIPI library. To do so, we will work with the following basic coupling exemple:

## A basic coupling

Let `code 1` and `code 2` be two solvers. We want to send a field of `code 1` defined on the nodes of the associated mesh to `code 2`. Let the coupling interface for each code be the following mesh:

![alt text](mesh.png)

In this particular case the coupling interface is the same mesh for both solvers but it is different most of the time in real case applications.

### Coupling definition

As seen in the text above, to define a coupling one has first to tell which codes are going to communicate together. In our case of this basic coupling, that gives us:

![alt text](schema_concept_coupling.svg)


### Coupling interface

If we want to couple those two solvers, it means that there is an area in which the physical phenomenon of `code 2` acts on `code 1` and vise versa. This area in a numerical simulation is defined as a portion of the mesh on which the solvers applies its numerical scheme. It is this portion of the mesh that we call coupling interface. Thus for each code participating in the coupling, we need to define this mesh portion:

![alt text](schema_concept_mesh.svg)


### Field definition

The point of the whole coupling is to let the physical phenomenon in the common area interact. This means that the field describing the phenomenon in `code 1` has to be shared with `code 2` (and vise versa but in this basic exemple we leave that appart). Each solver has its own numerical method thus location (cell centered, node centered or other) of the field on its mesh. To set what will be exchanged we define for each code its fields. Here `code 1` will send a field to `code 2` and `code 2` will create a buffer to receive what `code 1` sends it (again and vise versa). This leaves us with:

![alt text](schema_concept_field.svg)

# Pseudo-code coupling algorithm

The basic coupling described above, reads in pseudo-code as:

```{prf:algorithm} basic couling algorithm

**Inputs** Given $code1$ with a mesh $m1$ on which a field that will be sent is defined $sf1$. Given $code2$ with a mesh $m2$ on which a field that will be received is defined $rf2$.

**Output** $rf2$, which is $sf1$ interpolated on $m2$

1. Initialize CWIPI
2. Set coupling between $code1$ and $code2$
3. Describe codes:
   1. $code1$ has a mesh $m1$ on which we define a field $sf1$
   2. $code2$ has a mesh $m2$ and a receive buffer field $rf2$
4. Operate solver iterations:
   1. $code1$ sends $sf1$
   2. $code2$ receives $rf2$
5. Finalize CWIPI
```

In exercise 1, we will focus on implementing this basic coupling.
In exercise 2, we will adapt that code for a coupling interface mesh that is deformed over time.

# Exercise 0

But first, let's understand how to work with a Jupyter Notebook.

cf. Notebook Exercise 0
