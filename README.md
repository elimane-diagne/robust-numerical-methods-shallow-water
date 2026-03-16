# Robust Numerical Methods for the Shallow Water Equations (Fortran)

## Introduction

This project implements **robust numerical methods for the one-dimensional Shallow Water Equations** using the **finite volume method**.

The objective is to develop stable and accurate schemes capable of preserving important physical properties such as:

* conservation of mass
* positivity of water height
* preservation of steady states (lake at rest)

Two numerical schemes are implemented:

* **First-order scheme**
* **Second-order scheme**

The project also includes several **test cases** used to validate the numerical methods.

---

## Mathematical Model

The **Shallow Water Equations** are given by

Continuity equation

∂t h + ∂x (hu) = 0

Momentum equation

∂t (hu) + ∂x (hu
