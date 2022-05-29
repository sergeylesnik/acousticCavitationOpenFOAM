# Modeling acoustic cavitation with OpenFOAM
A collection of solvers and cases within the scope of OpenFOAM technology
(foam-extend) to model acoustic cavitation.

## Solvers
### coupledHelmholtzMUMPSFoam
solves the wave equation in frequency domain (Helmholtz equation). It utilizes
fvBlockMatrix class from foam-extend and an external multi-frontal (direct)
solver named MUMPS.

### acousticCavitationCloudFoam
The solver includes acoustics via the Helmholtz equation, discrete cavitation
bubbles and URANS modeling of the surrounding liquid. The effect of the
oscillating bubbles on the acoustics is achieved by introducing the attenuation
of the acoustic field due to losses during the bubble oscillations. The latter
are computed using 2D interpolation tables obtained from a bubble radial
dynamics solver (see corresponding repository "cavitationBubbleModeling"). A
similar approach is used to depict the effect of the acoustic waves on the
bubble motion (primary Bjerknes force). The coupling between the bubbles and the
liquid is treated with the standard OpenFOAM routines. The dynamic load
balancing provided within foam-extend is utilized in order to increase
performance which might be low if bubbles begin to cluster at few locations in
the domain.

## Libraries
The are several libraries in src folder which have only few changes compared to
the foam-extend release. This is due to bug fixes that are not included in the
official release yet.

## Installation
### Known to work with
* Ubuntu 16.04
* foam-extend 4.1
* MUMPS 4.10.0

### Steps
* Clone or download this repository to your machine.
* Install foam-extend
  (https://openfoamwiki.net/index.php/Installation/Linux/foam-extend-4.1).
* Install MUMPS
    - from the package repository on Ubuntu:
      ```sudo apt install mumps-test libmumps libmumps-dev```
    - or compile it on your own (http://mumps.enseeiht.fr).
* Compile a solver from solvers/ using wmake or an Allwmake script if provided

## Usage
The test cases are placed in run/ including Allrun and Allclean scripts.

## Feedback
My contact: sergey.lesnik@tu-clausthal.de






















