# acousticCavitationOpenFOAM
A collection of solvers and cases within the scope of OpenFOAM technology (foam-extend) to model acoustic cavitation.
The only solver available right now is coupledHelmholtzMUMPSFoam which solves the wave equation in frequency domain (Helmholtz equation). It utilizies fvBlockMatrix class from foam-extend and an external multi-frontal (direct) solver named MUMPS.

## Installation
### Known to work with
* Ubuntu 16.04
* foam-extend 4.1
* MUMPS 4.10.0

### Steps
* Clone or download this repository to your machine.
* Install foam-extend (https://openfoamwiki.net/index.php/Installation/Linux/foam-extend-4.1).
* Install MUMPS 
    - from the package repository: ```sudo apt install mumps-test libmumps libmumps-dev```
    - or compile it on your own (http://mumps.enseeiht.fr)
* Compile a solver from solvers/ using wmake.

## Usage
The test cases are placed in run/ including Allrun and Allclean scripts.

## Feedback
My contact: sergey.lesnik@tu-clausthal.de






















