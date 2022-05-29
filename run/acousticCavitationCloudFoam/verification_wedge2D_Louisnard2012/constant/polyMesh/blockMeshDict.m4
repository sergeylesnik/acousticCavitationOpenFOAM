/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | foam-extend: Open Source CFD                    |
|  \\    /   O peration     | Version:     4.1                                |
|   \\  /    A nd           | Web:         http://www.foam-extend.org         |
|    \\/     M anipulation  | For copyright notice see file Copyright         |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ------!!!!!!   run "m4 blockMeshDict.m4 > blockMeshDict"  !!!!!!---------//
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// General m4 macros
changecom(//)changequote([,]) dnl>
define(calc, [esyscmd(perl -e 'use Math::Trig; print ($1)')]) dnl>
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])
define(hex2D, hex (b$1 b$2 b$3 b$4 f$1 f$2 f$3 f$4))
define(quad2D, (b$1 b$2 f$2 f$1))
define(frontQuad, (f$1 f$2 f$3 f$4))
define(backQuad, (b$1 b$4 b$3 b$2))
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// User-defined parameters

convertToMeters 1;

define(R, 0.3)
define(H, 0.4)
define(r, 0.06)
define(k, 0.03)
define(h, calc(H-k))
define(l, calc(R-r))
define(a, 5) // a straddling angle in °

define(refine, 8)
define(nh, calc(refine*50))
define(nr, calc(refine*8))
define(nl, calc(refine*30))
define(nk, calc(refine*10))

define(pi, 3.14159265358979323844)
define(rad, [calc($1*pi/180.0)])
define(arad, rad(a)) // arad angle in rad with the degToRad-function

define(xr, calc(r*cos(arad/2.0)))
define(zr, calc(r*sin(arad/2.0)))
define(xR, calc(R*cos(arad/2.0)))
define(zR, calc(R*sin(arad/2.0)))
define(nzr, calc(-zr))
define(nzR, calc(-zR))

vertices
(
    (0   0  0) 	         // Vertex Nr. 0
    (xr  0  nzr)	 // Vertex Nr. 1
    (xr  h  nzr) 	 // Vertex Nr. 2
    (0   h  0)  	 // Vertex Nr. 3
    (xR  0  nzR) 	 // Vertex Nr. 4
    (xR  h  nzR) 	 // Vertex Nr. 5
    (xr  H  nzr) 	 // Vertex Nr. 6
    (xR  H  nzR) 	 // Vertex Nr. 7
    (xr  0  zr) 	 // Vertex Nr. 8
    (xr  h  zr) 	 // Vertex Nr. 9
    (xR  0  zR) 	 // Vertex Nr. 10
    (xR  h  zR) 	 // Vertex Nr. 11
    (xr  H  zr) 	 // Vertex Nr. 12
    (xR  H  zR) 	 // Vertex Nr. 13
);

edges
(
);

blocks
(
    hex (0 1 2 3 0  8  9 3)  (nr nh 1) simpleGrading (0.2 0.5  1)
    hex (1 4 5 2 8  10 11 9) (nl nh 1) simpleGrading (4   0.5  1)
    hex (2 5 7 6 9 11 13 12) (nl nk 1) simpleGrading (4   2.5  1)
);

boundary
(
    Transducer
    {
        type wall;
        faces
        (
            (2 3 3 9)
        );
    }

    TransducerLateralWall
    {
        type wall;
        faces
        (
            (6 2 9 12)
        );
    }

    FreeSurface
    {
        type patch;
        faces
        (
            (7 6 12 13)
        );
    }

    Walls
    {
        type wall;
        faces
        (
            (1 0 0 8)
            (4 1 8 10)
            (4 5 11 10)
            (5 7 13 11)
        );
    }

    Axis
    {
        type empty;
        faces
        (
            (0 3 3 0)
        );
    }

    Front
    {
        type wedge;
        faces
        (
            (8 9 3 0)
            (8 10 11 9)
            (11 13 12 9)
        );
    }

    Back
    {
        type wedge;
        faces
        (
            (0 3 2 1)
            (1 2 5 4)
            (2 6 7 5)
        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //
