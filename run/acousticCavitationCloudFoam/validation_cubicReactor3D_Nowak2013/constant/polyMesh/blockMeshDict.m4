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
// ------!!!!!!   run "m4 blockMeshDict.m4 > blockMeshDict"  !!!!!!--------- //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// General m4 macros
changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'use Math::Trig; use POSIX; print ($1)')])
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])

define(hex2D, hex ($1b $2b $3b $4b $1t $2t $3t $4t))
define(btQuad, ($1b $2b $2t $1t))
define(topQuad, ($1u $4u $3u $2u))
define(bottomQuad, ($1b $2b $3b $4b))

// Special treatment for the upper square
define(hex2Du, hex ($1t $2t $3t $4t $1u $2u $3u $4u))
define(tuQuad, ($1t $2t $2u $1u))
define(sonoBottomQuad, ($1t $4t $3t $2t))
//define(sonoLateralQuad, ($1t $2t $2u $1u))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 0.001;

// Global refinement parameter
define(ref, 2.5)

// Inner square side half
define(s, 2)

// Inner square side curvature
define(sc, 2.5)

// cylinder radius
define(r, 5)

// Cuvette height
define(Zc, 40)

// Sonotrode immersed height
define(ZSono, 8)

// Height of cylinder
define(z, calc(Zc-ZSono))

// Base z
define(Zb, 0)

// Outlet z
define(Zt, calc(Zb + z))

//// Cells
// Number of cells at inner square
define(Ns, calc(ceil(ref*4)))

// Number of cells between inner square and circle
define(Ni, calc(ceil(ref*2)))

// Number of cells in the cylinder height
define(Nz, calc(ceil(ref*10)))


//// Additionally: outer square region und upper block next to the sonotrode

// Outer square side half
define(o, 20)

// Number of cells of outer square in radial direction
define(No, calc(ceil(ref*8)))

// Number of cells of upper outer square in height
define(Nu, calc(ceil(ref*4)))

//// Simple grading

// Simple grading parameter for lower z
define(gzl, 0.15)

// Simple grading parameter for upper z
define(gzu, 3)

// Simple grading parameter for cirlce radial direction
define(grr, 1)

// Simple grading parameter for outer square radial deriction
define(gro, 3)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

define(vert, (x$1$2 y$1$2 $3))
define(evert, (ex$1$2 ey$1$2 $3))

// 45 degree points angle
define(a0, -45)
define(a1, -135)
define(a2, 135)
define(a3, 45)

// Half of 45 degree points angle
define(ea0, 0)
define(ea1, -90)
define(ea2, 180)
define(ea3, 90)

define(ca0, calc(cos((pi/180)*a0)))
define(ca1, calc(cos((pi/180)*a1)))
define(ca2, calc(cos((pi/180)*a2)))
define(ca3, calc(cos((pi/180)*a3)))

define(sa0, calc(sin((pi/180)*a0)))
define(sa1, calc(sin((pi/180)*a1)))
define(sa2, calc(sin((pi/180)*a2)))
define(sa3, calc(sin((pi/180)*a3)))

define(cea0, calc(cos((pi/180)*ea0)))
define(cea1, calc(cos((pi/180)*ea1)))
define(cea2, calc(cos((pi/180)*ea2)))
define(cea3, calc(cos((pi/180)*ea3)))

define(sea0, calc(sin((pi/180)*ea0)))
define(sea1, calc(sin((pi/180)*ea1)))
define(sea2, calc(sin((pi/180)*ea2)))
define(sea3, calc(sin((pi/180)*ea3)))

// Inner square x and y position

// x
define(x00, s)
define(x01, calc(-1.0*s))
define(x02, calc(-1.0*s))
define(x03, s)

// y
define(y00, calc(-1.0*s))
define(y01, calc(-1.0*s))
define(y02, s)
define(y03, s)

// Circle x and y positions

// x
define(x10, calc(r*ca0))
define(x11, calc(r*ca1))
define(x12, calc(r*ca2))
define(x13, calc(r*ca3))

// y
define(y10, calc(r*sa0))
define(y11, calc(r*sa1))
define(y12, calc(r*sa2))
define(y13, calc(r*sa3))

// Inner square x and y position middle curvatures

// x
define(ex00, sc)
define(ex01, 0)
define(ex02, calc(-1.0*sc))
define(ex03, 0)

// y
define(ey00, 0)
define(ey01, calc(-1.0*sc))
define(ey02, 0)
define(ey03, sc)

// Circle x and y positions middle curvatures

// x
define(ex10, calc(r*cea0))
define(ex11, calc(r*cea1))
define(ex12, calc(r*cea2))
define(ex13, calc(r*cea3))

// y
define(ey10, calc(r*sea0))
define(ey11, calc(r*sea1))
define(ey12, calc(r*sea2))
define(ey13, calc(r*sea3))

// Outer squre x and y positions

//x
define(x20, o)
define(x21, calc(-1.0*o))
define(x22, calc(-1.0*o))
define(x23, o)

// y
define(y20, calc(-1.0*o))
define(y21, calc(-1.0*o))
define(y22, o)
define(y23, o)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

vertices
(
    vert(0, 0, Zb) vlabel(s0b)
    vert(0, 1, Zb) vlabel(s1b)
    vert(0, 2, Zb) vlabel(s2b)
    vert(0, 3, Zb) vlabel(s3b)

    vert(1, 0, Zb) vlabel(r0b)
    vert(1, 1, Zb) vlabel(r1b)
    vert(1, 2, Zb) vlabel(r2b)
    vert(1, 3, Zb) vlabel(r3b)

    vert(0, 0, Zt) vlabel(s0t)
    vert(0, 1, Zt) vlabel(s1t)
    vert(0, 2, Zt) vlabel(s2t)
    vert(0, 3, Zt) vlabel(s3t)

    vert(1, 0, Zt) vlabel(r0t)
    vert(1, 1, Zt) vlabel(r1t)
    vert(1, 2, Zt) vlabel(r2t)
    vert(1, 3, Zt) vlabel(r3t)

    // Outer square lower
    vert(2, 0, Zb) vlabel(o0b)
    vert(2, 1, Zb) vlabel(o1b)
    vert(2, 2, Zb) vlabel(o2b)
    vert(2, 3, Zb) vlabel(o3b)

    vert(2, 0, Zt) vlabel(o0t)
    vert(2, 1, Zt) vlabel(o1t)
    vert(2, 2, Zt) vlabel(o2t)
    vert(2, 3, Zt) vlabel(o3t)

    // Outer square upper
    vert(1, 0, Zc) vlabel(r0u)
    vert(1, 1, Zc) vlabel(r1u)
    vert(1, 2, Zc) vlabel(r2u)
    vert(1, 3, Zc) vlabel(r3u)

    vert(2, 0, Zc) vlabel(o0u)
    vert(2, 1, Zc) vlabel(o1u)
    vert(2, 2, Zc) vlabel(o2u)
    vert(2, 3, Zc) vlabel(o3u)
);

blocks
(
    //block0
    hex2D(s1, s0, s3, s2) square (Ns Ns Nz) simpleGrading (1 1 gzl)

    //block1
    hex2D(s0, r0, r3, s3) innerCircle (Ni Ns Nz) simpleGrading (grr 1 gzl)

    //block2
    hex2D(s3, r3, r2, s2) innerCircle (Ni Ns Nz) simpleGrading (grr 1 gzl)

    //block3
    hex2D(s2, r2, r1, s1) innerCircle (Ni Ns Nz) simpleGrading (grr 1 gzl)

    //block4
    hex2D(s1, r1, r0, s0) innerCircle (Ni Ns Nz) simpleGrading (grr 1 gzl)

    //// Outer square lower
    //block5
    hex2D(r0, o0, o3, r3) outerSquare (No Ns Nz) simpleGrading (gro 1 gzl)

    //block6
    hex2D(r3, o3, o2, r2) outerSquare (No Ns Nz) simpleGrading (gro 1 gzl)

    //block7
    hex2D(r2, o2, o1, r1) outerSquare (No Ns Nz) simpleGrading (gro 1 gzl)

    //block8
    hex2D(r1, o1, o0, r0) outerSquare (No Ns Nz) simpleGrading (gro 1 gzl)

    //// Outer square upper
    //block9
    hex2Du(r0, o0, o3, r3) outerSquareUp (No Ns Nu) simpleGrading (gro 1 gzu)

    //block10
    hex2Du(r3, o3, o2, r2) outerSquareUp (No Ns Nu) simpleGrading (gro 1 gzu)

    //block11
    hex2Du(r2, o2, o1, r1) outerSquareUp (No Ns Nu) simpleGrading (gro 1 gzu)

    //block12
    hex2Du(r1, o1, o0, r0) outerSquareUp (No Ns Nu) simpleGrading (gro 1 gzu)
);

edges
(
    //Circle edges
    arc r3b r0b evert(1, 0, Zb)
    arc r0b r1b evert(1, 1, Zb)
    arc r1b r2b evert(1, 2, Zb)
    arc r2b r3b evert(1, 3, Zb)

    arc r3t r0t evert(1, 0, Zt)
    arc r0t r1t evert(1, 1, Zt)
    arc r1t r2t evert(1, 2, Zt)
    arc r2t r3t evert(1, 3, Zt)

    arc s3b s0b evert(0, 0, Zb)
    arc s0b s1b evert(0, 1, Zb)
    arc s1b s2b evert(0, 2, Zb)
    arc s2b s3b evert(0, 3, Zb)

    arc s3t s0t evert(0, 0, Zt)
    arc s0t s1t evert(0, 1, Zt)
    arc s1t s2t evert(0, 2, Zt)
    arc s2t s3t evert(0, 3, Zt)

    // Upper
    arc r3u r0u evert(1, 0, Zc)
    arc r0u r1u evert(1, 1, Zc)
    arc r1u r2u evert(1, 2, Zc)
    arc r2u r3u evert(1, 3, Zc)
);

patches
(
    wall metallicWall
    (
        btQuad(o2, o1)
        tuQuad(o2, o1)
    )

    wall acrylicWalls
    (
        btQuad(o0, o3)
        btQuad(o1, o0)
        btQuad(o3, o2)
        tuQuad(o0, o3)
        tuQuad(o1, o0)
        tuQuad(o3, o2)
        bottomQuad(s3, s0, s1, s2)
        bottomQuad(s3, r3, r0, s0)
        bottomQuad(s2, r2, r3, s3)
        bottomQuad(s1, r1, r2, s2)
        bottomQuad(s0, r0, r1, s1)
        bottomQuad(r3, o3, o0, r0)
        bottomQuad(r2, o2, o3, r3)
        bottomQuad(r1, o1, o2, r2)
        bottomQuad(r0, o0, o1, r1)
        topQuad(r3, o3, o0, r0)
        topQuad(r2, o2, o3, r3)
        topQuad(r1, o1, o2, r2)
        topQuad(r0, o0, o1, r1)
    )

    wall sonotrodeBottom
    (
        sonoBottomQuad(s3, s0, s1, s2)
        sonoBottomQuad(s3, r3, r0, s0)
        sonoBottomQuad(s2, r2, r3, s3)
        sonoBottomQuad(s1, r1, r2, s2)
        sonoBottomQuad(s0, r0, r1, s1)
    )

    wall sonotrodeLateralWall
    (
        tuQuad(r3, r0)
        tuQuad(r2, r3)
        tuQuad(r1, r2)
        tuQuad(r0, r1)
    )
);

mergePatchPairs
(
);
