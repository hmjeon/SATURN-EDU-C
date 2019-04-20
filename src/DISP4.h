/*
!
! =============================================================================
!
! Module - PlaneStress for the 4-node element
! Last Updated : 12/25/2018, by Hyungmin Jun (hyungminjun@outlook.com)
!
! =============================================================================
!
! This is part of SATURN-DISP4-C software that is the finite element analysis
! package for the plane stress condition with the 4-node finite element.
! Copyright 2018 Hyungmin Jun. All rights reserved.
!
! License - GPL version 3
! PERDIX is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or any later version.
! PERDIX is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU General Public License
! for more details.
! You should have received a copy of the GNU General Public License along with
! this program. If not, see <http://www.gnu.org/licenses/>.
!
! -----------------------------------------------------------------------------
!
*/
#pragma once

#include <iostream>
#include <cmath>
#include "Def.h"

using namespace std;

void InitArray(double Array[], int num, const double initvalue);
template <typename T, int bnum> void InitArray(T Array[][bnum], int anum, const T initvalue);
template <typename T, int A_b, int B_b> void Mulmat(T C[][B_b], T A[][A_b], T B[][B_b], int A_a, int B_a);

class CPlaneStress_DISP4R
{
public:
    CPlaneStress_DISP4R(void);
    ~CPlaneStress_DISP4R(void);

    // calculate stiffness matrix of plane stress element DOF vector (u1, u2, u3, u4, v1, v2, v3, v4)
    void Plane_Stiffness(double Young, double Poisson, double thick, double node[][nDIM + 1], double Ke[][nNPE * nDPN + 1]);

    // equivalent nodal loads
    void Plane_Load(double node[][nDIM + 1], double q[], double nodal_load[]);

    // calculate stress (Sx, Sy, Sxy)
    void Plane_Stress(double Young, double Poisson, double node[][nDIM + 1], double Displace[], double Stress[][3 + 1]);

private:
    // material law for plane stress condition
    void Material_law(double Young, double Poisson, double C[][3 + 1]);

    // Gauss integration point (2*2) and weight factor
    void Gauss22_Point(int i, int j, double* r, double* s, double* weight);

    // B-matrix (strain-displacement matrix)
    void Strain_displacement(double r, double s, double node[][nDIM + 1], double* det_j, double B[][nNPE * nDPN + 1]);

    // dHxy matrix. dHxy[1][:] = dH/dx, dHxy[2][:] = dH/dy
    void dHxy_Matrix(double r, double s, double node[][nDIM + 1], double* det_j, double dHxy[][nNPE + 1]);

    // derivatives of shape functions. dHrs[1][:] = dH/dr, dHrs[2][:] = dH/ds
    void dHrs_Matrix(double r, double s, double dHrs[][nNPE + 1]);

    // shape function (H matrix), determinant of Jacobian
    void H_Matrix(double r, double s, double node[][nDIM + 1], double H[], double* det_j);

    // calculate shape functions
    void Shape_Function(double r, double s, double H[]);
};