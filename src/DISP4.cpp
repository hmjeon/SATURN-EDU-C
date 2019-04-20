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

#include "DISP4.h"

CPlaneStress_DISP4R::CPlaneStress_DISP4R(void)
{
}

CPlaneStress_DISP4R::~CPlaneStress_DISP4R(void)
{
}

template <typename T>
void InitArray(T Amat[], int num, T val)
{
    for (int i = 0; i < num; i++)
        Amat[i] = val;
}

template <typename T, int A_b>
void InitArray(T Amat[][A_b], int A_a, T val)
{
    for (int i = 0; i < A_a; i++)
        for (int j = 0; j < A_b; j++)
            Amat[i][j] = val;
}

template <typename T, int A_b, int B_b>
void Mulmat(T Cmat[][B_b], T Amat[][A_b], T Bmat[][B_b], int A_a, int B_a)
{
    double sum = 0.0;

    for (int i = 1; i < A_a; i++)
    {
        for (int j = 1; j < B_b; j++)
        {
            sum = 0.0;

            for (int k = 1; k < A_b; k++)
                sum += Amat[i][k] * Bmat[k][j];

            Cmat[i][j] = sum;
        }
    }
}

// calculate stiffness matrix of plane stress element DOF vector (u1, u2, u3, u4, v1, v2, v3, v4)
void CPlaneStress_DISP4R::Plane_Stiffness(double Young, double Poisson, double thick, double node[][nDIM + 1], double Ke[][nNPE * nDPN + 1])
{
    // Young, Poisson              , Young's modulus, Poison's ratio
    // thick                       , thickness
    // node[nNPE][2]               , nodal position of element
    // Ke[nDPN * nNPE][nDPN * nNPE], stiffness matrix

    double r, s, weight;                          // integration point, weight factor
    double det_j;                                 // determinant of Jacobian
    double B[3 + 1][nNPE * nDPN + 1];             // B-matrix (strain-displacement matrix)
    double BT[nNPE * nDPN + 1][3 + 1];            // transpose B-matrix
    double BTC[nNPE * nDPN + 1][3 + 1];
    double BTCB[nNPE * nDPN + 1][nNPE * nDPN + 1];
    double C[3 + 1][3 + 1];                       // material matrix (material law)

    // calculate a material matrix
    Material_law(Young, Poisson, C);

    // numerical integration
    InitArray<double, nNPE * nDPN + 1>(Ke, nNPE * nDPN + 1, 0.0);

    // for quadrilateral element
    for (int i = 1; i < 2 + 1; i++)
    {
        for (int j = 1; j < 2 + 1; j++)
        {
            // Gauss integration point and weight factor
            Gauss22_Point(i, j, &r, &s, &weight);

            // calculate determinant of Jacobian and B-matrix
            Strain_displacement(r, s, node, &det_j, B);

            // integration
            for (int ii = 1; ii < 3 + 1; ii++)
                for (int jj = 1; jj < nNPE * nDPN + 1; jj++)
                    BT[jj][ii] = B[ii][jj];

            // Ke = Ke + weight * thick * BTCB * det_j 
            Mulmat<double, 3 + 1, 3 + 1>(BTC, BT, C, nNPE * nDPN + 1, 3 + 1);
            Mulmat<double, 3 + 1, nNPE * nDPN + 1>(BTCB, BTC, B, nNPE * nDPN + 1, 3 + 1);

            for (int ii = 1; ii < nNPE * nDPN + 1; ii++)
                for (int jj = 1; jj < nNPE * nDPN + 1; jj++)
                    Ke[ii][jj] = Ke[ii][jj] + weight * thick * BTCB[ii][jj] * abs(det_j);
        }
    }
}

// material law for plane stress condition
void CPlaneStress_DISP4R::Material_law(double Young, double Poisson, double C[][3 + 1])
{
    // Young, Poisson, material constants
    // C[3][3]       , material matrix

    InitArray<double, 3 + 1>(C, 3 + 1, 0.0);

    C[1][1] = Young / (1.0 - pow(Poisson, 2.0));
    C[2][2] = C[1][1];
    C[1][2] = Poisson * C[1][1];
    C[2][1] = C[1][2];
    C[3][3] = 0.5 * Young / (1.0 + Poisson);
}

// Gauss integration point (2*2) and weight factor
void CPlaneStress_DISP4R::Gauss22_Point(int i, int j, double* r, double* s, double* weight)
{
    double GaussPoint[2 + 1], w[2 + 1];

    GaussPoint[1] = -0.577350269189626;
    GaussPoint[2] = 0.577350269189626;
    w[1] = 1.000000000000000;
    w[2] = 1.000000000000000;

    *weight = w[i] * w[j];
    *r = GaussPoint[i];
    *s = GaussPoint[j];
}

// B-matrix (strain-displacement matrix)
void CPlaneStress_DISP4R::Strain_displacement(double r, double s, double node[][nDIM + 1], double* det_j, double B[][nNPE * nDPN + 1])
{
    double dHxy[nDIM + 1][nNPE + 1];

    // calculate dHxy. dHxy[1][:] = dH/dx, dH[2][:] = dHxy/dy
    dHxy_Matrix(r, s, node, det_j, dHxy);

    // B-matrix
    InitArray<double, nNPE * nDPN + 1>(B, 3 + 1, 0.0);

    for (int j = 1; j < nNPE + 1; j++)
    {
        B[1][j] = dHxy[1][j];
        B[2][j + nNPE] = dHxy[2][j];
        B[3][j] = dHxy[2][j];
        B[3][j + nNPE] = dHxy[1][j];
    }
}

// dHxy matrix. dHxy[1][:] = dH/dx, dHxy[2][:] = dH/dy
void CPlaneStress_DISP4R::dHxy_Matrix(double r, double s, double node[][nDIM + 1], double* det_j, double dHxy[][nNPE + 1])
{
    double buf, dHrs[nDIM + 1][nNPE + 1], Jacob[2 + 1][2 + 1];

    // dHrs = derivative of shape function (dH/dr, dH/ds)
    dHrs_Matrix(r, s, dHrs);

    // Jacob = Jacobian Matrix, Jacob = dHrs * node
    Mulmat<double, nNPE + 1, nDIM + 1>(Jacob, dHrs, node, nDIM + 1, nNPE + 1);

    // Jacob => inverse of jacobian Matrix
    *det_j = Jacob[1][1] * Jacob[2][2] - Jacob[1][2] * Jacob[2][1];
    Jacob[1][2] = -Jacob[1][2];
    Jacob[2][1] = -Jacob[2][1];
    buf = Jacob[1][1];
    Jacob[1][1] = Jacob[2][2];
    Jacob[2][2] = buf;

    for (int i = 1; i < 2 + 1; i++)
        for (int j = 1; j < 2 + 1; j++)
            Jacob[i][j] = Jacob[i][j] / (*det_j);

    // dHxy[1][:] = dH/dx, dHxy[2][:] = dH/dy, dHxy = Jacob * dHrs
    Mulmat<double, 2 + 1, nNPE + 1>(dHxy, Jacob, dHrs, 2 + 1, nDIM + 1);
}

// derivatives of shape functions. dHrs[1][:] = dH/dr, dHrs[2][:] = dH/ds
void CPlaneStress_DISP4R::dHrs_Matrix(double r, double s, double dHrs[][nNPE + 1])
{
    dHrs[1][1] = 0.25 * (1.0 + s);
    dHrs[1][2] = -0.25 * (1.0 + s);
    dHrs[1][3] = -0.25 * (1.0 - s);
    dHrs[1][4] = 0.25 * (1.0 - s);

    dHrs[2][1] = 0.25 * (1.0 + r);
    dHrs[2][2] = 0.25 * (1.0 - r);
    dHrs[2][3] = -0.25 * (1.0 - r);
    dHrs[2][4] = -0.25 * (1.0 + r);
}

// equivalent nodal loads
void CPlaneStress_DISP4R::Plane_Load(double node[][nDIM + 1], double q[], double nodal_load[])
{
    // node[nNPE][2]          , node position
    // q[2]                   , body forces in x- and y-directions
    // nodal_load[nDPN * nNPE], element load vector

    double H[nNPE + 1];     // shape functions
    double r, s, weight;    // Gauss point, weight factor
    double det_j;           // determinant of Jacobian

    InitArray(nodal_load, nNPE * nDPN + 1, 0.0);

    // numerical integration
    for (int i = 1; i < 2 + 1; i++)
    {
        for (int j = 1; j < 2 + 1; j++)
        {
            // Gauss points
            Gauss22_Point(i, j, &r, &s, &weight);

            // determinant of Jacobian and shape function (H)
            H_Matrix(r, s, node, H, &det_j);

            // equivalent nodal load vector
            for (int k = 1; k < nNPE + 1; k++)
            {
                nodal_load[k] = nodal_load[k] + weight * abs(det_j) * H[k] * q[1];
                nodal_load[k + nNPE] = nodal_load[k + nNPE] + weight * abs(det_j) * H[k] * q[2];
            }
        }
    }
}

// shape function (H matrix), determinant of Jacobian
void CPlaneStress_DISP4R::H_Matrix(double r, double s, double node[][nDIM + 1], double H[], double* det_j)
{
    double dHrs[nDIM + 1][nNPE + 1], Jacob[2 + 1][2 + 1];

    // shape function (H)
    Shape_Function(r, s, H);

    // dHrs = derivative of shape functions (dH/dr, dH/ds)
    dHrs_Matrix(r, s, dHrs);

    // Jacob = Jacobian matrix, Jacob = dHrs * node
    Mulmat<double, nNPE + 1, nDIM + 1>(Jacob, dHrs, node, nDIM + 1, nNPE + 1);

    // determinant of Jacobian matrix
    *det_j = Jacob[1][1] * Jacob[2][2] - Jacob[1][2] * Jacob[2][1];
}

// calculate shape functions
void CPlaneStress_DISP4R::Shape_Function(double r, double s, double H[])
{
    // r, s   , natural coordinate
    // H[nNPE], shape functions

    H[1] = 0.25 * (1.0 + r) * (1.0 + s);
    H[2] = 0.25 * (1.0 - r) * (1.0 + s);
    H[3] = 0.25 * (1.0 - r) * (1.0 - s);
    H[4] = 0.25 * (1.0 + r) * (1.0 - s);
}

// calculate stress (Sx, Sy, Sxy)
void CPlaneStress_DISP4R::Plane_Stress(double Young, double Poisson, double Node[][nDIM + 1], double Displace[], double Stress[][3 + 1])
{
    // Young, Poisson       , Young's modulus, Poisson's ratio
    // node[nNPE][2]        , nodal position of element
    // displace[nDPN * nNPE], displacement vector
    // Stress[nNPE][3]      , stress (Sx, Sy, Sxy)

    double C[3 + 1][3 + 1], B[3 + 1][nNPE * nDPN + 1];
    double CB[3 + 1][nNPE * nDPN + 1];
    double CBDisp[3 + 1];
    double buf_det;
    double r, s, weight;

    double r1[4 + 1] = { 0.0, 1.0, 0.0, 0.0, 1.0 };
    double s1[4 + 1] = { 0.0, 1.0, 1.0, 0.0, 0.0 };

    // calculate a material matrix
    Material_law(Young, Poisson, C);

    for (int i = 1; i < 2 + 1; i++)
    {
        for (int j = 1; j < 2 + 1; j++)
        {
            // set positions where stresses are out
            Gauss22_Point(i, j, &r, &s, &weight);   // at Gauss points

            //r = r1[i * 2 + j - 2];
            //s = s1[i * 2 + j - 2];

            // calculate B-matrix
            Strain_displacement(r, s, Node, &buf_det, B);

            Mulmat<double, 3 + 1, nNPE * nDPN + 1>(CB, C, B, 3 + 1, 3 + 1); // CB matrix

            for (int k = 1; k < 3 + 1; k++)                                 // CB * displacement
            {
                CBDisp[k] = 0.0;

                for (int kk = 1; kk < nNPE * nDPN + 1; kk++)
                    CBDisp[k] = CBDisp[k] + CB[k][kk] * Displace[kk];
            }

            // calculate stresses
            for (int k = 1; k < 3 + 1; k++)
                Stress[i * 2 + j - 2][k] = CBDisp[k];
        }
    }
}