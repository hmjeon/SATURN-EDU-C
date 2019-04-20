#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include "Def.h"

using namespace std;

// define a type for node
struct NodeType
{
    double x[nDIM + 1];        // nodal position (x, y)
    double pm[nDPN + 1];       // nodal force (Px, Py)
    int    bc[nDPN + 1];       // displacement BC (u, v) (1=fixed, 0=free)
    int    eq_n[nDPN + 1];     // equation number (u, v)
};

// define a type for element 
struct ElementType
{
    int    cn[nNPE + 1];       // connectivity
    double thickness;          // thickness
    double q[nDIM + 1];        // distributed load
};

// define a type for material
struct MaterialType
{
    double Young;              // Young's modulus
    double Poisson;            // Poison's ratio
};

void InitArray(double Array[], int num, const double initvalue);

template <typename T, int bnum> void InitArray(T Array[][bnum], int anum, const T initvalue);

template <typename T> T* AllocateDynamicArray(int num);

template <typename T> void InitDynamicArray(T* vector, int num, const T initvalue);

template <typename T> void FreeDynamicArray(T* vector);

double Dot_product(double A[], double B[], int size);

void Print_Info(void);

void Print_TimeConsuming(double m_end, double m_start, double a_end, double a_start, double s_end, double s_start);

void Read_Input_File(ifstream& fin);

// initialization of NodeType structure
void InitNode(NodeType* node);

// initialization of ElementType structure
void InitElement(ElementType* element);

// set rectangular domain with triangular element
void Set2DQuadEleRectDomain(void);

// calculate # of total DOF, # of free DOF, # of fixed DOF and assign equation numbers to DOFs
void Set_DOF_Number(int* tn, int* fn, int* cn, ofstream& fout3);

// calculate column index for skyline solver
void Calculate_Index(ofstream& fout3);

// assemble total stiffness matrix by using equation numbers
void Assemble_Kt(ofstream& fout2);

// assemble load vector
void Assemble_Load(void);

// love linear equations
void Solve_Equations();

// calculate stress and print solutions
void Displacement_Stress(ofstream& fout_r, ofstream& fout_t, ofstream& fout_m);

// print (nNPE * nDPN) x (nNPE * nDPN) matrix
void Print_Matrix(double M[][nNPE * nDPN + 1], ofstream& fout2);

// deallocate memory
void Deallocate(NodeType* Node, ElementType* Element, int* c_index, double* Kt, double* U, double* R);

// close files
void CloseFile(ofstream& fout1, ofstream& fout2, ofstream& fout3, ofstream& fout4);