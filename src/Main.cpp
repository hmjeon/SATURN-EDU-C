/*
!
! =============================================================================
!
! SATURN-DISP4-pub  v1.0
! Last Updated : 04/19/2019, by Hyungmin Jun (hyungminjun@outlook.com)
!
! =============================================================================
!
! SATURN-DISP4-pub is an open-source for the finite element analysis
! package for the planestress condition with the 4-node element.
! Copyright 2018 Hyungmin Jun. All rights reserved.
!
! License - GPL version 3
! This program is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or any later version.
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU General Public License
! for more details.
! You should have received a copy of the GNU General Public License along with
! this program. If not, see <http://www.gnu.org/licenses/>.
!
! -----------------------------------------------------------------------------
!
*/

#include "Funcs.h"
#include "Solver_Skyline.h"
#include "DISP4.h"
#include "Def.h"

// Define global variables
NodeType*               Node;           // Node
ElementType*            Element;        // Element
MaterialType            Material;       // Material
CPlaneStress_DISP4R     PlaneStress;    // Plane stress

// # of nodes, # of elements
int node_n, element_n;

// # of DOFs(total), # of DOFs (free), # of DOFs (fixed)
int total_dof_n, free_dof_n, fixed_dof_n;

// Stiffness, displacement, load
double *Kt ,*U, *R;
int    *c_index;        // Column index

int main(void)
{
    clock_t t_main_start, t_main_end;
    clock_t t_assemble_start, t_assemble_end;
    clock_t t_solve_start, t_solve_end;

    t_main_start = clock();

    // Set domain and properties
#ifdef _FILEINPUT
    string filename = "../inputs/ex1";
    ifstream fin(filename + ".txt");
    Read_Input_File(fin);
    fin.close();
#endif

#ifndef _FILEINPUT
    string filename = "ps_prob";
    Set2DQuadEleRectDomain();
#endif

    // Open files
    ofstream fout1(filename + "_out.txt");
    ofstream fout2(filename + "_res.txt");
    ofstream fout3(filename + "_tec.dat");        // For Tecplot plotting
    ofstream fout4(filename + "_pos.txt");        // For Matlab plotting

    // Calculate # of total DOF, # of free DOF, # of fixed DOF and assign equation numbers
    Set_DOF_Number(&total_dof_n, &free_dof_n, &fixed_dof_n, fout2);

    // Calculate column index
    Calculate_Index(fout2);

    // Print information
    Print_Info();

    // Allocate memory
    Kt = AllocateDynamicArray<double>(c_index[free_dof_n + 1]); // Total stiffness vector (Kt)
    U = AllocateDynamicArray<double>(total_dof_n + 1);          // Displacement vector
    R = AllocateDynamicArray<double>(total_dof_n + 1);          // Load vector

    // Assemble total stiffness matrix
    cout << " 1/ Assembling Stiffness and Load" << endl;
    t_assemble_start = clock();
    Assemble_Kt(fout1);

    // Assemble load vector
    Assemble_Load();
    t_assemble_end = t_solve_start = clock();

    // Solve linear system (find U in K*U=R)
    cout << " 2/ Solving Linear System" << endl;
    Solve_Equations();
    t_solve_end = clock();

    // Calculate stress and print solutions
    cout << " 3/ Printing Output Files" << endl;
    Displacement_Stress(fout2, fout3, fout4);

    // Close files
    CloseFile(fout1, fout2, fout3, fout4);
    cout << " 4/ Completed !" << endl << endl;

#ifndef _FILEINPUT
    cout << " [" << setw(3) << N_I_NODE - 1 << " X " << setw(3) << N_J_NODE - 1 << "]";
    cout << setw(3) << ", DOFs : " << setw(5) << total_dof_n << "[" << setw(5) << free_dof_n << "],";
#endif

#ifdef _FILEINPUT
    cout << " [" << setw(3) << filename << "]";
#endif

    cout.setf(ios::scientific); cout.precision(5); cout << "  Strain energy = " << 0.5 * Dot_product(R, U, free_dof_n + 1) << endl;

    // Deallocate memory
    Deallocate(Node, Element, c_index, Kt, U, R);
    t_main_end = clock();

    // Print time consuming
    Print_TimeConsuming(t_main_end, t_main_start, t_assemble_end, t_assemble_start, t_solve_end, t_solve_start);
 
    system("pause");
    return 0;
}

void InitArray(double Array[], int num, const double initvalue)
{
    for (int i = 0; i < num; i++)
        Array[i] = initvalue;
}

template <typename T, int bnum> void InitArray(T Array[][bnum], int anum, const T initvalue)
{
    for (int i = 0; i < anum; i++)
        for (int j = 0; j < bnum; j++)
            Array[i][j] = initvalue;
}

// Allocation of the dynamic array
template <typename T> T* AllocateDynamicArray(int num)
{
    T* vector;
    vector = (T*)malloc(sizeof(T) * num);
    return vector;
}

template <typename T> void InitDynamicArray(T* vector, int num, const T initvalue)
{
    for (int i = 0; i < num; i++)
        vector[i] = initvalue;
}

template <typename T> void FreeDynamicArray(T* vector)
{
    free(vector);
}

double Dot_product(double A[], double B[], int size)
{
    double dotproduct = 0;

    for (int i = 1; i < size + 1; i++)
        dotproduct = dotproduct + A[i] * B[i];

    return dotproduct;
}

void Print_Info(void)
{
#ifndef _DEBUG
    cout << "[RELEASE mode]" << endl;
#endif
#ifdef _DEBUG
    cout << "[DEBUG mode]" << endl;
#endif

    cout << " [SATURN-DISP4-C]" << endl << endl;
    cout << " --------------------------------------------------------------------" << endl;
    cout << " Young's modulus : " << Material.Young << ", Poisson's ratio : " << Material.Poisson << endl;

#ifdef _FILEINPUT
    cout << " thickness : " << Element[1].thickness << endl;
#endif

#ifndef _FILEINPUT
    cout << " Length x : " << X_WIDTH << ", and y : " << Y_WIDTH << ", Thickness : " << THICK << endl;
    cout << " # element in x-dir: " << N_I_NODE - 1 << ",  # element in y-dir: " << N_J_NODE - 1 << ",   [" << N_I_NODE - 1 << " X " << N_J_NODE - 1 << "]" << endl;
#endif

    cout << " # total elements : " << element_n << ",  # total nodes : " << node_n << endl;
    cout << " # total DOFs : " << total_dof_n << ",  # free DOFs : " << free_dof_n << ",  # fix DOFs : " << fixed_dof_n << endl;
    cout << " Required memory = " << c_index[free_dof_n + 1] * 8.0 / 1024.0 / 1024.0 << " MB" << endl;
    cout << " --------------------------------------------------------------------" << endl;
    cout << endl;
}

void Print_TimeConsuming(double m_end, double m_start, double a_end, double a_start, double s_end, double s_start)
{
    double    consuming_m = (double)(m_end - m_start) / CLOCKS_PER_SEC;
    double    consuming_a = (double)(a_end - a_start) / CLOCKS_PER_SEC;
    double    consuming_s = (double)(s_end - s_start) / CLOCKS_PER_SEC;

    cout.precision(3);
    cout << endl;
    cout << setw(20) << " time consuming : "  << consuming_m << " [sec], " << consuming_m / 60.0 << " [min], " << consuming_m / 60.0 / 60.0 << " [hr]" << endl;
    cout << setw(20) << " assembling time : " << consuming_a << " [sec], " << consuming_a / 60.0 << " [min], " << consuming_a / 60.0 / 60.0 << " [hr]" << endl;
    cout << setw(20) << " solving time : "    << consuming_s << " [sec], " << consuming_s / 60.0 << " [min], " << consuming_s / 60.0 / 60.0 << " [hr]\n" << endl;
}

void Read_Input_File(ifstream& fin)
{
    char   strs[256];
    int    bufi;
    double buff;

    // Read nodal information
    fin.getline(strs, 256);
    fin >> node_n;
    fin.getline(strs, 256);

    Node = AllocateDynamicArray<NodeType>(node_n + 1);      // Allocate a node array
    InitNode(Node);

    for (int i = 1; i < node_n + 1; i++)
    {
        fin >> bufi;

        fin >> Node[i].x[1] >> Node[i].x[2];
        fin >> Node[i].bc[1] >> Node[i].bc[2];
        fin >> Node[i].pm[1] >> Node[i].pm[2];
    }

    // Read elemental information
    fin.getline(strs, 256);    fin.getline(strs, 256);
    fin >> element_n;

    Element = AllocateDynamicArray<ElementType>(element_n + 1);    // Allocate a element array 
    InitElement(Element);

    for (int i = 1; i < element_n + 1; i++)
    {
        fin >> bufi;
        fin >> Element[i].cn[1] >> Element[i].cn[2] >> Element[i].cn[3] >> Element[i].cn[4];
        fin >> Element[i].q[1] >> Element[i].q[2];
    }

    // Read properties
    fin.getline(strs, 256); fin.getline(strs, 256);
    fin >> Material.Young;

    fin.getline(strs, 256); fin.getline(strs, 256);
    fin >> buff;

    for (int i = 1; i < element_n + 1; i++)
        Element[i].thickness = buff;

    fin.getline(strs, 256); fin.getline(strs, 256);
    fin >> Material.Poisson;
}

// Initialization of NodeType structure
void InitNode(NodeType* node)
{
    for (int i = 1; i < node_n + 1; i++)
    {
        node[i].x[1] = 0.0;
        node[i].x[2] = 0.0;

        for (int j = 1; j < nDPN + 1; j++)
        {
            node[i].pm[j] = 0.0;
            node[i].bc[j] = 0;
            node[i].eq_n[j] = 0;
        }
    }
}

// Initialization of ElementType structure
void InitElement(ElementType* element)
{
    for (int i = 1; i < element_n + 1; i++)
    {
        element[i].thickness = 0.0;
        element[i].q[1] = 0.0;
        element[i].q[2] = 0.0;

        for (int j = 1; j < nNPE + 1; j++)
            element[i].cn[j] = 0;
    }
}

// Calculate # of total DOF, # of free DOF, # of fixed DOF and assign equation numbers to DOFs
void Set_DOF_Number(int* tn, int* fn, int* cn, ofstream& fout)
{
    // tn, fn, cn # of total DOF, # of free DOF, # of fixed DOF
    fout << " EQUATION NUMBER" << endl;
    fout << " ---------------" << endl;
    fout << "     node   dof    eqn" << endl;

    *tn = node_n * nDPN;    // # of total DOF
    *fn = 0;
    *cn = 0;

    for (int i = 1; i < node_n + 1; i++)
    {
        for (int j = 1; j < nDPN + 1; j++)
        {
            if (Node[i].bc[j] == 0)
            {
                *fn = *fn + 1;
                Node[i].eq_n[j] = *fn;
                fout << setw(7) << i << setw(7) << j << setw(7) << *fn << endl;
            }
            else
            {
                *cn = *cn + 1;
                Node[i].eq_n[j] = *tn - (*cn) + 1;
            }
        }
    }

    fout << endl;
}

// Calculate column index for skyline solver
void Calculate_Index(ofstream& fout)
{
    int* column_h = AllocateDynamicArray<int>(total_dof_n + 1);     // Column height
    InitDynamicArray<int>(column_h, total_dof_n + 1, 0);

    int a_index[nDPN * nNPE + 1];                                   // Index for assemblage
    int buf_sum;

    // Allocate c_index array
    c_index = AllocateDynamicArray<int>(free_dof_n + 1 + 1);
    InitDynamicArray<int>(c_index, free_dof_n + 1 + 1, 0);

    // Column height
    for (int i = 1; i < element_n + 1; i++)
    {
        // Assemblage index
        for (int j = 1; j < nDPN + 1; j++)
            for (int k = 1; k < nNPE + 1; k++)
                a_index[nNPE * j + k - nNPE] = Node[Element[i].cn[k]].eq_n[j];

        // Column height
        for (int k = 1; k < nDPN * nNPE + 1; k++)
            for (int j = 1; j < nDPN * nNPE + 1; j++)
                if (a_index[j] <= free_dof_n && a_index[k] <= free_dof_n)
                    if (a_index[j] < a_index[k])
                        if (a_index[k] - a_index[j] > column_h[a_index[k]])
                            column_h[a_index[k]] = a_index[k] - a_index[j];
    }

    // c_index array
    buf_sum = 1;
    for (int i = 1; i < free_dof_n + 1; i++)
    {
        c_index[i] = buf_sum;
        buf_sum = buf_sum + column_h[i] + 1;
    }

    c_index[free_dof_n + 1] = buf_sum;
    fout << " REQUIRED MEMORY = " << buf_sum * 8.0 / 1000000.0 << " MB" << endl << endl;

    FreeDynamicArray(column_h);
}

// Assemble total stiffness matrix by using equation numbers
void Assemble_Kt(ofstream& fout2)
{
    double eNode[nNPE + 1][nDIM + 1];               // nodal position of 4-node element(x,y)
    double Ke[nNPE * nDPN + 1][nNPE * nDPN + 1];    // stiffness matrix of element
    int    a_index[nNPE * nDPN + 1];                // assemblage index
    int    address;

    InitDynamicArray(Kt, c_index[free_dof_n + 1] - 1 + 1, 0.0);

    for (int i = 1; i < element_n + 1; i++)
    {
        // nodal position of element
        for (int j = 1; j < nNPE + 1; j++)
            for (int k = 1; k < nDIM + 1; k++)
                eNode[j][k] = Node[Element[i].cn[j]].x[k];

        // calculate stiffness matrix of element
        PlaneStress.Plane_Stiffness(Material.Young, Material.Poisson, Element[i].thickness, eNode, Ke);
        //fout2<<" STIFFNESS of ELEMENT : "<<setw(2)<<i<<endl;
        //Print_Matrix(Ke, fout2); 

        // assemblage index
        for (int j = 1; j < nDPN + 1; j++)
            for (int k = 1; k < nNPE + 1; k++)
                a_index[nNPE * j + k - nNPE] = Node[Element[i].cn[k]].eq_n[j];

        // assemble total stiffness matrix
        for (int j = 1; j < nNPE * nDPN + 1; j++)
            for (int k = 1; k < nNPE * nDPN + 1; k++)
                if (a_index[j] <= free_dof_n && a_index[k] <= free_dof_n)
                    if (a_index[j] <= a_index[k])
                    {
                        address = c_index[a_index[k]] + a_index[k] - a_index[j];
                        Kt[address] = Kt[address] + Ke[j][k];
                    }
    }
}

// assemble load vector
void Assemble_Load(void)
{
    double eNode[nNPE + 1][nDIM + 1];
    double NodalLoad[nNPE * nDPN + 1];        // equivalent nodal load

    InitDynamicArray(R, total_dof_n + 1, 0.0);

    // assemble load vector for nodal load
    for (int i = 1; i < node_n + 1; i++)
        for (int j = 1; j < nDPN + 1; j++)
            R[Node[i].eq_n[j]] = Node[i].pm[j];

    // assemble load vector for body force
    for (int i = 1; i < element_n + 1; i++)
    {
        // nodal position of element
        for (int j = 1; j < nNPE + 1; j++)
            for (int k = 1; k < nDIM + 1; k++)
                eNode[j][k] = Node[Element[i].cn[j]].x[k];

        // calculate equivalent nodal load from body force
        PlaneStress.Plane_Load(eNode, Element[i].q, NodalLoad);

        // assemble load vector
        for (int j = 1; j < nDPN + 1; j++)
            for (int k = 1; k < nNPE + 1; k++)
                R[Node[Element[i].cn[k]].eq_n[j]] = R[Node[Element[i].cn[k]].eq_n[j]] + NodalLoad[nNPE * j + k - nNPE];
    }
}

// love linear equations
void Solve_Equations()
{
    InitDynamicArray(U, total_dof_n + 1, 0.0);

    for (int j = 1; j < free_dof_n + 1; j++)
        U[j] = R[j];

    COLSOL(Kt, U, c_index, free_dof_n, 1);
    COLSOL(Kt, U, c_index, free_dof_n, 2);
}

// calculate stress and print solutions
void Displacement_Stress(ofstream& fout_r, ofstream& fout_t, ofstream& fout_m)
{
    const int n_gp = 4;                  // the number of gauss points
    double eNode[nNPE + 1][nDIM + 1];    // nodal position of element
    double displace[nNPE * nDPN + 1];    // nodal displacement vector of element
    double Stress[n_gp + 1][3 + 1];      // Sxx, Syy, Sxy in Gauss points or nodal positions (2*2)

    // print strain energy
    fout_r << " STRAIN ENERGY = " << 0.5 * Dot_product(R, U, free_dof_n + 1) << endl;

    // print nodal displacement
    fout_r << endl;
    fout_r << " DISPLACEMENT " << endl;
    fout_r << " ----------------------------------------" << endl;
    fout_r << "      Node              u               v" << endl;

    for (int i = 1; i < node_n + 1; i++)
    {
        fout_r.precision(6);
        fout_r << setw(9) << i;
        fout_r << setw(16) << U[Node[i].eq_n[1]] << setw(16) << U[Node[i].eq_n[2]] << endl;
    }
    fout_r << endl;

    // set scaling factor for the plotting
    double scale_factor;
    double max_pos = 0.0;
    double max_disp = 0.0;

    for (int i = 1; i < node_n + 1; i++)
    {
        // find maximum value of the displacement and positions
        if (max_disp < abs(U[Node[i].eq_n[1]]))
        {
            max_disp = abs(U[Node[i].eq_n[1]]);
            max_pos = sqrt(pow(Node[i].x[1], 2.0) + pow(Node[i].x[2], 2.0));
        }

        if (max_disp < abs(U[Node[i].eq_n[2]]))
        {
            max_disp = abs(U[Node[i].eq_n[2]]);
            max_pos = sqrt(pow(Node[i].x[1], 2.0) + pow(Node[i].x[2], 2.0));
        }
    }

    scale_factor = (1.2 * max_pos - max_pos) / max_disp;    // 1.2 * max_pos = (scale_factor * max_disp + max_pos)

    // write scale factor into m file
    fout_m << element_n << "\t" << scale_factor << endl;

    double max_stress_vm = 1.0e-12;

    // for Tecplot output
    fout_t << "title = \"CFEM_planestress_DISP4\"" << endl;
    fout_t << "variables = \"x\", \"y\", \"Tvm\"" << endl;
    fout_t << "ZONE F=FEPOINT, N=" << element_n * nNPE << ", E=" << element_n << ", ET=QUADRILATERAL" << endl << endl;

    for (int i = 1; i < element_n + 1; i++)
    {
        // nodal position of element
        for (int j = 1; j < nNPE + 1; j++)
            for (int k = 1; k < nDIM + 1; k++)
                eNode[j][k] = Node[Element[i].cn[j]].x[k];

        // displacement vector of element
        for (int j = 1; j < nDPN + 1; j++)
            for (int k = 1; k < nNPE + 1; k++)
                displace[nNPE * j + k - nNPE] = U[Node[Element[i].cn[k]].eq_n[j]];

        // calculate stress of element
        PlaneStress.Plane_Stress(Material.Young, Material.Poisson, eNode, displace, Stress);

        // print stress
        fout_r << "  STRESS of ELEMENT : " << i << endl;
        fout_r << " -----------------------------------------------------------------" << endl;
        fout_r << "  Position        Txx           Tyy           Txy           Tvm" << endl;

        double VonMises[n_gp + 1];
        for (int j = 1; j < n_gp + 1; j++)            // gauss point
        {
            fout_r.precision(4);
            fout_r << setw(7) << j << setw(14) << Stress[j][1] << setw(14) << Stress[j][2] << setw(14) << Stress[j][3];

            // Von Mises stress
            VonMises[j] = sqrt((pow(Stress[j][1], 2.0) - (Stress[j][1] * Stress[j][2])
                + pow(Stress[j][2], 2.0) + 3.0 * (pow(Stress[j][3], 2.0))));

            fout_r << setw(14) << VonMises[j] << endl;
        }
        fout_r << endl;

        // print deformed shape and stress for post processing by using MATLAB
        fout_m.precision(4);
        for (int j = 1; j < nNPE + 1; j++)
        {
            for (int k = 1; k < nDIM + 1; k++)
                fout_m << setw(12) << eNode[j][k];

            fout_m << setw(12) << displace[j] << setw(12) << displace[j + nNPE];

            for (int k = 1; k < 3 + 1; k++)
                fout_m << setw(12) << Stress[j][k];

            fout_m << setw(12) << VonMises[j];

            fout_t << setw(20) << eNode[j][1] + scale_factor * displace[j];
            fout_t << setw(20) << eNode[j][2] + scale_factor * displace[j + nNPE];
            fout_t << setw(20) << VonMises[j] << endl;

            max_stress_vm = MAX(VonMises[j], max_stress_vm);
        }
        fout_m << endl;
    }
    fout_t << endl;

    // set connectivities into Tecplot file
    for (int i = 1; i < element_n + 1; i++)
    {
        for (int j = 1; j < nNPE + 1; j++)
            fout_t << nNPE * (i - 1) + j << "\t";
        fout_t << endl;
    }

    fout_r << "Maximum effective stress = " << max_stress_vm << endl;
    fout_r << "Strain energy = " << 0.5 * Dot_product(R, U, free_dof_n + 1) << endl;
}

// print (nNPE * nDPN) x (nNPE * nDPN) matrix
void Print_Matrix(double M[][nNPE * nDPN + 1], ofstream& fout2)
{
    fout2 << " ---------------------------" << endl;
    for (int i = 1; i < nNPE * nDPN + 1; i++)
    {
        for (int j = 1; j < nNPE * nDPN + 1; j++)
        {
            fout2.precision(3);
            fout2 << setw(12) << M[i][j];
        }
        fout2 << endl;
    }
    fout2 << endl;
}

// deallocate memory
void Deallocate(NodeType* Node, ElementType* Element, int* c_index, double* Kt, double* U, double* R)
{
    FreeDynamicArray(Node);
    FreeDynamicArray(Element);
    FreeDynamicArray(c_index);
    FreeDynamicArray(Kt);
    FreeDynamicArray(U);
    FreeDynamicArray(R);
}

// close files
void CloseFile(ofstream& fout1, ofstream& fout2, ofstream& fout3, ofstream& fout4)
{
    fout1.close();
    fout2.close();
    fout3.close();
    fout4.close();
}

// set rectangular domain with triangular element
void Set2DQuadEleRectDomain(void)
{
    // given parameter
    const double young    = YOUNG;          // Young's modulus
    const double possion  = POSSION;        // Poisson's ratio
    const double bxforce  = BXFORCE;        // x-direction body force
    const double byforce  = BYFORCE;        // y-direction body force
    const double thick    = THICK;          // thickness
    const double x_width  = X_WIDTH;        // x length
    const double y_width  = Y_WIDTH;        // y length
    const int    n_i_node = N_I_NODE;       // the number of nodes in the x-direction
    const int    n_j_node = N_J_NODE;       // the number of nodes in the y-direction

    // parameter for setting domain
    int    n_i_element = n_i_node - 1;
    int    n_j_element = n_j_node - 1;

    double del_x = x_width / double(n_i_node - 1);
    double del_y = y_width / double(n_j_node - 1);

    node_n = n_i_node * n_j_node;
    element_n = n_i_element * n_j_element;

    Node = AllocateDynamicArray<NodeType>(node_n + 1);          // allocate a node array 
    Element = AllocateDynamicArray<ElementType>(element_n + 1); // allocate a element array

    InitNode(Node);
    InitElement(Element);

    // set nodal position
    for (int j = 1; j < n_j_node + 1; j++)
    {
        for (int i = 1; i < n_i_node + 1; i++)
        {
            int numbering = n_i_node * (j - 1) + i;             // int numbering = n_i_node * j + i - n_i_node;

            Node[numbering].x[1] = del_x * double(i - 1);
            Node[numbering].x[2] = del_y * double(j - 1);
        }
    }

    // set connectivity
    for (int j = 1; j < n_j_element + 1; j++)
    {
        for (int i = 1; i < n_i_element + 1; i++)
        {
            int numbering = n_i_element * (j - 1) + i;

            Element[numbering].cn[1] = n_i_node * (j - 1) + i;
            Element[numbering].cn[2] = n_i_node * (j - 1) + i + 1;
            Element[numbering].cn[3] = n_i_node * j + i + 1;
            Element[numbering].cn[4] = n_i_node * j + i;
        }
    }

    // imposing boundary condition
    int x_fix_surf = 1;     // left side B.C. : u, v = 1
    for (int j = 1; j < n_j_node + 1; j++)
    {
        int numbering = n_i_node * (j - 1) + x_fix_surf;

        Node[numbering].bc[1] = 1;
        Node[numbering].bc[2] = 1;
    }

    // set body force
    for (int i = 1; i < element_n + 1; i++)
    {
        Element[i].q[1] = bxforce;
        Element[i].q[2] = byforce;
    }

    // set thickness
    for (int i = 1; i < element_n + 1; i++)
        Element[i].thickness = thick;

    // set properties
    Material.Young = young;
    Material.Poisson = possion;
}