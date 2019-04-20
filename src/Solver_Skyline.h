#pragma once

/*
!
! =============================================================================
!
! Module - Solver - Skyline storage
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

#include <iostream>
#include <fstream>
#include <cmath>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

void COLSOL(double A[],double V[], int MAXA[],int NN, int KKK)
{
    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    // .                                                                   .
    // .   P R O G R A M                                                   .
    // .        TO SOLVE FINITE ELEMENT STATIC EQUILIBRIUM EQUATIONS IN    .
    // .        CORE, USING COMPACTED STORAGE AND COLUMN REDUCTION SCHEME  .
    // .                                                                   .
    // .  - - INPUT VARIABLES - -                                          .
    // .        A(NWK)    = STIFFNESS MATRIX STORED IN COMPACTED FORM      .
    // .        V(NN)     = RIGHT-HAND-SIDE LOAD VECTOR                    .
    // .        MAXA(NNM) = VECTOR CONTAINING ADDRESSES OF DIAGONAL        .
    // .                    ELEMENTS OF STIFFNESS MATRIX IN A              .
    // .        NN        = NUMBER OF EQUATIONS                            .
    // .        NWK       = NUMBER OF ELEMENTS BELOW SKYLINE OF MATRIX     .
    // .        NNM       = NN + 1                                         .
    // .        KKK       = INPUT FLAG                                     .
    // .            EQ. 1   TRIANGULARIZATION OF STIFFNESS MATRIX          .
    // .            EQ. 2   REDUCTION AND BACK-SUBSTITUTION OF LOAD VECTOR .
    // .        IOUT      = UNIT NUMBER USED FOR OUTPUT                    .
    // .                                                                   .
    // .  - - OUTPUT - -                                                   .
    // .        A(NWK)    = D AND L - FACTORS OF STIFFNESS MATRIX          .
    // .        V(NN)     = DISPLACEMENT VECTOR                            .
    // .                                                                   .
    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    //
    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    // .   THIS PROGRAM IS USED IN SINGLE PRECISION ARITHMETIC ON CRAY     .
    // .   EQUIPMENT AND DOUBLE PRECISION ARITHMETIC ON IBM MACHINES,      .
    // .   ENGINEERING WORKSTATIONS AND PCS. DEACTIVATE ABOVE LINE FOR     .
    // .   SINGLE PRECISION ARITHMETIC.                                    .
    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    //
    //
    //     PERFORM L*D*L(T) FACTORIZATION OF STIFFNESS MATRIX
    //

    int        N, KN, KL, KU, KH, K, IC, KLT, J, KI, ND, KK, L;
    double    C, B;

    if(KKK - 2 < 0)
        goto here_40;
    else if(KKK - 2 == 0)
        goto here_150;
    else if(KKK - 2 > 0)
        goto here_150;

here_40:
    for(N = 1; N < NN + 1; N++){
        KN    = MAXA[N];
        KL    = KN + 1;
        KU    = MAXA[N + 1] - 1;
        KH    = KU - KL;

        if(KH < 0) 
            goto here_110;
        else if(KH == 0)
            goto here_90;
        else if(KH > 0)
            goto here_50;

here_50:
        K    = N - KH; 
        IC    = 0;
        KLT    = KU;
        for(J = 1; J < KH + 1; J++){
            IC    = IC + 1;
            KLT    = KLT - 1;
            KI    = MAXA[K];
            ND    = MAXA[K + 1] - KI - 1;

            if(ND < 0)
                goto here_80;
            else if(ND == 0)
                goto here_80;
            else if(ND > 0)
                goto here_60;

here_60:
            KK    = MIN(IC,ND);
            C    = 0.0;

            for(L = 1; L < KK + 1; L++)
                C = C + A[KI + L] * A[KLT + L];

            A[KLT]    = A[KLT] - C;

here_80:
            K    = K + 1;
        }

here_90:
        K    = N;
        B    = 0.0;
        for(KK = KL; KK < KU + 1; KK++){
            K        = K - 1;
            KI        = MAXA[K];
            C        = A[KK] / A[KI];
            B        = B + C * A[KK];
            A[KK]    = C;
        }

        A[KN] = A[KN] - B;

here_110:
        if(A[KN] < 0)
            goto here_120;
        else if(A[KN] == 0)
            goto here_120;
        else if(A[KN] > 0)
            goto here_140;

here_120:
        cout<<" STOP - STIFFNESS MATRIX NOT POSITIVE DEFINITE ,NONPOSITIVE PIVOT FOR EQUATION ";
        cout<<N<<"\t"<<" PIVOT = "<<A[KN]<<endl;
        goto here_800;

here_140:
        continue;
    }
    goto here_900;

    //
    // REDUCE RIGHT-HAND-SIDE LOAD VECTOR
    //

here_150:
    for(N = 1; N < NN + 1; N++){
        KL    = MAXA[N] + 1;
        KU    = MAXA[N + 1] - 1;

        if(KU - KL < 0)
            goto here_180;
        else if(KU - KL == 0)
            goto here_160;
        else if(KU - KL > 0)
            goto here_160;

here_160:
        K    = N;
        C    = 0.0;

        for(KK = KL; KK < KU + 1; KK++){
            K    = K - 1;
            C    = C + A[KK] * V[K];
        }

        V[N] = V[N] - C;

here_180:
        continue;
    }

    //
    //  BACK-SUBSTITUTE
    //

    for(N = 1; N < NN + 1; N++){
        K        = MAXA[N];
        V[N]    = V[N] / A[K];
    }

    if(NN == 1) 
        goto here_900;

    N  = NN;

    for(L = 2; L < NN + 1; L++){
        KL    = MAXA[N] + 1;
        KU    = MAXA[N + 1] - 1;

        if(KU - KL < 0)
            goto here_230;
        else if(KU - KL == 0)
            goto here_210;
        else if(KU - KL > 0)
            goto here_210;

here_210:
        K    = N;

        for(KK = KL; KK < KU + 1; KK++){
            K        = K - 1;
            V[K]    = V[K] - A[KK] * V[N];
        }

here_230:
        N    = N - 1;
    }

    goto here_900;

here_800:
    exit(0);

here_900:
    return;
}