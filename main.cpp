#include "LA.h"
#include <typeinfo>
#include <chrono>

int main(){
    //verificacion de los algoritmos, cambiar dim a 3
    MatrixShell A;
    // A = defineMatrix(93);
    A.matrix[0][0] = 0;
    A.matrix[0][1] = -51;
    A.matrix[0][2] = 4;
    A.matrix[1][0] = 0;
    A.matrix[1][1] = 51;
    A.matrix[1][2] = -68;
    A.matrix[2][0] = 0;
    A.matrix[2][1] = 7;
    A.matrix[2][2] = 9;

    // cout << "A = " << endl;
    // printMatrix(A);
    // QRShell qrA_GS = QRDecompositionGramSchmidt(A);
    // MatrixShell Q = qrA_GS.Q;
    // MatrixShell R = qrA_GS.R;
    // cout << "Q = " << endl;
    // printMatrix(Q);
    // cout << "R = " << endl;
    // printMatrix(R);
    // MatrixShell A_prime = Q*R;
    // cout << "Q*R = " << endl;
    // printMatrix(A_prime);
    // MatrixShell Q_t = transposeMatrix(Q);
    // cout << "Q*Q_t = " << endl;
    // printMatrix(Q*Q_t);

    cout << "-------" << endl;

    cout << "A = " << endl;
    printMatrix(A);
    QRShell qrA_h = QRDecompositionHouseholder(A);
    MatrixShell Q = qrA_h.Q;
    MatrixShell R = qrA_h.R;
    cout << "Q = " << endl;
    printMatrix(Q);
    cout << "R = " << endl;
    printMatrix(R);
    MatrixShell A_prime = Q*R;
    cout << "Q*R = " << endl;
    printMatrix(A_prime);
    MatrixShell Q_t = transposeMatrix(Q);
    cout << "Q*Q_t = " << endl;
    printMatrix(Q*Q_t);


    // cambiar dim a 50

    //comparacion de la duracion Gram-Schmidt

    // MatrixShell A = defineMatrix(93);
    // // cout << "A = " << endl;
    // // printMatrix(A);
    // auto start_GS = chrono::high_resolution_clock::now();
    // QRShell qrA_GS = QRDecompositionGramSchmidt(A);
    // auto stop_GS = chrono::high_resolution_clock::now();
    // std::chrono::duration<double, std::milli> interval_GS = stop_GS - start_GS;
    // cout << "Duracion QR Gram-Schmidt:" << interval_GS.count() << endl;
    
    // // comparacion de la duracion Householder

    // // cout << "A = " << endl;
    // // printMatrix(A);
    // auto start_H = chrono::high_resolution_clock::now();
    // QRShell qrA_H = QRDecompositionHouseholder(A);
    // auto stop_H = chrono::high_resolution_clock::now();
    // std::chrono::duration<double, std::milli> interval_H = stop_H - start_H;
    // cout << "Duracion QR Householder:" << interval_H.count() << endl;

    // se puede probar que la cantidad de operaciones de punto flotante es 2*n*m^2 - (2/3)*n^3, para matrices cuadradas esto es (4/3)n^3

    return 0;
}

