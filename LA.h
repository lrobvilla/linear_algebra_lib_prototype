#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <string.h>
#include <sstream>
#include <iomanip>
using namespace std;

// struct para trabajar con vectores **
struct VectorShell
{
    static constexpr size_t tamano {3};
    double vect[tamano];

    void setVector(double newVect[tamano]){
        for (size_t i = 0; i < tamano; i++)
        {
            vect[i] = newVect[i];
        }
    }
};

//devuelve el vector e_j de la base canonica de R^n
VectorShell canonical_e_j(size_t j){
    VectorShell e_j;
    for (size_t i = 0; i < e_j.tamano; i++)
    {
        if (i==j)
        {
            e_j.vect[i] = 1;
        } else {
            e_j.vect[i] = 0;
        }
    }
    return e_j;
}

//producto interno de vectores
double operator*(VectorShell vector1, VectorShell vector2){
    size_t tamano = vector1.tamano;
    double innerP = 0;
    for (size_t i = 0; i < tamano; i++)
    {
        innerP += (vector1.vect[i])*(vector2.vect[i]);
    }
    return innerP;
}

//suma de vectores por coordenadas
VectorShell operator+(VectorShell v1, VectorShell v2){
    VectorShell suma;
    for (size_t i = 0; i < v1.tamano; i++)
    {
        suma.vect[i] = 0;
        suma.vect[i] = v1.vect[i] + v2.vect[i];
    }
    return suma;
}

//resta de vectores por coordenadas
VectorShell operator-(VectorShell v1, VectorShell v2){
    VectorShell resta;
    for (size_t i = 0; i < v1.tamano; i++)
    {
        resta.vect[i] = 0;
        resta.vect[i] = v1.vect[i] - v2.vect[i];
    }
    return resta;
}

//multiplicacion de escalar para vectores
VectorShell operator*(double escalar, VectorShell v){
    VectorShell productoEscalar;
    for (size_t i = 0; i < v.tamano; i++)
    {
        productoEscalar.vect[i] = double(escalar)*double(v.vect[i]);
    }
    return productoEscalar;
}

//obtener norma2 de vectores
double norma2(VectorShell v){
    double norm2 = v*v;
    norm2 = double(pow(norm2, 1/double(2)));
    return norm2;
}

//devuelve el vector pasado como parametro normalizado
VectorShell normalizeVect(VectorShell v){
    double invNorm = double(1)/norma2(v);
    return invNorm*v;
}

//funcion de signo para un double
int sign(double x){
    double signo = 0;
    if (abs(x)<pow(10,-16))
    {
        signo = 1;
    } else {
        signo = double(x)/double(abs(x));
    }
    return signo;
}

//struct para trabajar con matrices**
struct MatrixShell
{
    static constexpr size_t squareDim {3};
    double matrix[squareDim][squareDim]; //cambiar esto luego

    void changeCol(VectorShell v, size_t col){
        for (size_t i = 0; i < v.tamano; i++)
        {
            matrix[i][col] = v.vect[i];
        }
    }

    void changeFila(VectorShell v, size_t fila){
        for (size_t i = 0; i < v.tamano; i++)
        {
            matrix[fila][i] = v.vect[i];
        }
    }

    void changeEntry(size_t fila, size_t col, double newEntry){
        matrix[fila][col] = newEntry;
    }
};

//struct para almacenar dos matrices, idealmente una descomposicion QR
struct QRShell
{
    MatrixShell Q;
    MatrixShell R;
    
    void setQ(MatrixShell newQ){
        size_t squareDim = newQ.squareDim;
        size_t dim = pow(squareDim, 2);
        for(size_t i = 0; i < dim; i++)
        {
            size_t col = i % squareDim;
            size_t fila = (i - col) / squareDim;
            Q.changeEntry(fila, col, newQ.matrix[fila][col]);
        }
    }

    void setR(MatrixShell newQ){
        size_t squareDim = newQ.squareDim;
        size_t dim = pow(squareDim, 2);
        for(size_t i = 0; i < dim; i++)
        {
            size_t col = i % squareDim;
            size_t fila = (i - col) / squareDim;
            R.changeEntry(fila, col, newQ.matrix[fila][col]);
        }
    }
};


//llena una matriz con numeros aleatorios, o bien la identidad para seed = 0, o puros ceros para seed = -1
MatrixShell defineMatrix(size_t seed){
    MatrixShell newMatrix;
    size_t squareDim = newMatrix.squareDim;
    size_t dim = pow(squareDim, 2);
    if (seed != -1)
    {
        srand(seed);
        for(size_t i = 0; i < dim; i++)
        {
            size_t col = i % squareDim;
            size_t fila = (i - col) / squareDim;
            if (seed == 0) // caso para la identidad
            {
                if (fila == col)
                {
                    newMatrix.matrix[fila][col] = 1;
                } else {
                    newMatrix.matrix[fila][col] = 0;
                }
            } else {
                newMatrix.matrix[fila][col] = double(rand())/double(RAND_MAX); 
            }          
            // newMatrix.matrix[fila][col] = double(rand());  //para pruebas con valores enteros    
        }   
    } else {
       for(size_t i = 0; i < dim; i++)
        {
            size_t col = i % squareDim;
            size_t fila = (i - col) / squareDim;
            newMatrix.matrix[fila][col] = 0;
        } 
    }
    return newMatrix;
}

//d
MatrixShell truncatedAsignMatrix(MatrixShell asign, MatrixShell toAsign, size_t fila, size_t col){
    MatrixShell newMatrix = asign;
    size_t squareDim = newMatrix.squareDim;
    for (size_t k = fila; k < squareDim; k++)
    {
        for (size_t l = col; l < squareDim; l++)
        {
            newMatrix.changeEntry(k, l, toAsign.matrix[k][l]);
        }
    }
    return newMatrix;
}

//devueve la columna deseada de la matriz pasada como parametro
VectorShell obtainCol(MatrixShell matrixShell, size_t col){
    VectorShell column;
    for (size_t i = 0; i < column.tamano; i++)
    {
        column.vect[i] = matrixShell.matrix[i][col];
    }
    return column;
}

//devuelve la fila deseada de la matriz
VectorShell obtainFila(MatrixShell matrixShell, size_t fil){
    VectorShell fila;
    for (size_t i = 0; i < fila.tamano; i++)
    {
        fila.vect[i] = matrixShell.matrix[fil][i];
    }
    return fila;
}

//funcion para imprimir un vector
void printVector(VectorShell vector){
    size_t tamano = vector.tamano;
    for (size_t i = 0; i < tamano; i++)
    {
        cout << setprecision(9) << vector.vect[i] << " ";
    }
    cout << endl;
}

//funcion para obtener la cantidad de digitos de una entra de una matriz
size_t entrySize(MatrixShell matrixShell, size_t fila, size_t col){
    ostringstream streamObj; //solamente convertir a string no preserva los digitos, con esto se arregla    
    streamObj << setprecision(9) << matrixShell.matrix[fila][col];
    string entry = streamObj.str();
    size_t size = entry.size();
    return size;
}

//esta funcion es para obtener la cantidad de espacios adecuada para desplegar las entradas de cada columna alineadas
int spaceMaster(MatrixShell matrixShell, size_t col){
    size_t squareDim = matrixShell.squareDim;
    size_t maxDigits = 0;
    for(size_t i = 0; i < squareDim; i++)
    {
        size_t fila = i;
        size_t size = entrySize(matrixShell, fila, col);
        if (size>maxDigits)
        {
            maxDigits = size;
        }
    }
    return maxDigits;
}

//funcion para imprimir una matriz
void printMatrix(MatrixShell matrixShell){
    size_t squareDim = matrixShell.squareDim;
    size_t dim = pow(squareDim,2);
    for(size_t i = 0; i < dim; i++)
    {
        size_t col = i % squareDim;
        size_t fila = (i - col) / squareDim;
        if (col % squareDim == squareDim-1)
        {
            cout << setprecision(9) << matrixShell.matrix[fila][col] << endl;
        } else {
            cout << setprecision(9) << matrixShell.matrix[fila][col] << " ";
            int blankSpaces = spaceMaster(matrixShell, col) - entrySize(matrixShell, fila, col);
            // int blankSpaces = 16 - entrySize(matrixShell, fila, col);
            for (int j = 0; j < blankSpaces; j++)
            {
                cout << " ";
            }  
        }
    }
    cout << endl;
}

//para multiplicacion de matrices de iguales dimensiones
MatrixShell operator*(MatrixShell matrix1, MatrixShell matrix2){
    size_t squareDim = matrix1.squareDim;
    size_t dim = squareDim*squareDim;
    MatrixShell mult;
    for(size_t i = 0; i < dim; i++)
    {
        size_t col = i % squareDim;
        size_t fila = (i - col) / squareDim;
        mult.matrix[fila][col] = 0;
        for(size_t k = 0; k < squareDim; k++)    
        {    
            mult.matrix[fila][col] += (matrix1.matrix[fila][k])*(matrix2.matrix[k][col]);    
        }       
    }
    return mult;
}

//multiplicacion de escalar por matriz
MatrixShell operator*(double c, MatrixShell matrix){
    size_t squareDim = matrix.squareDim;
    size_t dim = squareDim*squareDim;
    MatrixShell mult;
    for(size_t i = 0; i < dim; i++)
    {
        size_t col = i % squareDim;
        size_t fila = (i - col) / squareDim;
        mult.matrix[fila][col] = 0;
        mult.matrix[fila][col] = c*(matrix.matrix[fila][col]);
    }
    return mult;
}

//para suma de matrices
MatrixShell operator+(MatrixShell matrix1, MatrixShell matrix2){
    size_t squareDim = matrix1.squareDim;
    size_t dim = squareDim*squareDim;
    MatrixShell mult;
    for(size_t i = 0; i < dim; i++)
    {
        size_t col = i % squareDim;
        size_t fila = (i - col) / squareDim;
        mult.matrix[fila][col] = 0;
        mult.matrix[fila][col] = matrix1.matrix[fila][col] + matrix2.matrix[fila][col];       
    }
    return mult;
}

//resta de matrices
MatrixShell operator-(MatrixShell matrix1, MatrixShell matrix2){
    size_t squareDim = matrix1.squareDim;
    size_t dim = squareDim*squareDim;
    MatrixShell mult = defineMatrix(-1);
    for(size_t i = 0; i < dim; i++)
    {
        size_t col = i % squareDim;
        size_t fila = (i - col) / squareDim;
        mult.matrix[fila][col] = matrix1.matrix[fila][col] - matrix2.matrix[fila][col];          
    }
    return mult;
}

//para obtener matriz producto de un vector con el transpuesto de otro vector, v1*v2^T, con v1 un nx1 y v2 un 1xn
MatrixShell operator%(VectorShell v1, VectorShell v2){
    MatrixShell newMatrix = defineMatrix(-1);
    size_t squareDim = newMatrix.squareDim;
    size_t dim = squareDim*squareDim;
    for(size_t i = 0; i < dim; i++)
    {
        size_t col = i % squareDim;
        size_t fila = (i - col) / squareDim;
        newMatrix.matrix[fila][col] = v1.vect[fila]*v2.vect[col];
    }
    return newMatrix;
}

//producto de matriz con vector, para utilizar la matrix como funcional lineal, devuelve mat*v
VectorShell operator*(MatrixShell mat,VectorShell v){
    VectorShell newVector;
    VectorShell fila;
    size_t tamano = v.tamano; //tal deberia inicializar el vector, para luego...
    for (size_t i = 0; i < tamano; i++)
    {
        fila.setVector(mat.matrix[i]);
        newVector.vect[i] = fila*v;
    }
    return newVector;
}

//devuelve un vector con ceros en las primeras (j-tamano) entradas  y el resto de entradas las del vector pasado como parametro
VectorShell zeros_to_jth_entry(VectorShell vect, size_t j){
    VectorShell newVect;
    for (size_t i = 0; i < vect.tamano; i++)
    {
        newVect.vect[i] = 0;
        if (i>=j)
        {
            newVect.vect[i] = vect.vect[i];
        }
    }
    return newVect;
}

//devuelve la transpuesta de una matriz
MatrixShell transposeMatrix(MatrixShell mat){
    MatrixShell newMatrix = defineMatrix(-1);
    size_t squareDim = newMatrix.squareDim;
    size_t dim = squareDim*squareDim;
    for(size_t i = 0; i < dim; i++)
    {
        size_t col = i % squareDim;
        size_t fila = (i - col) / squareDim;
        newMatrix.matrix[fila][col] = mat.matrix[col][fila];
    }
    return newMatrix;
}

//funcion para aplicar un reflector de householder con respecto a un vector v
VectorShell HouseholderReflector(VectorShell v, VectorShell x){
    VectorShell reflection;
    for (size_t i = 0; i < reflection.tamano; i++)
    {
        reflection.vect[i] = 0;
    }
    reflection = x - (2/(v*v))*(v%v)*x;
    return reflection;
}

//descomposicion QR con GramSchmidt
QRShell QRDecompositionGramSchmidt(MatrixShell mat){
    QRShell matQR;
    MatrixShell matQ = defineMatrix(-1);
    MatrixShell matR = defineMatrix(-1);
    VectorShell u_n;
    VectorShell e_n;
    for (size_t i = 0; i < mat.squareDim; i++)
    {
        VectorShell a_n = obtainCol(mat, i);
        u_n = a_n;
        for (size_t j = 0; j < i; j++)
        {
            u_n = u_n - (a_n*obtainCol(matQ,j))*obtainCol(matQ,j);
        }
        e_n = normalizeVect(u_n);
        matQ.changeCol(e_n, i);
        for (size_t k = 0; k <= i; k++)
        {
            VectorShell e_k = obtainCol(matQ,k);
            matR.changeEntry(k, i, a_n*e_k);
        }
    }
    matQR.Q = matQ;
    matQR.R = matR;
    return matQR;
}

QRShell QRDecompositionHouseholder(MatrixShell mat){
    QRShell matQR;
    MatrixShell matR = mat;
    MatrixShell id = defineMatrix(0);
    MatrixShell matQ = id;
    VectorShell x_j;
    VectorShell v_j;
    for (size_t j = 0; j < mat.squareDim-1; j++)
    {
        x_j = obtainCol(matR, j);
        x_j = zeros_to_jth_entry(x_j, j);
        // cout << "x_" << j << "=" << endl;
        // printVector(x_j);
        v_j = x_j + sign(x_j.vect[j])*norma2(x_j)*canonical_e_j(j);
        // cout << "u_" << j << "=" << endl;
        // printVector(v_j);
        v_j = normalizeVect(v_j);
        // cout << "v_" << j << "=" << endl;
        // printVector(v_j);
        // cout << "-------" << endl;
        MatrixShell H = id;
        H = truncatedAsignMatrix(id, H - 2*(v_j%v_j), j, j); //aqui esta el error
        // cout << "Q_" << j << "=" << endl;
        // printMatrix(H);
        matR = H*matR; //*****
        // cout << "QA_" << j << "=" << endl;
        // printMatrix(matR);
        // cout << "-------" << endl;
        matQ = matQ*H;
        // cout << "Q_" << j << "=" << endl;
        // printMatrix(matQ);
    }
    matQR.Q = matQ;
    matQR.R = matR;
    return matQR;
}