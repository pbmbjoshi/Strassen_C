/*
 * Application to demonstrate Strassens multiplication.
 * https://en.wikipedia.org/wiki/Strassen_algorithm
 *
 *
 * Author : Prasanna Joshi
 * bmbjoshi@gmail.com
 *
 * Please write to me for feedback/questions/improvements or general comment.
 * This application is built in gcc compiler (gcc version 7.3.0) in ubuntu 18.04.
 *
 * Next, C++ version of the strassens multiplication.
*/

#include <stdio.h>
#include <stdlib.h>
#include<time.h>

/*
 * Enum to identify portion of a matrix whose rows and columns are even numbers
 */
enum Esection_of_matrix{top_left, top_right, bottom_left, bottom_right};

/*
 * Matrix structure
 */
struct Matrix
{
    int m_nRow;
    int m_nColumn;
    int **MatrixArray;
};

// Function creates a matrix for valid nRow and valid nColumn values (validity is > 0)
struct Matrix CreateMatrix(const int nRow, const int nColumn);

// Releases matrix and sets the row and column to 0.
void ReleaseMatrix(struct Matrix* A);

// Function populates Matrix 'matrix' with random number for upperlimt of 'AllowedTill'
void GenerateMatrix(struct Matrix matrix, const int AllowedTill);

// Helper function to check whether or not given integer is a power of 2
int isPowerOfTwo(int n);

// Displays matrix in row and clumn fashion with a heading.
void DisplayMatrix(struct Matrix matrix, const char* heading);

// Returns a resultant matrix of addition of given matrices.
struct Matrix AddMatrices(const struct Matrix A, const struct Matrix B);

// Returns a resultant matrix of subracting 'A' from 'B'.
struct Matrix SubtractMatrices(const struct Matrix A, const struct Matrix B);

// Returns portion of defined by 'esection_of_matrix' from a given matrix 'A'
struct Matrix GetMatrixPortion(const struct Matrix A, const enum Esection_of_matrix esection_of_matrix);

// Returns Assembled matrix from given matrices.
struct Matrix AssembleMatrix(   const struct Matrix top_left,
                                const struct Matrix top_right,
                                const struct Matrix bottom_left,
                                const struct Matrix bottom_right);

// Function for calling StrassenMult.
void MatrixMult(const int nOrder, const int nLimitOfRandomNumber);

// Main function for Strassens multiplication. Carries out multiplication ONLY if orders of both matrices are same and they are square matrices.
struct Matrix StrassenMult(const struct Matrix A, const struct Matrix B);

int main()
{
    MatrixMult(2,10);
    MatrixMult(4,10);
    return 0;
}


void MatrixMult(const int nOrder, const int nLimitOfRandomNumber)
{
    struct Matrix A = CreateMatrix(nOrder, nOrder);
    GenerateMatrix(A,nLimitOfRandomNumber);
    DisplayMatrix(A,"Matrix A:");
    struct Matrix B = CreateMatrix(nOrder, nOrder);
    GenerateMatrix(B,nLimitOfRandomNumber);
    DisplayMatrix(B,"Matrix B:");
    printf("\n");
    struct Matrix CStrassMult = StrassenMult(A,B);
    DisplayMatrix(CStrassMult,"Strassen Multiplication of A and B :");
    printf("\n");
    ReleaseMatrix(&A);
    ReleaseMatrix(&B);
    ReleaseMatrix(&CStrassMult);
}

struct Matrix CreateMatrix(const int nRow, const int nColumn)
{
    struct Matrix matrix;
    if (nRow > 0 && nColumn > 0)
    {
        matrix.MatrixArray = (int**)malloc(sizeof(int*)*nRow);
        matrix.m_nRow = nRow;
        matrix.m_nColumn = nColumn;
        for (int i=0;i<nRow;++i)
        {
            matrix.MatrixArray[i] = (int*)malloc(sizeof(int)*nColumn);
            for (int j=0;j<nColumn;++j)
            {
                matrix.MatrixArray[i][j] = 0;
            }
        }
    }
    else {
        matrix.m_nRow = 0;
        matrix.m_nColumn = 0;
        matrix.MatrixArray = NULL;
    }
    return matrix;
}

void ReleaseMatrix(struct Matrix* A)
{
    if (A != NULL)
    {
        for (int i=0;i<A->m_nRow;++i)
        {
            free(A->MatrixArray[i]);
            A->MatrixArray[i] = NULL;
        }
        free(A->MatrixArray);
        A->MatrixArray = NULL;
    }
}

struct Matrix AddMatrices(const struct Matrix A, const struct Matrix B)
{
    struct Matrix C = CreateMatrix(0,0);
    if (A.m_nRow == B.m_nRow && A.m_nColumn == B.m_nColumn && A.m_nRow > 0 && A.m_nColumn > 0)
    {
        C = CreateMatrix(A.m_nRow, A.m_nColumn);
        for (int i=0;i<C.m_nRow;++i)
        {
            for (int j=0;j<C.m_nColumn;++j)
            {
                C.MatrixArray[i][j] = A.MatrixArray[i][j] + B.MatrixArray[i][j];
            }
        }
    }
    return C;
}

struct Matrix SubtractMatrices(const struct Matrix A, const struct Matrix B)
{
    struct Matrix C = CreateMatrix(0,0);
    if (A.m_nRow == B.m_nRow && A.m_nColumn == B.m_nColumn && A.m_nRow > 0 && A.m_nColumn > 0)
    {
        C = CreateMatrix(A.m_nRow, A.m_nColumn);
        for (int i=0;i<C.m_nRow;++i)
        {
            for (int j=0;j<C.m_nColumn;++j)
            {
                C.MatrixArray[i][j] = A.MatrixArray[i][j] - B.MatrixArray[i][j];
            }
        }
    }
    return C;
}


void DisplayMatrix(struct Matrix matrix,const char* heading)
{
    if (heading != NULL)
    {
        printf("%s\n",heading);
    }
    for (int i=0;i<matrix.m_nRow;++i)
    {
        for (int j=0;j<matrix.m_nColumn;++j)
        {
            printf("%d\t",matrix.MatrixArray[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void GenerateMatrix(struct Matrix matrix, const int AllowedTill)
{
    srand(time(0));
    for (int i=0;i<matrix.m_nRow;++i)
    {
        for (int j=0;j<matrix.m_nColumn;++j)
        {
            matrix.MatrixArray[i][j]= rand()%AllowedTill;
        }
    }
}

struct Matrix StrassenMult(const struct Matrix A, const struct Matrix B)
{
    struct Matrix C = CreateMatrix(0,0);
    if (A.m_nRow == A.m_nColumn && A.m_nRow == B.m_nRow && B.m_nRow == B.m_nColumn && isPowerOfTwo(B.m_nColumn))
    {
        int nOrder = B.m_nColumn;

        if (nOrder == 2)
        {
            C = CreateMatrix(nOrder,nOrder);
            /*
                M1= (A11+A22) (B11+B22)
                M2= (A21+A22)B11
                M3=A11(B12−B22)
                M4=A22(B21−B11)
                M5= (A11+A12)B22
                M6= (A21−A11)(B11+B12)
                M7= (A12−A22)(B21+B22)

                C11=M1+M4−M5+M7
                C12=M3+M5
                C21=M2+M4
                C22=M1−M2+M3 +M6
             */

            // Reached leaf. Apply Strassens operation directly on values. indices start with 0, hence A11 in mathematical matrix is equal to A[0][0] of C-matrix.
            int M1 = (A.MatrixArray[0][0]+A.MatrixArray[1][1])*(B.MatrixArray[0][0]+B.MatrixArray[1][1]);
            int M2 = (A.MatrixArray[1][0]+A.MatrixArray[1][1])*B.MatrixArray[0][0];
            int M3 = A.MatrixArray[0][0]*(B.MatrixArray[0][1]-B.MatrixArray[1][1]);
            int M4 = A.MatrixArray[1][1]*(B.MatrixArray[1][0]-B.MatrixArray[0][0]);
            int M5 = (A.MatrixArray[0][0]+A.MatrixArray[0][1])*B.MatrixArray[1][1];
            int M6 = (A.MatrixArray[1][0]-A.MatrixArray[0][0])*(B.MatrixArray[0][0]+B.MatrixArray[0][1]);
            int M7 = (A.MatrixArray[0][1]-A.MatrixArray[1][1])*(B.MatrixArray[1][0]+B.MatrixArray[1][1]);

            C.MatrixArray[0][0] = M1+M4-M5+M7;
            C.MatrixArray[0][1] = M3+M5;
            C.MatrixArray[1][0] = M2+M4;
            C.MatrixArray[1][1] = M1-M2+M3+M6;
        }
        else
        {
            /*
                M1= (A11+A22) (B11+B22)
                M2= (A21+A22)B11
                M3=A11(B12−B22)
                M4=A22(B21−B11)
                M5= (A11+A12)B22
                M6= (A21−A11)(B11+B12)
                M7= (A12−A22)(B21+B22)

                C11=M1+M4−M5+M7
                C12=M3+M5
                C21=M2+M4
                C22=M1−M2+M3 +M6
             */
            // Divide and conquer. Apply Strassens operation directly on divided matrices.

            struct Matrix A11 = GetMatrixPortion(A,top_left);
            struct Matrix A12 = GetMatrixPortion(A,top_right);
            struct Matrix A21 = GetMatrixPortion(A,bottom_left);
            struct Matrix A22 = GetMatrixPortion(A,bottom_right);

            struct Matrix B11 = GetMatrixPortion(B,top_left);
            struct Matrix B12 = GetMatrixPortion(B,top_right);
            struct Matrix B21 = GetMatrixPortion(B,bottom_left);
            struct Matrix B22 = GetMatrixPortion(B,bottom_right);

            struct Matrix A11_PLUS_A22 = AddMatrices(A11,A22);
            struct Matrix B11_PLUS_B22 = AddMatrices(B11,B22);
            struct Matrix A21_PLUS_A22 = AddMatrices(A21,A22);
            struct Matrix B12_MINUS_B22 = SubtractMatrices(B12,B22);
            struct Matrix B21_MINUS_B11 = SubtractMatrices(B21,B11);
            struct Matrix A11_PLUS_A12 = AddMatrices(A11,A12);
            struct Matrix A21_MINUS_A11 = SubtractMatrices(A21,A11);
            struct Matrix B11_PLUS_B12 = AddMatrices(B11,B12);
            struct Matrix A12_MINUS_A22 = SubtractMatrices(A12,A22);
            struct Matrix B21_PLUS_B22 = AddMatrices(B21,B22);

            struct Matrix M1 = StrassenMult(A11_PLUS_A22,B11_PLUS_B22);
            struct Matrix M2 = StrassenMult(A21_PLUS_A22,B11);
            struct Matrix M3 = StrassenMult(A11,B12_MINUS_B22);
            struct Matrix M4 = StrassenMult(A22,B21_MINUS_B11);
            struct Matrix M5 = StrassenMult(A11_PLUS_A12,B22);
            struct Matrix M6 = StrassenMult(A21_MINUS_A11,B11_PLUS_B12);
            struct Matrix M7 = StrassenMult(A12_MINUS_A22,B21_PLUS_B22);

            struct Matrix M1_PLUS_M4 = AddMatrices(M1,M4);
            struct Matrix M1_PLUS_M4_MINUS_M5 = SubtractMatrices(M1_PLUS_M4,M5);
            struct Matrix M1_MINUS_M2 = SubtractMatrices(M1,M2);
            struct Matrix M1_MINUS_M2_PLUS_M3 = AddMatrices(M1_MINUS_M2,M3);

            struct Matrix C11 = AddMatrices(M1_PLUS_M4_MINUS_M5,M7);
            struct Matrix C12 = AddMatrices(M3,M5);
            struct Matrix C21 = AddMatrices(M2,M4);
            struct Matrix C22 = AddMatrices(M1_MINUS_M2_PLUS_M3,M6);

            C = AssembleMatrix(C11,C12,C21,C22);

            // free the memory by calling release
            ReleaseMatrix(&A11);
            ReleaseMatrix(&A12);
            ReleaseMatrix(&A21);
            ReleaseMatrix(&A22);

            ReleaseMatrix(&B11);
            ReleaseMatrix(&B12);
            ReleaseMatrix(&B21);
            ReleaseMatrix(&B22);

            ReleaseMatrix(&A11_PLUS_A22);
            ReleaseMatrix(&B11_PLUS_B22);
            ReleaseMatrix(&A21_PLUS_A22);
            ReleaseMatrix(&B12_MINUS_B22);
            ReleaseMatrix(&B21_MINUS_B11);
            ReleaseMatrix(&A11_PLUS_A12);
            ReleaseMatrix(&A21_MINUS_A11);
            ReleaseMatrix(&B11_PLUS_B12);
            ReleaseMatrix(&A12_MINUS_A22);
            ReleaseMatrix(&B21_PLUS_B22);

            ReleaseMatrix(&M1);
            ReleaseMatrix(&M2);
            ReleaseMatrix(&M3);
            ReleaseMatrix(&M4);
            ReleaseMatrix(&M5);
            ReleaseMatrix(&M6);
            ReleaseMatrix(&M7);

            ReleaseMatrix(&M1_PLUS_M4);
            ReleaseMatrix(&M1_PLUS_M4_MINUS_M5);
            ReleaseMatrix(&M1_MINUS_M2);
            ReleaseMatrix(&M1_MINUS_M2_PLUS_M3);

            ReleaseMatrix(&C11);
            ReleaseMatrix(&C12);
            ReleaseMatrix(&C21);
            ReleaseMatrix(&C22);
        }
    }
    return C;
}

struct Matrix GetMatrixPortion(const struct Matrix matrix, const enum Esection_of_matrix esection_of_matrix)
{
    struct Matrix portionMatrix = CreateMatrix(0,0);
    if (matrix.m_nRow % 2 == 0 && matrix.m_nColumn % 2 == 0)
    {
        int nRow = matrix.m_nRow / 2;
        int nColumn = matrix.m_nColumn / 2;
        portionMatrix = CreateMatrix(nRow,nColumn);

        int AddRowIndex = 0;
        int AddColumnIndex = 0;

        switch (esection_of_matrix)
        {
            case top_left :
                AddRowIndex = 0;
                AddColumnIndex = 0;
            break;

            case top_right :
                AddRowIndex = 0;
                AddColumnIndex = nColumn;
            break;

            case bottom_left:
                AddRowIndex = nRow;
                AddColumnIndex = 0;
            break;

            case bottom_right:
                AddRowIndex = nRow;
                AddColumnIndex = nColumn;
            break;
        }

        for (int i=0;i<nRow;++i)
        {
            for (int j=0;j<nColumn;++j)
            {
                portionMatrix.MatrixArray[i][j]= matrix.MatrixArray[i+AddRowIndex][j+AddColumnIndex];
            }
        }

    }
    return portionMatrix;
}

struct Matrix AssembleMatrix(   const struct Matrix top_left,
                                const struct Matrix top_right,
                                const struct Matrix bottom_left,
                                const struct Matrix bottom_right)
{
    struct Matrix assembled_matrix = CreateMatrix(0,0);

    if (    top_left.m_nRow == top_right.m_nRow &&
            top_right.m_nRow == bottom_left.m_nRow &&
            bottom_left.m_nRow == bottom_right.m_nRow &&
            bottom_right.m_nRow > 0 &&
            top_left.m_nColumn == top_right.m_nColumn &&
            top_right.m_nColumn == bottom_left.m_nColumn &&
            bottom_left.m_nColumn == bottom_right.m_nColumn &&
            bottom_right.m_nColumn > 0)
    {
        int nPortionRow = top_left.m_nRow;
        int nPortionColumn = top_left.m_nColumn;
        int nRow = nPortionRow+nPortionRow;
        int nColumn = nPortionColumn+nPortionColumn;

        assembled_matrix = CreateMatrix(nRow,nColumn);

        for (int i=0;i<nPortionRow;++i)
        {
            for (int j=0;j<nPortionColumn;++j)
            {
                assembled_matrix.MatrixArray[i][j] = top_left.MatrixArray[i][j];
                assembled_matrix.MatrixArray[i][j+nPortionColumn] = top_right.MatrixArray[i][j];
                assembled_matrix.MatrixArray[i+nPortionRow][j] = bottom_left.MatrixArray[i][j];
                assembled_matrix.MatrixArray[i+nPortionRow][j+nPortionColumn] = bottom_right.MatrixArray[i][j];
            }
        }
    }
    return assembled_matrix;
}

int isPowerOfTwo(int n)
{
    if (n == 0)
        return 0;
    while (n != 1)
    {
        if (n%2 != 0)
            return 0;
        n = n/2;
    }
    return 1;
}
