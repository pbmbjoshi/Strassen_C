#ifndef STRASSENS_MATRIX_MULTIPLICATION_H
#define STRASSENS_MATRIX_MULTIPLICATION_H

#include <vector>
#include <iostream>
#include <ostream>
#include <algorithm>
#include <cstdlib>
#include <ctime>

void strassens_matrix_multiplication();

template <typename T>
class CMatrix
{
public:
    CMatrix();
    CMatrix(const std::size_t& row, const std::size_t& column);
    CMatrix(const std::size_t& order); // square matrix
    CMatrix(const CMatrix& rhs);
    CMatrix& operator = (const CMatrix& rhs);

    CMatrix(CMatrix&& rhs) noexcept;
    CMatrix& operator = (CMatrix&& rhs) noexcept;

    constexpr std::size_t row() const noexcept;
    constexpr std::size_t column() const noexcept;

    T& operator()(const std::size_t i, const std::size_t j);
    const T& operator()(const std::size_t i, const std::size_t j) const;

    CMatrix operator +(const CMatrix& rhs) const;
    CMatrix operator -(const CMatrix& rhs) const;
    CMatrix operator *(const CMatrix& rhs) const;

    friend constexpr std::ostream& operator << (std::ostream& o, const CMatrix& rhs) noexcept
    {
        std::for_each(rhs.matrix.begin(),rhs.matrix.end(),[](const typename CMatrix<T>::OneDArray& aRow)
        {
            for_each(aRow.begin(),aRow.end(),[](const T& anElement) {
                std::cout << anElement << "  ";
            });
            std::cout << std::endl;
        });

        return o;
    }

private:
    using OneDArray = std::vector<T>;
    using TwoDArray = std::vector<OneDArray>;

    TwoDArray matrix;
};

template <typename T>
CMatrix<T>::CMatrix()
{

}

template <typename T>
CMatrix<T>::CMatrix(const std::size_t& row, const std::size_t& column)
{
    if (row > 0 && column > 0)
    {
        OneDArray aRow(column,0);
        matrix.insert(matrix.begin(),row,aRow);
    }
}

template <typename T>
CMatrix<T>::CMatrix(const std::size_t& order) : CMatrix(order,order) // square matrix
{

}

template <typename T>
CMatrix<T>::CMatrix(const CMatrix& rhs)
{
    *this = rhs;
}

template <typename T>
CMatrix<T>& CMatrix<T>::operator = (const CMatrix& rhs)
{
    if (this != &rhs)
    {
        matrix.clear();
        matrix = rhs.matrix;
    }
    return *this;
}

template <typename T>
constexpr std::size_t CMatrix<T>::row() const noexcept
{
    return matrix.size();
}

template <typename T>
constexpr std::size_t CMatrix<T>::column() const noexcept
{
    return matrix.size() == 0 ? 0 : matrix.at(0).size();
}

template <typename T>
CMatrix<T>::CMatrix(CMatrix&& rhs) noexcept
{
    *this = std::move(rhs);
}

template <typename T>
CMatrix<T>& CMatrix<T>::operator = (CMatrix&& rhs) noexcept
{
    if (this != &rhs)
    {
        matrix = std::move(rhs.matrix);
    }
    return *this;
}

template <typename T>
T& CMatrix<T>::operator()(const std::size_t i, const std::size_t j)
{
    if (i >= 0 && i < row() && j >= 0 && j < column())
    {
        OneDArray& aRow = matrix[i];
        T& anElement = aRow[j];
        return anElement;
    }
    throw "out of index exception";
}

template <typename T>
const T& CMatrix<T>::operator()(const std::size_t i, const std::size_t j) const
{
    if (i >= 0 && i < row() && j >= 0 && j < column())
    {
        const OneDArray& aRow = matrix[i];
        const T& anElement = aRow[j];
        return anElement;
    }
    throw "out of index exception";
}

template<typename T>
CMatrix<T> CMatrix<T>::operator +(const CMatrix& rhs) const
{
    std::size_t m = row();
    std::size_t n = column();


    if (m == rhs.row() && n == rhs.column())
    {
        CMatrix Temp(m,n);
        for (std::size_t i = 0; i < m; ++i)
        {
            for (std::size_t j = 0; j < n; ++j)
            {
                Temp(i,j) = this->operator()(i,j) + rhs(i,j);
            }
        }
        return Temp;
    }
    throw "order mismatch for operation";
}

template<typename T>
CMatrix<T> CMatrix<T>::operator -(const CMatrix& rhs) const
{
    std::size_t m = row();
    std::size_t n = column();


    if (m == rhs.row() && n == rhs.column())
    {
        CMatrix Temp(m,n);
        for (std::size_t i = 0; i < m; ++i)
        {
            for (std::size_t j = 0; j < n; ++j)
            {
                Temp(i,j) = this->operator()(i,j) - rhs(i,j);
            }
        }
        return Temp;
    }
    throw "order mismatch for operation";
}

template<typename T>
CMatrix<T> CMatrix<T>::operator *(const CMatrix& rhs) const
{
    std::size_t m = row();
    std::size_t n = column();
    std::size_t q = rhs.column();


    if (m == q && n == rhs.row())
    {
        CMatrix Temp(m,q); // in mxn and pxq matrix multiplication, resultaant is mxq matrix
        for (std::size_t i = 0; i < m; ++i)
        {
            for (std::size_t j = 0; j < q; ++j)
            {
                Temp(i,j) = static_cast<T>(0);
                for (std::size_t k = 0; k < n; ++k)
                {
                    Temp(i,j) = Temp(i,j) + this->operator()(i,k) * rhs(k,j);
                }
            }
        }
        return Temp;
    }
    throw "order mismatch for operation";
}

template<typename T>
CMatrix<T> GetMatrix(const std::size_t& row, const std::size_t& column, const int& MaxOfRandomValue)
{
    CMatrix<T> matrix(row,column);

    if (row > 0 && column > 0)
    {

        for (std::size_t i = 0; i < row; ++i)
        {
            for (std::size_t j = 0; j < column; ++j)
            {
                matrix(i,j) = static_cast<T>(rand()%MaxOfRandomValue);
            }
        }
    }
    return  matrix;
}

template <typename T>
class StrassenMult
{
public:
    static CMatrix<T> StrasMult(const CMatrix<T>& A, const CMatrix<T>& B);
private:
    static CMatrix<T> d_and_q_mult(const CMatrix<T>& A, const CMatrix<T>& B); // divide and conquer
    static bool partition_matrix(const CMatrix<T>& A,CMatrix<T>& A00, CMatrix<T>& A01, CMatrix<T>& A10, CMatrix<T>& A11);
    static bool merge_matrix(CMatrix<T>& A,const CMatrix<T>& A00, const CMatrix<T>& A01, const CMatrix<T>& A10, const CMatrix<T>& A11);
    static bool isPowerOfTwo(std::size_t n);
};

template <typename T>
CMatrix<T> StrassenMult<T>::StrasMult(const CMatrix<T>& A, const CMatrix<T>& B)
{
    if (A.row() == A.column() && A.row() == B.column() && B.row() == B.column())
    {
        std::size_t n = A.row();
        if (isPowerOfTwo(n) )
        {
            CMatrix<T> Result = d_and_q_mult(A, B);
            return Result;
        }
        throw "square matrices required";
    }
    throw "order mismatch for operation";
}

template <typename T>
CMatrix<T> StrassenMult<T>::d_and_q_mult(const CMatrix<T>& A, const CMatrix<T>& B)
{
    std::size_t n = A.row();
    CMatrix<T> C(n);
    if ( n == 2 )
    {
        T M1 = (A(0,0) + A(1,1))*(B(0,0)+B(1,1));
        T M2 = (A(1,0) + A(1,1))*B(0,0);
        T M3 = A(0,0)*(B(0,1)-B(1,1));
        T M4 = A(1,1)*(B(1,0)-B(0,0));
        T M5 = (A(0,0)+A(0,1))*B(1,1);
        T M6 = (A(1,0) - A(0,0))*(B(0,0)+B(0,1));
        T M7 = (A(0,1) - A(1,1))*(B(1,0)+B(1,1));

        // Assemble
        C(0,0) = M1+M4-M5+M7;
        C(0,1) = M3+M5;
        C(1,0) = M2+M4;
        C(1,1) = M1-M2+M3+M6;
    }
    else
    {
        std::size_t order = n/2;
        CMatrix<T> A00(order), A01(order), A10(order), A11(order), B00(order), B01(order), B10(order), B11(order);
        partition_matrix(A, A00,A01,A10,A11);
        partition_matrix(B, B00,B01,B10,B11);

        CMatrix<T> M1 = d_and_q_mult((A00+A11),(B00+B11));
        CMatrix<T> M2 = d_and_q_mult((A10 + A11),B00);
        CMatrix<T> M3 = d_and_q_mult(A00,(B01-B11));
        CMatrix<T> M4 = d_and_q_mult(A11,(B10-B00));
        CMatrix<T> M5 = d_and_q_mult((A00+A01),B11);
        CMatrix<T> M6 = d_and_q_mult((A10 - A00),(B00+B01));
        CMatrix<T> M7 = d_and_q_mult((A01 - A11),(B10+B11));

        // Assemble
        merge_matrix(C,M1+M4-M5+M7,M3+M5,M2+M4,M1-M2+M3+M6);
    }
    return C;
}

template <typename T>
bool StrassenMult<T>::partition_matrix(const CMatrix<T>& A,CMatrix<T>& A00, CMatrix<T>& A01, CMatrix<T>& A10, CMatrix<T>& A11)
{
    std::size_t order = A.row() / 2;
    for(std::size_t i = 0; i < order; ++i)
    {
        for(std::size_t j = 0; j < order; ++j)
        {
            A00(i,j) = A(i,j);
            A01(i,j) = A(i,j+order);
            A10(i,j) = A(i+order,j);
            A11(i,j) = A(i+order,j+order);
        }
    }
    return true;
}

template <typename T>
bool StrassenMult<T>::merge_matrix(CMatrix<T>& A,const CMatrix<T>& A00, const CMatrix<T>& A01, const CMatrix<T>& A10, const CMatrix<T>& A11)
{
    std::size_t order = A00.row();
    for(std::size_t i = 0; i < order; ++i)
    {
        for(std::size_t j = 0; j < order; ++j)
        {
            A(i,j) = A00(i,j);
            A(i,j+order) = A01(i,j);
            A(i+order,j) = A10(i,j);
            A(i+order,j+order) = A11(i,j);
        }
    }
    return true;
}

template <typename T>
bool StrassenMult<T>::isPowerOfTwo(std::size_t n)
{
    if (n <= 0)
        return false;
    while (n != 1)
    {
        if (n%2 != 0)
            return false;
        n = n/2;
    }
    return true;
}


void strassens_matrix_multiplication()
{
    std::size_t order = 4;
    CMatrix<double> A = GetMatrix<double>(order,order,5);
    CMatrix<double> B = GetMatrix<double>(order,order,5);
    std::cout << "A : " << std::endl << A << std::endl;
    std::cout << "B : " << std::endl << B << std::endl;
    CMatrix<double> C = StrassenMult<double>::StrasMult(A, B);

    std::cout << "C : " << std::endl << C << std::endl;
    std::cout << "A*B : " << std::endl << A*B << std::endl;
}

#endif // STRASSENS_MATRIX_MULTIPLICATION_H
