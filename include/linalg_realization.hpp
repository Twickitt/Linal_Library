#include <iostream>
#include <initializer_list>


namespace linalg{

using std::size_t;
using std::initializer_list;

class Matrix {
private:
    //sources
    double* m_ptr = nullptr; //отметим его нулевым указателем ради избежания возможных ошибок(UB) в дальнейшем(например при вызове дистрактора)
    size_t m_rows, m_columns, m_capacity;

public:
    //metods 
    bool matr_empty() const {return m_rows == 0 || m_columns == 0;}
    size_t matr_rows() const {return m_rows;}
    size_t matr_columns() const {return m_columns;}
    size_t matr_capacity() const {return m_capacity;}
    size_t matr_size() const {return m_columns * m_rows;}
    const double* matr_begin() const {return m_ptr;} //указатель на 1 элемент матрицы  
    const double* matr_end () const {return m_ptr ? m_ptr + matr_size(): m_ptr;} //если m_ptr не nullptr, то вернёт m_ptr + size(). а если m_ptr =nullptr, тогда вернет nullptr; ?-тернарный оператор(сокращенный вар. if ... else)
    
    //matrix modification methods
    void reshape(size_t rows, size_t cols);
    void reserve(size_t n);
    void clear() noexcept; // очистка матрицы не дает исключений и ошибок, вместимость прежняя 
    void shrink_to_fit();
    void swap(Matrix& other) noexcept; //смена элементов тоже не дает исключений и ошибок
    friend void swap(Matrix& first, Matrix& second) noexcept {first.swap(second);} //глобальная функция обмена между двумя матрицами 
    
    //matrix constructors 
    Matrix() noexcept : m_rows{0}, m_columns{0}, m_capacity{0}, m_ptr(nullptr){} //default constr, m_capacity = 0 because of 0 size of the matrix 
    Matrix(const Matrix& other); //copying constr, no noexcept потому, что не можем гарантировать отсутствие исключений(например при выделении памяти)
    Matrix(Matrix&& other) noexcept; //moving constr
    Matrix(int rows, int columns = 1); //добавил конструктор с параметрами
    //unified initial constructors
    Matrix(initializer_list<initializer_list<double>> list);
    Matrix(initializer_list<double> list, size_t columns = 1);
    
    ~Matrix() {delete[] m_ptr;} //destructor(~ обозначает деструктор)
    //copying and assignment opperators(используем operator= для перегрузки)
    Matrix& operator=(const Matrix& other); // оператор копирующего присваивания 
    Matrix& operator=(Matrix&& other) noexcept; // оператор перемещающего присваивания(не const так как забирает русурсы у объекта при копировании)
    
    double& operator()(size_t row, size_t column);
    const double& operator()(size_t row, size_t column) const; 
    
    //matrix computation operators
    //subtraction  
    Matrix operator -() const;
    friend Matrix operator -(const Matrix& first, const Matrix& second);
    Matrix& operator -=(const Matrix& other);
    
    //smmation
    Matrix operator +() const;
    friend Matrix operator +(const Matrix& first, const Matrix& second);
    Matrix& operator +=(const Matrix& other);
    
    //checking equalities and inequalities
    bool operator ==(const Matrix& other) const;
    bool operator !=(const Matrix& other) const;
    
    //multiplication
    Matrix operator *(const Matrix& other) const; //матрица на матрицу 
    Matrix operator *(double x) const;
    friend Matrix operator *(double x, const Matrix& other); //число на матрицу  
    Matrix& operator *=(const Matrix& other); //умножение и присваивание результата операции
    Matrix& operator *=(double x); //поэлементарное умножение и присваивание результата 

    //computational operations
    double norm() const noexcept;
    double trace() const noexcept;
    double det() const;
    int rank() const noexcept;
    //friend is using to access internal resources of matrix
    Matrix& gauss_forward(); //Gaussian method forward 
    Matrix& gauss_backward(); //Gaussian method backward 
    friend Matrix concatenate(const Matrix& first, const Matrix& second); //union of two matrices
    friend Matrix transpose(const Matrix& matr); //transposing our matrix 
    friend Matrix invert(const Matrix& matr); //inverting rows and cols of uor matrix
    friend Matrix power(const Matrix& matr, int digit); //raising our matrix to a power
    friend Matrix solve(const Matrix& matr, const Matrix& vector); //solving matrix equation

    //helping methods
    friend inline bool size_check(const Matrix& first, const Matrix& second){return first.m_rows == second.m_rows && first.m_columns ==second.m_columns;};
    Matrix upper_triang(int& expression_sign) const;
    Matrix lover_triang(int& expression_sign) const;
};
}
