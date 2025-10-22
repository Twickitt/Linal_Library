#include <stdio.h>
#include "linalg_realization.hpp"

using namespace linalg;
using std::size_t;
using std::initializer_list;
using std::swap;


//constructors
Matrix::Matrix(const Matrix& other) : m_rows(other.m_rows), m_columns(other.m_columns){ //copying constructor 
    m_capacity = m_rows * m_columns;
    m_ptr = new double[m_capacity];
    
    for(size_t i = 0; i < m_capacity; ++i)
        m_ptr[i] = other.m_ptr[i];
}


Matrix::Matrix(Matrix&& other) noexcept : m_ptr(other.m_ptr), m_rows(other.m_rows), m_columns(other.m_columns), m_capacity(other.m_capacity){ //moving constructor
    other.m_ptr = nullptr;
    other.m_rows = other.m_columns = other.m_capacity = 0; //we have moved all values from target matrix to other matrix and imitializing old matrix with nulls
} 


Matrix::Matrix(int rows, int columns) : Matrix(){ //constructor with parameters using default constructor
    //we are using default constructor btw
    m_rows = rows;
    m_columns = columns;
    m_capacity = m_rows*m_columns;
    delete[] m_ptr;
    m_ptr = new double[m_capacity]{0}; //allocating memmory with nulls 
}

//unified initial constructor
Matrix::Matrix(initializer_list<initializer_list<double>> list){
    //realization
}

//operators
Matrix& Matrix::operator=(const Matrix& other){ //copy assignment operator
    if(this == &other)//checking if matrix is the same to avoid copying same matrices
        return *this;
    
    size_t needed_capacity = other.m_rows * other.m_columns;
    
    if(m_capacity < needed_capacity){//comparing matrices capacity to avoid useless memmory allocation
        delete[] m_ptr; //deliting old results of previous memmory allocations(in case they were complited)
        m_capacity = needed_capacity;
        m_ptr = new double[m_capacity];
    }

    for(size_t i = 0; i < m_capacity; ++i){
        m_ptr[i] = other.m_ptr[i];
    }
    
    m_rows = other.m_rows;
    m_columns = other.m_columns;
    
    return *this; 
} 

Matrix& Matrix::operator=(Matrix&& other) noexcept{
    if(this == &other)//checking if matrix is the same to avoid copying same matrices
        return *this;
    
    delete[] m_ptr;
    m_ptr = other.m_ptr; //copying new array as lvalue 
    m_rows = other.m_rows; 
    m_columns = other.m_columns;
    m_capacity = other.m_capacity;
    
    other.m_ptr = nullptr;
    other.m_rows = 0;
    other.m_columns = 0;
    other.m_capacity = 0;
    return *this;
}

//matrix modifications
void Matrix::reshape(size_t rows, size_t cols){
    if(rows *cols != m_capacity)
        throw std::invalid_argument("Matrix capacity must be the same");
    
    m_rows = rows;
    m_columns = cols;
}

void Matrix::reserve(size_t n){
    if(n < m_capacity)
        return;

    size_t new_size = m_rows * m_columns;
    double* new_ptr = new double[n];
    for(size_t i = 0; i < new_size; ++i){
        new_ptr[i] = m_ptr[i];
    }

    delete[] m_ptr;
    m_ptr = new_ptr;
    m_capacity = n;
}

void Matrix::clear() noexcept{
    m_rows = 0;
    m_columns = 0;
}

void Matrix::shrink_to_fit(){
    size_t new_capacity = m_rows * m_columns;
    if(new_capacity == m_capacity)
        return;

    double* new_ptr = nullptr;
    if(new_capacity == 0){
        delete[] m_ptr;
        m_ptr = nullptr;
        m_capacity = 0;
        return;
    }
    
    if(new_capacity > 0){
        new_ptr = new double[new_capacity];
        for(size_t i = 0; i < new_capacity; ++i)
            new_ptr[i] = m_ptr[i];
    }

    delete[] m_ptr;
    m_ptr = new_ptr;
    m_capacity = new_capacity;
}

void Matrix::swap(Matrix& other) noexcept{
    ::swap(m_rows, other.m_rows);
    ::swap(m_columns, other.m_columns);
    ::swap(m_capacity, other.m_capacity);
    ::swap(m_ptr, other.m_ptr);
}

//helping methods 
Matrix Matrix::upper_triang() const{
    
}

//Matrix substraction
Matrix Matrix::operator -() const{

}

