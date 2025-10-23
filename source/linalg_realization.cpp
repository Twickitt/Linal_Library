#include <stdio.h>
#include "linalg_realization.hpp"
#include <cstdlib>
#include <cmath> 

using namespace linalg;
using std::cout;
using std::cin;
using std::size_t;
using std::initializer_list;
using std::swap;
using std::cerr;
using std::abs;
using std::fabs;


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
    m_ptr = new double[m_capacity]{0}; //allocating matrix with nulls 
}

//unified initial constructors
Matrix::Matrix(initializer_list<initializer_list<double>> list){
    m_rows = list.size();
    m_columns = list.begin()->size(); //computing size of sublist to know how many columns we need
    m_capacity = m_rows*m_columns;
    m_ptr = new double[m_capacity];
    size_t k = 0;
    
    for(const auto& row : list){//range-based cycle for each sublist, const because of list(list is always const)
        if(row.size() != m_columns)
            cout<<"This matrix is not rectangular"; //checking if m_columns equal to number of elements in sublist(<...>, <...>) 
        
        for(double val : row)//assigning elements from each sublist to matrix, increasing k on each itaration(starating from k = 0)
            m_ptr[k++] = val;
    }
}

Matrix::Matrix(initializer_list<double> list, size_t columns){
    m_rows = list.size()/columns;
    m_columns = columns;
    m_capacity = m_rows * m_columns;

    if(list.size() % columns != 0)
        cout<<"Matrix with such patameters is invalid"; 
    
    size_t i = 0;
    m_ptr = new double[m_capacity];
    for(double val : list)
        m_ptr[i++] = val;
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
        cerr<<"Matrix capacity must be the same";
    
    m_rows = rows;
    m_columns = cols;
}

void Matrix::reserve(size_t n){
    if(n <= m_capacity)
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

//matrix element changing methods
double& Matrix::operator()(size_t row, size_t column){
    if(row >= m_rows || column >= m_columns)
        cerr<<"This indexes aren't in range. The range is:" <<' '<<"[0;" <<m_rows - 1<<"," <<"]"<<"[;"<<m_columns - 1<<"]";
    
    return m_ptr[m_columns * row + column];
    
}

const double& Matrix::operator()(size_t row, size_t column) const{
    if(row >= m_rows || column >= m_columns)
        cerr<<"This indexes aren't in range. The range is:" <<' '<<"[0;" <<m_rows - 1<<"," <<"]"<<"[;"<<m_columns - 1<<"]";
    
    return m_ptr[m_columns * row + column];
}

Matrix Matrix::operator -() const{
    Matrix res = *this;

    for(size_t i = 0; i < m_capacity; ++i)
        m_ptr[i] = -m_ptr[i];

    return res; 
}

Matrix Matrix::operator +() const{
    return *this;
}

Matrix linalg::operator -(const Matrix& first, const Matrix& second){
    if(!size_check(first, second)){
        cerr <<"Operation can't be complited because of size diff. Lets return default matrix";
        return Matrix();
    }
    
    Matrix res = first;
    res -= second;
    return res;
}

Matrix linalg::operator +(const Matrix& first, const Matrix& second){
    if(!size_check(first, second)){
        cerr <<"Operation can't be complited because of size diff. Lets return default matrix";
        return Matrix();
    }

    Matrix res = first;
    res += second;
    return res;
}

Matrix& Matrix::operator -=(const Matrix& other){
    if(!size_check(*this, other)){
        cerr <<"Operation can't be complited because of size diff. Lets return first matrix";
        return *this;
    }

    for(size_t i = 0; i < m_capacity; ++i)
        m_ptr[i] -= other.m_ptr[i];
    
    return *this;
}

Matrix& Matrix::operator +=(const Matrix& other){
    if(!size_check(*this, other)){
        cerr <<"Operation can't be complited because of size diff. Lets return first matrix";
        return *this;
    }

    for(size_t i = 0; i < m_capacity; ++i)
        m_ptr[i] += other.m_ptr[i];
    
    return *this;
}

bool Matrix::operator ==(const Matrix& other) const{
    if(!size_check(*this, other)){
        cerr <<"Operation can't be complited because of size diff. Lets return false";
        return false;
    }

    for(size_t i = 0; i < m_capacity; ++i){
        if(!are_equal(m_ptr[i], other.m_ptr[i]))
            return false;
    }

    return true;
}

bool Matrix::operator !=(const Matrix& other) const{
    if(!size_check(*this, other)){
        cerr<<"Operation can't be complited because of size diff. Lets return false";
        return false;
    }

    return !(*this == other);  
}


Matrix Matrix::operator *(const Matrix& other) const{
    if(!size_check(*this, other)){
        cerr<<"Operation can't be complited because of size diff. Lets return first matrix";
        return *this;
    }

    Matrix res(*this);
    res *= other;
    return res; 
} 

Matrix Matrix::operator *(double x) const{
    Matrix res(m_rows, m_columns);
    for(size_t i = 0; i < m_capacity; ++i)
        res.m_ptr[i] = m_ptr[i] * x;
    
    return res; 
}

Matrix linalg::operator *(double x, const Matrix& other){
    return other*x; 
}   

Matrix& Matrix::operator *=(const Matrix& other){
    if(m_columns != other.m_rows){
        cerr<<"Operation can't be complited because of multiplication rules of matrices. Lets return first matrix";
        return *this;
    }

    Matrix res(m_rows, other.m_columns);
    for(size_t i = 0; i < m_rows; ++i){
        
        for(size_t j = 0; j < other.m_columns; ++j){
            double sum = 0;
            
            for(size_t k = 0; k < m_columns; ++k)
                sum += m_ptr[i*m_columns + k] * other.m_ptr[k*other.m_columns + j];
            res.m_ptr[i*other.m_columns + j] = sum;
            
        }

    }
    
    *this = res;
    return *this;
    
} 

Matrix& Matrix::operator*=(double x) {
    for (size_t i = 0; i < m_capacity; ++i) {
        m_ptr[i] *= x;
    }
    
    return *this;
}

//helping methods 
bool are_equal(double x, double y, double eps = exp(-12)) {
    return abs(x - y) < eps;
}