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

    for(size_t i = 0; i < needed_capacity; ++i){
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
    if(rows *cols != matr_size())
        throw std::invalid_argument("Matrix capacity must be the same");

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
        throw std::out_of_range("This indexes aren't in range");
    
    return m_ptr[m_columns * row + column];
    
}

const double& Matrix::operator()(size_t row, size_t column) const{
    if(row >= m_rows || column >= m_columns)
        throw std::out_of_range("This indexes aren't in range");;

    return m_ptr[m_columns * row + column];
}

Matrix Matrix::operator -() const{
    Matrix res = *this;

    for(size_t i = 0; i < matr_size(); ++i)
        m_ptr[i] = -m_ptr[i];

    return res; 
}

Matrix Matrix::operator +() const{
    return *this;
}

Matrix linalg::operator -(const Matrix& first, const Matrix& second){
    if(!size_check(first, second))
        throw std::invalid_argument("Operation can't be complited because of size diff");
    
    Matrix res = first;
    res -= second;
    return res;
}

Matrix linalg::operator +(const Matrix& first, const Matrix& second){
    if(!size_check(first, second))
        throw std::invalid_argument("Operation can't be complited because of size diff");

    Matrix res = first;
    res += second;
    return res;
}

Matrix& Matrix::operator -=(const Matrix& other){
    if(!size_check(*this, other))
        throw std::invalid_argument("Operation can't be complited because of size diff");

    for(size_t i = 0; i < matr_size(); ++i)
        m_ptr[i] -= other.m_ptr[i];
    
    return *this;
}

Matrix& Matrix::operator +=(const Matrix& other){
    if(!size_check(*this, other))
        throw std::invalid_argument("Operation can't be complited because of size diff");

    for(size_t i = 0; i < matr_size(); ++i)
        m_ptr[i] += other.m_ptr[i];
    
    return *this;
}

bool Matrix::operator ==(const Matrix& other) const{
    if(!size_check(*this, other))
        return false;

    for(size_t i = 0; i < matr_size(); ++i){
        if(!are_equal(m_ptr[i], other.m_ptr[i]))
            return false;
    }

    return true;
}

bool Matrix::operator !=(const Matrix& other) const{
    return !(*this == other);  
}


Matrix Matrix::operator *(const Matrix& other) const{
    if(m_columns != other.m_rows)
        throw std::invalid_argument("Operation can't be complited because of because of multiplication rules of matrices");
    
    Matrix res(*this);
    res *= other;
    return res; 
} 

Matrix Matrix::operator *(double x) const{
    Matrix res(m_rows, m_columns);
    for(size_t i = 0; i < matr_size(); ++i)
        res.m_ptr[i] = m_ptr[i] * x;
    
    return res; 
}

Matrix linalg::operator *(double x, const Matrix& other){
    return other*x; 
}   

Matrix& Matrix::operator *=(const Matrix& other){
    if(m_columns != other.m_rows)
        throw std::invalid_argument("Operation can't be complited because of multiplication rules of matrices");

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
    for (size_t i = 0; i < matr_size(); ++i) {
        m_ptr[i] *= x;
    }
    
    return *this;
}

double Matrix::norm() const noexcept{
    double sum{0};
    for (const double* el = matr_begin(); el != matr_end(); ++el)
        sum += (*el) * (*el);
    return std::sqrt(sum);
}

double Matrix::trace() const noexcept{
    if(m_rows != m_columns)
        throw std::invalid_argument("Operation can't be complited because finding trace is only possible for square matrices");
    
    double sum = 0.0;
    for(size_t i = 0; i < matr_size(); i+=(m_columns + 1)){
        sum+= m_ptr[i];
    }

    return sum;
}

double Matrix::det() const{
    if(m_rows != m_columns)
        throw std::invalid_argument("Operation can't be complited because finding det is only possible for square matrices");

    int expression_sign = 0;
    double determinant = 1.0;
    Matrix Upper_Triangular = this->upper_triang(expression_sign);
    for(size_t i = 0; i < m_rows; ++i)
        determinant*=Upper_Triangular(i, i);

    if(expression_sign % 2 != 0)
        determinant= -determinant;

    return determinant;
    
}

int Matrix::rank() const noexcept{
    int rank = 0;
    int useless = 0;
    Matrix Upper_Triangular = this->upper_triang(useless);
    for(size_t i = 0; i < m_rows; ++i){
        
        for(size_t j = 0; i < m_columns; ++j){
            
            if(!are_equal(Upper_Triangular(i, j), 0.0)){
                rank++;
                break;//collapsing cycle for j 
            }

        }

    }
    return rank;
}

Matrix& Matrix::gauss_forward(){
    int useless = 0; 
    Matrix Upper = upper_triang(useless);
    for(size_t i = 0; i < m_rows; ++i){

        for(size_t j = 0; i < m_columns; ++j)
            (*this)(i, j) = Upper(i, j);         
    }
    return *this;
}

Matrix& Matrix::gauss_backward(){
    int useless = 0; 
    Matrix Lower = lover_triang(useless);
    for(size_t i = 0; i < m_rows; ++i){

        for(size_t j = 0; i < m_columns; ++j)
            (*this)(i, j) = Lower(i, j);         
    }
    return *this;
}

//helping methods 
bool are_equal(double x, double y, double eps = exp(-12)) {
    return abs(x - y) < eps;
}

Matrix Matrix::upper_triang(int& expression_sign) const{
    Matrix copy(*this);
    expression_sign = 0;
    for(size_t i = 0; i < m_rows; ++i){
        
        if(are_equal(copy(i, i), 0.0)){
            size_t next_row = i + 1;
            
            while(next_row < m_rows && are_equal(copy(next_row, i), 0.0))
                ++next_row;
            
            if(next_row == m_rows)
                continue;
            
            for(size_t k = 0; k < m_columns; ++k)
                ::swap(copy(i,k), copy(next_row,k)); 

            expression_sign++;
        }
        // making nulls under diag 
        for(size_t j = (i+1); j < m_rows; ++j){
            double factor = copy(j,i) / copy(i,i);
            
            for(size_t k = i; k < m_columns; ++k)
                copy(j,k) -= factor * copy(i, k);
        }
    
    }

    return copy;
}

Matrix Matrix::lover_triang(int& expression_sign) const{
    Matrix copy(*this);
    expression_sign = 0;
    for (size_t i = 0; i < m_rows; ++i) {
        
        if (are_equal(copy(i, i), 0.0)) {
            size_t next_row = i + 1;

            while (next_row < m_rows && are_equal(copy(next_row, i), 0.0))
                ++next_row;

            if (next_row == m_rows)
                continue; 

            for (size_t k = 0; k < m_columns; ++k)
                ::swap(copy(i, k), copy(next_row, k));

            expression_sign++;
        }
        // making nulls above diag 
        for (size_t j = 0; j < i; ++j) {
            double factor = copy(j, i) / copy(i, i);
            for (size_t k = i; k < m_columns; ++k)
                copy(j, k) -= factor * copy(i, k);
        }
    }

    return copy;
}