#include "linalg_realization.hpp"

using namespace linalg;
using std::cout;
using std::cin;
using std::size_t;
using std::initializer_list;
using std::swap;
using std::abs;
using std::fabs;

constexpr double EPS = 1e-12;


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

Matrix Matrix::uno(size_t size) {
    Matrix res(size, size);
    for (size_t i = 0; i < size; ++i)
        res(i, i) = 1.0;  
    return res;
    }

linalg::Matrix::Matrix(size_t rows, size_t cols, double value) : m_rows(rows), m_columns(cols), m_capacity(rows*cols){
    m_ptr = new double[m_capacity];
    std::fill(m_ptr, m_ptr + m_capacity, value);//start, finish, chosen elem
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
    if (rows == 0 || cols == 0){
            m_rows = 0; 
            m_columns = 0;
        }
        if (m_rows == rows && m_columns == cols)
            return; 

        reserve(rows * cols);
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

    for(size_t i = 0; i < size(); ++i)
        res.m_ptr[i] = -m_ptr[i];

    return res; 
}

Matrix Matrix::operator +() const{
    return *this;
}

Matrix linalg::operator -(const Matrix& first, const Matrix& second){
    if(!size_check(first, second))
        throw std::runtime_error("Operation can't be complited because of size diff");
    
    Matrix res = first;
    res -= second;
    return res;
}

Matrix linalg::operator +(const Matrix& first, const Matrix& second){
    if(!size_check(first, second))
        throw std::runtime_error("Operation can't be complited because of size diff");

    Matrix res = first;
    res += second;
    return res;
}

Matrix& Matrix::operator -=(const Matrix& other){
    if(!size_check(*this, other))
        throw std::runtime_error ("Operation can't be complited because of size diff");

    for(size_t i = 0; i < size(); ++i)
        m_ptr[i] -= other.m_ptr[i];
    
    return *this;
}

Matrix& Matrix::operator +=(const Matrix& other){
    if(!size_check(*this, other))
        throw std::runtime_error("Operation can't be complited because of size diff");

    for(size_t i = 0; i < size(); ++i)
        m_ptr[i] += other.m_ptr[i];
    
    return *this;
}

bool Matrix::operator ==(const Matrix& other) const{
    if(!size_check(*this, other))
        return false;

    for(size_t i = 0; i < size(); ++i){
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
        throw std::runtime_error("Operation can't be complited because of because of multiplication rules of matrices");
    
    Matrix res(*this);
    res *= other;
    return res; 
} 

Matrix Matrix::operator *(double x) const{
    Matrix res(m_rows, m_columns);
    for(size_t i = 0; i < size(); ++i)
        res.m_ptr[i] = m_ptr[i] * x;
    
    return res; 
}

Matrix linalg::operator *(double x, const Matrix& other){
    return other*x; 
}   

Matrix& Matrix::operator *=(const Matrix& other){
    if(m_columns != other.m_rows)
        throw std::runtime_error("Operation can't be complited because of multiplication rules of matrices");

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
    for (size_t i = 0; i < size(); ++i) {
        m_ptr[i] *= x;
    }
    
    return *this;
}

double Matrix::norm() const noexcept{
    double sum{0};
    for (const double* el = begin(); el != end(); ++el)
        sum += (*el) * (*el);
    return std::sqrt(sum);
}

double Matrix::trace() const{
    if(m_rows != m_columns)
        throw std::runtime_error("Operation can't be complited because finding trace is only possible for square matrices");
    
    double sum = 0.0;
    for(size_t i = 0; i < size(); i+=(m_columns + 1)){
        sum+= m_ptr[i];
    }

    return sum;
}

double Matrix::det() const {
    if (m_rows != m_columns || m_columns == 0 || m_rows == 0)
        throw std::runtime_error("Operation can't be completed: determinant is only defined for square matrices");

    int expression_sign = 0;
    Matrix upper = linalg::upper_triang(*this, expression_sign);

    double determinant = 1.0;
    for (size_t i = 0; i < m_rows; ++i)
        determinant *= upper(i, i);

    if (expression_sign % 2 != 0)
        determinant = -determinant;


    return determinant;
}

int Matrix::rank() const {
    if(this->rows()==0 || this->columns() == 0)
        throw std::runtime_error("Can't find rank of inproper matrix");
    int rank = 0;
    int useless = 0;
    Matrix Upper_Triangular = linalg::upper_triang(*this, useless);
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m_columns; ++j) {
            if (!are_equal(Upper_Triangular(i, j), 0.0)) {
                rank++;
                break; // переходим к следующей строке
            }
        }
    }
    return rank;
}

Matrix& Matrix::gauss_forward() {
    int expression_sign = 0;
    *this = linalg::upper_triang(*this, expression_sign);
    return *this;
}

Matrix& Matrix::gauss_backward() {
    int expression_sign = 0;
    *this = linalg::lover_triang(*this, expression_sign);
    return *this;
}

Matrix linalg::concatenate(const Matrix& first, const Matrix& second) {
    if (first.rows() != second.rows() || first.rows() ==0 || first.columns() == 0 || second.rows() == 0 || second.columns() == 0)
        throw std::runtime_error("Matrices must have the same amount of rows to proceed");

    Matrix res(first.rows(), first.columns() + second.columns());
    for (size_t i = 0; i < first.rows(); ++i) {
        for (size_t j = 0; j < first.columns(); ++j)
            res(i, j) = first(i, j);
        for (size_t j = 0; j < second.columns(); ++j)
            res(i, j + first.columns()) = second(i, j);
    }
    return res;
}


Matrix linalg::transpose(const Matrix& matr) {
    if(matr.rows() == 0 || matr.columns() == 0)
        throw std::runtime_error("Matrix can't be transposed");
    Matrix res(matr.columns(), matr.rows());
    for (size_t i = 0; i < matr.rows(); ++i)
        for (size_t j = 0; j < matr.columns(); ++j)
            res(j, i) = matr(i, j);
    return res;
}

Matrix linalg::invert(const Matrix& matr) {
    if (matr.rows() != matr.columns())
        throw std::runtime_error("Matrix must be square to invert");
    if(matr.empty())
        throw std::runtime_error("Empty matrix isn't valid for this action");
    if(std::abs(matr.det()) < EPS)
        throw std::runtime_error("Matrices with zero det can't have opposite matrix");
    size_t size = matr.rows();
    Matrix hz = concatenate(matr, Matrix::uno(size));
    hz.gauss_forward();
    hz.gauss_backward();

    for (size_t i = 0; i < size; ++i) {
        double diag = hz(i, i);
        for (size_t j = 0; j < 2 * size; ++j)
            hz(i, j) /= diag;
    }

    Matrix inv(size, size);
    for (size_t i = 0; i < size; ++i)
        for (size_t j = 0; j < size; ++j)
            inv(i, j) = hz(i, j + size);

    return inv;
}

Matrix linalg::power(const Matrix& matr, int digit) {
    if (digit < 0)
        throw std::runtime_error("Only positive powers allowed");
    if (digit == 0)
        return Matrix::uno(matr.rows());
    if (matr.rows() != matr.columns())
        throw std::runtime_error("Matrix must be square");
    if(digit == 1)
        return matr;
    Matrix res = Matrix::uno(matr.rows());
    Matrix main = matr;
    while (digit > 1) {
        if (digit % 2 == 1)
            res *= main;
        main *= main;
        digit /= 2;
    }
    return res;
}

Matrix linalg::solve(const Matrix& matr, const Matrix& vector) {
    if (matr.rows() != matr.columns() || matr.rows() == 0 || matr.columns() == 0)
        throw std::runtime_error("This solve is only available for square matrices");
    if (vector.rows() != matr.rows() || vector.columns() != 1)
        throw std::runtime_error("Vector must have the same number of rows as matrix");
    if(matr.empty() && vector.empty())
        throw std::runtime_error("Infinity ammount of sollutions");
    if((matr.empty() && !vector.empty()) || (!matr.empty() && vector.empty()))
        throw std::runtime_error("One of this matrices is empty, so no solution");

    Matrix united = concatenate(matr, vector);
    united.gauss_forward();
    united.gauss_backward();

    Matrix solution(vector.rows(), 1);
    for (size_t i = 0; i < vector.rows(); ++i)
        solution(i, 0) = united(i, matr.columns());

    return solution;
}


//helping methods 
Matrix linalg::upper_triang(const Matrix& matr, int& expression_sign) {
    Matrix copy(matr);
    expression_sign = 0;
    for (size_t i = 0; i < matr.rows(); ++i) {
        if (are_equal(copy(i, i), 0.0)) {
            size_t next_row = i + 1;
            while (next_row < matr.rows() && are_equal(copy(next_row, i), 0.0))
                ++next_row;
            if (next_row == matr.rows())
                continue;
            for (size_t k = 0; k < matr.columns(); ++k)
                std::swap(copy(i, k), copy(next_row, k));
            expression_sign++;
        }
        for (size_t j = i + 1; j < matr.rows(); ++j) {
            double factor = copy(j, i) / copy(i, i);
            for (size_t k = i; k < matr.columns(); ++k)
                copy(j, k) -= factor * copy(i, k);
        }
    }
    return copy;
}

Matrix linalg::lover_triang(const Matrix& matr, int& expression_sign) {
    Matrix copy(matr);
    expression_sign = 0;
    for (size_t i = 0; i < matr.rows(); ++i) {
        if (are_equal(copy(i, i), 0.0)) {
            size_t next_row = i + 1;
            while (next_row < matr.rows() && are_equal(copy(next_row, i), 0.0))
                ++next_row;
            if (next_row == matr.rows())
                continue;
            for (size_t k = 0; k < matr.columns(); ++k)
                std::swap(copy(i, k), copy(next_row, k));
            expression_sign++;
        }
        for (size_t j = 0; j < i; ++j) {
            double factor = copy(j, i) / copy(i, i);
            for (size_t k = i; k < matr.columns(); ++k)
                copy(j, k) -= factor * copy(i, k);
        }
    }
    return copy;
}

//output
std::ostream& linalg::operator<<(std::ostream& out, const Matrix& m) {
    if (m.empty()) {
        out << "||";
        return out;
    }

    // Находим максимальную ширину в каждом столбце
    std::vector<size_t> col_widths(m.columns(), 0);
    for (size_t j = 0; j < m.columns(); ++j) {
        for (size_t i = 0; i < m.rows(); ++i) {
            std::ostringstream oss;
            oss << m(i, j);
            col_widths[j] = std::max(col_widths[j], oss.str().length());
        }
    }

    // Выводим матрицу с выравниванием
    for (size_t i = 0; i < m.rows(); ++i) {
        out << "|";
        for (size_t j = 0; j < m.columns(); ++j) {
            if (j > 0) out << " ";
            out << std::setw(col_widths[j]) << m(i, j);
        }
        out << "|\n"; // \n для каждой строки, как ожидает тест
    }

    return out;
}
