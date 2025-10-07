#include <iostream>

namespace linalg
{
    std::size_t;
    
    class Matrix {
    private:
        //sources
        double* m_ptr = nullptr; //отметим его нулевым указателем ради избежания возможных ошибок(UB) в дальнейшем(например при вызове дистрактора)
        size_t m_rows;
        size_t m_columns;
        size_t m_capacity;

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
        void shrink_to_fit(size_t size);
        void swap(Matrix& other) noexcept; //смена элементов тоже не дает исключений и ошибок
        void swap(Matrix& first, Matrix& second) noexcept {first.swap(second);} //глобальная функция обмена между двумя матрицами 

        //matrix constructors 
        Matrix() noexcept {}; //default constr
        Matrix(const Matrix& other); //copying constr, no noexcept потому, что не можем гарантировать отсутствие исключений(например при выделении памяти)
        Matrix(Matrix&& other) noexcept; //moving constr
        Matrix(int rows, int columns); //добавил конструктор с параметрами
        ~Matrix() {delete[] m_ptr;} //destructor(~ обозначает деструктор)

        //copying and assignment opperators(используем operator= для перегрузки)
        Matrix& operator=(const Matrix& other); // оператор копирующего присваивания 
        Matrix& operator=(Matrix&& other) noexcept; // оператор перемещающего присваивания(не const так как забирает русурсы у объекта при копировании)


        //matrix Operators
        //subtraction  
        Matrix operator -() const;
        Matrix operator -(const Matrix& other) const;
        Matrix operator -=(const Matrix& other);
        //smmation
        Matrix operator +() const;
        Matrix operator +(const Matrix& other) const;
        Matrix operator +=(const Matrix& other);
        //checking equalities and inequalities
        bool operator ==(const Matrix&) const;
        bool operator !=(const Matrix&) const;
        //multiplication
        Matrix operator *(const Matrix& other); //матрица на матрицу 
        friend Matrix operator *(double x, const Matrix& other); //число на матрицу  
        Matrix operator *(double x) const;
        Matrix operator *=(const Matrix& other); //умножение и присваивание результата операции
        Matrix operator *=(double x); //поэлементарное умножение и присваивание результата 
    
    };
}
