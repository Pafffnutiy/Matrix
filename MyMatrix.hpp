#ifndef Matrix_hpp
#define Matrix_hpp

#include <ostream>
#include <istream>

struct Matrix
{
    size_t rows; //строки
    size_t cols; //столбцы
    double **data; //элементы
    
    Matrix();//конструктор умолчания
   ~Matrix();//деструктор
    Matrix( const Matrix & other );//конструктор копирования
    Matrix(size_t rows, size_t cols);//конструктор прямоугольной матрицы
    Matrix(size_t order);//конструктор квалратной матрицы

    void SwapRows(size_t r1, size_t r2); //
    void Triangle();//приведение к треугольному виду
    long double CalcDeterminant() const;//Вычисление определителя
    void FillMagickSE();//заполнение матрицы
    void FillRandomInt();//заполнение случайными целыми
    void FillRandomDouble();//заполнение случайными double    
   
    Matrix operator+ (const Matrix & right) const;
    Matrix operator- (const Matrix & right) const;
    
};

std::ostream & operator<< (std::ostream & left, const Matrix & right);
std::istream & operator>> (std::istream & in, Matrix & right);

Matrix operator* (double left, const Matrix & right);
#endif //Matrix_h
