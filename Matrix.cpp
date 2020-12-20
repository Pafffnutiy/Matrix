#include "Matrix.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

Matrix::Matrix()
{
    rows = 0;
    cols = 0;
    
    data=nullptr;
}

Matrix::~Matrix()
{
    if( data!=nullptr )
    {
        for (size_t i=0; i<rows; ++i) delete[] data[i];
        delete[] data;
        //освоободить памяти
        data=nullptr;
    }
}

Matrix::Matrix( const Matrix & other )
{
    rows=other.rows;
    cols=other.cols;
    
    try
    {
        data = new double * [rows];
        
        for(size_t r = 0; r<rows; ++r)
            data[r] = new double [cols];
        
        for(size_t r = 0; r<rows; ++r)
            for(size_t c = 0; c<cols; ++c)
            data[r][c] = other.data[r][c];
    }
    catch(...)
    {
        std::cerr << "Error allocating memory for matrix!" <<  std::endl;
        exit(1);
    }
}

Matrix::Matrix(size_t rows, size_t cols)
{
    this->rows=rows;
    this->cols=cols;
    
    try
    {
        data = new double * [rows];
        
        for(size_t r = 0; r<rows; ++r)
            data[r] = new double [cols];
        
        for(size_t r = 0; r<rows; ++r)
            for(size_t c = 0; c<cols; ++c)
            data[r][c] = 0.0;
    }
    catch(...)
    {
        std::cerr << "Error allocating memory for matrix!" <<  std::endl;
        exit(1);
    }
}

Matrix::Matrix(size_t order) 
{
    this->rows=order;
    this->cols=order;
    
    try
    {
        data = new double * [order];
        
        for(size_t r = 0; r<order; ++r)
            data[r] = new double [order];
        
        for(size_t r = 0; r<order; ++r)
            for(size_t c = 0; c<order; ++c)
            data[r][c] = 0.0;
    }
    catch(...)
    {
        std::cerr << "Error allocating memory for matrix!" <<  std::endl;
        exit(1);
    }
}

void Matrix::SwapRows(size_t r1, size_t r2){
        for(size_t i = 0;i<this->cols;++i){
            //std::swap(this->data[r1][i],this->data[r2][i]);
            this->data[r1][i]+=this->data[r2][i];
            this->data[r2][i]=this->data[r1][i]-this->data[r2][i];
            this->data[r1][i]-=this->data[r2][i];
        }
        
    }
    
void Matrix::Triangle() {
        for (size_t l = 0; l < this->rows; ++l) {
            size_t k{ 1 };
            for (size_t r = l + 1; r < this->rows; ++r) {
                double x = -(this->data[r][l] / this->data[l][l]);
                for (size_t c =0; c < this->rows; ++c) {
                    this->data[r][c] = this->data[r][c] + this->data[r - k][c] * x;
                }
                k++;
            }
        }
    }
    
long double Matrix::CalcDeterminant() const {
    Matrix M(this->rows,this->cols);
        
    for(size_t i=0; i<M.rows; ++i)
        for(size_t j=0; j<M.cols; ++j) {
            M.data[i][j]=this->data[i][j];
        }
        
    assert(M.rows==M.cols);
        
    //if(this->rows!=this->cols) {std::cerr<<"The determinant can't be calculated for non-square matrices"<<'\n'; return 1;}
    double det{1};
    size_t i{1};
    while(M.data[0][0]==0){
        M.SwapRows(0,i);
        ++i;
        if (i==M.rows-1){ return (det=0); }
    }
    if (i!=1) {det*=(-1);}
    M.Triangle();
    std::cout<<'\n';
    for (size_t i=0; i < M.rows; ++i) {
        det *= M.data[i][i];
    }
        
    if ((fabs(det)==0) or (std::isnan(fabs(det)))) {return (det=0);} else {return det;}
}

void Matrix::FillMagickSE(){
    unsigned cnt = 1;
    
    for (size_t i = 0; i < this->rows; ++i)
        for (size_t j = 0; j < this->cols; ++j)
            if (i +j > this->cols-1){
                this->data[i][j] = cnt;
                ++cnt;
            }
}

void Matrix::FillRandomInt(){
    for (size_t i = 0; i < this->rows; ++i)
        for (size_t j = 0; j < this->cols; ++j) {
            this->data[i][j]=rand() % 20 - 10 + 1;
        }
}    
    
void Matrix::FillRandomDouble(){
    for (size_t i = 0; i < this->rows; ++i)
        for (size_t j = 0; j < this->cols; ++j) {
            this->data[i][j] = (double)rand() / (double)RAND_MAX * (20 - -12) + -12;
    }
}

Matrix Matrix::operator+(const Matrix & right) const
{
    if( (rows==right.rows) && (cols==right.cols) )
    {
        Matrix result(right); //дописать сложение матриц
        for (size_t i = 0; i < this->rows; ++i)
            for (size_t j = 0; j < this->cols; ++j)
            {
                result.data[i][j]=data[i][j]+right.data[i][j];
            }
        return result;
    }
    else{
        std::cerr<<"Error: Adding matrices of different sizes"<<'\n';
        exit(1);
    }
    
    return Matrix();
}

Matrix Matrix::operator- (const Matrix & right) const
{
    if( (rows==right.rows) && (cols==right.cols) )
    {
        Matrix result(right); //дописать сложение матриц
        for (size_t i = 0; i < this->rows; ++i)
            for (size_t j = 0; j < this->cols; ++j)
            {
                result.data[i][j]=data[i][j]-right.data[i][j];
            }
        return result;
    }
    else {
        std::cerr<<"Error: Adding matrices of different sizes"<<'\n';
        exit(1);
    }
    
    return Matrix();
}

Matrix Matrix::operator* (const Matrix & right) const
{
    if(this->cols!=right.rows){
        std::cerr<<"Error: The number of columns of the matrix on the left must be equal to the number of rows of the matrix on the right"<<'\n';
        exit(1);
    }
    Matrix result(this->rows,right.cols);
    for(size_t r=0;r<this->rows;++r)
        for(size_t c=0;c<right.cols;++c)
            for(size_t s=0;s<this->cols;++s){
                result.data[r][c]+=(this->data[r][s]*right.data[s][c]);
            }
    return result;
}

std::ostream & operator<< (std::ostream & left, const Matrix & right)
{
    left<<right.rows<<' '<<right.cols<<std::endl;
    
    for(size_t r=0; r<right.rows;++r)
    {
        for(size_t c=0;c<right.cols;++c)
            left<<right.data[r][c]<<' ';
        
        left<<std::endl;
    }
    return left;
}

std::istream & operator>> (std::istream & in, Matrix & right)
{
    
    if (right.data != nullptr) {
        for (size_t i=0; i<right.rows; ++i) delete[] right.data[i];
        delete[] right.data;
        right.data = nullptr;
    };
    
    in >> right.rows >> right.cols;
    
    try
    {
        right.data = new double * [right.rows];
        
        for(size_t r = 0; r<right.rows; ++r)
            right.data[r] = new double [right.cols];
        
        for(size_t r = 0; r<right.rows; ++r)
            for(size_t c = 0; c<right.cols; ++c)
                in >> right.data[r][c]; 
    }
    catch(...)    
    {
        std::cerr<<"Error allocating memory for matrix!"<<std::endl;
        exit(1);
    }
    return in;
}

Matrix operator* (double left, const Matrix & right)
{
    Matrix result(right.rows, right.cols);
    for(size_t r=0; r<right.rows;++r)
    {
        for(size_t c=0;c<right.cols;++c)
            result.data[r][c]=left*right.data[r][c];
    }
    return result;
}

