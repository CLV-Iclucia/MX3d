#ifndef MX3D_MX3D_MATRIX_H
#define MX3D_MX3D_MATRIX_H
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <initializer_list> 
#include <cstring>
#include "mx3d_Vector.h"

namespace mx3d
{
    template<typename T> class Mat3x3;
    template<typename T>
    Mat3x3<T> operator*(T v, const Mat3x3<T>& A);
    // TODO: implement matrices of more shapes
    template<typename T>
    class Mat3x3
    {
            static_assert(std::is_arithmetic_v<T>, "ERROR: in instantiation of Vector: the type must be arithmetic type.");
            template<typename Any>
            friend Mat3x3 operator*<Any>(Any, const Mat3x3&);
        public:
            static Mat3x3 Identity()
            {
                Mat3x3 ret;
                ret.value[0] = ret.value[4] = ret.value[8] = static_cast<T>(1);
                return ret;
            }
            Mat3x3() { for(uint i = 0; i < 9; i++) value[i] = 0; }
            Mat3x3(std::initializer_list<T> lst)
            {
                auto it = lst.begin();
                for (int i = 0; i < 9; i++, it++)
                    value[i] = *it;
            }
            Mat3x3(const Mat3x3& M) = default;
            Mat3x3(const Vec3& A, const Vec3& B, const Vec3& C)
            {
                value[0] = A[0];
                value[1] = A[1];
                value[2] = A[2];
                value[3] = B[0];
                value[4] = B[1];
                value[5] = B[2];
                value[6] = C[0];
                value[7] = C[1];
                value[8] = C[2];
            }
            T& operator()(int i, int j) { return value[3 * i + j]; }
            const T& operator()(int i, int j) const { return value[3 * i + j]; }
            T& operator[](int i) { return value + 3 * i; }
            const T& operator[](int i) const { return value + 3 * i; }
            Vec3 operator*(const Vec3& A) const
            {
                Vec3 V;
                V[0] = value[0] * A[0] + value[3] * A[1] + value[6] * A[2];
                V[1] = value[1] * A[0] + value[4] * A[1] + value[7] * A[2];
                V[2] = value[2] * A[0] + value[5] * A[1] + value[8] * A[2];
                return V;
            }
            Mat3x3 operator*(const Mat3x3& A) const
            {
                Mat3x3 M;
                M.value[0] = value[0] * A.value[0] + value[3] * A.value[1] + value[6] * A.value[2];
                M.value[3] = value[0] * A.value[3] + value[3] * A.value[4] + value[6] * A.value[5];
                M.value[6] = value[0] * A.value[6] + value[3] * A.value[7] + value[6] * A.value[8];
                M.value[1] = value[1] * A.value[0] + value[4] * A.value[1] + value[7] * A.value[2];
                M.value[4] = value[1] * A.value[3] + value[4] * A.value[4] + value[7] * A.value[5];
                M.value[7] = value[1] * A.value[6] + value[4] * A.value[7] + value[7] * A.value[8];
                M.value[2] = value[2] * A.value[0] + value[5] * A.value[1] + value[8] * A.value[2];
                M.value[5] = value[2] * A.value[3] + value[5] * A.value[4] + value[8] * A.value[5];
                M.value[8] = value[2] * A.value[6] + value[5] * A.value[7] + value[8] * A.value[8];
                return M;
            }
            Mat3x3 operator+(const Mat3x3& A) const
            {
                Mat3x3 M;
                for (int i = 0; i < 9; i++)
                    M.value[i] = value[i] + A.value[i];
                return M;
            }
            Mat3x3 operator/(const Mat3x3& A) const
            {
                for (int i = 0; i < 9; i++)
                {
                    try
                    {
                        if (A.value[i] == 0)throw std::runtime_error("Division by zero!");
                    }
                    catch (const std::exception& e)
                    {
                        std::cerr << e.what() << '\n';
                        exit(-1);
                    }
                }
                Mat3x3 M;
                for (int i = 0; i < 9; i++)
                    M.value[i] = value[i] / A.value[i];
                return M;
            }
            Mat3x3 operator-(const Mat3x3& A) const
            {
                Mat3x3 M;
                for (int i = 0; i < 9; i++)
                    M.value[i] = value[i] - A.value[i];
                return M;
            }
            Mat3x3 operator-() const
            {
                Mat3x3 M;
                for (int i = 0; i < 9; i++)
                    M.value[i] = -value[i];
                return M;
            }
            template<typename Any>
            Mat3x3 operator/(Any v) const
            {
                try
                {
                    if (isZero(v))throw std::runtime_error("Division by zero!");
                }
                catch (const std::exception& e)
                {
                    std::cerr << e.what() << '\n';
                    exit(-1);
                }
                Mat3x3 M;
                for (int i = 0; i < 9; i++)
                    M.value[i] = value[i] / v;
                return M;
            }
            template<typename Any>
            Mat3x3 operator*(Any v) const
            {
                static_assert(std::is_arithmetic_v<Any>);
                Mat3x3 M;
                for (int i = 0; i < 9; i++)
                    M.value[i] = v * value[i];
                return M;
            }
            Mat3x3& operator*=(const Mat3x3& A)
            {
                for (int i = 0; i < 9; i++)
                    value[i] *= A.value[i];
                return *this;
            }
            Mat3x3& operator+=(const Mat3x3& A)
            {
                for (int i = 0; i < 9; i++)
                    value[i] += A.value[i];
                return *this;
            }
            Mat3x3& operator/=(const Mat3x3& A)
            {
                for (int i = 0; i < 9; i++)
                {
                    try
                    {
                        if (isZero(A.value[i]) == 0)throw std::runtime_error("Division by zero!");
                        else value[i] /= A.value[i];
                    }
                    catch (const std::exception& e)
                    {
                        std::cerr << e.what() << '\n';
                        exit(-1);
                    }
                }
                return *this;
            }
            template<typename Any>
            Mat3x3& operator/=(Any v)
            {
                static_assert(std::is_arithmetic_v<Any>, "ERROR: divider must be arithmetic type");
                try
                {
                    if (v == 0)throw std::runtime_error("Division by zero!");
                }
                catch (const std::exception& e)
                {
                    std::cerr << e.what() << '\n';
                    exit(-1);
                }
                for (int i = 0; i < 9; i++)
                    value[i] /= v;
                return *this;
            }
            Mat3x3 transpose() const
            {
                Mat3x3 ret;
                for(int i = 0; i < 3; i++)
                    for(int j = 0; j < 3; j++)
                        ret.value[3 * i + j] = *(value + 3 * i + j);
                return ret;
            }
            Mat3x3& operator-=(const Mat3x3& A)
            {
                for (int i = 0; i < 9; i++)
                    value[i] -= A.value[i];
                return *this;
            }
            template<typename Any>
            Mat3x3& operator*=(Any v)
            {
                static_assert(std::is_arithmetic_v<Any>, "ERROR: must be arithmetic type.");
                for (int i = 0; i < 9; i++)
                    value[i] *= v;
                return *this;
            }
            ~Mat3x3() = default;
            Real det() const
            {
                return value[0] * value[4] * value[8] + value[3] * value[7] * value[2] + value[1] * value[5] * value[7]
                        - value[2] * value[4] * value[6] - value[0] * value[5] * value[7] - value[1] * value[3] * value[8];
            }
            Mat3x3<Real> inv() const
            {
                Mat3x3 mat;
                Real detv = det();
                mat.value[0] = (value[4] * value[8] - value[7] * value[5]) / detv;
                mat.value[1] = (value[6] * value[5] - value[3] * value[8]) / detv;
                mat.value[2] = (value[3] * value[7] - value[4] * value[6]) / detv;
                mat.value[3] = (value[2] * value[7] - value[1] * value[8]) / detv;
                mat.value[4] = (value[0] * value[8] - value[2] * value[6]) / detv;
                mat.value[5] = (value[1] * value[6] - value[0] * value[7]) / detv;
                mat.value[6] = (value[1] * value[5] - value[2] * value[4]) / detv;
                mat.value[7] = (value[3] * value[4] - value[0] * value[5]) / detv;
                mat.value[8] = (value[0] * value[4] - value[1] * value[3]) / detv;
                return mat;
            }
        private:
            T value[9];//Column-Major Order
    };
    template<typename Any, typename T>
    Mat3x3<T> operator*(Any v, const Mat3x3<T>& A) { return A * v; }
    using Mat3 = Mat3x3<Real>;
}
#endif
