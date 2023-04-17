//
// Created by creeper on 22-10-26.
//
/**
 * \brief Vector class for Rendering
 */
#ifndef MX3D_MX3D_VECTOR_H
#define MX3D_MX3D_VECTOR_H

#include <mx3d_utils.h>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <initializer_list>

namespace mx3d
{
    /**
     * @tparam T type of the vector
     * @tparam N dimension of the vector.
     *           This vector usually has 2, 3 or 4 dims.
     *           That suits the need in graphics.
     * @class Vector is specifically designed for graphics and image processing.
     */
    template<typename T, unsigned N> class Vector;
    template<typename T, unsigned N>
    Vector<T, N> operator*(T v, const Vector<T, N>& A);
    template<typename T, unsigned N>
    Vector<T, N> operator/(T v, const Vector<T, N>& A);

    template<typename T, unsigned N>
    class Vector
    {
        static_assert(std::is_arithmetic_v<T>, "ERROR: in instantiation of Vector: the type must be arithmetic type.")
        public:
            friend Vector operator*<T, N>(const T&, const Vector&);
            friend Vector operator/<T, N>(const T&, const Vector&);
            Vector() {for (int i = 0; i < N; i++)value[i] = 0;}
            Vector(std::initializer_list<T> lst)
            {
                auto it = lst.begin();
                for (int i = 0; i < N; i++, it++)
                    value[i] = *it;
            }
            Vector(Vector&& V) noexcept
            {
                for (int i = 0; i < N; i++)
                    value[i] = std::move(V.value[i]);
            }
            Vector& operator=(Vector&& V) noexcept
            {
                for (int i = 0; i < N; i++)
                    value[i] = std::move(V.value[i]);
                return *this;
            }
            Vector& operator=(const Vector& V)
            {
                if(&V == this) return *this;
                for (int i = 0; i < N; i++)
                    value[i] = V.value[i];
                return *this;
            }
            Vector(const Vector& V)
            {
                for (int i = 0; i < N; i++)value[i] = V.value[i];
            }
            explicit Vector(T v)
            {
                for (int i = 0; i < N; i++)value[i] = v;
            }
            Real norm() const
            {
                Real sum = 0.0;
                for (int i = 0; i < N; i++)
                    sum += value[i] * value[i];
                return std::sqrt(sum);
            }
            Vector& normalize()
            {
                Real length = this->norm();
                if (isZero(length))return *this;
                else
                {
                    for (int i = 0; i < N; i++)
                        value[i] /= length;
                    return *this;
                }
            }
            Vector normalized() const
            {
                Real length = this->norm();
                if (isZero(length)) return Vector();
                else return (*this) / length;
            }
            T dot(const Vector& A) const
            {
                T sum = 0;
                for (int i = 0; i < N; i++)
                    sum += this->value[i] * A.value[i];
                return sum;
            }
            T& operator[](int i)
            {
                return value[i];
            }
            const T& operator[](int i) const
            {
                return value[i];
            }
            Vector operator*(const Vector& A) const
            {
                Vector V;
                for (int i = 0; i < N; i++)
                    V.value[i] = value[i] * A.value[i];
                return V;
            }
            Vector operator+(const Vector& A) const
            {
                Vector V;
                for (int i = 0; i < N; i++)
                    V.value[i] = value[i] + A.value[i];
                return V;
            }
            Vector operator-()
            {
                Vector V;
                for (int i = 0; i < N; i++)
                    V.value[i] = -value[i];
                return V;
            }
            Vector operator/(const Vector& A) const
            {
                Vector V;
                for (int i = 0; i < N; i++)
                {
                    assert(!isZero(A.value[i]));
                    V.value[i] = value[i] / A.value[i];
                }
                return V;
            }
            Vector operator-(const Vector& A) const
            {
                Vector V;
                for (int i = 0; i < N; i++)
                    V.value[i] = value[i] - A.value[i];
                return V;
            }
            Vector operator/(T v) const
            {
                assert(isZero(v));
                Vector V;
                for (int i = 0; i < N; i++)
                    V.value[i] = value[i] / v;
                return V;
            }
            Vector operator*(T v) const
            {
                Vector V;
                for (int i = 0; i < N; i++)
                    V.value[i] = v * value[i];
                return V;
            }
            Vector& operator*=(const Vector& A)
            {
                for (int i = 0; i < N; i++)
                    value[i] *= A.value[i];
                return *this;
            }
            Vector& operator+=(const Vector& A)
            {
                for (int i = 0; i < N; i++)
                    value[i] += A.value[i];
                return *this;
            }
            Vector& operator/=(const Vector& A)
            {
                for (int i = 0; i < N; i++)
                {
                    try
                    {
                        if (A.value[i] == 0)throw std::runtime_error("Division by zero!");
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
            Vector& operator/=(T v)
            {
                try
                {
                    if (v == 0)throw std::runtime_error("Division by zero!");
                }
                catch (const std::exception& e)
                {
                    std::cerr << e.what() << '\n';
                    exit(-1);
                }
                for (int i = 0; i < N; i++)
                    value[i] /= v;
                return *this;
            }
            Vector& operator-=(const Vector& A)
            {
                for (int i = 0; i < N; i++)
                    value[i] -= A.value[i];
                return *this;
            }
            Vector& operator*=(T v)
            {
                for (int i = 0; i < N; i++)
                    value[i] *= v;
                return *this;
            }
            ~Vector() = default;
        private:
            T value[N];
    };
    template<typename T, unsigned N>
    Vector<T, N> operator*(T v, const Vector<T, N>& A)
    {
        return A * v;
    }
    template<typename T, unsigned N>
    Vector<T, N> operator/(T v, const Vector<T, N>& A)
    {
        return Vector<T, N>(v) / A;
    }

    template<typename T>
    class Vector<T, 2u>
    {
            static_assert(std::is_arithmetic_v<T>, "ERROR: must be arithmetic type.");
        public:
            T x = static_cast<T>(0), y = static_cast<T>(0);
            friend Vector operator*<T>(const T&, const Vector&);
            friend Vector operator/<T>(const T&, const Vector&);
            Vector(T _x, T _y) : x(_x), y(_y) {}
            Vector(Vector&& V) noexcept : x(std::move(V.x)), y(std::move(V.y)) {}
            Vector& operator=(Vector&& V) noexcept
            {
                x = std::move(V.x);
                y = std::move(V.y);
                return *this;
            }
            Vector& operator=(const Vector& V)
            {
                if(&V == this) return *this;
                x = V.x;
                y = V.y;
                return *this;
            }
            Vector(const Vector& V) : x(V.x), y(V.y) {}
            explicit Vector(T v) : x(v), y(v) {}
            Real norm() const { return std::sqrt(x * x + y * y); }
            Vector& normalize()
            {
                Real length = this->norm();
                if (isZero(length))return *this;
                else
                {
                    x /= length;
                    y /= length;
                    return *this;
                }
            }
            Vector normalized() const
            {
                Real length = this->norm();
                if (isZero(length))return Vector();
                else return (*this) / length;
            }
            T dot(const Vector& A) const { return x * A.x + y * A.y; }
            T& operator[](int i)
            {
                if(i == 0) return x;
                else if(i == 1) return y;
            }
            T operator[](int i) const
            {
                if(i == 0) return x;
                else if(i == 1) return y;
            }
            Vector operator*(const Vector& A) const { return {x * A.x, y * A.y}; }
            Vector operator+(const Vector& A) const { return {x + A.x, y + A.y}; }
            Vector operator-() const { return {-x, -y}; }
            Vector operator/(const Vector& A) const
            {
                try
                {
                    if (A.x == 0 || A.y == 0)throw std::runtime_error("Division by zero!");
                    return {x / A.x, y / A.y};
                }
                catch (const std::exception& e)
                {
                    std::cerr << e.what() << '\n';
                    exit(-1);
                }
            }
            Vector operator-(const Vector& A) const { return {x - A.x, y - A.y}; }
            Vector operator/(T v) const
            {
                try
                {
                    if (isZero(v))throw std::runtime_error("Division by zero!");
                    return {v.x, v.y};
                }
                catch (const std::exception& e)
                {
                    std::cerr << e.what() << '\n';
                    exit(-1);
                }
            }
            Vector operator*(T v) const
            {
                return {v.x * v, v.y * v};
            }
            Vector& operator*=(const Vector& A)
            {
                x *= A.x;
                y *= A.y;
                return *this;
            }
            Vector& operator+=(const Vector& A)
            {
                x += A.x;
                y += A.y;
                return *this;
            }
            Vector& operator/=(const Vector& A)
            {
                try
                {
                    if (isZero(A.x) || isZero(A.y))throw std::runtime_error("Division by zero!");
                    else
                    {
                        x /= A.x;
                        y /= A.y;
                    }
                    return *this;
                }
                catch (const std::exception& e)
                {
                    std::cerr << e.what() << '\n';
                    exit(-1);
                }
            }
            Vector& operator/=(T v)
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
                v.x /= v;
                v.y /= v;
                return *this;
            }
            Vector& operator-=(const Vector& A)
            {
                x -= A.x;
                y -= A.y;
                return *this;
            }
            Vector& operator*=(T v)
            {
                x *= v;
                y *= v;
                return *this;
            }
            ~Vector() = default;
    };

    template<typename T>
    class Vector<T, 3u>
    {
            static_assert(std::is_arithmetic_v<T>, "ERROR: elements of Vector must be arithmetic type.")
        public:
            T x = static_cast<T>(0), y = static_cast<T>(0), z = static_cast<T>(0);
            friend Vector operator*<T>(const T&, const Vector&);
            friend Vector operator/<T>(const T&, const Vector&);
            Vector(T _x, T _y, T _z) : x(_x), y(_y), z(_z) {}
            Vector(Vector&& V) noexcept : x(std::move(V.x)), y(std::move(V.y)), z(std::move(V.z)) {}
            Vector& operator=(Vector&& V) noexcept
            {
                x = std::move(V.x);
                y = std::move(V.y);
                z = std::move(V.z);
                return *this;
            }
            Vector& operator=(const Vector& V)
            {
                if(&V == this) return *this;
                x = V.x;
                y = V.y;
                z = V.z;
                return *this;
            }
            Vector(const Vector& V) : x(V.x), y(V.y), z(V.z) {}
            explicit Vector(T v) : x(v), y(v), z(v) {}
            Real norm() const { return std::sqrt(x * x + y * y + z * z); }
            Vector& normalize()
            {
                Real length = this->norm();
                if (isZero(length))return *this;
                else
                {
                    x /= length;
                    y /= length;
                    z /= length;
                    return *this;
                }
            }
            Vector normalized() const
            {
                Real length = this->norm();
                if (isZero(length))return Vector();
                else return (*this) / length;
            }
            T dot(const Vector& A) const { return x * A.x + y * A.y + z * A.z; }
            T& operator[](int i)
            {
                if(i == 0) return x;
                else if(i == 1) return y;
                else if(i == 2) return z;
            }
            T operator[](int i) const
            {
                if(i == 0) return x;
                else if(i == 1) return y;
                else if(i == 2) return z;
            }
            Vector operator*(const Vector& A) const { return {x * A.x, y * A.y, z * A.z}; }
            Vector operator+(const Vector& A) const { return {x + A.x, y + A.y, z + A.z}; }
            Vector operator-() const { return {-x, -y, -z}; }
            Vector operator/(const Vector& A) const
            {
                try
                {
                    if (isZero(A.x) || isZero(A.y) || isZero(A.z) )throw std::runtime_error("Division by zero!");
                    return {x / A.x, y / A.y, z / A.z};
                }
                catch (const std::exception& e)
                {
                    std::cerr << e.what() << '\n';
                    exit(-1);
                }
            }
            Vector operator-(const Vector& A) const { return {x - A.x, y - A.y, z - A.z}; }
            Vector operator/(T v) const
            {
                try
                {
                    if (isZero(v))throw std::runtime_error("Division by zero!");
                    return {v.x, v.y, v.z};
                }
                catch (const std::exception& e)
                {
                    std::cerr << e.what() << '\n';
                    exit(-1);
                }
            }
            Vector operator*(T v) const { return {v.x * v, v.y * v, v.z * v}; }
            Vector& operator*=(const Vector& A)
            {
                x *= A.x;
                y *= A.y;
                z *= A.z;
                return *this;
            }
            Vector& operator+=(const Vector& A)
            {
                x += A.x;
                y += A.y;
                z += A.z;
                return *this;
            }
            Vector& operator/=(const Vector& A)
            {
                try
                {
                    if (A.x == 0 || A.y == 0 || A.z == 0)throw std::runtime_error("Division by zero!");
                    else
                    {
                        x /= A.x;
                        y /= A.y;
                        z /= A.z;
                    }
                    return *this;
                }
                catch (const std::exception& e)
                {
                    std::cerr << e.what() << '\n';
                    exit(-1);
                }
            }
            Vector& operator/=(T v)
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
                v.x /= v;
                v.y /= v;
                v.z /= v;
                return *this;
            }
            Vector& operator-=(const Vector& A)
            {
                x -= A.x;
                y -= A.y;
                z -= A.z;
                return *this;
            }
            Vector& operator*=(T v)
            {
                x *= v;
                y *= v;
                z *= v;
                return *this;
            }
            Vector<T, 2> xy() const { return {x, y}; }
            Vector<T, 2> yx() const { return {y, x}; }
            Vector<T, 2> xz() const { return {x, z}; }
            Vector<T, 2> zx() const { return {z, x}; }
            Vector<T, 2> yz() const { return {y, z}; }
            Vector<T, 2> zy() const { return {z, y}; }
            ~Vector() = default;
    };

    using Vec2f = Vector<float, 2u>;
    using Vec3f = Vector<float, 3u>;
    using Vec3b = Vector<bool, 3u>;
    using Vec2 = Vector<Real, 2u>;
    using Vec3 = Vector<Real, 3u>;
    using Vec2i = Vector<int, 2u>;
    using Vec3i = Vector<int, 3u>;
    using Vec2u = Vector<uint, 2u>;
    using Vec3u = Vector<uint, 3u>;
}
#endif
