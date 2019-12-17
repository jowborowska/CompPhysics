#ifndef VEC3_H
#define VEC3_H
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
using namespace std;

class vec3
{
public:
    vec3();
    vec3(double x, double y, double z);
    double lengthSquared();
    double length();

    // Actions
    void zeros();
    vec3 cross(vec3 otherVector);
    double dot(vec3 otherVector);
    void normalize();
    vec3 normalized();

    // Getters and setters
    double x() const { return components[0]; }
    double y() const { return components[1]; }
    double z() const { return components[2]; }
    void setX(double x) { components[0] = x; }
    void setY(double y) { components[1] = y; }
    void setZ(double z) { components[2] = z; }

    // Convenience functions
    void print();
    void print(std::string name);
    friend std::ostream& operator<<(std::ostream& os, const vec3& myVector); // Allows cout << myVector << endl;

    // Operators
    double &operator()(int index) { return components[index]; } // Allows access like myVector(0)
    double &operator[](int index) { return components[index]; } // Allows access like myVector[0]
    vec3 &operator+=(double rhs); // Componentwise addition with scalar
    vec3 &operator+=(vec3 rhs);   // Componentwise addition with vector
    vec3 &operator*=(double rhs); // Componentwise multiplication with scalar
    vec3 &operator*=(vec3 rhs);   // Componentwise multiplicationwith vector
    vec3 &operator-=(double rhs); // Componentwise subtraction with scalar
    vec3 &operator-=(vec3 rhs);   // Componentwise subtraction with vector
    vec3 &operator/=(double rhs); // Componentwise division with scalar
    vec3 &operator/=(vec3 rhs);   // Componentwise division with vector
private:
    double components[3];
};

inline vec3 operator+(vec3 lhs, double rhs) {
    lhs += rhs;
    return lhs;
}

inline vec3 operator+(double lhs, vec3 rhs) {
    rhs += lhs;
    return rhs;
}

inline vec3 operator+(vec3 lhs, vec3 rhs) {
    lhs += rhs;
    return lhs;
}

inline vec3 operator-(vec3 lhs, double rhs) {
    lhs -= rhs;
    return lhs;
}

inline vec3 operator-(double lhs, vec3 rhs) {
    rhs -= lhs;
    return rhs;
}

inline vec3 operator-(vec3 lhs, vec3 rhs) {
    lhs -= rhs;
    return lhs;
}


inline vec3 operator*(vec3 lhs, double rhs) {
    lhs *= rhs;
    return lhs;
}

inline vec3 operator*(double lhs, vec3 rhs) {
    rhs *= lhs;
    return rhs;
}

inline vec3 operator*(vec3 lhs, vec3 rhs) {
    lhs *= rhs;
    return lhs;
}


inline vec3 operator/(vec3 lhs, double rhs) {
    lhs /= rhs;
    return lhs;
}

inline vec3 operator/(double lhs, vec3 rhs) {
    rhs /= lhs;
    return rhs;
}

inline vec3 operator/(vec3 lhs, vec3 rhs) {
    lhs /= rhs;
    return lhs;
}

#endif // VEC3_H


vec3::vec3()
{
    zeros();
}

vec3::vec3(double x, double y, double z)
{
    components[0] = x;
    components[1] = y;
    components[2] = z;
}

void vec3::print()
{
    // Will print matlab syntax vector. Output will be like: [2.09, 5.3, 9.1];
    cout << "[" << components[0] << ", " << components[1] << ", " << components[2] << "]" << endl;
}

void vec3::print(string name)
{
    // Will print matlab syntax vector with a name. Output will be like: A = [2.09, 5.3, 9.1];
    cout << name << " = ";
    print();
}

vec3 vec3::cross(vec3 otherVector)
{
    return vec3(y()*otherVector.z()-z()*otherVector.y(), z()*otherVector.x()-x()*otherVector.z(), x()*otherVector.y()-y()*otherVector.x());
}

double vec3::dot(vec3 otherVector)
{
    return otherVector[0]*components[0] + otherVector[1]*components[1] + otherVector[2]*components[2];
}

void vec3::normalize()
{
    double length = this->length();
    if(length > 0) {
        components[0] /= length;
        components[1] /= length;
        components[2] /= length;
    }
}

vec3 vec3::normalized()
{
    vec3 newVector = *this;
    newVector.normalize();
    return newVector;
}

double vec3::lengthSquared()
{
    // Returns the square of the length (or norm) of the vector
    return components[0]*components[0]+components[1]*components[1]+components[2]*components[2];
}

double vec3::length()
{
    // Returns the length (or norm) of the vector
    return sqrt(lengthSquared());
}

void vec3::zeros()
{
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
}

vec3 &vec3::operator+=(double rhs)
{
    components[0] += rhs;
    components[1] += rhs;
    components[2] += rhs;
    return *this;
}

vec3 &vec3::operator+=(vec3 rhs)
{
    components[0] += rhs[0];
    components[1] += rhs[1];
    components[2] += rhs[2];
    return *this;
}

vec3 &vec3::operator*=(double rhs)
{
    components[0] *= rhs;
    components[1] *= rhs;
    components[2] *= rhs;
    return *this;
}

vec3 &vec3::operator*=(vec3 rhs)
{
    components[0] *= rhs[0];
    components[1] *= rhs[1];
    components[2] *= rhs[2];
    return *this;
}

vec3 &vec3::operator-=(double rhs)
{
    components[0] -= rhs;
    components[1] -= rhs;
    components[2] -= rhs;
    return *this;
}

vec3 &vec3::operator-=(vec3 rhs)
{
    components[0] -= rhs[0];
    components[1] -= rhs[1];
    components[2] -= rhs[2];
    return *this;
}

vec3 &vec3::operator/=(double rhs)
{
    components[0] /= rhs;
    components[1] /= rhs;
    components[2] /= rhs;
    return *this;
}

vec3 &vec3::operator/=(vec3 rhs)
{
    components[0] /= rhs[0];
    components[1] /= rhs[1];
    components[2] /= rhs[2];
    return *this;
}

std::ostream &operator<<(std::ostream &os, const vec3 &myVector) // Allows cout << myVector << endl;
{
    os << "[" << myVector.x() << ", " << myVector.y() << ", " << myVector.z() << "];";
    return os;
}
