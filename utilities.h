#ifndef UTILITIES_H
#define UTILITIES_H

#include <iostream>
#include <fstream>
using namespace std;

#include <math.h>
#include <stdlib.h>

// Matrix Model in 0 dim
// Compile with g++ -O3 *.cpp -llapack -lblas -lg2c
// Reads parameters from file parameters

const int L = 1;
const int T = 1;
const int D = 1;
const int NSCALAR = 4; // Change this according to model
const int SITES = L * T;
const int NCOLOR = 3;
const int RANK = NCOLOR * NCOLOR - 1;
const double GAUGETOL = 0.000001;
const double TRACETOL = 1e-8;       // Will be used to test norm()
const double NORM = 0.5 / sqrt(D);
const double PBC = 1.0;

extern double BETA,DT,MASS,TIME,MU;
extern int SWEEPS,GAP,START,THERM,READIN,SEED;
extern int TRAJECTORY_LENGTH;

class Complex
{
private:
  double re, im;
public:
  Complex();
  Complex(double, double);
  double real() const;
  double imag() const;
  double norm();
  void print() const;
  friend ostream& operator<<(ostream&,Complex);
  friend istream& operator>>(istream&,Complex &);
};

inline Complex conjug(const Complex &o1)
{
  return(Complex(o1.real(),-o1.imag()));
}

inline Complex operator +(const Complex &o1, const Complex &o2)
{
  return(Complex(o1.real()+o2.real(),o1.imag()+o2.imag()));
}

inline Complex operator -(const Complex &o1, const Complex &o2)
{
  return(Complex(o1.real()-o2.real(),o1.imag()-o2.imag()));
}

inline Complex operator *(const Complex &o1, const Complex &o2)
{
  return(Complex(o1.real()*o2.real()-o1.imag()*o2.imag(),
                 o1.real()*o2.imag()+o1.imag()*o2.real()));
}

inline Complex operator *(const Complex &o1, const double o2)
{
  return(Complex(o1.real()*o2,o1.imag()*o2));
}

inline Complex operator *(const double o1, const Complex &o2)
{
  return(Complex(o2.real()*o1,o2.imag()*o1));
}

Complex operator /(const Complex &, const Complex &);
Complex pow(const Complex &, const int);

class Umatrix
{
private:
  Complex mat[NCOLOR][NCOLOR];
public:
  Umatrix();
  Umatrix(int);
  Umatrix(Complex [NCOLOR][NCOLOR]);
  Complex get(int,int) const;
  void set(int,int,const Complex);
  void print();
  friend ostream& operator<<(ostream &, Umatrix);
  friend istream& operator>>(istream &, Umatrix &);
};

Umatrix operator +(const Umatrix &o1, const Umatrix &o2);
Umatrix operator -(const Umatrix &o1, const Umatrix &o2);
Umatrix operator *(const Umatrix &, const Umatrix &);
Umatrix operator *(const Umatrix &, const Complex &);
Umatrix operator *(const Complex &, const Umatrix &);
Umatrix operator *(const Umatrix &, const double);
Umatrix operator *(const double, const Umatrix &);
Umatrix comm(const Umatrix &, const Umatrix &);
Umatrix exp(const Umatrix &u);
Umatrix Adj(const Umatrix &u);
Complex Tr(const Umatrix &);
Umatrix real_gaussian_Umatrix();

class Lattice_Vector
{
private:
  int coords[D];
public:
  Lattice_Vector();
  Lattice_Vector(int);
  void set(int, int);
  int get(int) const;
  void print() const;
};

Lattice_Vector operator +(const Lattice_Vector &x, const Lattice_Vector &y);
Lattice_Vector operator -(const Lattice_Vector &x, const Lattice_Vector &y);
Lattice_Vector operator -(const Lattice_Vector &x);
double BC(const Lattice_Vector &x, const Lattice_Vector &y);
int loop_over_lattice(Lattice_Vector &, int &);

class Site_Field
{
private:
  Umatrix points[SITES];
public:
  Site_Field();
  Site_Field(int);
  Umatrix get(const Lattice_Vector &) const;
  void set(const Lattice_Vector &, const Umatrix &);
  void print();
};

Site_Field Adj(const Site_Field &);

Site_Field operator +(const Site_Field &, const Site_Field &);
Site_Field operator -(const Site_Field &, const Site_Field &);
Site_Field operator *(const double, const Site_Field &);
Site_Field operator *(const Complex &, const Site_Field &);
Umatrix operator *(const Site_Field &, const Site_Field &);
Site_Field comm(const Site_Field &, const Site_Field &);

extern Umatrix Lambda[RANK];

double gasdev();
#endif
