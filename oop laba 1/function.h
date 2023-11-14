#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <numeric>
#include <algorithm>
#include <string>
using namespace std;
struct HuberDistribution
{
	double v;
	double K;
	double scale;
	double shift;
};
HuberDistribution* init_huber_distribution(double v, double K, double scale = 1, double shift = 0);
double Huber(double x, HuberDistribution* HD);
double phi(double x);
double phi_lower(double x);
double Mksi_huber(HuberDistribution* HD);
double Dksi_huber(HuberDistribution* HD);
double asymmetry_huber(HuberDistribution* HD);
double kurtosis_huber(HuberDistribution* HD);
double P(HuberDistribution* HD);
double K(double v);
double algorithm(HuberDistribution* HD);


struct Mixture
{
	double p;
	HuberDistribution* HD1;
	HuberDistribution* HD2;
};
Mixture* init_mixture(double p, double v1, double v2, double scale1 = 1, double scale2 = 1, double shift1 = 0, double shift2 = 0);
double mixture_ditribution(double x, Mixture* MD);
double 𝑀ksi_mixture(Mixture* MD);
double Dksi_mixture(Mixture* MD);
double asymmetry_mixture(Mixture* MD);
double kurtosis_mixture(Mixture* MD);



double Mn(int n, vector<double> x_s);
double Dn(int n, vector<double> x_s);
double asymmetry_empirical(int n, vector<double> x_selection);
double kurtosis_empirical(int n, vector<double> x_selection);
double empirical(int n, double x, vector<double> x_selection);
vector<double> generate_sequence(int n, HuberDistribution* HD);