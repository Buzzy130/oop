#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "function.h"
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



void file_output_distribution(int n, vector<double>& x_selection, HuberDistribution* HD)
{
	ofstream xs;
	ofstream fs_theoretical;
	ofstream fs_empirical;
	xs.open("xs.txt");
	fs_theoretical.open("fs_theoretical.txt");
	fs_empirical.open("fs_empirical.txt");
	sort(x_selection.begin(), x_selection.end());
	for (double& x : x_selection)
	{
		xs << x << endl;
		double f_theoretical = Huber(x, HD);
		fs_theoretical << f_theoretical << endl;
		double f_empirical = empirical(n, x, x_selection);
		fs_empirical << f_empirical << endl;
	}
	xs.close();
	fs_theoretical.close();
	fs_empirical.close();
}

void file_output_mixture(int n, vector<double>& x_selection, Mixture* MD)
{
	ofstream xs;
	ofstream fs_theoretical;
	ofstream fs_empirical;
	xs.open("xs.txt");
	fs_theoretical.open("fs_theoretical.txt");
	fs_empirical.open("fs_empirical.txt");
	sort(x_selection.begin(), x_selection.end());
	for (double& x : x_selection)
	{
		xs << x << endl;
		double f_theoretical = mixture_ditribution(x, MD);
		fs_theoretical << f_theoretical << endl;
		double f_empirical = empirical(n, x, x_selection);
		fs_empirical << f_empirical << endl;
	}
	xs.close();
	fs_theoretical.close();
	fs_empirical.close();
}


void general_distribution()
{
	int n = 1000;
	int distribution_params_option;
	int file_option;
	double v, scale, shift, x;
	vector<double> x_selection;
	HuberDistribution* HD;
	cout << "Параметры распределения" << endl << "1) Стандартное распределение" << endl << "2) Другие" << endl;
	cin >> distribution_params_option;
	if (distribution_params_option == 1)
	{
		cout << "Введите v: ";
		cin >> v;
		HD = init_huber_distribution(v, K(v));
		x = 0;
	}
	else
	{
		cout << "Введите v, lambda, mu, x: ";
		cin >> v >> scale >> shift >> x;
		HD = init_huber_distribution(v, K(v), scale, shift);
	}
	x_selection = generate_sequence(n, HD);
	cout << "Для основного распределения:" << endl;
	cout << "Математическое ожидание: " << Mksi_huber(HD) << endl;
	cout << "Дисперсия: " << Dksi_huber(HD) << endl;
	cout << "Коэффициент асимметрии: " << asymmetry_huber(HD) << endl;
	cout << "Коэффициент эксцесса: " << kurtosis_huber(HD) << endl;
	cout << "Значение плотности в точке " << x << ": " << Huber(x, HD) << endl << endl;


	cout << "Для эмпирического распределения:" << endl;
	cout << "Математическое ожидание: " << Mn(n, x_selection) << endl;
	cout << "Дисперсия: " << Dn(n, x_selection) << endl;
	cout << "Коэффициент асимметрии: " << asymmetry_empirical(n, x_selection) << endl;
	cout << "Коэффициент эксцесса: " << kurtosis_empirical(n, x_selection) << endl;
	cout << "Значение плотности в точке " << x << ": " << empirical(n, x, x_selection) << endl;
	file_output_distribution(n, x_selection, HD);
}

void mixture_distribution1()
{
	int n = 1000;
	int file_option;
	double v1, scale1, shift1, v2, scale2, shift2, x, p;
	vector<double> x_selection;
	Mixture* MD;
	cout << "Введите v1, lambda1, mu1, v2, lambda2, mu2, p, x: ";
	cin >> v1 >> scale1 >> shift1 >> v2 >> scale2 >> shift2 >> p >> x;
	MD = init_mixture(p, v1, v2, scale1, scale2, shift1, shift2);
	x_selection = generate_sequence(n, MD->HD1);
	cout << "Для смеси распределений: " << endl;
	cout << "Математическое ожидание: " << 𝑀ksi_mixture(MD) << endl;
	cout << "Дисперсия: " << Dksi_mixture(MD) << endl;
	cout << "Коэффициент асимметрии: " << asymmetry_mixture(MD) << endl;
	cout << "Коэффициент эксцесса: " << kurtosis_mixture(MD) << endl;
	cout << "Значение плотности в точке " << x << ": " << mixture_ditribution(x, MD) << endl << endl;


	cout << "Для эмпирического распределения:" << endl;
	cout << "Математическое ожидание: " << Mn(n, x_selection) << endl;
	cout << "Дисперсия: " << Dn(n, x_selection) << endl;
	cout << "Коэффициент асимметрии: " << asymmetry_empirical(n, x_selection) << endl;
	cout << "Коэффициент эксцесса: " << kurtosis_empirical(n, x_selection) << endl;
	cout << "Значение плотности в точке " << x << ": " << empirical(n, x, x_selection) << endl;
	file_output_mixture(n, x_selection, MD);
}




int main(int argc, char** argv)
{
	setlocale(LC_ALL, "ru");
	int n = 30000;
	int result = Catch::Session().run(argc, argv);
	return result;
	int distribution_option;
	int file_option;
	cout << "Выберите распределение" << endl << "1) Основное распределение" << endl << "2) Смесь распределений" << endl << "3) Выйти" << endl;
	cin >> distribution_option;
	switch (distribution_option)
	{
	case 1:
		general_distribution();
		break;
	case 2:
		mixture_distribution1();
		break;
	}
}