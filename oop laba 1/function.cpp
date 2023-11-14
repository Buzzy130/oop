#include "function.h"

//--------------------------------------------------------------------------huber_distribution--------------------------------------------------------------------------------//
HuberDistribution* init_huber_distribution(double v, double K, double scale, double shift)
{
	HuberDistribution* HD = new HuberDistribution();
	HD->v = v;
	HD->K = K;
	HD->scale = scale;
	HD->shift = shift;
	return HD;
}
double Huber(double x, HuberDistribution* HD)
{
	if (abs((x - HD->shift) / HD->scale) <= HD->v)
	{
		return (1. / (sqrt(2. * M_PI) * HD->K) * exp(-pow((x - HD->shift) / HD->scale, 2.) / 2.)) / HD->scale;
	}
	if (abs((x - HD->shift) / HD->scale) > HD->v)
	{
		return (1. / (sqrt(2. * M_PI) * HD->K) * exp(pow(HD->v, 2.) / 2. - HD->v * abs((x - HD->shift) / HD->scale))) / HD->scale;
	}
}
double phi(double x)//Ф(x)
{
	return 0.5 * (1. + erf(x / sqrt(2.)));
}
double phi_lower(double x)//ф(x)
{
	return 1. / sqrt(2. * M_PI) * exp(-1. / 2. * pow(x, 2.));
}
double Mksi_huber(HuberDistribution* HD)//мат ожидание
{
	return HD->shift;
}
double Dksi_huber(HuberDistribution* HD)//дисперсия
{
	return 1. + 2. * phi_lower(HD->v) * (pow(HD->v, 2.) + 2.) / (pow(HD->v, 3.) * HD->K);
}
double asymmetry_huber(HuberDistribution* HD)//ассиметрия
{
	return 0.;
}
double kurtosis_huber(HuberDistribution* HD)//коэфф эксцесса
{
	return (3. * (2. * phi(HD->v) - 1.) + 2. * phi_lower(HD->v) * (24. / pow(HD->v, 5.) + 24. / pow(HD->v, 3.) + 12. /
		HD->v + HD->v)) / (pow(Dksi_huber(HD), 2.) * HD->K) - 3.;
}
double P(HuberDistribution* HD)//вероятности попадания в центральный интервал 
{
	return (2. * phi(HD->v) - 1.) / HD->K;
}
double K(double v)//значение K зная V
{
	return 2. / v * phi_lower(v) + 2. * phi(v) - 1.;
}
//-----------------------------------------------------//
double algorithm(HuberDistribution* HD)
{
	std::random_device rd;
	std::default_random_engine gen(rd());
	std::uniform_real_distribution<> d(0, 1);
	//шаг 1
	double r1 = d(gen);
	if (r1 <= P(HD))
	{
		//шаг 2
		double r2, r3, x1;
		do {
			r2 = d(gen);
			r3 = d(gen);
			x1 = sqrt(-2 * log(r2)) * cos(2 * M_PI * r3);
			//double x1 = sqrt(-2 * log(r2)) * sin(2 * M_PI * r3)
		} while (!(-HD->v <= x1 && x1 <= HD->v)); //шаг 3
		return x1 * HD->scale + HD->shift;
	}
	else
	{
		//шаг 4
		double r4 = d(gen);
		double x2 = HD->v - log(r4) / HD->v;
		//шаг 5
		return r1 < (1 + P(HD)) / 2 ? x2 * HD->scale + HD->shift : -x2 * HD->scale + HD->shift;
	}
}
//--------------------------------------------------------------------------mixture_distribution--------------------------------------------------------------------------------//
Mixture* init_mixture(double p, double v1, double v2, double scale1, double scale2, double shift1, double shift2)
{
	Mixture* MD = new Mixture();
	HuberDistribution* HD1 = new HuberDistribution();
	HuberDistribution* HD2 = new HuberDistribution();
	double K1 = K(v1);
	double K2 = K(v2);
	HD1->v = v1;
	HD1->K = K1;
	HD1->scale = scale1;
	HD1->shift = shift1;
	HD2->v = v2;
	HD2->K = K2;
	HD2->scale = scale2;
	HD2->shift = shift2;
	MD->p = p;
	MD->HD1 = HD1;
	MD->HD2 = HD2;
	return MD;
}

double  mixture_ditribution(double x, Mixture* MD)//функции плотности
{
	return (1 - MD->p) * Huber(x, MD->HD1) + MD->p * Huber(x, MD->HD2);
}
double 𝑀ksi_mixture(Mixture* MD)
{
	return (1 - MD->p) * Mksi_huber(MD->HD1) + MD->p * Mksi_huber(MD->HD2);
}
double Dksi_mixture(Mixture* MD)//шибка
{
	return (1 - MD->p) * (pow(Mksi_huber(MD->HD1), 2) + Dksi_huber(MD->HD1)) + MD->p *
		(pow(Mksi_huber(MD->HD2), 2) + Dksi_huber(MD->HD2)) - pow(𝑀ksi_mixture(MD), 2);
}
double asymmetry_mixture(Mixture* MD)
{
	return (1 / pow(Dksi_mixture(MD), 3 / 2)) *
		((1 - MD->p) *
			(pow(Mksi_huber(MD->HD1) - 𝑀ksi_mixture(MD), 3) + 3 *
				(Mksi_huber(MD->HD1) - 𝑀ksi_mixture(MD)) * Dksi_huber(MD->HD1) +
				pow(Dksi_huber(MD->HD1), 3 / 2) * asymmetry_huber(MD->HD1)) +
			(MD->p) *
			(pow(Mksi_huber(MD->HD2) - 𝑀ksi_mixture(MD), 3) + 3 *
				(Mksi_huber(MD->HD2) - 𝑀ksi_mixture(MD)) * Dksi_huber(MD->HD2) +
				pow(Dksi_huber(MD->HD2), 3 / 2) * asymmetry_huber(MD->HD2)
				));
}
double kurtosis_mixture(Mixture* MD)
{
	return (1 / pow(Dksi_mixture(MD), 2)) * (
		(1 - MD->p) * (pow(Mksi_huber(MD->HD1) - 𝑀ksi_mixture(MD), 4) +
			6 * pow(Mksi_huber(MD->HD1) - 𝑀ksi_mixture(MD), 2) * Dksi_huber(MD->HD1) +
			4 * (Mksi_huber(MD->HD1) - 𝑀ksi_mixture(MD)) * pow(Dksi_huber(MD->HD1), 3 / 2) * asymmetry_huber(MD->HD1) +
			pow(Dksi_huber(MD->HD1), 2) * kurtosis_huber(MD->HD1)) +
		(MD->p) * (pow(Mksi_huber(MD->HD2) - 𝑀ksi_mixture(MD), 4) +
			6 * pow(Mksi_huber(MD->HD2) - 𝑀ksi_mixture(MD), 2) * Dksi_huber(MD->HD2) +
			4 * (Mksi_huber(MD->HD2) - 𝑀ksi_mixture(MD)) * pow(Dksi_huber(MD->HD2), 3 / 2) * asymmetry_huber(MD->HD2) +
			pow(Dksi_huber(MD->HD2), 2) * kurtosis_huber(MD->HD2))) - 3;
}
//--------------------------------------------------------------------------empirical_distribution--------------------------------------------------------------------------------//
double Mn(int n, vector<double> x_selection)//мат ожидание
{
	double sum = 0;
	for (double& x : x_selection)
		sum += x;
	return sum / n;
}

double Dn(int n, vector<double> x_selection)
{
	double Mn_ksi = Mn(n, x_selection);
	double sum = 0;
	for (double& x : x_selection)
		sum += pow((x - Mn_ksi), 2);
	return sum / n;
}

double asymmetry_empirical(int n, vector<double> x_selection)
{
	double Dn_ksi = Dn(n, x_selection);
	double Mn_ksi = Mn(n, x_selection);
	double sum = 0;
	for (double& x : x_selection)
		sum += pow((x - Mn_ksi), 3);
	return sum / (n * pow(Dn_ksi, 3 / 2));
}

double kurtosis_empirical(int n, vector<double> x_selection)
{
	double Dn_ksi = Dn(n, x_selection);
	double Mn_ksi = Mn(n, x_selection);
	double sum = 0;
	for (double& x : x_selection)
		sum += pow((x - Mn_ksi), 4);
	return (sum / (n * pow(Dn_ksi, 2))) - 3;
}


double empirical(int n, double x, vector<double> x_selection)
{
	vector<double> f_s;
	int k = (int)trunc(log2((double)n)) + 1;//кол-во промежутков
	double x_min = *min_element(x_selection.begin(), x_selection.end());
	double x_max = *max_element(x_selection.begin(), x_selection.end());
	double delta = (1 / (double)k) * (x_max - x_min);//ширина каждого интервала
	for (int i = 0; i < k; i++)
	{
		if (x_min + delta * i <= x && x < x_min + delta * (i + 1))
		{
			int n_i = count_if(x_selection.begin(), x_selection.end(), //считаем кол во элементов в x_selection
				[i, k, x_min, x_max, delta](double x) //Лямбда-выражение проверяет принадлежит ли x определенному интервалу
				{ return i == k - 1 ? x_min +
				delta * (double)i <= x && x <= x_min + delta * (double)(i + 1) : x_min + delta * (double)i <= x && x < x_min +
				delta * (double)(i + 1); });
			return n_i / (n * delta);
		}
	}
}

vector<double> generate_sequence(int n, HuberDistribution* HD)
{
	vector<double> x_selection;
	for (int i = 0; i < n; i++)
	{
		double x = algorithm(HD);
		x_selection.push_back(x);
	}
	return x_selection;
}