#include<iostream>
#include<cmath>
#include<iomanip>
#include<random>
#include<map>
#include<fstream>
#include<stdio.h>
#include<math.h>
#include<vector>

const double cosPiOver8 = 0.923879532511286;
const double pi = 3.141592653589793238;

std::map<double, int> m;
std::vector<double> f_values;
std::vector<double> density_values;

/*
double generate(int L) {
  std::random_device rd;
  std::mt19937 gen(rd());  
  std::uniform_real_distribution<> dis(0.0, 1.0);
  double x = dis(gen);
  //std::cout << x << std::endl;
  
  double y = pow(8*L*cosPiOver8/(3*pi*x)-1, 0.375);
  return y;
}
*/

//probability power law distribution
double p(double L, double a) {
  return 1/(1+pow(L*a, 8.0/3.0));
}

//calculates variance given a vector and average
double calculate_variance(std::vector<double> v, double a) {
  double res = 0, res_c = 0, res_t, res_y;
  for(int i=0; i<v.size(); i++) {
    res_y = (v[i]-a)*(v[i]-a) - res_c;
    res_t = res + res_y;
    res_c = (res_t - res) - res_y;
    res = res_t;
  }

  if(v.size() > 1) {
    res/=(v.size()-1);
  }
  return res;
}

//returns delta f value for given x, y
double sum(double x, double y, int steps, double variance, double L) {
  std::random_device rd;
  std::mt19937 gen(rd());  
  std::uniform_real_distribution<> dis(0.0, 1.0);

  double k = 0.05*L;
  double res = 0;
  double ratio = pow(100, 1.0/steps);
  double delta_k = k;
  double prev_k = k;
  //return mult;
  double c, temp = 0;

  //find normalization constant
  for(int i=0; i<steps; i++) {
    temp += (2 * pi * k * delta_k * p(L, k));
    prev_k = k;
    k *= ratio;
    delta_k = k - prev_k;
  }

  c = variance/temp;
  //return c;  
  
  k = 0.05*L;
  delta_k = k;
  prev_k = k;

  //sum over wave modes
  for(int i=0; i<steps; i++) {
    //generate random angles between 0 and 2pi	  
    double theta = dis(gen);
    theta *= (2*pi);
    double phi = dis(gen);
    phi *= (2*pi);
    
    //calculate
    double f1 = sqrt(c * 4 * pi * k * delta_k * p(L, k));
    double f2 = cos(k * x * cos(theta) + k * y * sin(theta) + phi);
    res += (f1 * f2);
    //f_values.push_back(f1*f2);
    prev_k = k;
    k *= ratio;
    delta_k = k - prev_k;
  }
  return res;  
}

int main()
{
  freopen("output.txt","w",stdout);
  //std::cout << sum(1.0, 1.0, 200, 0.25403, 3); return 0;
  double f_tot = 0, f_c = 0, f_y, f_t, f_i;
  double density_tot = 0, d_c = 0, d_y, d_t, d_i; //variables for Kahan summation
  double L = 4;
  for(double i=-5*L; i<=5*L; i+=0.25*L) {
    for(double j=-5*L; j<=5*L; j+=0.25*L) {
	double add = sum(i, j, 100, 0.25403, L);
//  	double f_var = calculate_variance(f_values, add/100);
//  	std::cout << "\naverage value and variance of f: "
//	       	<< add/100 << " " << f_var << "\n";
	std::cout << i << "\t" << j << "\t" << exp(-0.11957+add) << "\n";
	
	f_y = add - f_c;
	f_t = f_tot + f_y;
	f_c = (f_t - f_tot) - f_y;
	f_tot = f_t;
	
	d_y = exp(-0.11957+add) - d_c;
	d_t = density_tot + d_y;
	d_c = (d_t - density_tot) - d_y;
	density_tot = d_t;

	density_values.push_back(exp(-0.11957+add));
	f_values.push_back(add);
	//	f_values.clear();
    }
  }
//  for(int i=0; i<f_values.size(); i++) {
//    std::cout << f_values[i] << "\n";
//  }
/*
  double f_avg = f_tot/(f_values.size());
  double f_var = calculate_variance(f_values, f_avg);
  std::cout << "\naverage value and variance of f: " << f_avg << " " << f_var;

  double density_avg = density_tot/(density_values.size());
  double density_var = calculate_variance(density_values, density_avg);
  std::cout << "\naverage value and variance of density: " << density_avg << " " << density_var;
  std::cout << "\n0.2*n0*n0: " << 0.2*density_avg*density_avg;
*/
  fclose(stdout);
  
  return 0;
}
