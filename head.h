#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <cstdlib> // for exit function
#include <random>  
#include <cmath> 
#include <numeric>  
#include <algorithm> 
#include <cstdio>
#include <chrono>

const int num_genotypes=6;

// Contains the total population numbers over all space (all patches) for different mosquito types
struct Pars {
	double phi,omegaH,omegaR,play,e,nj,mate; 
	double in[num_genotypes][num_genotypes][num_genotypes];
	double omega[num_genotypes];
	int gen_size,num_generations,num_reps;
};	

int random_poisson(double lambda);
int random_binomial(int n, double p);
std::vector<int> random_multinomial(int n, const std::vector<double>& probs);

void OneGen(std::array<int, 2*num_genotypes> &state); 

