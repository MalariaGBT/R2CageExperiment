#include "head.h"

// random number seed
std::random_device rd;
std::mt19937 twister(rd());

// global structs
Pars pa; 
// output files
std::ostringstream os1;
std::ofstream local_data; // population values of males with each genotype in each patch for each simulated day
std::ofstream global_data; // total population values (in all patches) of males with each genotype for each simulated day
std::ofstream coord_list; // coordinate list (x, y) of the patches


int main()
{	
	// input parameters
	std::cin >> pa.phi;
	std::cin >> pa.omegaH;
	std::cin >> pa.omegaR;
	std::cin >> pa.play; 
	std::cin >> pa.e; 
	std::cin >> pa.nj; 
	std::cin >> pa.mate; 
	std::cin >> pa.gen_size; 
	std::cin >> pa.num_generations; 
	std::cin >> pa.num_reps; 
	std::array<int,2*num_genotypes> init;
	std::array<int,2*num_genotypes> state;
	for (int i=0; i < 2*num_genotypes; ++i)std::cin >> init[i]; 
	

	for (int i=0; i < num_genotypes; ++i) {pa.omega[i]=0;}
	pa.omega[0]=1;
	pa.omega[1]=pa.omegaH;
	pa.omega[3]=pa.omegaR;
	double eplus=(1+pa.e)*0.5;
	double eminus=(1-pa.e)*0.5;
	double gamma=pa.nj*(1-pa.e)*0.5;
	double theta=(1-pa.nj)*(1-pa.e)*0.5;

	double inherit[6][6][6] ={{{1,0,0,0,0,0},{theta,eplus,0,gamma,0,0}, {0,1,0,0,0,0},{0.5,0,0,0.5,0,0}, {0,0,0,1,0,0},{0,0.5,0,0.5,0,0}}, {{theta,eplus,0,gamma,0,0},{theta*theta,2*theta*eplus,eplus*eplus,2*theta*gamma,gamma*gamma,2*eplus*gamma}, {0,theta,eplus,0,0,gamma},{ theta*0.5,eplus*0.5,0,( theta*0.5+(gamma*0.5)),gamma*0.5,eplus*0.5}, {0,0,0,theta,gamma,eplus},{0, theta*0.5,eplus*0.5, theta*0.5,gamma*0.5,(eplus*0.5+gamma*0.5)}}, {{0,1,0,0,0,0},{0,eminus,eplus,0,0,0}, {0,0,1,0,0,0},{0,0.5,0,0,0,0.5}, {0,0,0,0,0,1},{0,0,0.5,0,0,0.5}}, {{0.5,0,0,0.5,0,0},{ theta*0.5,eplus*0.5,0,(gamma*0.5+ theta*0.5),gamma*0.5,eplus*0.5}, {0,0.5,0,0,0,0.5},{0.25,0,0,0.5,0.25,0},{0,0,0,0.5,0.5,0}, {0,0.25,0,0.25,0.25,0.25}}, {{0,0,0,1,0,0},{0,0,0,eminus,0,eplus}, {0,0,0,0,0,1},{0,0,0,0.5,0.5,0}, {0,0,0,0,1,0},{0,0,0,0,0.5,0.5}}, {{0,0.5,0,0.5,0,0},{0,eminus*0.5,eplus*0.5,eminus*0.5,0,eplus*0.5}, {0,0,0.5,0,0,0.5},{0,0.25,0,0.25,0.25,0.25}, {0,0,0,0,0.5,0.5},{0,0,0.25,0,0.25,0.5}}};

	
	for (int i=0; i < num_genotypes; ++i) 
	{
		for (int j=0; j < num_genotypes; ++j) 
		{
		for (int k=0; k < num_genotypes; ++k) { pa.in[i][j][k]=inherit[i][j][k];
	       
	
//		std::cout<<" in0 "<<i<<"    "<<j<<"   "<<k<<"   "<<pa.in[i][j][k]<<std::endl;
//		std::cout<<pa.in[i][j][k]<<std::endl;
	
	
		}
		}
	}


	for(int rep=0;rep<pa.num_reps;++rep)
	{
		for (int i=0; i < 2*num_genotypes; ++i)state[i]= init[i]; 
		for (int i=0; i < 2*num_genotypes; ++i)std::cout << state[i]<<"    "; 
		std::cout<<std::endl;

		for(int i=0;i<pa.num_generations;++i)
		{
		OneGen(state);
		for (int i=0; i < 2*num_genotypes; ++i)std::cout << state[i]<<"    "; 
		std::cout<<std::endl;
		};
	};

return 0;}



void OneGen(std::array<int,2*num_genotypes> &state) {

	int totV=0;int totM=0;int totE=0;
	int FLay;
	for (int i=0;i<num_genotypes;++i){totV+=state[i];totM+=state[i+6];};
	if (totV==0 ||totM==0) return;
	std::array<int, num_genotypes> v;
	std::array<int, num_genotypes> eggs;
	std::array<int, num_genotypes> Feggs;
	std::array<int, num_genotypes> Meggs;
	std::array<std::array<int, num_genotypes>, num_genotypes> F;
	std::vector<double> m(num_genotypes);
	std::vector<double> in_probs(num_genotypes);
	std::vector<double> eggprobs(2*num_genotypes);
	std::vector<int> v_c;
	std::vector<int> eg;
	std::vector<int> take_eg;
	for (int i=0; i < num_genotypes; ++i) 
		{ 
		v[i]=0;
		eggs[i]=0;
		Feggs[i]=0;
		Meggs[i]=0;
		for (int j=0; j < num_genotypes; ++j) { F[i][j]=0; }
		};
	for (int i=0; i < num_genotypes; ++i) 
	{
		m[i] = state[i+num_genotypes]/(1.0*totM);
	}
//	std::cout<<"male probs  "<<std::endl;
//	for (int i=0; i < num_genotypes; ++i)std::cout << m[i]<<"    "; 
//	std::cout<<std::endl;
	for (int i=0; i < num_genotypes; ++i) 
	{
		v[i] = random_binomial(state[i],pa.mate); // how many V will mate
		if (v[i] > 0) 
		{
		v_c = random_multinomial(v[i], m); // how many V with given genotype will carry each of the male genotypes
		for (int j=0; j < num_genotypes; j++) { F[i][j] = v_c[j]; }
		};
	}
	int rb;
	for (int i=0; i < num_genotypes; ++i) 
	{
		for (int j=0; j < num_genotypes; ++j) 
		{
		rb=random_binomial(F[i][j],pa.play);//how many Fij females will lay eggs
		FLay=random_poisson(pa.phi*pa.omega[i]*rb);//how many eggs will be layed from Fij females
		if(FLay>0)
		{
		for (int k=0; k < num_genotypes; ++k){ in_probs[k]=pa.in[i][j][k];
		}
		eg = random_multinomial(FLay,in_probs);
		for (int k=0; k < num_genotypes; ++k){ eggs[k]+=eg[k];totE+=eg[k];}
		}
		}
	}
	if(totE==0) return;
	for (int i=0; i < num_genotypes; ++i) 
	{
	Feggs[i]=random_binomial(eggs[i],0.5);Meggs[i]=eggs[i]-Feggs[i];
	eggprobs[i]=Feggs[i]/(1.0*totE);
	eggprobs[i+num_genotypes]=Meggs[i]/(1.0*totE);
	};	
	take_eg = random_multinomial(pa.gen_size,eggprobs);
	for (int i=0; i < 2*num_genotypes; ++i) {state[i]=take_eg[i];}
	return;
	}









// Returns a random draw (non-negative integer) from the Poisson distribution with mean lambda (using normal distribution approximation when lambda > 30)
 int random_poisson(double lambda) 
{
	 int result;
	if (lambda < 1e-5) { 
		result = 0;
	}
	else if (lambda > 30) {
		// use normal approximation	
		std::normal_distribution<> dist(lambda, std::sqrt(lambda)); // distribution(mean, standard deviation)
		int x = std::round(dist(twister));
		result = std::max(0, x);
	}
	else {
		// sample poisson directly
		std::poisson_distribution<> dist(lambda); // distribution(mean)
		result = dist(twister);
	}

	return result;
}

// Returns a random draw (non-negative integer) from the Binomial distribution B(N,p)
// Uses Normal and Poisson distribution approximations for large N
 int random_binomial( int n, double p) 
{
	 int result;
	if (n*p > 10 && n*(1 - p) > 10) {
		// use normal approximation
		std::normal_distribution<> dist(n*p, std::sqrt(n*p*(1 - p))); // distribution(mean, standard deviation)
		 int x = std::round(dist(twister));
		if (x<0) x=0;
		if (x>n) x=n;
		result = x;
	}
	else if ((n > 20 && p < 0.05) || (n > 100 && n*p < 10)) {
		// use Poisson approximation
		result = random_poisson(n*p);
	}
	else if ((n > 20 && p > 0.95) || (n > 100 && n*(1-p) < 10)) {
		// use Poisson approximation
		result = n - random_poisson(n*(1 - p));
	}
	else {
		// use binomial distribution directly
		std::binomial_distribution<> dist(n, p);
		result = dist(twister);
	}

	return result;
}

// Returns a vector of outcomes from a random draw of the Multinomial distribution with N trials where each trial has a vector of probabilities <probs>
std::vector<int> random_multinomial(int n, const std::vector<double>& probs) 
{
	int num_outcomes = probs.size();
	double sum_p = 0.0;
	for (int i = 0; i < num_outcomes; ++i) {
		sum_p += probs[i];
	}

	 int n_used = n;
	std::vector< int> result(num_outcomes, 0);
	for(int i=0; i < num_outcomes; ++i) {
		if (n_used > 0) {
			result[i] = random_binomial(n_used, probs[i] / sum_p);
			sum_p -= probs[i];
			n_used -= result[i];
		}
		else {
			result[i] = 0;
		}
	}

	return result;
}
