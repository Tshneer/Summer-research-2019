// Frozen.cpp : This file contains the 'main' function. Program execution begins and ends there.
// THE GOAL OF THIS FILE IS TO HAVE INITIAL CONDITIONS FOLLOWING A BIMODAL DISTRUBTION. Also calculating the order parameter
// Title:Many dispesal fractions for one r and many radio value
//This gets you the snapshots for ratio 0.5 and df greater than 0.5
//

#include "pch.h"

// All information and variables names are based on those in Doi 10.1038/ncomms7664
//Emergent long-range synchronization of oscillating
//ecological populations without external forcing
//described by Ising universality
//Independent study in fall 20// Frozen.cpp : This file contains the 'main' function. Program execution begins and ends there.






#include <iostream>
//#include < omp.h > //http://www.bowdoin.edu/~ltoma/teaching/cs3225-GIS/fall17/Lectures/openmp.html
#include <time.h>
#include <fstream>
#include <sstream>
//#include <ppl.h>
#include <algorithm>                    
#include <cmath> 
#include <math.h> 
#include <cstdlib> 
#include <stdio.h>
//#include <armadillo>
#include <fstream>
#include <time.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <deque>
#include <time.h>
#include <random>
#include <valarray>
using namespace std;
//using namespace arma;


#define L 4  /* system size */
#define N L*L  /* number of spins */
#define D 4  /* number of neighbors */
int nn[N][D]; //neighbor table


int const skip = 0;
double not_chosen;
int const test = 100;
long int const burn = 90; //No data collected over the first burn# of steps
long int const collect = test - burn;
double lattice[N]; //each spot holds a population density between 0 and 1. Not including 0 and 1. 
//long double last_last_lattice[N];
double last_lattice[N];
double last_last_lattice[N];
double neighbors[N]; //Is not directly related to the function neighor on line 82. 
double one = 1;
const  double r = 2.2; //the r for the logistc map: rx(1-x)
//const double radio = 0.5;
const double df_step = 0.1;
const double radio_step = 0;
double ls = 0;
const int di = 12;
const int dr = 1;
double df[di];
double radio = 0.5;
//const int snap_shots = di;upper_limit_of_high_value
//Total number of interations for each local noise strengh

int it;
int itt;



double high_value_in_ricker_two_cycle = 1.5029;
double low_value_in_ricker_two_cycle = 0.4971;

double tolerance_to_be_considered_within_cycle;
double upper_limit_of_high_value;
double lower_limit_of_high_value;
double upper_limit_of_low_value;
double lower_limit_of_low_value;
double mjt[N]; //As named in cited paper in line 1. In the paper it is called the two-cycle amplitude
double mt = 0; //instantaneous synchronization. It is the sum of all mjt over the lattice for one discrete step in time
//double plus_end=N/2; //Amount of patches that start around x+
//double minus_start=N-plus_end;//Amount of patches that start around x-
double mt_store[di][test];//instantenous syncronizations

double mL[di];//Order parameter
double extinct;
double extinct_table[di];

double lattice_store[di][test];

double choice;

void neighbor(void) {   //two dimensional neighbor ... interesting note to remember. In c++ division of ints result in ints. Not floats! That's why this modular arith works. We need this because we have periodic boundary conditions.
	int disp;
	for (int k = 0; k < N; k++) {
		disp = L;
		nn[k][0] = (k / disp)*disp + (k + 1) % disp;  //+x
		nn[k][1] = (k / disp)*disp + (k - 1 + disp) % disp;  //-x
		disp = L * L;
		nn[k][2] = (k / disp)*disp + (k + L) % disp; //+y
		nn[k][3] = (k / disp)*disp + (k - L + disp) % disp; //-y
	}
}
int indices_for_horizontal_stripe_test[L];

int* create_indices_for_horizontal_stripe_test_for_given_size_length(int length_of_side_of_lattice) {
	for (int i = 0; i <= (length_of_side_of_lattice -1); i++) {
		indices_for_horizontal_stripe_test[i] = i;
	}
	return indices_for_horizontal_stripe_test;
}


int indices_for_vertical_stripe_test[L];

int* create_indices_for_vertical_stripe_test_for_given_size_length(int length_of_side_of_lattice) {
	for (int i = 0; i <= (length_of_side_of_lattice -1); i++) {
		indices_for_vertical_stripe_test[i] = 1+i* length_of_side_of_lattice;
	}
	return indices_for_vertical_stripe_test;
}

int indices_for_descending_diagonal_stripe_test[L];

int* create_indices_for_descending_diagonal_stripe_test_for_given_size_length(int length_of_side_of_lattice) {
	for (int i = 0; i <= (length_of_side_of_lattice - 1); i++) {
		indices_for_descending_diagonal_stripe_test[i] = length_of_side_of_lattice - 1 + i * (length_of_side_of_lattice -1);
	}
	return indices_for_descending_diagonal_stripe_test;
}

int converted_test_points[L];

int* convert_test_points_of_Ricker_map_to_integers_for_any_stripe_test(double original_lattice[L],int indices_for_given_stripe_test[L],int length_of_side_of_lattice) {
	
	double test_points_to_be_converted[L];
	
	int index;
	for (int i = 0; i <= (length_of_side_of_lattice - 1); i++) {
		index = indices_for_given_stripe_test[i];
		test_points_to_be_converted[i] = original_lattice[index];
		
	}
	
	for (int i = 0; i <= (length_of_side_of_lattice - 1); i++) {
	
		if (test_points_to_be_converted[i] < lower_limit_of_low_value || test_points_to_be_converted[i] > upper_limit_of_high_value)
		{
			converted_test_points[i] = 0;//Possible values on border of domains
		}
		if (test_points_to_be_converted[i] > upper_limit_of_low_value &&  test_points_to_be_converted[i] < lower_limit_of_high_value)
		{
			converted_test_points[i] = 0;//Possible values on border of domains
		}
		if (test_points_to_be_converted[i] > lower_limit_of_high_value)
		{
			converted_test_points[i] = 1; //Equivalent to spin up in Ising model
		}
		if (test_points_to_be_converted[i] > lower_limit_of_low_value && test_points_to_be_converted[i] < upper_limit_of_low_value)
		{
			converted_test_points[i] = -1;//Equivalent to spin down in Ising model
		}
	
	}
	return converted_test_points;
}
int is_this_a_horizontal_stripe = 0;
int horizontal_stripe_test(int horizontal_stripe_with_converted_points, int length_of_side_of_lattice) {
	int sum = 0;
	for (int i = 0; i <= (length_of_side_of_lattice - 1); i++) {
		sum = sum + horizontal_stripe_with_converted_points[i];
	}
	if ( abs(sum) == length_of_side_of_lattice)
	{
		is_this_a_horizontal_stripe = 1;
	}
	return is_this_a_horizontal_stripe;
}

int is_this_a_vertical_stripe = 0;
int vertical_stripe_test(int vertical_stripe_with_converted_points, int length_of_side_of_lattice) {
	int sum = 0;
	for (int i = 0; i <= (length_of_side_of_lattice - 1); i++) {
		sum = sum + vertical_stripe_with_converted_points[i];
	}
	if (abs(sum) == length_of_side_of_lattice)
	{
		is_this_a_vertical_stripe = 1;
	}
	return is_this_a_vertical_stripe;
}

int is_this_a_descending_diagonal_stripe = 0;
int descending_diagonal_stripe_test(int descending_diagonal_stripe_with_converted_points, int length_of_side_of_lattice) {
	int sum = 0;
	for (int i = 0; i <= (length_of_side_of_lattice - 1); i++) {
		sum = sum + descending_diagonal_stripe_with_converted_points[i];
	}
	if (abs(sum) == length_of_side_of_lattice)
	{
		is_this_a_descending_diagonal_stripe = 1;
	}
	return is_this_a_descending_diagonal_stripe;
}

int is_this_an_ascending_diagonal_stripe = 0;



//int snap = 0; //index of snapshot
//int snap_last = 0;
int main() {

	create_indices_for_horizontal_stripe_test_for_given_size_length(L);
	create_indices_for_vertical_stripe_test_for_given_size_length(L);
	create_indices_for_descending_diagonal_stripe_test_for_given_size_length(L);



	df[0] = 0;
	for (int d = 1; d < di; d++) {
		df[d] = df[d - 1] + df_step;
	}


	clock_t start = clock();
	


	tolerance_to_be_considered_within_cycle = 0.05;

	upper_limit_of_high_value = high_value_in_ricker_two_cycle + tolerance_to_be_considered_within_cycle; //upper lim of x +

	lower_limit_of_high_value = high_value_in_ricker_two_cycle - tolerance_to_be_considered_within_cycle; //lower lim of x +

	upper_limit_of_low_value = low_value_in_ricker_two_cycle + tolerance_to_be_considered_within_cycle; //upper lim of x -

	lower_limit_of_low_value = low_value_in_ricker_two_cycle - tolerance_to_be_considered_within_cycle; //lower lim of x -

	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution_plus(lower_limit_of_high_value, upper_limit_of_high_value);
	std::uniform_real_distribution<double> distribution_minus(lower_limit_of_low_value, upper_limit_of_low_value);
	std::uniform_real_distribution<double> distribution_decide(0, 1);
	neighbor(); //create neighbor table



	for (int d = 0; d < di; d++) //dispersal (df) loop
	{


		//snap = 0;
		//snap_last = 0;
		for (int i = 0; i < N; i++) { //first lattice generation
			choice = distribution_decide(generator);
			if (choice >= radio) {
				lattice[i] = distribution_plus(generator);
			}
			if (choice < radio) {
				lattice[i] = distribution_minus(generator);
			}

		}

		

		convert_test_points_of_Ricker_map_to_integers_for_any_stripe_test(lattice,indices_for_horizontal_stripe_test,L);
		
		for (int i = 0; i < L; i++) {

			cout << converted_test_points[i] << " ";
		}
		cout << "\n";
		convert_test_points_of_Ricker_map_to_integers_for_any_stripe_test(lattice, indices_for_vertical_stripe_test, L);


		for (int i = 0; i < L; i++) {

			cout << converted_test_points[i] << " ";
		}
		cout << "\n";
		convert_test_points_of_Ricker_map_to_integers_for_any_stripe_test(lattice, indices_for_descending_diagonal_stripe_test, L);


		for (int i = 0; i < L; i++) {

			cout << converted_test_points[i] << " ";
		}




		for (long int year = 0; year < test; year++) {


			for (int i = 0; i < (N); i++) {//remembering the last lattice
				last_lattice[i] = lattice[i];

			}


			if (year == test - 1)
			{
				for (int i = 0; i < N; i++)
				{
					lattice_store[d][i] = lattice[i];

				}

			}



			for (int i = 0; i < (N); i++) {
				lattice[i] = last_lattice[i] * exp(r*(1 - last_lattice[i]));

			}

			for (int i = 0; i < (N); i++) {
				neighbors[i] = lattice[i];
			}

			for (int i = 0; i < (N); i++) {//adding neighbors 

				lattice[i] = (1 - df[d]) * neighbors[i] + 0.25*df[d] * (neighbors[nn[i][0]] + neighbors[nn[i][1]] + neighbors[nn[i][2]] + neighbors[nn[i][3]]);//no leaking of df (dispersal \\fraction). 1-.25*4=0
				if (lattice[i] <= 0)
				{
					extinct = 1;
					//add lower bound check
				}
				if (lattice[i] > 0)
				{
					extinct = 0;
					//add lower bound check
				}

				for (int i = 0; i < N; i++) {
					mjt[i] = .5 *  pow(-1, year) * (lattice[i] - last_lattice[i]);
				}

				for (int i = 0; i < N; i++) {
					mt = mt + mjt[i];

				}
				mt_store[d][year] = mt;

				mt = 0;
			}









		}
		extinct_table[d] = extinct;

		mL[d] = 0;
		for (long int year = burn; year < test; year++) {

			mL[d] = mL[d] + mt_store[d][year];
		}

	}







	ofstream myfile("notes.csv");

	myfile.is_open();//Recording testing conditions
	if (myfile.is_open())
	{

		myfile << L << " ";
		myfile << collect << " ";
		myfile << di << " ";
		myfile << r << " ";
		myfile << radio << " ";
		myfile << df[0] << " ";
		myfile << test << " ";
		myfile << skip << " ";
		myfile << ls << endl;
		myfile.close();
	}
	else cout << "Unable to open file";





	myfile.open("lattice_store.csv");//All the files are sent to matlab for visualization. In the matlab file I do variance = sus1 - sus2.^2 
	if (myfile.is_open())
	{
		for (int d = 0; d < di; d++) {

			for (int i = 0; i < N; i++) {
				myfile << lattice_store[d][i] << " ";
			}
		}
		myfile.close();
	}
	else cout << "Unable to open file";



	myfile.open("df_store.csv");
	if (myfile.is_open())
	{
		for (int d = 0; d < di; d++) {


			myfile << df[d] << " ";


		}
		myfile.close();
	}
	else cout << "Unable to open file";

	myfile.open("radio_store.csv");
	if (myfile.is_open())
	{



		myfile << radio << " ";



		myfile.close();
	}
	else cout << "Unable to open file";

	myfile.open("mL_store.csv");
	if (myfile.is_open())
	{

		for (int d = 0; d < di; d++) {


			myfile << mL[d] << " ";


		}


		myfile.close();
	}
	else cout << "Unable to open file";
	myfile.open("extinct_table.csv");
	if (myfile.is_open())
	{

		for (int d = 0; d < di; d++) {


			myfile << extinct_table[d] << " ";


		}


		myfile.close();
	}
	else cout << "Unable to open file";



	myfile.open("mt_store.csv");//All the files are sent to matlab for visualization. In the matlab file I do variance = sus1 - sus2.^2 
	if (myfile.is_open())
	{
		for (int d = 0; d < di; d++) {
			for (int year = 0; year < collect - 1; year++) {
				myfile << mt_store[d][year] << " ";
			}
		}
		myfile.close();
	}
	else cout << "Unable to open file";

	return 0;

}




