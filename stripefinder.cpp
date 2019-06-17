// All information and variables names are based on those in Doi 10.1038/ncomms7664

//#include "pch.h"
#include <iostream>
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <gsl/gsl_randist.h>

#include <time.h>
#include <fstream>
#include <sstream>

#include <algorithm>                    
#include <cmath> 
#include <math.h> 
#include <cstdlib> 
#include <stdio.h>
#include <fstream>
#include <time.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <deque>
#include <time.h>

#include <valarray>
using namespace std;

#define L 64  /* system size */
#define N L*L  /* number of spins */
#define D 4  /* number of neighbors */
int nearest_neighbors[N][D]; //neighbor table


const int length_of_side_of_lattice = L;
const int total_number_of_patches_on_lattice = N;
const int number_of_neighbors = D;

int const test = 20000;
long int const burn = 9990; 
long int const collect = test - burn;


double lattice[total_number_of_patches_on_lattice]; 

double last_lattice[total_number_of_patches_on_lattice];

double lattice_store[total_number_of_patches_on_lattice];

double neighbors[total_number_of_patches_on_lattice]; 
const  double reproductive_rate = 2.2;
//const  double reproductive_rate = 4;

double dispersal_fraction;
double radio = 0.5;

double high_value_in_ricker_two_cycle = 1.5029;
double low_value_in_ricker_two_cycle = 0.4971;

double tolerance_to_be_considered_within_cycle = 0.05;

double upper_limit_of_high_value ; //upper lim of x +

double lower_limit_of_high_value ; //lower lim of x +

double upper_limit_of_low_value ; //upper lim of x -

double lower_limit_of_low_value ; //lower lim of x -

void neighbor(void) {   //two dimensional neighbor ... interesting note to remember. In c++ division of ints result in ints. Not floats! That's why this modular arith works. We need this because we have periodic boundary conditions.
	int disp;
	for (int k = 0; k < total_number_of_patches_on_lattice; k++) {
		disp = length_of_side_of_lattice;
		nearest_neighbors[k][0] = (k / disp)*disp + (k + 1) % disp;  //+x
		nearest_neighbors[k][1] = (k / disp)*disp + (k - 1 + disp) % disp;  //-x
		disp = length_of_side_of_lattice * length_of_side_of_lattice;
		nearest_neighbors[k][2] = (k / disp)*disp + (k + length_of_side_of_lattice) % disp; //+y
		nearest_neighbors[k][3] = (k / disp)*disp + (k - length_of_side_of_lattice + disp) % disp; //-y
	}
}
double* update_map_then_nearest_neighbor_interaction(double lattice[], double last_lattice[], double dispersal_fraction)
{
	double neighbors[total_number_of_patches_on_lattice];

	for (int i = 0; i < total_number_of_patches_on_lattice; i++) {
		lattice[i] = last_lattice[i] * exp(reproductive_rate*(1 - last_lattice[i]));
	}

	for (int i = 0; i < total_number_of_patches_on_lattice; i++) {
		neighbors[i] = lattice[i];
	}

	for (int i = 0; i < (N); i++) {//adding neighbors 

		lattice[i] = (1 - dispersal_fraction) * neighbors[i] + 0.25*dispersal_fraction * (neighbors[nearest_neighbors[i][0]] + neighbors[nearest_neighbors[i][1]] + neighbors[nearest_neighbors[i][2]] + neighbors[nearest_neighbors[i][3]]);//no leaking of dispersal_fraction (dispersal \\fraction). 1-.25*4=0

	}

	return lattice;
}

double* remembering_the_last_lattice(double lattice[]) {
	for (int i = 0; i < total_number_of_patches_on_lattice; i++) {//remembering the last lattice
		last_lattice[i] = lattice[i];
	}

	return last_lattice;
}

double* storing_a_final_snap_shot(double lattice[], int year)
{
	if (year == test - 1)
	{
		for (int i = 0; i < total_number_of_patches_on_lattice; i++)
		{
			lattice_store[i] = lattice[i];

		}
	}

	return lattice_store;
}


double instantaneous_synchronization_across_time=0;

double unscaled_synchronization_parameter;

double scaled_syncronization_parameter;

double two_cycle_amplitude[length_of_side_of_lattice];
double instantaneous_synchronization;
double calculate_instantaneous_synchronization(double lattice[],double last_lattice[], const int length_of_side_of_lattice, int year)
{
	
	instantaneous_synchronization = 0;
	
	for (int i = 0; i < total_number_of_patches_on_lattice; i++) {
		two_cycle_amplitude[i] = .5 *  pow(-1, year) * (lattice[i] - last_lattice[i]);
	}

	for (int i = 0; i < total_number_of_patches_on_lattice; i++) {
		instantaneous_synchronization = instantaneous_synchronization + two_cycle_amplitude[i];
	}
	return instantaneous_synchronization;
}
double b_ricker = 2.028; //as defined in 10.1038/ncomms7664

double choice;

int converted_lattice_points_at_time_t[total_number_of_patches_on_lattice];
int converted_lattice_points_at_time_t_plus_1[total_number_of_patches_on_lattice];

int time_t_minus_two[total_number_of_patches_on_lattice];

int converted_lattice_points[total_number_of_patches_on_lattice];
int* convert_all_lattice_points(double test_points_to_be_converted[total_number_of_patches_on_lattice])
{
	for (int i = 0; i < total_number_of_patches_on_lattice; i++) {

		if (test_points_to_be_converted[i] < lower_limit_of_low_value || test_points_to_be_converted[i] > upper_limit_of_high_value)
		{
			converted_lattice_points[i] = 0;//Possible values on border of domains
		}
		if (test_points_to_be_converted[i] > upper_limit_of_low_value &&  test_points_to_be_converted[i] < lower_limit_of_high_value)
		{
			converted_lattice_points[i] = 0;//Possible values on border of domains
		}
		if (test_points_to_be_converted[i] > lower_limit_of_high_value)
		{
			converted_lattice_points[i] = 1; //Equivalent to spin up in Ising model
		}
		if (test_points_to_be_converted[i] > lower_limit_of_low_value && test_points_to_be_converted[i] < upper_limit_of_low_value)
		{
			converted_lattice_points[i] = -1;//Equivalent to spin down in Ising model
		}

	}
	return converted_lattice_points;
}


int counter = 0;
double percent_of_spins_that_flip[test];
double find_percent_of_spins_that_flip(int converted_lattice_points_at_time_t[total_number_of_patches_on_lattice], int converted_lattice_points_at_time_t_plus_1[total_number_of_patches_on_lattice],int counter) {
	
	double difference_between_time_t_plus_1_and_time_t=0;
	//cout << difference_between_time_t_plus_1_and_time_t << endl;

	for (int i = 0; i < total_number_of_patches_on_lattice; i++)
	{
		difference_between_time_t_plus_1_and_time_t = difference_between_time_t_plus_1_and_time_t + abs(converted_lattice_points_at_time_t_plus_1[i] - converted_lattice_points_at_time_t[i]);
	}
	//cout << difference_between_time_t_plus_1_and_time_t << endl;
	percent_of_spins_that_flip[counter] = difference_between_time_t_plus_1_and_time_t/ total_number_of_patches_on_lattice/2*100;
	
	return percent_of_spins_that_flip[counter];
}



double mean_energy_of_a_patch;
double find_mean_energy_of_each_patch(int converted_lattice_points[], int nearest_neighbors[total_number_of_patches_on_lattice][number_of_neighbors]) {
	//The motivation for this definition of energy can be found at : https://iopscience.iop.org/article/10.1088/0253-6102/51/4/18/pdf
	int energy_of_each_spin;
	double energy_across_lattice = 0;
	int energy_from_this_neighbor_interaction;
	int index_nearest_neighbor;

	for (int i = 0; i < total_number_of_patches_on_lattice; i++) {
		energy_of_each_spin = 0;
		for (int j = 0; j < number_of_neighbors; j++) {
			index_nearest_neighbor = nearest_neighbors[i][j];
			if (converted_lattice_points[i] == converted_lattice_points[index_nearest_neighbor]) {
				energy_from_this_neighbor_interaction = 0;
			}
			if (converted_lattice_points[i] != converted_lattice_points[index_nearest_neighbor]) {
				energy_from_this_neighbor_interaction = 1;
			}
			energy_of_each_spin = energy_of_each_spin + energy_from_this_neighbor_interaction;
		}
		energy_across_lattice = energy_across_lattice + energy_of_each_spin;
	}
	mean_energy_of_a_patch = energy_across_lattice / total_number_of_patches_on_lattice / 4;
	//The zero/boundary line is has two nearest neighbors that don't match, but really there should be only 1. 
	//The other way to do this would be to give interactions between 0 and 1 or 0 and -1 
	//energies of 1/2 instead of 1. Since there is only a 'half' difference in spins instead of a full one.

	return mean_energy_of_a_patch;
}
double number_of_stripes = 0;
double find_number_of_stripes(double mean_energy_of_a_patch) {

	number_of_stripes = length_of_side_of_lattice * mean_energy_of_a_patch;

	return number_of_stripes;
}


int indices_for_horizontal_stripe_test[length_of_side_of_lattice];

int* create_indices_for_horizontal_stripe_test_for_given_size_length() {
	for (int i = 0; i < length_of_side_of_lattice ; i++) {
		indices_for_horizontal_stripe_test[i] = i;
	}
	return indices_for_horizontal_stripe_test;
}


int indices_for_vertical_stripe_test[length_of_side_of_lattice];

int* create_indices_for_vertical_stripe_test_for_given_size_length() {
	for (int i = 0; i < length_of_side_of_lattice ; i++) {
		indices_for_vertical_stripe_test[i] = 1+i* length_of_side_of_lattice;
	}
	return indices_for_vertical_stripe_test;
}

int indices_for_descending_diagonal_stripe_test[length_of_side_of_lattice];

int* create_indices_for_descending_diagonal_stripe_test_for_given_size_length() {
	for (int i = 0; i < length_of_side_of_lattice ; i++) {
		indices_for_descending_diagonal_stripe_test[i] = length_of_side_of_lattice - 1 + i * (length_of_side_of_lattice -1);
	}
	return indices_for_descending_diagonal_stripe_test;
}

int indices_for_ascending_diagonal_stripe_test[length_of_side_of_lattice];

int* create_indices_for_ascending_diagonal_stripe_test_for_given_size_length() {
	for (int i = 0; i < length_of_side_of_lattice; i++) {
		indices_for_ascending_diagonal_stripe_test[i] = i * (length_of_side_of_lattice + 1);
	}
	return indices_for_ascending_diagonal_stripe_test;
}

int converted_test_points[length_of_side_of_lattice];

int* convert_test_points_of_Ricker_map_to_integers_for_any_stripe_test(double original_lattice[],int indices_for_given_stripe_test[]) {
	
	double test_points_to_be_converted[length_of_side_of_lattice];
	
	int index;

	for (int i = 0; i < length_of_side_of_lattice ; i++) {
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
int horizontal_stripe_test(int horizontal_stripe_with_converted_points[]) {
	is_this_a_horizontal_stripe = 0;
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
int vertical_stripe_test(int vertical_stripe_with_converted_points[]) {
	is_this_a_vertical_stripe = 0;
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
int descending_diagonal_stripe_test(int descending_diagonal_stripe_with_converted_points[]) {
	is_this_a_descending_diagonal_stripe = 0;
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
int ascending_diagonal_stripe_test(int ascending_diagonal_stripe_with_converted_points[]) {
	is_this_an_ascending_diagonal_stripe = 0;
	int sum = 0;
	for (int i = 0; i <= (length_of_side_of_lattice - 1); i++) {
		sum = sum + ascending_diagonal_stripe_with_converted_points[i];
	}

	if (abs(sum) == length_of_side_of_lattice)
	{
		is_this_an_ascending_diagonal_stripe = 1;
	}
	return is_this_an_ascending_diagonal_stripe;
}

int stripe_test_result ;
int what_kind_of_stripe_is_it(int is_this_an_ascending_diagonal_stripe, int is_this_a_descending_diagonal_stripe, int is_this_a_vertical_stripe, int is_this_a_horizontal_stripe) {
	
	if(is_this_an_ascending_diagonal_stripe==1)
	{
		stripe_test_result = 1;
	}

	if (is_this_a_horizontal_stripe == 1)
	{
		stripe_test_result = 2;
	}

	if (is_this_a_vertical_stripe == 1)
	{
		stripe_test_result = 3;
	}

	if (is_this_a_descending_diagonal_stripe == 1)
	{
		stripe_test_result = 4;
	}

	int sum_of_tests;
	sum_of_tests = is_this_an_ascending_diagonal_stripe + is_this_a_descending_diagonal_stripe + is_this_a_vertical_stripe + is_this_a_horizontal_stripe;
	if ((sum_of_tests > 1 && sum_of_tests!=4 )|| (sum_of_tests > 1 && number_of_stripes > 0) || (sum_of_tests == 0 && number_of_stripes > 0))
	{
		stripe_test_result = 5;
	}

	if (sum_of_tests == 4 )
	{
		stripe_test_result = 6;
	}
	return stripe_test_result;
}

int main(int argc, char** argv) {
	const gsl_rng_type * T;
	gsl_rng * r;


	int man_seed = atol(argv[1]);
	int i, n = 1;


	gsl_rng_env_setup();

	//void gsl_rng_set(const gsl_rng * r, unsigned long int '121');


	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, man_seed);


	dispersal_fraction = 0.5;

	double tolerance_to_be_considered_within_cycle = 0.05;

	upper_limit_of_high_value = high_value_in_ricker_two_cycle + tolerance_to_be_considered_within_cycle; //upper lim of x +

	lower_limit_of_high_value = high_value_in_ricker_two_cycle - tolerance_to_be_considered_within_cycle; //lower lim of x +

	upper_limit_of_low_value = low_value_in_ricker_two_cycle + tolerance_to_be_considered_within_cycle; //upper lim of x -

	lower_limit_of_low_value = low_value_in_ricker_two_cycle - tolerance_to_be_considered_within_cycle; //lower lim of x -

	neighbor(); //create neighbor table

	for (int i = 0; i < total_number_of_patches_on_lattice; i++) { //first lattice generation
		choice = gsl_ran_flat(r, 0, 1);
		if (choice >= radio) {
			lattice[i] = gsl_ran_flat(r, lower_limit_of_high_value, upper_limit_of_high_value);
		}
		if (choice < radio) {
			lattice[i] = gsl_ran_flat(r, lower_limit_of_low_value, upper_limit_of_low_value);;
		}

	}

	for (long int year = 0; year < test; year++) {

		remembering_the_last_lattice(lattice);

		storing_a_final_snap_shot(lattice, year);

		update_map_then_nearest_neighbor_interaction(lattice, last_lattice, dispersal_fraction);
		
		
		if (year > burn) {
			instantaneous_synchronization_across_time = instantaneous_synchronization_across_time + calculate_instantaneous_synchronization(lattice, last_lattice, length_of_side_of_lattice, year);
		}

		convert_all_lattice_points(lattice);

		if ((year % 2 == 0) && (year % 4 == 0))
		{
			for (int i = 0; i < total_number_of_patches_on_lattice; i++)
			{
				time_t_minus_two[i] = converted_lattice_points[i];				
			}
		}
		if (year % 2 == 0 && (year % 4 == 2))
		{
			counter++;
			find_percent_of_spins_that_flip(converted_lattice_points,time_t_minus_two,counter);
		}

	}
	
	unscaled_synchronization_parameter = instantaneous_synchronization_across_time;
	
	scaled_syncronization_parameter = b_ricker * unscaled_synchronization_parameter * pow(length_of_side_of_lattice,1/8) / total_number_of_patches_on_lattice / collect;

	convert_all_lattice_points(lattice_store);

	find_mean_energy_of_each_patch(converted_lattice_points, nearest_neighbors);

	find_number_of_stripes(mean_energy_of_a_patch);

	
		create_indices_for_horizontal_stripe_test_for_given_size_length();
		create_indices_for_vertical_stripe_test_for_given_size_length();
		create_indices_for_descending_diagonal_stripe_test_for_given_size_length();
		create_indices_for_ascending_diagonal_stripe_test_for_given_size_length();

		convert_test_points_of_Ricker_map_to_integers_for_any_stripe_test(lattice_store, indices_for_horizontal_stripe_test);

		horizontal_stripe_test(converted_test_points);

		convert_test_points_of_Ricker_map_to_integers_for_any_stripe_test(lattice_store, indices_for_vertical_stripe_test);

		vertical_stripe_test(converted_test_points);

		convert_test_points_of_Ricker_map_to_integers_for_any_stripe_test(lattice_store, indices_for_descending_diagonal_stripe_test);

		descending_diagonal_stripe_test(converted_test_points);

		convert_test_points_of_Ricker_map_to_integers_for_any_stripe_test(lattice_store, indices_for_ascending_diagonal_stripe_test);

		ascending_diagonal_stripe_test(converted_test_points);

		what_kind_of_stripe_is_it(is_this_an_ascending_diagonal_stripe, is_this_a_descending_diagonal_stripe, is_this_a_vertical_stripe, is_this_a_horizontal_stripe);

	ofstream myfile("stripe_conditions.txt");
	if (myfile.is_open())
	{
		myfile << "1 means this is an ascending diagonal stripe.\n";
		myfile << "2 means this is a horizontal stripe.\n";
		myfile << "3 means this is a vertical stripe.\n";
		myfile << "4 means this is an descending diagonal stripe.\n";
		myfile << "5 means this is neither a striped nor a coherent state.\n";
		myfile << "6 means this is a coherent state.\n";
		myfile.close();
	}
	else cout << "Unable to open file";

	myfile.open("notes.txt");//Recording testing conditions
	if (myfile.is_open())
	{
		myfile << length_of_side_of_lattice << " ";
		myfile << collect << " ";
		myfile << scaled_syncronization_parameter << " ";
		myfile << reproductive_rate << " ";
		myfile << radio << " ";
		myfile << dispersal_fraction << " ";
		myfile << test << " ";
		myfile << burn << " ";
		myfile << mean_energy_of_a_patch << " ";
		myfile << number_of_stripes << " ";
		myfile << percent_of_spins_that_flip << " ";
		myfile << man_seed << " ";
		myfile << stripe_test_result << endl;
		myfile.close();
	}
	else cout << "Unable to open file";


	myfile.open("notes_in_words.txt");//Recording testing conditions
	if (myfile.is_open())
	{
		myfile << "length_of_side_of_lattice.\n" ;
		myfile << "collect.\n ";
		myfile << "scaled_syncronization_parameter.\n  ";
		myfile << "r.\n  ";
		myfile << "radio.\n  ";
		myfile << "dispersal_fraction.\n ";
		myfile << "test.\n ";
		myfile << "burn.\n ";
		myfile << "mean_energy_of_a_patch.\n ";
		myfile << "number_of_stripes.\n ";
		myfile << "percent_of_spins_that_flip.\n ";
		myfile << "seed.\n";
		myfile << "stripe_test_result" << endl;
		myfile.close();
	}
	else cout << "Unable to open file";


	myfile.open("converted_lattice_points.csv");//All the files are sent to matlab for visualization. In the matlab file I do variance = sus1 - sus2.^2 
	if (myfile.is_open())
	{
		for (int i = 0; i < total_number_of_patches_on_lattice; i++) {
			myfile << converted_lattice_points[i] << " ";
		}

		myfile.close();
	}
	else cout << "Unable to open file";

	myfile.open("lattice_store.csv");//All the files are sent to matlab for visualization. In the matlab file I do variance = sus1 - sus2.^2 
	if (myfile.is_open())
	{
		for (int i = 0; i < total_number_of_patches_on_lattice; i++) {
			myfile << lattice_store[i] << " ";
		}

		myfile.close();
	}
	else cout << "Unable to open file";


	myfile.open("percent_of_spins_that_flip.csv");//All the files are sent to matlab for visualization. In the matlab file I do variance = sus1 - sus2.^2 
	if (myfile.is_open())
	{
		for (int i = 0; i < counter; i++) {
			myfile << percent_of_spins_that_flip[i] << " ";
		}

		myfile.close();
	}
	else cout << "Unable to open file";

	return 0;

}




