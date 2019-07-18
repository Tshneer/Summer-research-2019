// All information and variables names are based on those in Doi 10.1038/ncomms7664
//Version 1 6/18/19
//#include "pch.h"
#include <iostream>
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <fstream>
#include <sstream>

#include <algorithm>                    
//#include <cmath> 
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


long int end_year = 0;

long int const test = 100000;
long int const burn = 0;
long int const collect = test - burn;

const  double reproductive_rate = 2.2;

double choice;

double upper_limit_of_high_value; //upper lim of x +

double lower_limit_of_high_value; //lower lim of x +

double upper_limit_of_low_value; //upper lim of x -

double lower_limit_of_low_value; //lower lim of x -

void neighbor(const int length_of_side_of_lattice, const int total_number_of_patches_on_lattice, const int number_of_neighbors, int nearest_neighbors[][4]) {   //two dimensional neighbor ... interesting note to remember. In c++ division of ints result in ints. Not floats! That's why this modular arith works. We need this because we have periodic boundary conditions.
	int disp;
	//int nearest_neighbors[total_number_of_patches_on_lattice][number_of_neighbors];
	for (int k = 0; k < total_number_of_patches_on_lattice; k++) {
		disp = length_of_side_of_lattice;
		nearest_neighbors[k][0] = (k / disp)*disp + (k + 1) % disp;  //+x
		nearest_neighbors[k][1] = (k / disp)*disp + (k - 1 + disp) % disp;  //-x
		disp = length_of_side_of_lattice * length_of_side_of_lattice;
		nearest_neighbors[k][2] = (k / disp)*disp + (k + length_of_side_of_lattice) % disp; //+y
		nearest_neighbors[k][3] = (k / disp)*disp + (k - length_of_side_of_lattice + disp) % disp; //-y
	}
}

double* update_map_then_nearest_neighbor_interaction(double lattice[], double last_lattice[], double dispersal_fraction, const int total_number_of_patches_on_lattice, int nearest_neighbors[][4], double neighbors[])
{


	for (int i = 0; i < total_number_of_patches_on_lattice; i++) {
		lattice[i] = last_lattice[i] * exp(reproductive_rate*(1 - last_lattice[i]));

	}

	for (int i = 0; i < total_number_of_patches_on_lattice; i++) {
		neighbors[i] = lattice[i];
	}

	for (int i = 0; i < total_number_of_patches_on_lattice; i++) {//adding neighbors 

		lattice[i] = (1 - dispersal_fraction) * neighbors[i] + 0.25*dispersal_fraction * (neighbors[nearest_neighbors[i][0]] + neighbors[nearest_neighbors[i][1]] + neighbors[nearest_neighbors[i][2]] + neighbors[nearest_neighbors[i][3]]);//no leaking of dispersal_fraction (dispersal \\fraction). 1-.25*4=0

	}

	return lattice;
}

double* remembering_the_last_lattice(double lattice[], int total_number_of_patches_on_lattice, double last_lattice[]) {
	for (int i = 0; i < total_number_of_patches_on_lattice; i++) {//remembering the last lattice
		last_lattice[i] = lattice[i];
	}

	return last_lattice;
}
int* lower_diagonal_neighbor_table(const int total_number_of_patches_on_lattice, int nearest_neighbors[][4], int descending_diagonal_neighbors[])
{
	int lower_neighbor;
	int descending_diagonal;
	for (int k = 0; k < total_number_of_patches_on_lattice; k++) {
		lower_neighbor = nearest_neighbors[k][3];
		descending_diagonal = nearest_neighbors[lower_neighbor][0];

		descending_diagonal_neighbors[k] = descending_diagonal;
	}

	return descending_diagonal_neighbors;
}

int* upper_diagonal_neighbor_table(const int total_number_of_patches_on_lattice, int nearest_neighbors[][4], int ascending_diagonal_neighbors[])
{
	int lower_neighbor;
	int ascending_diagonal;

	for (int k = 0; k < total_number_of_patches_on_lattice; k++) {
		lower_neighbor = nearest_neighbors[k][2];
		ascending_diagonal = nearest_neighbors[lower_neighbor][0];

		ascending_diagonal_neighbors[k] = ascending_diagonal;
	}

	return ascending_diagonal_neighbors;
}


int* create_row_indicies_for_horizontal_stripe_test_for_given_size_length(int j, const int length_of_side_of_lattice, int row_indicies_for_horizontal_stripe_test[], int indicies_for_start_of_horizontal_stripe[]) {
	


	int k;

		for (int i = 0; i < length_of_side_of_lattice; i++)
		{
			k = indicies_for_start_of_horizontal_stripe[j] + i ;
			row_indicies_for_horizontal_stripe_test[i] = k;
		}
	return row_indicies_for_horizontal_stripe_test;
}



int* create_row_indicies_for_vertical_stripe_test_for_given_size_length(int j, const int length_of_side_of_lattice, int row_indicies_for_vertical_stripe_test[], int indicies_for_start_of_vertical_stripe[]) {

	int k;


		for (int i = 0; i < length_of_side_of_lattice; i++) {
			k = indicies_for_start_of_vertical_stripe[j] + i * length_of_side_of_lattice;

			row_indicies_for_vertical_stripe_test[i] = k;
		}
	
	return row_indicies_for_vertical_stripe_test;
}


int* create_row_indicies_for_descending_diagonal_stripe_test_for_given_size_length(int j, const int length_of_side_of_lattice, int row_indicies_for_descending_diagonal_stripe_test[], int descending_diagonal_neighbors[], int indicies_for_start_of_descending_diagonal_stripe[]) {


	
	int start;
		start = indicies_for_start_of_descending_diagonal_stripe[j];


		for (int i = 0; i < length_of_side_of_lattice; i++) {
			row_indicies_for_descending_diagonal_stripe_test[i] = descending_diagonal_neighbors[start];
			start = descending_diagonal_neighbors[start];

		}
	
	return row_indicies_for_descending_diagonal_stripe_test;
}

int* create_row_indicies_for_ascending_diagonal_stripe_test_for_given_size_length(int j, const int length_of_side_of_lattice, int row_indicies_for_ascending_diagonal_stripe_test[], int ascending_diagonal_neighbors[],int indicies_for_start_of_ascending_diagonal_stripe[]) {

	

	int start;
		start = indicies_for_start_of_ascending_diagonal_stripe[j];


		for (int i = 0; i < length_of_side_of_lattice; i++) {
			row_indicies_for_ascending_diagonal_stripe_test[i] = ascending_diagonal_neighbors[start];
			start = ascending_diagonal_neighbors[start];

		}
	
	return row_indicies_for_ascending_diagonal_stripe_test;
}



int is_this_a_vertical_stripe = 0;
int vertical_stripe_test(int vertical_stripe_with_converted_points[], const int length_of_side_of_lattice) {
	is_this_a_vertical_stripe = 0;
	int sum = 0;
	int sum_is_not_zero = 0;
	for (int i = 0; i < length_of_side_of_lattice; i++) {
		sum = sum + vertical_stripe_with_converted_points[i];
		if (sum != 0)
		{
			sum_is_not_zero = 1;
		}
	}
	if (abs(sum) == length_of_side_of_lattice || (sum == 0 && sum_is_not_zero == 0))
	{
		is_this_a_vertical_stripe = 1;
	}
	return is_this_a_vertical_stripe;
}

int is_this_a_horizontal_stripe = 0;
int horizontal_stripe_test(int horizontal_stripe_with_converted_points[], const int length_of_side_of_lattice) {
	is_this_a_horizontal_stripe = 0;
	int sum_is_not_zero = 0;
	int sum = 0;
	for (int i = 0; i < length_of_side_of_lattice; i++) {
		sum = sum + horizontal_stripe_with_converted_points[i];
		if (sum != 0)
		{
			sum_is_not_zero = 1;
		}
	}

	if (abs(sum) == length_of_side_of_lattice || (sum == 0 && sum_is_not_zero == 0))
	{
		is_this_a_horizontal_stripe = 1;
	}
	return is_this_a_horizontal_stripe;
}

int is_this_a_descending_diagonal_stripe = 0;
int descending_diagonal_stripe_test(int descending_diagonal_stripe_with_converted_points[], const int length_of_side_of_lattice) {
	is_this_a_descending_diagonal_stripe = 0;
	int sum = 0;
	int sum_is_not_zero = 0;
	int sum_is_32_and_is_a_diagonal = 0;
	for (int i = 0; i < length_of_side_of_lattice; i++) {
		sum = sum + descending_diagonal_stripe_with_converted_points[i];
		if (sum != 0)
		{
			sum_is_not_zero = 1;
		}
	}

	if (abs(sum) == length_of_side_of_lattice / 2)
	{
		double mean_of_even_test_indicies;
		double sum_for_even_test_indicies = 0;
		double variance_for_even_test_indicies = 0;

		for (int i = 0; i < length_of_side_of_lattice; i = i + 2)
		{
			sum_for_even_test_indicies = sum_for_even_test_indicies + descending_diagonal_stripe_with_converted_points[i];
		}

		mean_of_even_test_indicies = sum_for_even_test_indicies / (2 * length_of_side_of_lattice);

		for (int i = 0; i < length_of_side_of_lattice; i = i + 2)
		{
			variance_for_even_test_indicies = variance_for_even_test_indicies + pow(descending_diagonal_stripe_with_converted_points[i] - mean_of_even_test_indicies, 2);
		}
		double mean_of_odd_test_indicies;
		double sum_for_odd_test_indicies = 0;
		double variance_for_odd_test_indicies = 0;

		for (int i = 0; i < length_of_side_of_lattice; i = i + 2)
		{
			sum_for_odd_test_indicies = sum_for_even_test_indicies + descending_diagonal_stripe_with_converted_points[i];
		}

		mean_of_odd_test_indicies = sum_for_odd_test_indicies / (2 * length_of_side_of_lattice);

		for (int i = 0; i < length_of_side_of_lattice; i = i + 2)
		{
			variance_for_odd_test_indicies = variance_for_odd_test_indicies + pow(descending_diagonal_stripe_with_converted_points[i] - mean_of_odd_test_indicies, 2);
		}

		if (variance_for_odd_test_indicies + variance_for_even_test_indicies == 0)
		{
			sum_is_32_and_is_a_diagonal = 1;
		}
	}


	if (abs(sum) == length_of_side_of_lattice || (sum == 0 && sum_is_not_zero == 0) || sum_is_32_and_is_a_diagonal == 1)
	{
		is_this_a_descending_diagonal_stripe = 1;

	}


	return is_this_a_descending_diagonal_stripe;
}

int is_this_an_ascending_diagonal_stripe = 0;
int ascending_diagonal_stripe_test(int ascending_diagonal_stripe_with_converted_points[], const int length_of_side_of_lattice) {
	is_this_an_ascending_diagonal_stripe = 0;
	int sum = 0;
	int sum_is_not_zero = 0;
	int sum_is_32_and_is_a_diagonal = 0;
	for (int i = 0; i < length_of_side_of_lattice; i++) {
		sum = sum + ascending_diagonal_stripe_with_converted_points[i];
		if (sum != 0)
		{
			sum_is_not_zero = 1;
		}
	}

	if (abs(sum) == length_of_side_of_lattice / 2)
	{
		double mean_of_even_test_indicies;
		double sum_for_even_test_indicies = 0;
		double variance_for_even_test_indicies = 0;
		for (int i = 0; i < length_of_side_of_lattice; i = i + 2)
		{
			sum_for_even_test_indicies = sum_for_even_test_indicies + ascending_diagonal_stripe_with_converted_points[i];
		}

		mean_of_even_test_indicies = sum_for_even_test_indicies / (2 * length_of_side_of_lattice);

		for (int i = 0; i < length_of_side_of_lattice; i = i + 2)
		{
			variance_for_even_test_indicies = variance_for_even_test_indicies + pow(ascending_diagonal_stripe_with_converted_points[i] - mean_of_even_test_indicies, 2);
		}
		double mean_of_odd_test_indicies;
		double sum_for_odd_test_indicies = 0;
		double variance_for_odd_test_indicies = 0;

		for (int i = 0; i < length_of_side_of_lattice; i = i + 2)
		{
			sum_for_odd_test_indicies = sum_for_even_test_indicies + ascending_diagonal_stripe_with_converted_points[i];
		}

		mean_of_odd_test_indicies = sum_for_odd_test_indicies / (2 * length_of_side_of_lattice);

		for (int i = 0; i < length_of_side_of_lattice; i = i + 2)
		{
			variance_for_odd_test_indicies = variance_for_odd_test_indicies + pow(ascending_diagonal_stripe_with_converted_points[i] - mean_of_odd_test_indicies, 2);
		}

		if (variance_for_odd_test_indicies + variance_for_even_test_indicies == 0)
		{
			sum_is_32_and_is_a_diagonal = 1;
		}
	}

	if (abs(sum) == length_of_side_of_lattice || (sum == 0 && sum_is_not_zero == 0) || sum_is_32_and_is_a_diagonal == 1)
	{
		is_this_an_ascending_diagonal_stripe = 1;
	}
	return is_this_an_ascending_diagonal_stripe;
}

int* convert_all_lattice_points(double test_points_to_be_converted_for_lattice[], int total_number_of_patches_on_lattice, int converted_lattice_points_for_lattice[])
{
	for (int i = 0; i < total_number_of_patches_on_lattice; i++) {

		if (test_points_to_be_converted_for_lattice[i] < lower_limit_of_low_value || test_points_to_be_converted_for_lattice[i] > upper_limit_of_high_value)
		{
			converted_lattice_points_for_lattice[i] = 0;//Possible values on border of domains
		}
		if (test_points_to_be_converted_for_lattice[i] > upper_limit_of_low_value &&  test_points_to_be_converted_for_lattice[i] < lower_limit_of_high_value)
		{
			converted_lattice_points_for_lattice[i] = 0;//Possible values on border of domains
		}
		if (test_points_to_be_converted_for_lattice[i] > lower_limit_of_high_value)
		{
			converted_lattice_points_for_lattice[i] = 1; //Equivalent to spin up in Ising model
		}
		if (test_points_to_be_converted_for_lattice[i] > lower_limit_of_low_value && test_points_to_be_converted_for_lattice[i] < upper_limit_of_low_value)
		{
			converted_lattice_points_for_lattice[i] = -1;//Equivalent to spin down in Ising model
		}
		//cout<< "converted_lattice_points " <<converted_lattice_points_for_lattice[i]<<endl;

	}
	return converted_lattice_points_for_lattice;
}

int* convert_test_points_of_Ricker_map_to_integers_for_any_stripe_test(double original_lattice[], int indicies_for_given_stripe_test[], const int length_of_side_of_lattice, int converted_test_points_for_length[], double test_points_to_be_converted_for_length[]) {



	int index;

	for (int i = 0; i < length_of_side_of_lattice; i++) {
		index = indicies_for_given_stripe_test[i];

		test_points_to_be_converted_for_length[i] = original_lattice[index];

	}

	for (int i = 0; i < length_of_side_of_lattice; i++) {

		if (test_points_to_be_converted_for_length[i] < lower_limit_of_low_value || test_points_to_be_converted_for_length[i] > upper_limit_of_high_value)
		{
			converted_test_points_for_length[i] = 0;//Possible values on border of domains
		}
		if (test_points_to_be_converted_for_length[i] > upper_limit_of_low_value &&  test_points_to_be_converted_for_length[i] < lower_limit_of_high_value)
		{
			converted_test_points_for_length[i] = 0;//Possible values on border of domains
		}
		if (test_points_to_be_converted_for_length[i] > lower_limit_of_high_value)
		{
			converted_test_points_for_length[i] = 1; //Equivalent to spin up in Ising model
		}
		if (test_points_to_be_converted_for_length[i] > lower_limit_of_low_value && test_points_to_be_converted_for_length[i] < upper_limit_of_low_value)
		{
			converted_test_points_for_length[i] = -1;//Equivalent to spin down in Ising model
		}

	}
	return converted_test_points_for_length;
}


int stripe_test_result = 0;
int what_kind_of_stripe_is_it(int is_this_an_ascending_diagonal_stripe, int is_this_a_descending_diagonal_stripe, int is_this_a_vertical_stripe, int is_this_a_horizontal_stripe) {

	if (is_this_an_ascending_diagonal_stripe == 1)
	{
		stripe_test_result = 1;
	}


	if (is_this_a_vertical_stripe == 1)
	{
		stripe_test_result = 3;
	}

	if (is_this_a_horizontal_stripe == 1)
	{
		stripe_test_result = 2;
	}

	if (is_this_a_descending_diagonal_stripe == 1)
	{
		stripe_test_result = 4;
	}

	int sum_of_tests;
	sum_of_tests = is_this_an_ascending_diagonal_stripe + is_this_a_descending_diagonal_stripe + is_this_a_vertical_stripe + is_this_a_horizontal_stripe;
	if ((sum_of_tests > 1 && sum_of_tests != 4) || (sum_of_tests == 0 ))
	{
		stripe_test_result = 5;
	}

	if (sum_of_tests == 4)
	{
		stripe_test_result = 6;
	}
	
	return stripe_test_result;
}


int wavy_stripe_test_result = 0;
int what_kind_of_wavy_stripe_is_it(int is_this_an_ascending_diagonal_stripe, int is_this_a_descending_diagonal_stripe, int is_this_a_vertical_stripe, int is_this_a_horizontal_stripe) {

	if (is_this_an_ascending_diagonal_stripe == 1)
	{
		wavy_stripe_test_result = 1;
	}


	if (is_this_a_vertical_stripe == 1)
	{
		wavy_stripe_test_result = 3;
	}

	if (is_this_a_horizontal_stripe == 1)
	{
		wavy_stripe_test_result = 2;
	}

	if (is_this_a_descending_diagonal_stripe == 1)
	{
		wavy_stripe_test_result = 4;
	}

	int sum_of_tests;
	sum_of_tests = is_this_an_ascending_diagonal_stripe + is_this_a_descending_diagonal_stripe + is_this_a_vertical_stripe + is_this_a_horizontal_stripe;
	if ((sum_of_tests > 1 && sum_of_tests != 4) || (sum_of_tests == 0 ))
	{
		wavy_stripe_test_result = 5;
	}

	if (sum_of_tests == 4)
	{
		wavy_stripe_test_result = 6;
	}
	
	return wavy_stripe_test_result;
}

int main(int argc, char** argv) {	

	double tolerance_to_be_considered_within_cycle = 0.01;

	double high_value_in_ricker_two_cycle = 1.5029;
	double low_value_in_ricker_two_cycle = 0.4971;
	double bb = atol(argv[4]);

	double dispersal_fraction=bb/10;
	double radio = 0.5;

	upper_limit_of_high_value = high_value_in_ricker_two_cycle + tolerance_to_be_considered_within_cycle; //upper lim of x +

	lower_limit_of_high_value = high_value_in_ricker_two_cycle - tolerance_to_be_considered_within_cycle; //lower lim of x +

	upper_limit_of_low_value = low_value_in_ricker_two_cycle + tolerance_to_be_considered_within_cycle; //upper lim of x -

	lower_limit_of_low_value = low_value_in_ricker_two_cycle - tolerance_to_be_considered_within_cycle; //lower lim of x -

	const int length_of_side_of_lattice = atol(argv[2]);
	const int total_number_of_patches_on_lattice = length_of_side_of_lattice * length_of_side_of_lattice;
	const int number_of_neighbors = 4;

	double number = total_number_of_patches_on_lattice;

	double tolerance_for_difference=0.1;


	int nearest_neighbors[total_number_of_patches_on_lattice][number_of_neighbors]; //neighbor table
	int descending_diagonal_neighbors[total_number_of_patches_on_lattice];
	int ascending_diagonal_neighbors[total_number_of_patches_on_lattice];

	double lattice[total_number_of_patches_on_lattice];
	double neighbors[total_number_of_patches_on_lattice];

	double last_lattice[total_number_of_patches_on_lattice];
	double last_last_lattice[total_number_of_patches_on_lattice];

	double lattice_time_t_minus_one[total_number_of_patches_on_lattice];

	int converted_lattice_points_for_lattice[total_number_of_patches_on_lattice];

	int final_lattice[total_number_of_patches_on_lattice];

	double test_points_to_be_converted_for_length[length_of_side_of_lattice];

	double test_points_to_be_converted_for_lattice[total_number_of_patches_on_lattice];

	int converted_test_points_for_length[length_of_side_of_lattice];

	double tolerance_for_kink = 0.5 / length_of_side_of_lattice;

	int indicies_for_horizontal_stripe_test[length_of_side_of_lattice][length_of_side_of_lattice];
	int indicies_for_vertical_stripe_test[length_of_side_of_lattice][length_of_side_of_lattice];
	int indicies_for_ascending_diagonal_stripe_test[length_of_side_of_lattice][length_of_side_of_lattice];
	int indicies_for_descending_diagonal_stripe_test[length_of_side_of_lattice][length_of_side_of_lattice];

	int indicies_for_start_of_horizontal_stripe[length_of_side_of_lattice];
	int indicies_for_start_of_vertical_stripe[length_of_side_of_lattice];
	int indicies_for_start_of_descending_diagonal_stripe[length_of_side_of_lattice];
	int indicies_for_start_of_ascending_diagonal_stripe[length_of_side_of_lattice];

	int row_indicies_for_horizontal_stripe_test[length_of_side_of_lattice];
	int row_indicies_for_vertical_stripe_test[length_of_side_of_lattice];
	int row_indicies_for_descending_diagonal_stripe_test[length_of_side_of_lattice];
	int row_indicies_for_ascending_diagonal_stripe_test[length_of_side_of_lattice];

	int single_vertical_fail;
	int single_horizontal_fail;
	int single_descending_diagonal_fail;
	int single_ascending_diagonal_fail;

	int vertical_stripe_test_failed;
	int horizontal_stripe_test_failed;
	int descending_diagonal_stripe_test_failed;
	int ascending_diagonal_stripe_test_failed;

	int vertical_stripe_test_counter;
	int horizontal_stripe_test_counter;
	int ascending_diagonal_stripe_test_counter;
	int descending_diagonal_stripe_test_counter;

	double sum_of_difference_between_time_t_plus_1_and_time_t;


	neighbor(length_of_side_of_lattice, total_number_of_patches_on_lattice, number_of_neighbors, nearest_neighbors);
	lower_diagonal_neighbor_table(total_number_of_patches_on_lattice, nearest_neighbors, descending_diagonal_neighbors);
	upper_diagonal_neighbor_table(total_number_of_patches_on_lattice, nearest_neighbors, ascending_diagonal_neighbors);

	int skip = 100;

	const gsl_rng_type * T;
	gsl_rng * r;


	int man_seed = atol(argv[1]);
	int i, n = 1;

	gsl_rng_env_setup();


	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, man_seed);

	for (int i = 0; i < total_number_of_patches_on_lattice; i++) { //first lattice generation
		choice = gsl_ran_flat(r, 0, 1);
		if (choice >= radio) {
			lattice[i] = gsl_ran_flat(r, lower_limit_of_high_value, upper_limit_of_high_value);
		}
		if (choice < radio) {
			lattice[i] = gsl_ran_flat(r, lower_limit_of_low_value, upper_limit_of_low_value);;
		}

	}

	convert_all_lattice_points(lattice, total_number_of_patches_on_lattice, converted_lattice_points_for_lattice);

	for (int i = 0; i < length_of_side_of_lattice; i++) {
		indicies_for_start_of_horizontal_stripe[i] = i * length_of_side_of_lattice;
	}

	for (int i = 0; i < length_of_side_of_lattice; i++) {
		indicies_for_start_of_vertical_stripe[i] = i;
	}

	for (int i = 0; i < length_of_side_of_lattice; i++) {
		indicies_for_start_of_descending_diagonal_stripe[i] = 1 + i * length_of_side_of_lattice;
	}

	for (int i = 0; i < length_of_side_of_lattice; i++) {
		indicies_for_start_of_ascending_diagonal_stripe[i] = 1 + i * length_of_side_of_lattice;
	}


	for (int j = 0; j < length_of_side_of_lattice; j++) {
		create_row_indicies_for_horizontal_stripe_test_for_given_size_length(j, length_of_side_of_lattice, row_indicies_for_horizontal_stripe_test, indicies_for_start_of_horizontal_stripe);
		create_row_indicies_for_vertical_stripe_test_for_given_size_length(j, length_of_side_of_lattice, row_indicies_for_vertical_stripe_test, indicies_for_start_of_vertical_stripe);
		create_row_indicies_for_descending_diagonal_stripe_test_for_given_size_length(j, length_of_side_of_lattice, row_indicies_for_descending_diagonal_stripe_test, descending_diagonal_neighbors, indicies_for_start_of_descending_diagonal_stripe);			
		create_row_indicies_for_ascending_diagonal_stripe_test_for_given_size_length(j ,length_of_side_of_lattice, row_indicies_for_ascending_diagonal_stripe_test, ascending_diagonal_neighbors, indicies_for_start_of_ascending_diagonal_stripe);

		for (int b = 0; b < length_of_side_of_lattice; b++) {
			indicies_for_horizontal_stripe_test[j][b] = row_indicies_for_horizontal_stripe_test[b];
			indicies_for_vertical_stripe_test[j][b] = row_indicies_for_vertical_stripe_test[b];
			indicies_for_ascending_diagonal_stripe_test[j][b] = row_indicies_for_ascending_diagonal_stripe_test[b];
			indicies_for_descending_diagonal_stripe_test[j][b] = row_indicies_for_descending_diagonal_stripe_test[b];

		}
	}


	sum_of_difference_between_time_t_plus_1_and_time_t=100000;

	for (int year = 0; year < test; year++)
	{
		vertical_stripe_test_failed=0;
		horizontal_stripe_test_failed=0;
		ascending_diagonal_stripe_test_failed=0;
		descending_diagonal_stripe_test_failed=0;

		single_descending_diagonal_fail = 0;
		single_ascending_diagonal_fail = 0;
		single_vertical_fail = 0;
		single_horizontal_fail = 0;

		vertical_stripe_test_counter  = 0;
		horizontal_stripe_test_counter = 0;
		ascending_diagonal_stripe_test_counter = 0;
		descending_diagonal_stripe_test_counter = 0;

		if (year % skip == 0)
		{
			single_horizontal_fail = 0;
			single_vertical_fail = 0;
			single_ascending_diagonal_fail = 0;
			single_descending_diagonal_fail = 0;

			for (int j = 0; j < length_of_side_of_lattice; j++)
			{
				for (int b = 0; b < length_of_side_of_lattice; b++)
				{
					if (horizontal_stripe_test_failed == 0)
					{
						row_indicies_for_horizontal_stripe_test[b] = indicies_for_horizontal_stripe_test[j][b];
					}
					if (vertical_stripe_test_failed == 0)
					{
						row_indicies_for_vertical_stripe_test[b] = indicies_for_vertical_stripe_test[j][b];
					}
					if (ascending_diagonal_stripe_test_failed == 0)
					{
						row_indicies_for_ascending_diagonal_stripe_test[b] = indicies_for_ascending_diagonal_stripe_test[j][b];
					}
					if (descending_diagonal_stripe_test_failed == 0)
					{
						row_indicies_for_descending_diagonal_stripe_test[b] = indicies_for_descending_diagonal_stripe_test[j][b];
					}
				}

				if (vertical_stripe_test_failed == 0)
				{
					
					convert_test_points_of_Ricker_map_to_integers_for_any_stripe_test(lattice, row_indicies_for_vertical_stripe_test, length_of_side_of_lattice, converted_test_points_for_length, test_points_to_be_converted_for_length);
					vertical_stripe_test(converted_test_points_for_length, length_of_side_of_lattice);
				if (is_this_a_vertical_stripe == 0 && single_vertical_fail == 1)
				{
					vertical_stripe_test_failed = 1;
				}
				if (is_this_a_vertical_stripe == 0)
				{
					single_vertical_fail = 1;
				}
				if (is_this_a_vertical_stripe == 1)
				{
					single_vertical_fail = 0;
				}
				}

				if (horizontal_stripe_test_failed == 0)
				{
					
					convert_test_points_of_Ricker_map_to_integers_for_any_stripe_test(lattice, row_indicies_for_horizontal_stripe_test, length_of_side_of_lattice, converted_test_points_for_length, test_points_to_be_converted_for_length);
					horizontal_stripe_test(converted_test_points_for_length, length_of_side_of_lattice);
					if (is_this_a_horizontal_stripe == 0 && single_horizontal_fail == 1)
					{
						horizontal_stripe_test_failed = 1;
					}
					if (is_this_a_horizontal_stripe == 0)
					{
						single_horizontal_fail = 1;
					}
					if (is_this_a_horizontal_stripe == 1)
					{
						single_horizontal_fail = 0;
					}
				}

				if (ascending_diagonal_stripe_test_failed == 0)
				{
					
					convert_test_points_of_Ricker_map_to_integers_for_any_stripe_test(lattice, row_indicies_for_ascending_diagonal_stripe_test, length_of_side_of_lattice, converted_test_points_for_length, test_points_to_be_converted_for_length);
					ascending_diagonal_stripe_test(converted_test_points_for_length, length_of_side_of_lattice);
					if (is_this_an_ascending_diagonal_stripe == 0 && single_ascending_diagonal_fail == 1)
					{
						ascending_diagonal_stripe_test_failed = 1;
					}
					if (is_this_an_ascending_diagonal_stripe == 0)
					{
						single_ascending_diagonal_fail = 1;
					}
					if (is_this_an_ascending_diagonal_stripe == 1)
					{
						single_ascending_diagonal_fail = 0;
					}
				}
				if (descending_diagonal_stripe_test_failed == 0)
				{
					
					convert_test_points_of_Ricker_map_to_integers_for_any_stripe_test(lattice, row_indicies_for_descending_diagonal_stripe_test, length_of_side_of_lattice, converted_test_points_for_length, test_points_to_be_converted_for_length);
					descending_diagonal_stripe_test(converted_test_points_for_length, length_of_side_of_lattice);
					if (is_this_a_descending_diagonal_stripe == 0 && single_descending_diagonal_fail == 1)
					{
						descending_diagonal_stripe_test_failed = 1;
					}
					if (is_this_a_descending_diagonal_stripe == 0)
					{
						single_descending_diagonal_fail = 1;
					}
					if (is_this_a_descending_diagonal_stripe == 1)
					{
						single_descending_diagonal_fail = 0;
					}
				}
			}


			what_kind_of_stripe_is_it(is_this_an_ascending_diagonal_stripe, is_this_a_descending_diagonal_stripe, is_this_a_vertical_stripe, is_this_a_horizontal_stripe);
		}

		if (stripe_test_result != 5)
			{
				end_year = year;

				break;
			}

		if (sum_of_difference_between_time_t_plus_1_and_time_t < tolerance_for_difference  && (year % skip == 0))
		{
			for (int j = 0; j < length_of_side_of_lattice; j++)
			{
				for (int b = 0; b < length_of_side_of_lattice; b++)
				{
					
						row_indicies_for_horizontal_stripe_test[b] = indicies_for_horizontal_stripe_test[j][b];
					
					
						row_indicies_for_vertical_stripe_test[b] = indicies_for_vertical_stripe_test[j][b];
					
					
						row_indicies_for_ascending_diagonal_stripe_test[b] = indicies_for_ascending_diagonal_stripe_test[j][b];
					
					
						row_indicies_for_descending_diagonal_stripe_test[b] = indicies_for_descending_diagonal_stripe_test[j][b];
					
				}

				
					
				convert_test_points_of_Ricker_map_to_integers_for_any_stripe_test(lattice, row_indicies_for_vertical_stripe_test, length_of_side_of_lattice, converted_test_points_for_length, test_points_to_be_converted_for_length);
				vertical_stripe_test(converted_test_points_for_length, length_of_side_of_lattice);
				if (is_this_a_vertical_stripe == 0 && single_vertical_fail == 1)
				{
					vertical_stripe_test_failed = 1;
				}
				if (is_this_a_vertical_stripe == 0)
				{
					single_vertical_fail = 1;
				}
				if (is_this_a_vertical_stripe == 1)
				{
					single_vertical_fail = 0;
				}
				
				if (single_vertical_fail == 0)
				{
				vertical_stripe_test_counter = vertical_stripe_test_counter + 1;
				}

				if (single_vertical_fail == 1)
				{
				vertical_stripe_test_counter = vertical_stripe_test_counter - 1;
				}

				convert_test_points_of_Ricker_map_to_integers_for_any_stripe_test(lattice, row_indicies_for_horizontal_stripe_test, length_of_side_of_lattice, converted_test_points_for_length, test_points_to_be_converted_for_length);
				horizontal_stripe_test(converted_test_points_for_length, length_of_side_of_lattice);
				if (is_this_a_horizontal_stripe == 0 && single_horizontal_fail == 1)
				{
					horizontal_stripe_test_failed = 1;
				}
				if (is_this_a_horizontal_stripe == 0)
				{
					single_horizontal_fail = 1;
				}
				if (is_this_a_horizontal_stripe == 1)
				{
					single_horizontal_fail = 0;
				}								

				if (single_horizontal_fail == 0)
				{
				horizontal_stripe_test_counter = horizontal_stripe_test_counter + 1;
				}

				if (single_horizontal_fail == 1)
				{
				horizontal_stripe_test_counter = horizontal_stripe_test_counter - 1;
				}					

				convert_test_points_of_Ricker_map_to_integers_for_any_stripe_test(lattice, row_indicies_for_ascending_diagonal_stripe_test, length_of_side_of_lattice, converted_test_points_for_length, test_points_to_be_converted_for_length);
				ascending_diagonal_stripe_test(converted_test_points_for_length, length_of_side_of_lattice);
				if (is_this_an_ascending_diagonal_stripe == 0 && single_ascending_diagonal_fail == 1)
				{
					ascending_diagonal_stripe_test_failed = 1;
				}
				if (is_this_an_ascending_diagonal_stripe == 0)
				{
					single_ascending_diagonal_fail = 1;
				}
				if (is_this_an_ascending_diagonal_stripe == 1)
				{
					single_ascending_diagonal_fail = 0;
				}

				if (single_ascending_diagonal_fail == 0)
				{
				ascending_diagonal_stripe_test_counter = ascending_diagonal_stripe_test_counter + 1;
				}

				if (single_ascending_diagonal_fail == 1)
				{
				ascending_diagonal_stripe_test_counter = ascending_diagonal_stripe_test_counter - 1;
				}

				convert_test_points_of_Ricker_map_to_integers_for_any_stripe_test(lattice, row_indicies_for_descending_diagonal_stripe_test, length_of_side_of_lattice, converted_test_points_for_length, test_points_to_be_converted_for_length);
				descending_diagonal_stripe_test(converted_test_points_for_length, length_of_side_of_lattice);
				if (is_this_a_descending_diagonal_stripe == 0 && single_descending_diagonal_fail == 1)
				{
					descending_diagonal_stripe_test_failed = 1;
				}
				if (is_this_a_descending_diagonal_stripe == 0)
				{
					single_descending_diagonal_fail = 1;
				}
				if (is_this_a_descending_diagonal_stripe == 1)
				{
					single_descending_diagonal_fail = 0;
				}

				if (single_descending_diagonal_fail == 0)
				{
				descending_diagonal_stripe_test_counter = descending_diagonal_stripe_test_counter + 1;
				}

				if (single_descending_diagonal_fail == 1)
				{
				descending_diagonal_stripe_test_counter = descending_diagonal_stripe_test_counter - 1;
				}
			}

			is_this_a_descending_diagonal_stripe = 0;
			if (descending_diagonal_stripe_test_counter > 0)
			{
				is_this_a_descending_diagonal_stripe = 1;
			}

			is_this_an_ascending_diagonal_stripe = 0;
			if (ascending_diagonal_stripe_test_counter > 0)
			{
				is_this_an_ascending_diagonal_stripe = 1;
			}

			is_this_a_horizontal_stripe = 0;
			if (horizontal_stripe_test_counter > 0)
			{
				is_this_a_horizontal_stripe = 1;
			}

			is_this_a_vertical_stripe = 0;
			if (vertical_stripe_test_counter > 0)
			{
				is_this_a_vertical_stripe = 1;
			}

			what_kind_of_wavy_stripe_is_it(is_this_an_ascending_diagonal_stripe, is_this_a_descending_diagonal_stripe, is_this_a_vertical_stripe, is_this_a_horizontal_stripe);
			
			end_year = year;

			break;
							
	
		}


		

	remembering_the_last_lattice(lattice, total_number_of_patches_on_lattice, last_lattice);

	update_map_then_nearest_neighbor_interaction(lattice, last_lattice, dispersal_fraction, total_number_of_patches_on_lattice, nearest_neighbors, neighbors);

	

		if ((year % skip == 0) && (year % (skip * 2) == 0))
		{
			convert_all_lattice_points(lattice, total_number_of_patches_on_lattice, converted_lattice_points_for_lattice);

			for (int i = 0; i < total_number_of_patches_on_lattice; i++)
			{
				last_last_lattice[i] = converted_lattice_points_for_lattice[i];
			}

		}
		if (year % skip == 0 && (year % (skip * 2) == skip))
		{
			convert_all_lattice_points(lattice, total_number_of_patches_on_lattice, converted_lattice_points_for_lattice);

			sum_of_difference_between_time_t_plus_1_and_time_t = 0;

			for (int i = 0; i < total_number_of_patches_on_lattice; i++)
			{
				sum_of_difference_between_time_t_plus_1_and_time_t = sum_of_difference_between_time_t_plus_1_and_time_t + abs(last_last_lattice[i] - converted_lattice_points_for_lattice[i]);

			}

			}



	}

	




	ofstream myfile("indicies_for_horizontal_stripe_test.txt");
	if (myfile.is_open())
	{
		for (int j = 0; j < length_of_side_of_lattice; j++) {
			for (int i = 0; i < length_of_side_of_lattice; i++) {
				myfile << indicies_for_horizontal_stripe_test[j][i] << " ";
			}
		}
		myfile.close();
	}
	else cout << "Unable to open file";

	myfile.open("indicies_for_vertical_stripe_test.txt");//Recording testing conditions
	if (myfile.is_open())

	{
		for (int j = 0; j < length_of_side_of_lattice; j++) {
			for (int i = 0; i < length_of_side_of_lattice; i++) {
				myfile << indicies_for_vertical_stripe_test[j][i] << " ";
			}
		}
		myfile.close();
	}
	else cout << "Unable to open file";

	myfile.open("indicies_for_descending_diagonal_stripe_test.txt");//Recording testing conditions
	if (myfile.is_open())

	{
		for (int j = 0; j < length_of_side_of_lattice; j++) {
			for (int i = 0; i < length_of_side_of_lattice; i++) {
				myfile << indicies_for_descending_diagonal_stripe_test[j][i] << " ";
			}
		}
		myfile.close();
	}
	else cout << "Unable to open file";

	myfile.open("indicies_for_ascending_diagonal_stripe_test.txt");//Recording testing conditions
	if (myfile.is_open())

	{
		for (int j = 0; j < length_of_side_of_lattice; j++) {
			for (int i = 0; i < length_of_side_of_lattice; i++) {
				myfile << indicies_for_ascending_diagonal_stripe_test[j][i] << " ";
			}
		}
		

		myfile.close();
	}
	else cout << "Unable to open file";

	myfile.open("indicies_for_start_of_horizontal_stripe.txt");//Recording testing conditions
	if (myfile.is_open())

	{
		
			for (int i = 0; i < length_of_side_of_lattice; i++) {
				myfile << indicies_for_start_of_horizontal_stripe[i] << " ";
			}
		
		

		myfile.close();
	}
	else cout << "Unable to open file";

	myfile.open("indicies_for_start_of_vertical_stripe.txt");//Recording testing conditions
	if (myfile.is_open())

	{
		
			for (int i = 0; i < length_of_side_of_lattice; i++) {
				myfile << indicies_for_start_of_vertical_stripe[i] << " ";
			}
		
		

		myfile.close();
	}
	else cout << "Unable to open file";

	myfile.open("indicies_for_start_of_ascending_diagonal_stripe.txt");//Recording testing conditions
	if (myfile.is_open())

	{
		
			for (int i = 0; i < length_of_side_of_lattice; i++) {
				myfile << indicies_for_start_of_ascending_diagonal_stripe[i] << " ";
			}
		
		

		myfile.close();
	}
	else cout << "Unable to open file";

	myfile.open("indicies_for_start_of_descending_diagonal_stripe.txt");//Recording testing conditions
	if (myfile.is_open())

	{
		
			for (int i = 0; i < length_of_side_of_lattice; i++) {
				myfile << indicies_for_start_of_descending_diagonal_stripe[i] << " ";
			}
		
		

		myfile.close();
	}
	else cout << "Unable to open file";

	string a=argv[3];

	myfile.open("notes"+ a +".txt");//Recording testing conditions
	if (myfile.is_open())
	{
		myfile << length_of_side_of_lattice << " ";
		myfile << reproductive_rate << " ";
		myfile << radio << " ";
		myfile << dispersal_fraction << " ";
		myfile << test << " ";
		myfile << man_seed << " ";
		myfile << sum_of_difference_between_time_t_plus_1_and_time_t << " ";
		myfile << end_year << " ";
		myfile << wavy_stripe_test_result << " ";
		myfile << stripe_test_result << endl;
		myfile.close();
	}
	else cout << "Unable to open file";


	myfile.open("notes_in_words.txt");//Recording testing conditions
	if (myfile.is_open())
	{
		myfile << "length_of_side_of_lattice.\n";
		myfile << "r.\n  ";
		myfile << "radio.\n  ";
		myfile << "dispersal_fraction.\n ";
		myfile << "test.\n ";
		myfile << "seed.\n";
		myfile << "sum_of_difference_between_time_t_plus_1_and_time_t\n";
		myfile << "year that final state is reached.\n";
		myfile << "wavy_stripe_test_result\n";
		myfile << "stripe_test_result" << endl;
		myfile.close();
	}
	else cout << "Unable to open file";

	myfile.open("stripe_conditions.txt");
		if (myfile.is_open())
	{
		myfile << "0 means this test was never triggered.\n";
		myfile << "1 means this is an ascending diagonal stripe.\n";
		myfile << "2 means this is a horizontal stripe.\n";
		myfile << "3 means this is a vertical stripe.\n";
		myfile << "4 means this is an descending diagonal stripe.\n";
		myfile << "5 means this is neither a striped nor a coherent state.\n";
		myfile << "6 means this is a coherent state.\n";
		myfile.close();
	}
	else cout << "Unable to open file";

	myfile.open("converted_lattice_points_for_lattice.csv");
	if (myfile.is_open())
	{
		for (int i = 0; i < total_number_of_patches_on_lattice; i++) {
			myfile << converted_lattice_points_for_lattice[i] << " ";
		}

		myfile.close();
	}
	else cout << "Unable to open file";

	return 0;

}




