#include<iostream>
#include<masa.h>
#include<vector>

std::vector<std::vector<double>> A_matrix(int N, int dimension, int fin_order,std::string loglevel);
std::vector<std::vector<double>> D1_2nd(int N);
std::vector<std::vector<double>> D1_4th(int N);
std::vector<std::vector<double>> D2_2nd(int N);
std::vector<std::vector<double>> D2_4th(int N);
void print_matrix(std::vector<std::vector<double>> matrix, int n);
