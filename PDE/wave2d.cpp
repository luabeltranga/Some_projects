#include <iostream>
#include <cmath>
#include <vector>

const double RHO = 0.01;
const double TENSION = 40;
const double MAX = 10000;
const double C = std::sqrt(TENSION/RHO);
const double C1 = 3*C ;
const double RATIO = C*C/(C1*C1);
const int SIZE = 101;


void initial_conditions(std::vector<double> &x);
void first_chance(std::vector<double> &x, int &n);
void advance(std::vector<double> &x);
void print_gnuplot(std::vector<double> &x ,int &n);
void start_gnuplot(void);
bool relax(std::vector<double> &x,std::vector<double> &x_tmp,int n);

int main(void){
  std::vector<double> x (SIZE*3);
  int graph = 0;
  start_gnuplot();
  initial_conditions(x);
  print_gnuplot(x,graph);
  graph = 1;
  first_chance(x,graph);
  print_gnuplot(x,graph);
  graph = 2;
  for(int jj = 0; jj < MAX; jj++){
    advance(x);
    first_chance(x,graph);
    print_gnuplot(x,graph);
  }
  return 0;
}

void initial_conditions(std::vector<double> &x){
  //for(int ii = 0; ii < 81; ii++) x[ii+0*SIZE] = 0.00125*ii;
  //for(int ii = 81; ii < 101; ii++) x[ii+0*SIZE] = 0.1-0.005*(ii-80);
  for(int ii = 0; ii < 101; ii++) x[ii+(0*SIZE)] = 0.001*std::sin(2*M_PI*ii/100);
}

void first_chance(std::vector<double> &x, int &n){
  std::vector<double> x_tmp = x ;
  int row =1;
  do{
    x_tmp = x;
    for(int ii = 1; ii < 100; ii++){
      x[ii+(n*SIZE)]= x[ii+((n-1)*SIZE)] + 0.5*RATIO*(x[ii+1+((n-1)*SIZE)]+x[ii-1+((n-1)*SIZE)]-2*x[ii+((n-1)*SIZE)]);
    }
  }while(relax(x,x_tmp,row));
}

void advance(std::vector<double> &x){
     for(int ii = 1; ii < 100; ii++){
       x[ii+(2*SIZE)]= 2*x[ii+(1*SIZE)]-x[ii+(0*SIZE)] + RATIO*(x[ii+1+(1*SIZE)]+x[ii-1+(1*SIZE)]-2*x[ii+(1*SIZE)]);    }
     
     for(int kk = 0; kk < 101; kk++){
       x[kk+(0*SIZE)] = x[kk+(1*SIZE)];
       x[kk+(1*SIZE)] = x[kk+(2*SIZE)];
     }    
}
     
void print_gnuplot(std::vector<double> &x ,int &n){
  std::cout << "splot   '-'  " << std::endl;
  for(int ii = 0; ii<101 ; ii++){
    std::cout <<ii*0.01 <<" " << ii*0.01<<" "<< x[n*101+ii] << std::endl;
  }
  std::cout << "e" << std::endl;
}

void start_gnuplot(void){
  std::cout << "set terminal gif animate" << std::endl;
  std::cout << "set out 'wave1d.gif' " << std::endl;
}

bool relax(std::vector<double> &x,std::vector<double> &x_tmp,int n){
  for(int ii = 1 ; ii < 100 ; ii++){
    if((std::fabs(x_tmp[ii+n*SIZE]-x[ii+n*SIZE])/x_tmp[ii+n*SIZE])>1.0e-7) return true;
  }
  return false;
}
