#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/LU> 
#include <unsupported/Eigen/MatrixFunctions>

#define error_maximo 0.00001

using namespace Eigen;
using namespace std;

int N_T;  //N
int N_0;  //valor inicial
int N_1;  //fin
int N_L;  //N future steps 
int S;

double T;
double z_c;
double g;

		
double Gpx;   //Gp control previo mejorado
double Gpy;

double error_estimado;
double P_k;
double P_km1;


/* devuelve "a - b" en segundos */
double timeval_diff(struct timeval *a, struct timeval *b)
{
  return
    (double)(a->tv_sec + (double)a->tv_usec/1000000) -
    (double)(b->tv_sec + (double)b->tv_usec/1000000);
}


int main(int argc, char *argv[])
{
  if (argc != 6){
    cout << "Uso: ZMP_CM N_T N N_0 N_1 T\n";
	return 0;
  }
  else
    {
      N_T = atoi(argv[1]);
      N_L = atoi(argv[2]);
	  N_0 = atoi(argv[3]);
	  N_1 = atoi(argv[4]);
	  T   = atof(argv[5]);
      cout << N_T << ":" << N_L << ":" << N_0 << ":" << N_1 << ":" << T <<endl;
	  
    }	
	
	
  struct timeval t_ini, t_fin;
  double secs;
  
  //T = 0.005;
  z_c = 0.814;
  g = 9.8;
  
  VectorXd x(3);
  VectorXd y(3);
  
  double p_x;      //error en el tamaño del vector
  double p_y;
  double u_x;
  double u_y;
  
  double sum_ex_i;
  double sum_ey_i;
  
  MatrixXd A(3,3);
  MatrixXd B(3,1);
  MatrixXd C(1,3);
  
  MatrixXd Ap(4,4);
  MatrixXd Bp(4,1);
  MatrixXd Cp(1,4);
  
  MatrixXd Ac(4,4);
  MatrixXd Ip(4,1);
  MatrixXd Gp(4,N_L);
  MatrixXd Xp(4,N_L);
  
  
  MatrixXd Q = MatrixXd::Identity(1,1);
  MatrixXd R = MatrixXd::Identity(1,1);
  Q(0,0) = 1;
  R(0,0) = 1.0e-6;
  
  MatrixXd Qp = MatrixXd::Identity(4,4); //solo para probar
  
  MatrixXd Pp = MatrixXd::Identity(4,4);
  
  MatrixXd Kp(4,4);
  
  VectorXd px_kref(N_T);
  VectorXd py_kref(N_T);
  
  VectorXd px_error(N_T);
  VectorXd py_error(N_T);
  
  VectorXd Fp(N_L);
  //VectorXd G(N_L);
  
  VectorXd px_kref_temp(N_L);
  VectorXd py_kref_temp(N_L);
  
  //MatrixXd s_xy(2,8);
  
  //***************************************************************************************************** 
  //lee archivo de pasos, las posiciones se guardan en s_x y s_y
  
  ifstream in("s_xy.dat"); 
  //s_xy << 0.0, 0.3, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0,
          //0.1, 0.2, 0.2, 0.2, 0.1, 0.0, 0.0, 0.0;

  vector<double> s_xs;
  vector<double> s_ys;
  double s_x, s_y;

  while(in >> s_x){
	s_xs.push_back(s_x);
        in >> s_y;  
	s_ys.push_back(s_y);	
  }
	  
  S = s_xs.size();
  cout << S << endl;
  //matriz temporal... prueba
  MatrixXd s_xy(2,S);
	  
  for(int i = 0; i < S; i++){
	  s_xy(0,i) = s_xs[i];
	  s_xy(1,i) = s_ys[i];
	  cout << s_xs[i] << "\t" << s_ys[i] << endl;
	 }
	  
  MatrixXd p_xy(2,S+1);
  
  p_xy(0,0)=0.0;
  p_xy(1,0)=0.0;
  
  //definición de la posición de referencia del ZMP
  
  for(int n=1; n<S+1; n=n+1){
	  p_xy(0,n)=p_xy(0,n-1)+s_xy(0,n-1);                 // genera p_x
	  p_xy(1,n)=p_xy(1,n-1)-(n%2?-1:1)*s_xy(1,n-1);      // genera p_y
     }
  cout << "ZMP(x,y) = " << endl << p_xy << endl;
  
  //genera los puntos pk_ref (ZMP)
  ofstream out("ZMPxy.dat");
  
  for(int k=0; k<N_T; k=k+1){
      px_kref(k) = p_xy(0,round(k/(N_T/S)));
	  py_kref(k) = p_xy(1,round(k/(N_T/S)));
      out << k+1 << "\t"<< px_kref(k) << "\t"<< py_kref(k) <<"\n";	  
	 }
	  
	  //cambie de k a k+1 en el archivo de salida, el pk cae sobre el pk_ref
  
  //***********************************************************************************************
  //matrices y vectores que definen la representacion en el espacio de estados
 
  //vectores estado x y y (posicion, velocidad, aceleración)
  x << 0, 0, 0;
  y << 0, 0, 0;
	   
  A << 1, T, T*T/2,
       0, 1, T,
	   0, 0, 1;
  
  B << T*T*T/6, 
       T*T/2,
       T;
	   
  C << 1, 0, -z_c/g;	  
  
  
  //ZMP mejorado
  
  Ap << 1, C*A,
        MatrixXd::Zero(3, 1), A;
		
  cout << "Ap: " << endl <<  Ap << endl;
	   
  Bp << C*B, B;
  
  cout << "Bp: " << endl <<  Bp << endl;

  Cp << 1, 0, 0, 0; 
  
  
	  
  Pp = Ap.transpose()*Pp*Ap + Cp.transpose()*Q*Cp - Ap.transpose()*Pp*Bp*(R + Bp.transpose()*Pp*Bp).inverse()*Bp.transpose()*Pp*Ap;	  
  P_km1 = Pp.diagonal().transpose()*Pp.diagonal();
	  
  for(int k = 2; ;k++){	  
	  Pp = Ap.transpose()*Pp*Ap + Cp.transpose()*Q*Cp - Ap.transpose()*Pp*Bp*(R + Bp.transpose()*Pp*Bp).inverse()*Bp.transpose()*Pp*Ap;	  
	  P_k = Pp.diagonal().transpose()*Pp.diagonal();      //producto punto de la diagonal para calcular el error estimado
	  
	  error_estimado = fabs(P_km1 - P_k);
	  if(error_estimado <= error_maximo){
		  cout <<  endl << k << " iteraciones" << endl;
		  break;
		  }
		  P_km1 = P_k;
	  }	   
	  
  cout << "Pp: " << endl <<  Pp << endl;
  
  Kp = (R + Bp.transpose()*Pp*Bp).inverse()*Bp.transpose()*Pp*Ap;  
	  
  cout << "Kp: " << endl <<  Kp << endl;
	  
  //calcula los pesos f_i
  for(int i=1; i<N_L+1; i++){	  
	  Fp(i-1) = ((R + Bp.transpose()*Pp*Bp).inverse()*Bp.transpose()*((Ap-Bp*Kp).transpose()).pow(i-1)*Cp.transpose()*Q).value();
	  }	  
	  
  //archivos para datos de salida, CoM y ZMP
	  
  ofstream out_COMx_imp("COMx_imp.dat"); 
  ofstream out_COMy_imp("COMy_imp.dat"); 
	  
  ofstream out_px_imp("Px_imp.dat"); 
  ofstream out_py_imp("Py_imp.dat");

  ofstream out_Gp("Gp.dat"); 
	
  //tiempo de inicio
  gettimeofday(&t_ini, NULL);

//Control Previo mejorado

  cout << "CONTROL PREVIO MEJORADO... " << endl;

  Ac = Ap - Bp*Kp;

  cout << "Ac: " << endl <<  Ac << endl;

  Ip << 1, 0 ,0, 0;

  cout << "Ip: " << endl << Ip << endl;

Gp.block<4,1>(0,0) = (-Kp.block<1,1>(0,0)).value()*Ip;   //Gp(1)
out_Gp << 1 << "\t"<< Gp.block<1,1>(0,0) <<"\n";

cout << "Gp(1): " << endl << Gp.block<4,1>(0,0) << endl;

Xp.block<4,1>(0,0) = -Ac.transpose()*Pp*Ip;				 //Xp(1)

cout << "Xp(1): " << endl << Xp.block<4,1>(0,0) << endl;

for(int l=2; l<=N_L; l=l+1){	
	Gp.block<1,1>(0,l-1) = (R + Bp.transpose()*Pp*Bp).inverse()*Bp.transpose()*Xp.block<4,1>(0,l-2);
	if((Gp.block<1,1>(0,l-1)).value() > 0)Gp.block<1,1>(0,l-1) = -Gp.block<1,1>(0,l-1);
	Xp.block<4,1>(0,l-1) = -Ac.transpose()*Xp.block<4,1>(0,l-2);
	out_Gp << l << "\t"<< Gp.block<1,1>(0,l-1) <<"\n";
	}
	
cout << "Xp: " << endl << Xp.block<4,5>(0,0) << endl;

cout << "Gp: " << endl << Gp.block<4,5>(0,0) << endl;


for(int k=N_0; k<N_1; k=k+1){
	
	p_x = (C*x).value();	  
	out_px_imp << k << "\t"<< p_x <<"\n";
	  
	p_y = (C*y).value();
	out_py_imp << k << "\t"<< p_y << "\n"; 
	  
	px_error(k-N_0) =  p_x - px_kref(k);
	py_error(k-N_0) =  p_y - py_kref(k);
	
    sum_ex_i = 0;
    sum_ey_i = 0;
	
	for(int i=0; i<k; i=i+1){
	    sum_ex_i = sum_ex_i + px_error(i);         //error acumulado 
		sum_ey_i = sum_ey_i + py_error(i);
		}
	
	Gpx = 0;
	Gpy = 0;
	
	for(int j=0; j<N_L; j=j+1){
		Gpx = Gpx + px_kref(j+k)*(Gp.block<1,1>(0,j)).value();
		Gpy = Gpy + py_kref(j+k)*(Gp.block<1,1>(0,j)).value();
		}  
	
	  
	  u_x = (-Kp.block<1,1>(0,0)*sum_ex_i-Kp.block<1,3>(0,1)*x).value() - Gpx;
	  //u_x = (-Kp.block<1,3>(0,1)*x).value() - Gpx;
	  
	  u_y = (-Kp.block<1,1>(0,0)*sum_ey_i-Kp.block<1,3>(0,1)*y).value() - Gpy;
	  //u_y = (-Kp.block<1,3>(0,1)*y).value() - Gpy;
	  
	  x = A*x + B*u_x;
	  out_COMx_imp << k << "\t"<< x(0) <<"\n";
	  
	  y = A*y + B*u_y;
	  out_COMy_imp << k << "\t"<< y(0) << "\n"; 
  }
  
  //*******************************************************************************************************
  //Calcula el tiempo de ejecución

  gettimeofday(&t_fin, NULL);

  secs = timeval_diff(&t_fin, &t_ini);
  printf("%.16g milliseconds\n", secs * 1000.0);
  //*******************************************************************************************************
}

