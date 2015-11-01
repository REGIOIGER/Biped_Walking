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

int N_T;  	// N 
int N_0;  	// extremo izquiero del intervalo 
int N_1;  	// extremo derecho del intervalo 
int N_L;  	// N future steps 
int S;

double T; 	// sampling time
double z_c;	// diastancia del piso al CoM
double g;   // gravedad

		
double Gpx; // Gp ganancia usando valores incrementales para control previo mejorado
double Gpy;

double error_estimado; // error estimado en el algoritmo de Punto Fijo
double P_k;            // variables de paso para calcular el error estimado 
double P_km1;

 
// estimacion de tiempo de ejecuciçon 
                                             
double timeval_diff(struct timeval *a, struct timeval *b)  // devuelve "a - b" en segundos
{
  return
    (double)(a->tv_sec + (double)a->tv_usec/1000000) -
    (double)(b->tv_sec + (double)b->tv_usec/1000000);
}


// comienza main
// el proposito es convertirlo en una clase que genere objetos trayectoria del CoM a partir de los parámetros
// de control que se encontraran en un archivo, incluyendo la definición de pasos (huellas) del robot bipedo.

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
  
  z_c = 0.814;
  g = 9.8;
  
  VectorXd x(3);	// vectores de estado
  VectorXd y(3);
  
  double p_x;		// salida del contorl para los puntos de referencia que originan las posiciones de los ZMP.
  double p_y;
  double u_x;		// jerk usado como entrada de control
  double u_y;
  
  double sum_ex_i;	// error acumulado de los puntos ZMP de referencia y los generados por el control
  double sum_ey_i;
  
  MatrixXd A(3,3);	// matrices de la ecuacion del ZMP representada en el espacio de estados 
  MatrixXd B(3,1);
  MatrixXd C(1,3);
  
  MatrixXd Ap(4,4);	// matrices de la ecuacion del ZMP en representación en el espacio de estados 
  MatrixXd Bp(4,1);	// en control previo mejorado
  MatrixXd Cp(1,4);
  
  MatrixXd Ac(4,4); //matrices auxiliares para control previo mejorado
  MatrixXd Ip(4,1);
  MatrixXd Gp(4,N_L);
  MatrixXd Xp(4,N_L);
  
  
  MatrixXd Q = MatrixXd::Identity(1,1); // Matrices peso usadas en control previo
  MatrixXd R = MatrixXd::Identity(1,1);
  Q(0,0) = 1;
  R(0,0) = 1.0e-6;
  
  MatrixXd Qp = MatrixXd::Identity(4,4); //solo para probar
  
  MatrixXd Pp = MatrixXd::Identity(4,4); // solución para Pp en la ecuación algebraica de Riccati 
										 // que origina el cálculo de ganancias para optimizar el control
  MatrixXd Kp(4,4);						 // matriz de ganancias para control previo mejorado 
  
  VectorXd px_kref(N_T);				 // puntos del ZMP
  VectorXd py_kref(N_T);
  
  VectorXd px_error(N_T);				 // errores de seguimiento del ZMP
  VectorXd py_error(N_T);
  
  //VectorXd Fp(N_L);					 // fi - pesos de control previo mejorado (no utilizados)
  
  VectorXd px_kref_temp(N_L);			 // puntos del ZMP de referencia auxiliares
  VectorXd py_kref_temp(N_L);
  
  
  //***************************************************************************************************** 
  //lee archivo de pasos (huellas), las posiciones se guardan en s_x y s_y
  
						
  ifstream in("s_xy.dat");	//archivo de texto de entrada

  vector<double> s_xs;		// vectores auxiliares
  vector<double> s_ys;
  
  double s_x, s_y;

  while(in >> s_x){			// lee mientras haya datos en el archivo de texto de entrada
	s_xs.push_back(s_x);
        in >> s_y;  
	s_ys.push_back(s_y);	
  }
	  
  S = s_xs.size();			// número de pasos (incluidos pasos con incremento cero
  
  cout << S << endl;	    // salida por pantalla

  MatrixXd s_xy(2,S);		// matriz cargada con incrementos de los ZMP
	  
  for(int i = 0; i < S; i++){
	  s_xy(0,i) = s_xs[i];
	  s_xy(1,i) = s_ys[i];
	  cout << s_xs[i] << "\t" << s_ys[i] << endl;	// salida por pantalla
	 }
	  
	  
  MatrixXd p_xy(2,S+1);		// matriz de las posiciones (x,y) de los ZMP
  
  // definición de la posición de referencia del ZMP
  
  p_xy(0,0)=0.0;	// posición inicial (en una forma mas general no tiene que ser (0,0), por ejemplo si inicia
  p_xy(1,0)=0.0;	// desde una posición con un pie en el aire

  
  for(int n=1; n<S+1; n=n+1){
	  p_xy(0,n)=p_xy(0,n-1)+s_xy(0,n-1);                  
	  p_xy(1,n)=p_xy(1,n-1)-(n%2?-1:1)*s_xy(1,n-1);      
     }
  cout << "ZMP(x,y) = " << endl << p_xy << endl;
  
  ofstream out("ZMPxy.dat");							
  
  for(int k=0; k<N_T; k=k+1){
      px_kref(k) = p_xy(0,round(k/(N_T/S)));						// interpolación de puntos para generar p_x
	  py_kref(k) = p_xy(1,round(k/(N_T/S)));						// y py de referencia
      out << k+1 << "\t"<< px_kref(k) << "\t"<< py_kref(k) <<"\n";	  
	 }
	  
  
  //cambie de k a k+1 en el archivo de salida, el pk cae sobre el pk_ref
  
  //***********************************************************************************************
  //matrices y vectores que definen la representacion en el espacio de estados
 
  
  x << 0, 0, 0;							// vectores estado x y y (posicion, velocidad, aceleración)
  y << 0, 0, 0;
	   
  A << 1, T, T*T/2,						// inicialización de A, B y C
       0, 1, T,
	   0, 0, 1;
  
  B << T*T*T/6, 
       T*T/2,
       T;
	   
  C << 1, 0, -z_c/g;	  
  
  
  
  
  Ap << 1, C*A,								// inicialización de matrices para control previo mejorado
        MatrixXd::Zero(3, 1), A;
		
  cout << "Ap: " << endl <<  Ap << endl;
	   
  Bp << C*B, B;
  
  cout << "Bp: " << endl <<  Bp << endl;

  Cp << 1, 0, 0, 0; 
  
  
  // algoritmo de Punto Fijo para resolver la ecuacion algebraica de Riccati para Pp
  
  // primera aproximación
	  
  Pp = Ap.transpose()*Pp*Ap + Cp.transpose()*Q*Cp - Ap.transpose()*Pp*Bp*(R + Bp.transpose()*Pp*Bp).inverse()*Bp.transpose()*Pp*Ap;	  
  P_km1 = Pp.diagonal().transpose()*Pp.diagonal();
  
  // inicia iteracciones, para cuando se alcance el error maximo permitido
  
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
	  
  cout << "Pp: " << endl <<  Pp << endl;	// salida por pantalla de Pp
  
  Kp = (R + Bp.transpose()*Pp*Bp).inverse()*Bp.transpose()*Pp*Ap;  
	  
  cout << "Kp: " << endl <<  Kp << endl;	// salida por pantalla de Kp (ganancias para el control del error
											// de seguimiento del ZMP y la retroalimentación del estado
  
  for(int i=1; i<N_L+1; i++){	  			// calcula los pesos Fp_i (*************NO UTILIZADOS*************)
	  Fp(i-1) = ((R + Bp.transpose()*Pp*Bp).inverse()*Bp.transpose()*((Ap-Bp*Kp).transpose()).pow(i-1)*Cp.transpose()*Q).value();
	  }	  
	  
	  
  ofstream out_COMx_imp("COMx_imp.dat"); 	//archivos para datos de salida, CoM y ZMP
  ofstream out_COMy_imp("COMy_imp.dat"); 
	  
  ofstream out_px_imp("Px_imp.dat"); 
  ofstream out_py_imp("Py_imp.dat");

  ofstream out_Gp("Gp.dat"); 
	
  gettimeofday(&t_ini, NULL);				// tiempo de inicio

  // *********************************************************************************************************
  // Control Previo mejorado

  Ac = Ap - Bp*Kp;											// calculo de Gp (ganancia previa)

  cout << "Ac: " << endl <<  Ac << endl;

  Ip << 1, 0 ,0, 0;

  cout << "Ip: " << endl << Ip << endl;
  
  // se calculan la ganancia previa por recursión

  Gp.block<4,1>(0,0) = (-Kp.block<1,1>(0,0)).value()*Ip;     // Gp(1) 
  out_Gp << 1 << "\t"<< Gp.block<1,1>(0,0) <<"\n";

  cout << "Gp(1): " << endl << Gp.block<4,1>(0,0) << endl;

  Xp.block<4,1>(0,0) = -Ac.transpose()*Pp*Ip;				 // Xp(1)

  cout << "Xp(1): " << endl << Xp.block<4,1>(0,0) << endl;

  for(int l=2; l<=N_L; l=l+1){	
	 Gp.block<1,1>(0,l-1) = (R + Bp.transpose()*Pp*Bp).inverse()*Bp.transpose()*Xp.block<4,1>(0,l-2);
	 
	 // AQUI FORCE EL SIGNO A SER NEGATIVO, AUN NO ENCUENTRO EL ERROR
	 if((Gp.block<1,1>(0,l-1)).value() > 0)Gp.block<1,1>(0,l-1) = -Gp.block<1,1>(0,l-1);
	 
	 Xp.block<4,1>(0,l-1) = -Ac.transpose()*Xp.block<4,1>(0,l-2);
	 out_Gp << l << "\t"<< Gp.block<1,1>(0,l-1) <<"\n";
	 }
	
  cout << "Xp: " << endl << Xp.block<4,5>(0,0) << endl;		//salida de los 5 primeros valores de Gp

  cout << "Gp: " << endl << Gp.block<4,5>(0,0) << endl;


  for(int k=N_0; k<N_1; k=k+1){
	
	p_x = (C*x).value();	  					// salida p_x del control
	out_px_imp << k << "\t"<< p_x <<"\n";
	  
	p_y = (C*y).value();						// salida p_y del control
	out_py_imp << k << "\t"<< p_y << "\n"; 
	  
	px_error(k-N_0) =  p_x - px_kref(k);		// error de seguimiento del ZMP
	py_error(k-N_0) =  p_y - py_kref(k);
	
    sum_ex_i = 0;								// se reinicia cuenta del error acumulado
    sum_ey_i = 0;
	
	for(int i=0; i<k; i=i+1){
	    sum_ex_i = sum_ex_i + px_error(i);      // error acumulado de seguimiento del ZMP
		sum_ey_i = sum_ey_i + py_error(i);
		}
	
	Gpx = 0;									// se reicnicia valor del termino de control previo en u(k)
	Gpy = 0;
	
	for(int j=0; j<N_L; j=j+1){
		Gpx = Gpx + px_kref(j+k)*(Gp.block<1,1>(0,j)).value();
		Gpy = Gpy + py_kref(j+k)*(Gp.block<1,1>(0,j)).value();
		}  
	
	  
	  u_x = (-Kp.block<1,1>(0,0)*sum_ex_i-Kp.block<1,3>(0,1)*x).value() - Gpx;		// u(k)
	  //u_x = (-Kp.block<1,3>(0,1)*x).value() - Gpx;
	  
	  u_y = (-Kp.block<1,1>(0,0)*sum_ey_i-Kp.block<1,3>(0,1)*y).value() - Gpy;
	  //u_y = (-Kp.block<1,3>(0,1)*y).value() - Gpy;
	  
	  x = A*x + B*u_x;																// actualización de los estados
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

