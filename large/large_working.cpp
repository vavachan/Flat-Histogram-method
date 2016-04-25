using namespace std;
#include<iostream>
#include<math.h>
#include<fstream>
#include<stdlib.h>
float EN,beta,h;
long double omega[10000][10000]={0},omega_m[10000]={0},Partition_function;
float f(float x,float beta,float h=2.0)
	{
	return 0.5*log(cosh(2*x)+cosh(2*beta*h));
	}
long double max_y(float x,float beta=beta,float h=h)
	{
	long double a=sinh(2*beta*h);
	long double b=cosh(2*beta*h);
	//cout<<"a\t"<<a<<"b\t"<<b<<"\n";
	if(x==1 or x==-1)	
		{
		if(x<0)	x=x+.001;
		if(x>0)	x=x-.001;
		return 0.5*asinh(x/(1-x*x)*(b+sqrt(1+x*x*a*a)));
		}
	return 0.5*asinh(x/(1-x*x)*(b+sqrt(1+x*x*a*a)));
	}
float R(float x,float beta=beta,float h=h)
	{
	long double y_max=max_y(x,beta,h);
//	cout<<y_max<<"y_max\n";
	//return 0;
	return x*y_max-f(y_max,beta,h)+f(0,beta,h);
	}
float E(int n)
	{
	return 3*EN-2*n;
	}
long double Part_func(float beta=beta,float h=h,float n_min=0,float n_max=3*EN+1)
	{
	long double P_fnc=0;float m=0;
	for(int n=n_min;n<n_max;n++)
		for(int j=0;j<EN+1;j++)
			if(omega[n][j]!=0)
			{
			m=(2*j-EN)/EN;	
			P_fnc=P_fnc+exp(-beta*(E(n)+3*EN))*omega[n][j]*exp(-EN*(R(m,beta,h)-R(m,beta,0)));//*pow(cosh(beta*h),EN);///omega_m[j];
			}
	return P_fnc;
	}
long double P(int N_plus,float beta=beta)
	{
	long double P_m=0;
	float m=(2*N_plus-EN)/EN;
	for(int i=0;i<3*EN+1;i++)
		if(omega[i][N_plus]!=0)
		{
		P_m=P_m+omega[i][N_plus]*exp(-beta*(E(i)+3*EN))*exp(-EN*(R(m,beta)-R(m,beta,0)));//*pow(cosh(beta*h),EN);///omega_m[N_plus];
	//	cout<<i<<"\t"<<P_m<<"\t"<<omega[i][N_plus]<<"\t"<<exp(-beta*(E(i)+3*EN))<<"\t"<<exp(-EN*(R(m,beta)-R(m,beta,0)))<<"\n";
		}
	//	cout<<m<<"\t"<<exp(-EN*(R(m,beta)-R(m,beta,0)))<<"\t"<<Partition_function<<"\n";
	return (P_m)/Partition_function;
	}	
void readfile(int L)
	{
	fstream infile;
	long double a,b,c,d,N_plus;
	  char buffer[32];
	    snprintf(buffer,sizeof(char)*32,"entropy_2para_%i.dat",L);
	    infile.open(buffer);

	for(int i=0;!infile.eof();i++)
		{
		infile>>a>>b>>c>>d;
//		cout<<a<<"\t"<<b<<"\t"<<d<<"\n";
		N_plus=(b+EN)/2;	
		omega[int(a)][int(N_plus)]=d;
		}	
	for(int n=0;n<EN+1;n++)
		for(int i=0;i<3*EN+1;i++)
			omega_m[n]=omega_m[n]+omega[i][n];			
	}
int main(int argc,char* argv[])
	{
	int L=atof(argv[1]);
	float m;
	beta=1/atof(argv[2]);
	h=atof(argv[3]);
	EN=L*L*L;
//	cout<<EN<<"\n";
	readfile(L);
	Partition_function=Part_func(beta);
//	cout<<spec_heat_2(beta)<<"\n";
//	for(float T=1.0;T<8.0;T=T+.02)
//		cout<<T<<"\t"<<spec_heat(1/T)<<"\n";	
//	cout<<P(.5,beta)<<"\n";
//	cout<<R(0.5)<<"\n";
	long double sum=0,sum1=0;	
//	cout<<P(108,beta);
      for(int i=0;i<EN+1;i++)
//      //	int i=215;
      	{
      	m=(2*i-EN)/EN;
//		sum=sum+P(i,beta);
//      //	cout<<m<<"\t"<<i<<"\n";
     //	cout<<m<<"\t"<<P(i,beta)<<"\n";
//		cout<<EN*sqrt(m*m)<<"\t"<<sum<<"\n";
       	sum=sum+sqrt(m*m)*P(i,beta);
	sum1=sum1+m*m*P(i,beta);
       	}		
	cout<<1/beta<<"\t"<<sum1-sum*sum<<"\n";//<<P_m/Partition_function<<"\n";
	//cout<<R(0.5,beta)<<"\n";
	//for(float x=-1+.01;x<1;x=x+.01)
	//	cout<<x<<"\t"<<R(x,beta)<<"\n";
	//for(int i=0;i<EN+1;i++)
	//	cout<<2*i-EN<<"\t"<<omega_m[i]<<"\n";
	
	}
