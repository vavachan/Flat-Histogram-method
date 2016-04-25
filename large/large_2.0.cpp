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
			P_fnc=P_fnc+exp(-beta*(E(n)+3*EN))*omega[n][j];//*exp(-EN*(R(m,beta,h)-R(m,beta,0)));//*pow(cosh(beta*h),EN);///omega_m[j];
//			cout<<n<<"\t"<<j<<"\t"<<omega[n][j]*exp(-beta*(E(n)+3*EN))<<"\n";//<<"\t"<<omega[n][j]/omega_m[j]<<"\t"<<exp(-EN*R(m,beta,h))<<"\n";
			}
	//cout<<P_fnc<<"\n";
	return P_fnc;
	}
long double P(int n,float beta=beta)
	{
	long double P_m=0;
	int N_plus;
	float m=0;
	for(int i=0;i<EN+1;i++)
		if(omega[n][i]!=0)
		{
		m=(2*i-EN)/EN;
		P_m=P_m+sqrt(m*m)*omega[n][i]*exp(-beta*(E(n)+3*EN));//*exp(-EN*(R(m,beta)-R(m,beta,0)));//*pow(cosh(beta*h),EN);///omega_m[N_plus];
//		cout<<exp(-beta*(E(i)+3*EN))*omega[i][N_plus]*exp(-EN*R(m,beta,h)-R(m,beta,0))<<"\t"<<exp(-beta*(E(i)+3*EN))*omega[i][N_plus]<<"\t"<<exp(-EN*R(m,beta,h)-R(m,beta,0))<<"\n";
//		cout<<n<<"\t"<<i<<"\t"<<P_m<<"\n";//omega[i][N_plus]*exp(-beta*E(i))<<"\t"<<exp(-EN*R(m,beta))/omega_m[N_plus]<<"\n";
		}
	return (P_m)/Partition_function;
	}	
////////float spec_heat(float beta=beta,float n_min=0,float n_max=3*EN+1)
////////	{
////////	float delta_beta=.001;
////////	float der_P1=0,der_P2=0,der_2P=0;
////////	der_P1=(log(Part_func(beta+delta_beta))-log(Part_func(beta)))/delta_beta;
////////	der_P2=(log(Part_func(beta+2*delta_beta))-log(Part_func(beta+delta_beta)))/delta_beta;
////////	der_2P=(der_P2-der_P1)/delta_beta;
////////	return der_2P*beta*beta;
////////	}
////////float spec_heat_2(float beta=beta,long double Part=0)	
////////	{
////////	long double P_fnc_1=0,P_fnc_2=0,P_fnc=0;
////////	for(int n=0;n<3*EN+1;n++)
////////		for(int j=0;j<EN+1;j++)
////////				if(omega[n][j]!=0)
////////				{
////////				P_fnc=P_fnc+exp(-beta*(E(n)+3*EN))*omega[n][j];
////////				P_fnc_1=P_fnc_1+exp(-beta*(E(n)+3*EN))*omega[n][j]*(E(n)+3*EN)*(E(n)+3*EN);
////////				P_fnc_2=P_fnc_2+exp(-beta*(E(n)+3*EN))*omega[n][j]*(E(n)+3*EN);
////////				}
////////		return beta*beta*((P_fnc_1/P_fnc)-P_fnc_2*P_fnc_2/(P_fnc*P_fnc));
////////}
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
	long double sum=0;	
        for(int i=0;i<3*EN+1;i++)
        	{
        	//cout<<m<<"\t"<<P(m,beta)<<"\n";
        	sum=sum+P(i,beta);
//		cout<<EN*sqrt(m*m)<<"\t"<<sum<<"\n";
      //  	sum=sum+P(m,beta);
        	}		
	cout<<1/beta<<"\t"<<sum<<"\t"<<Partition_function<<"\n";
	//cout<<R(0.5,beta)<<"\n";
	//for(float x=-1+.01;x<1;x=x+.01)
	//	cout<<x<<"\t"<<R(x,beta)<<"\n";
	//for(int i=0;i<EN+1;i++)
	//	cout<<2*i-EN<<"\t"<<omega_m[i]<<"\n";
	
	}
