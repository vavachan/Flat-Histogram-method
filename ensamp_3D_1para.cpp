using namespace std;
#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<fstream>
int L=4,m,n,N_plus,N_minus,m1;
int ISING[100][100][100];
//DEFINING THE INITIAL APPROXIMATION FOR OMEGA(n,m,L)
double entropy[18000];
double entropy_mp[18000][6000];
double Histogram_mp[18000][6000];
double Histogram[18000];
double Mag[18000];
float avg,sum;
float EN;
//int KB=1;
double s(float a)
{
    return -1*log(2)*((2.0/3)*a-1)*((2.0/3)*a-1);
}
int ini_dist()
{
    for(int i=0; i<3*EN+1; i++)
        entropy[i]=0;
}
int initial( int k)
{
    float p=0;
    if (k==0)
        for(int i=0; i<L; i++)
            for(int j=0; j<L; j++)
                for(int k=0; k<L; k++)
                    ISING[i][j][k]=-1;

    if(k==1)
        for(int i=0; i<L; i++)
            for(int j=0; j<L; j++)
                for(int k=0; k<L; k++)
                    ISING[i][j][k]=1;

    if(k==2)
        for(int i=0; i<L; i++)
            for(int j=0; j<L; j++)
                for(int k=0; k<L; k++)
                {
                    ISING[i][j][k]=pow(-1,i+j+k);
                }
    if(k==3)
        for(int i=0; i<L; i++)
            for(int j=0; j<L; j++)
                for(int k=0; k<L; k++)
                {
                    ISING[i][j][k]=pow(-1,i+j+k+1);
                }
    if(k==4 or k==5)
    {
        for(int i=0; i<L; i++)
            for(int j=0; j<L; j++)
                for(int k=0; k<L; k++)
                {
                    p=float(rand())/RAND_MAX;
                    if (p>0.5)
                    {
                        ISING[i][j][k]=1;
                    }
                    else
                    {
                        ISING[i][j][k]=-1;
                    }
                }
    }
    if (k==6)
    {
        for(int i=0; i<L; i++)
            for(int j=0; j<L; j++)
                for(int k=0; k<L; k++)
                {
                    ISING[i][j][k]=-1;
                }
        ISING[1][1][1]=1;
        ISING[1][2][1]=1;
    }
    if(k==7)
    {
        for(int i=0; i<L; i++)
            for(int j=0; j<L; j++)
                for(int k=0; k<L; k++)
                    ISING[i][j][k]=1;
        ISING[1][1][1]=-1;
        ISING[1][2][1]=-1;
    }
    if(k==8)
    {
        for(int i=0; i<L; i++)
            for(int j=0; j<L; j++)
                for(int k=0; k<L; k++)
                {
                    ISING[i][j][k]=pow(-1,i+j+k);
                }

        ISING[1][1][1]=-1*ISING[1][1][1];
        ISING[1][2][1]=-1*ISING[1][2][1];

    }
    if(k==9)
    {
        for(int i=0; i<L; i++)
            for(int j=0; j<L; j++)
                for(int k=0; k<L; k++)
                {
                    ISING[i][j][k]=pow(-1,i+j+k+1);
                }
        ISING[1][1][1]=-1*ISING[1][1][1];
        ISING[1][2][1]=-1*ISING[1][2][1];

    }
}
int magnetization()
{
    int sum=0;
    for(int i=0; i<L; i++)
        for(int j=0; j<L; j++)
            for(int k=0; k<L; k++)
                sum=sum+ISING[i][j][k];
    return sum;
}
int energy_n()
{
    float Energy=0;
    int Ni_plus=0,Ni_minus=0,Nj_plus=0,Nj_minus=0,Nk_plus=0,Nk_minus=0;
    for(int i=0; i<L; i++)
        for(int j=0; j<L; j++)
            for(int k=0; k<L; k++)
            {
                Ni_plus=i+1;
                Ni_minus=i-1;
                Nj_plus=j+1;
                Nj_minus=j-1;
                Nk_plus=k+1;
                Nk_minus=k-1;
                if (i==0)
                    Ni_minus=L-1;
                if (i==L-1)
                    Ni_plus=0;
                if (j==0)
                    Nj_minus=L-1;
                if (j==L-1)
                    Nj_plus=0;
                if (k==0)
                    Nk_minus=L-1;
                if (k==L-1)
                    Nk_plus=0;

                Energy=Energy+-1*(ISING[i][j][k]*ISING[i][Nj_plus][k]+ISING[i][j][k]*ISING[i][Nj_minus][k]+ISING[i][j][k]*ISING[Ni_plus][j][k]+ISING[i][j][k]*ISING[Ni_minus][j][k]+ISING[i][j][k]*ISING[i][j][Nk_plus]+ISING[i][j][k]*ISING[i][j][Nk_minus]);
            }

    //cout<<(3*EN-Energy/2)/2<<"\n";
    return (3*EN-Energy/2)/2;
}

int flip(int E)
{
    //cout<<Histogram_mp[15][0]<<"\t"<<"before"<<"\n";
    int Ni_plus=0,Ni_minus=0,Nj_plus=0,Nj_minus=0,Nk_plus=0,Nk_minus=0;
    int X,Y,Z,n_1,n_2,delta_n,m2,i,j,delta_m;
    int EN=L*L*L;
    double alpha;
    float p;
    X=rand()%L;
    Y=rand()%L;
    Z=rand()%L;
    n_1=E;
    Ni_plus=X+1;
    Ni_minus=X-1;
    Nj_plus=Y+1;
    Nj_minus=Y-1;
    Nk_plus=Z+1;
    Nk_minus=Z-1;
    if (X==0)
        Ni_minus=L-1;
    if (X==L-1)
        Ni_plus=0;
    if (Y==0)
        Nj_minus=L-1;
    if (Y==L-1)
        Nj_plus=0;
    if (Z==L-1)
        Nk_plus=0;
    if(Z==0)
        Nk_minus=L-1;
    //m1=magnetization();
    if(ISING[X][Y][Z]==-1)
    {
        delta_m=2;
    }
    else
    {
        delta_m=-2;
    }
    i=(m1+EN)/2;
    m2=m1+delta_m;
    j=(m2+EN)/2;
    delta_n=-1*(ISING[X][Y][Z]*ISING[X][Nj_plus][Z]+ISING[X][Y][Z]*ISING[X][Nj_minus][Z]+ISING[X][Y][Z]*ISING[Ni_plus][Y][Z]+ISING[X][Y][Z]*ISING[Ni_minus][Y][Z]+ISING[X][Y][Z]*ISING[X][Y][Nk_plus]+ISING[X][Y][Z]*ISING[X][Y][Nk_minus]);
    n_2=n_1+delta_n;
    alpha=exp(EN*(entropy[n_1]-entropy[n_2]));
//    cout<<n_1<<"\t"<<n_2<<"\t"<<i<<"\t"<<j<<"\t"<<alpha<<"\n";
    if (alpha<1)
    {
        p=float(rand())/RAND_MAX;
        if (p>alpha)
        {
            Histogram_mp[n_1][i]++;//=Histogram_mp[n_1][i]+1;
            return n_1;
        }
        else
        {
            Histogram_mp[n_2][j]++;//=Histogram_mp[n_2][j]+1;
            m1=m2;
            ISING[X][Y][Z]=-1*ISING[X][Y][Z];
            return n_2;
        }
    }
    else
    {
        ISING[X][Y][Z]=-1*ISING[X][Y][Z];
        Histogram_mp[n_2][j]++;//=Histogram_mp[n_2][j]+1;
        m1=m2;
        return n_2;
    }
}
float average()
{
    sum=0;
    float avg=0,count=0,count_2=0;
    for(int i=0; i<3*EN+1; i++)
    {
        if(Histogram[i]!=0)
            count_2++;
        //	cout<<Histogram[i]<<"\n";
        sum=sum+Histogram[i];
        for(int j=0; j<EN+1; j++)
            if(Histogram_mp[i][j]!=0)
            {
                count++;
                avg=avg+Histogram_mp[i][j];
            }
    }
    sum=sum/count_2;
    return avg/(count);
}
////////float his()
////////{
////////for(int i=0;i<EN+1;i++)
////////	for(int j=0;j<EN+1;j++)
////////		Histogram[i]=Histogram[i]+Histogram_mp[i][j];
////////return 0;
////////}
float mag(float t)
{
    double P_fnc=0,numer=0;
    for(int i=0; i<EN+1; i++)
        if(i!=1 and i!=99)
        {
            P_fnc=P_fnc+exp(EN*entropy[i])*exp(-1/t*2*(EN-2*i));
            numer=numer+(Mag[i]/Histogram[i])*exp(EN*entropy[i])*exp(-1/t*2*(EN-2*i));
        }
    return numer/(P_fnc);
}
int main(int argc,char* argv[])
{
    L=atof(argv[1]);
    EN=L*L*L;
    ini_dist();
    int N,E,E1;
    float a,b;
    srand(time(NULL));
    ofstream en,cc,his,en_1;
    ifstream en_ini;
    char buffer[32];
    snprintf(buffer,sizeof(char)*32,"entropy_2para_%i.dat",L);
    en.open(buffer);
    snprintf(buffer,sizeof(char)*32,"Histogram_%i.dat",L);
    his.open(buffer);
    snprintf(buffer,sizeof(char)*32,"entropy_1para_%i.dat",L);
    en_1.open(buffer);
    en_ini.open("entropy_1para_ini.dat");
    N=EN*1000000;
    for(int i=0; i<3*EN+1; i+=2)
    {
        en_ini>>a>>b;
        entropy[int(a)]=b;
        //cout<<a<<"\t"<<b<<"\n";
    }
    for(int j=0; j<10; j++)
    {
        for(int i=0; i<3*EN+1; i++)
            for(int j=0; j<EN+1; j++)
            {
                Histogram[i]=0;
                Histogram_mp[i][j]=0;
                Mag[i]=0;
            }
        for(int k=0; k<10; k++)
        {

            initial(k);
            E=energy_n();
            m1=magnetization();
            cout<<E<<"\t"<<m1<<"\n";
            for(int i=0; i<N; i++)
            {
                E1=flip(E);
                E=E1;
            }
        }
        for(int i=0; i<3*EN+1; i++)
            for(int j=0; j<EN+1; j++)
            {
                Histogram[i]=Histogram[i]+Histogram_mp[i][j];
                //cout<<Histogram_mp[i][j]<<"\t"<<j<<"\t"<<i<<"\n";
            }
//	    his();
        double corr_his=0;
        int count=0;
        for(int i=0; i<3*EN+1; i++)
        {
            if(Histogram[i]!=0)
            {
                corr_his=corr_his+((Histogram[i]-sum)/sum)*((Histogram[i]-sum)/sum);
                count++;
            }
        }
	
        cout<<float(corr_his/count)<<"\n";
        if(float(corr_his/count) < .005)
        {
            break;
        }

        avg=average();
        cout<<sum<<"\n";
        for(int i=0; i<3*EN+1; i++)
        {
            if(Histogram[i]!=0)
                entropy[i]=log(Histogram[i]/sum)/EN+entropy[i];
        }
    }
    for(int i=0; i<3*EN+1; i++)
        for(int j=0; j<EN+1; j++)
        {
            if (Histogram_mp[i][j]!=0)
                entropy_mp[i][j]=entropy[i]+log(Histogram_mp[i][j]/Histogram[i])/EN;
        }
    for(int i=0; i<3*EN+1; i++)
    {
        if(Histogram[i]!=0)
        {
            his<<i<<"\t"<<(Histogram[i]-sum)/sum<<"\n";
        }
    }
    for(int i=0; i<3*EN+1; i++)
    {
        for(int n=0; n<EN+1; n+=1)
            if(Histogram_mp[i][n]!=0)
            {
                en<<i<<"\t"<<float(2*n-EN)<<"\t"<<EN*(entropy_mp[i][n]-entropy_mp[3*int(EN)][int(EN)])<<"\t"<<exp(EN*(entropy_mp[i][n]-entropy_mp[3*int(EN)][int(EN)]))<<"\n";
            }
        if(Histogram[i]!=0)
            en_1<<i<<"\t"<<EN*(entropy[i]-entropy[3*int(EN)])<<"\t"<<exp(EN*(entropy[i]-entropy[3*int(EN)]))<<"\n";
    }

//	for(float t=1.0;t<4.0;t+=0.03)
//	cout<<t<<"\t"<<mag(t)<<"\n";
    return 0;
}
