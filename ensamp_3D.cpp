using namespace std;
#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<fstream>
#include<ctime>
//Defining the variables requred
int L,m,n,N_plus,N_minus,m1;
int ISING[100][100][100];
double entropy[18000];
double entropy_mp[18000][6000];
double Histogram_mp[18000][6000];
double Histogram[18000];
double Mag[18000];
float avg,sum;
float EN;
//int KB=1;
//The different initial conditions 15 of them indexed by variable k from 0-14
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
                    ISING[i][j][k]=pow(-1,i+j);
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
        for(int i=0; i<L/2; i++)
            for(int j=0; j<L/2; j++)
                for(int k=0; k<L/2; k++)
                {
                    ISING[i][j][k]=pow(-1,i);
                }

        //  ISING[1][1][1]=-1*ISING[1][1][1];
        //  ISING[1][2][1]=-1*ISING[1][2][1];

    }
    if(k==9)
    {
        for(int i=0; i<L; i++)
            for(int j=0; j<L; j++)
                for(int k=0; k<L; k++)
                    ISING[i][j][k]=-1;
        int count=0;
        for(int i=0; i<L; i++)
            for(int j=0; j<L; j++)
                for(int k=0; k<L; k++)
                {
                    if (count==108)
                        break;
                    ISING[i][j][k]=1;
                    count=count+1;
                }
    }
    if(k==10)
    {
        for(int i=0; i<L; i++)
            for(int j=0; j<L; j++)
                for(int k=0; k<L; k++)
                    ISING[i][j][k]=1;
        for(int i=0; i<L; i++)
            ISING[1][1][i]=-1;
        for(int i=0; i<L; i++)
            ISING[1][3][i]=-1;
        for(int i=0; i<L; i++)
            ISING[1][5][i]=-1;
        for(int i=0; i<L; i++)
            ISING[3][1][i]=-1;
//	for(int i=0;i<L;i++)
//		ISING[3][3][i]=-1;
        for(int i=0; i<L; i++)
            ISING[3][5][i]=-1;
        for(int i=0; i<L; i++)
            ISING[5][1][i]=-1;
        for(int i=0; i<L; i++)
            ISING[5][3][i]=-1;


    }
    if(k==11)
    {
        for(int i=0; i<L; i++)
            for(int j=0; j<L; j++)
                for(int k=0; k<L; k++)
                    ISING[i][j][k]=-1;
        for(int i=0; i<L; i++)
            ISING[1][1][i]=1;
        for(int i=0; i<L; i++)
            ISING[1][3][i]=1;
        for(int i=0; i<L; i++)
            ISING[1][5][i]=1;
        for(int i=0; i<L; i++)
            ISING[3][1][i]=1;
//	for(int i=0;i<L;i++)
//		ISING[3][3][i]=1;
        for(int i=0; i<L; i++)
            ISING[3][5][i]=1;
        for(int i=0; i<L; i++)
            ISING[5][1][i]=1;
        for(int i=0; i<L; i++)
            ISING[5][3][i]=1;


    }
    if(k==12)
    {
        for(int i=0; i<L; i++)
            for(int j=0; j<L; j++)
                for(int k=0; k<L; k++)
                    ISING[i][j][k]=pow(-1,i+j+k+1);
        for(int i=0; i<L; i++)
            ISING[1][1][i]=pow(-1,i);
        for(int i=0; i<L; i++)
            ISING[1][3][i]=pow(-1,i);
        for(int i=0; i<L; i++)
            ISING[1][5][i]=pow(-1,i);
        for(int i=0; i<L; i++)
            ISING[3][1][i]=pow(-1,i);
        for(int i=0; i<L; i++)
            ISING[3][3][i]=pow(-1,i);
        for(int i=0; i<L; i++)
            ISING[3][5][i]=pow(-1,i);

    }
    if(k==13)
    {
        for(int i=0; i<L; i++)
            for(int j=0; j<L; j++)
                for(int k=0; k<L; k++)
                    ISING[i][j][k]=-1;
        int count=0;
        for(int i=0; i<L; i++)
            for(int j=0; j<L; j++)
                for(int k=0; k<L; k++)
                {
                    if (count==58)
                        break;
                    ISING[i][j][k]=1;
                    count=count+1;
                }
    }
    if(k==14)
    {
        for(int i=0; i<L; i++)
            for(int j=0; j<L; j++)
                for(int k=0; k<L; k++)
                    ISING[i][j][k]=-1;
        int count=0;
        for(int i=0; i<L; i++)
            for(int j=0; j<L; j++)
                for(int k=0; k<L; k++)
                {
                    if (count==158)
                        break;
                    ISING[i][j][k]=1;
                    count=count+1;
                }
    }


}
//Function to calculate the magnetisation which is just a simple sum of the ISING matrx
int magnetization()
{
    int sum=0;
    for(int i=0; i<L; i++)
        for(int j=0; j<L; j++)
            for(int k=0; k<L; k++)
                sum=sum+ISING[i][j][k];
    return sum;
}
//Functon to calculate the number of like bonds 'n' (energy is a function of n)
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
                //The following impose Periodic boundary conditions
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
                //Hamiltonian is given by H=-1/2 sum_{nearest neighbours} S_i x S_j
                Energy=Energy+-1*(ISING[i][j][k]*ISING[i][Nj_plus][k]+ISING[i][j][k]*ISING[i][Nj_minus][k]+ISING[i][j][k]*ISING[Ni_plus][j][k]+ISING[i][j][k]*ISING[Ni_minus][j][k]+ISING[i][j][k]*ISING[i][j][Nk_plus]+ISING[i][j][k]*ISING[i][j][Nk_minus]);
            }
    //H and # of like bonds(n) are related by (for 3d)  E=6*L^3-n

    return (3*EN-Energy/2)/2;
}

int flip(int E)
{
    int Ni_plus=0,Ni_minus=0,Nj_plus=0,Nj_minus=0,Nk_plus=0,Nk_minus=0;
    int X,Y,Z,n_1,n_2,delta_n,m2,i,j,delta_m;
    double alpha;
    float p;
    X=rand()%L;
    Y=rand()%L;
    Z=rand()%L;
    n_1=E;
    // the following are again used for PBC implemetation.
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
    //Change in m if you flip a spin
    if(ISING[X][Y][Z]==-1)
    {
        delta_m=2;
    }
    else
    {
        delta_m=-2;
    }
    //i and j are the number of up spins before and after  flipping
    i=(m1+EN)/2;
    m2=m1+delta_m;
    j=(m2+EN)/2;
    //The change in number of like bonds when you flip the spin. n_1 and n_2 are like bonds before and after flipping.
    n_2=n_1+-1*(ISING[X][Y][Z]*ISING[X][Nj_plus][Z]+ISING[X][Y][Z]*ISING[X][Nj_minus][Z]+ISING[X][Y][Z]*ISING[Ni_plus][Y][Z]+ISING[X][Y][Z]*ISING[Ni_minus][Y][Z]+ISING[X][Y][Z]*ISING[X][Y][Nk_plus]+ISING[X][Y][Z]*ISING[X][Y][Nk_minus]);//delta_n;
    alpha=exp(EN*(entropy_mp[n_1][i]-entropy_mp[n_2][j]));
    //alpha is the acceptance probability of the flip.
    //cout<<n_1<<"\t"<<m1<<"\t"<<n_2<<"\t"<<m2<<"\t"<<alpha<<"\t"<<Histogram_mp[n_1][i]<<"\n";
    if (alpha<1)
    {
        p=float(rand())/RAND_MAX;
        if(p>alpha)
        {
            //No flip
            Histogram_mp[n_1][i]++;//=Histogram_mp[n_1][i]+1;
            return n_1;
        }
        else
        {
            //flip
            Histogram_mp[n_2][j]++;//=Histogram_mp[n_2][j]+1;
            m1=m2;
            ISING[X][Y][Z]=-1*ISING[X][Y][Z];
            return n_2;
        }

    }
    else
    {
        //flip
        ISING[X][Y][Z]=-1*ISING[X][Y][Z];
        Histogram_mp[n_2][j]++;//=Histogram_mp[n_2][j]+1;
        m1=m2;
        return n_2;
    }
}
// the average of Histogram , here the average is calulated as the average over all config classes visited , as opposed to all config classes possib
//e.
//Note : Histogram_mp[i][j]!=0 means config (i,j) has been visited.
float average()
{
    sum=0;
    float avg=0,count=0,count_2=0;
    for(int i=0; i<3*EN+1; i++)
        for(int j=0; j<EN+1; j++)
            if(Histogram_mp[i][j]!=0)
            {
                avg=avg+Histogram_mp[i][j];
                count++;
            }
    cout<<avg<<"\t"<<count<<"\n";
    return avg/(count);
}
//The main program
int main(int argc,char* argv[])
{
    // Various varibles required , L (linear size of the lattice is given as an command line argument to the program
    L=atof(argv[1]);
    EN=L*L*L;
    int E,E1,count=0;
    long double N;
    float a,b,c,flag;
    long double d;
    srand(time(NULL));
    ofstream en,cc,his,en_1,diff;
    ifstream en_ini,infile;
    char buffer[32];
    //These two file streams are for writing the 1_d entropy and histogram , ie as a funciton of just n
    snprintf(buffer,sizeof(char)*32,"entropy_1para_%i.dat",L);
    en_1.open(buffer);
    en_ini.open("en.dat");
    // The number of mc moves
    N=EN*1000000;
    cout<<"N="<<N<<"\n";
    //The initial configuration is opened from a file
    infile.open("entropy_8_ini.dat");
    for(int i=0; !infile.eof(); i++)
    {
        infile>>a>>b>>c;
        N_plus=float(a+EN)/2;
        //Initial config is written to the entropy array
        entropy_mp[int(b)][int(N_plus)]=float(c)/EN;
    }
    infile.close();

    //main loop, the no-of time the entropy is updated.
    for(int j=0; j<15; j++)
    {
        //setting all data to zero
        for(int i=0; i<3*EN+1; i++)
            for(int j=0; j<EN+1; j++)
            {
                Histogram[i]=0;
                Histogram_mp[i][j]=0;
            }
        //loop to run the simulations with all inital condition.
        cout<<"n\tm"<<"\n";
        for(int k=3; k<15; k++)
        {
		if(k==6 or k==7)
			continue;
            //initializing with the k'th initial-config.
            initial(k);
            E=energy_n();
            //energy of the ising model-given as no of unlike bonds.
            m1=magnetization();
            cout<<E<<"\t"<<m1<<"\n";
            //this loop has N MC moves.
            clock_t begin =clock();
            for(long double i=0; i<N; i++)
            {
                E1=flip(E);
                E=E1;
            }
            clock_t end = clock();
            double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
            //	cout<<elapsed_secs<<"\n";


        }
        //After going through all the initial conditions , the averge of the Histogram is calulated
        avg=average();

        cout<<"average H="<<avg<<"\n";
        //The 2 paramenter entropy array is updated and sampling is again done with the modified entropy
	snprintf(buffer,sizeof(char)*32,"diff_mp_%i.dat",j);
        diff.open(buffer);
        for(int i=0; i<3*EN+1; i++)
            for(int j=0; j<EN+1; j++)
            {
                if (Histogram_mp[i][j]!=0)
		    {
                    	entropy_mp[i][j]=entropy_mp[i][j]+log(Histogram_mp[i][j]/avg)/EN;
		    	diff<<i<<"\t"<<j<<"\t"<<log(Histogram_mp[i][j]/avg)<<"\n";
		    }
		else 
			diff<<i<<"\t"<<j<<"\t"<<" "<<"\n";
            }
	diff.close();
        snprintf(buffer,sizeof(char)*32,"Histogram_mp_%i.dat",j);
        his.open(buffer);
        count=0;
        flag=0;
        //The Histogram is written after each sampling
        for(int i=0; i<3*EN+1; i++)
            for(int j=0; j<EN+1; j++)
            {
                if(Histogram_mp[i][j]!=0)
                {
                    his<<i<<"\t"<<2*j-EN<<"\t"<<Histogram_mp[i][j]<<"\t"<<(Histogram_mp[i][j]-avg)/avg<<"\n";
                    //The flatness of the histogram is calculated
                    count=count+((Histogram_mp[i][j]-avg)/avg)*((Histogram_mp[i][j]-avg)/avg);
                    flag=flag+1;
                }
            }
        cout<<"Flatness="<<float(count)/flag<<"\t"<<j<<"th iteration\n";
        his.close();
//	entropy is written after every sampling to a file called 'entropy_mp_i.dat'
        snprintf(buffer,sizeof(char)*32,"entropy_mp_%i.dat",j);
        en.open(buffer);
        for(int i=0; i<3*EN+1; i++)
        {
            for(int n=0; n<EN+1; n+=1)
                if(Histogram_mp[i][n]!=0)
                {
                    en<<i<<"\t"<<float(2*n-EN)<<"\t"<<EN*(entropy_mp[i][n])-EN*(entropy_mp[3*int(EN)][int(EN)])<<"\t"<<exp(EN*(entropy_mp[i][n])-EN*(entropy_mp[3*int(EN)][int(EN)]))<<"\n";
                }
        }
        en.close();
	long double dos;
	for(int i=0; i<3*EN+1; i++)
	    {
		dos=0;
            for(int j=0; j<EN+1; j++)
            {
                if (Histogram_mp[i][j]!=0)
			dos=dos+exp(EN*(entropy_mp[i][j]-entropy_mp[3*int(EN)][int(EN)]));
		}
		if(dos!=0)
			en_1<<i<<"\t"<<log(dos)<<"\n";
	}
	en_1.close();

        //IF the histogram is flat enough break
        if(float(count)/flag<.005)
            break;


    }

    return 0;
}
