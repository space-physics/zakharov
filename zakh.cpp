//1D Langmuir turbulence simulation according to Guio and Forme 2006
// Last update: 9/5/2013
// For Kappa on: 01/05/2016

//Michael Hirsch Oct 2013 -- updated vbeam,tetabeam to use C++ vector format

#include <boost/filesystem.hpp>
#include "boost/program_options.hpp" 
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <ctime>
#include <fstream>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const double pi=3.141592653589793;
const double me=9.109e-31;
const double electroncharge=1.602e-19;
const double mi=16*1.66e-27; // atomic oxygen
const double Kb=1.38e-23;    // Boltzmann cte
const double eV=1.6e-19;
const double epsilon0=8.854e-12;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int Z=1;
double Te=3000.0;
double Ti=1000.0;
double nuic=1.0;   // ion collision freq
double nuec=100.0; // electron collision freq
double n0=5.0e11;  // background density
std::vector<double> nbeam {(60.0e-7)*n0};
int Nnbeam=nbeam.size();
std::vector<double> vbeam_ev {500.0};
int Nvbeam=vbeam_ev.size();
double power_n_cte=7.7735e6/sqrt(53.5);
double power_E_cte=0.5033*0.5033/sqrt(3.0);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////    Kappa parameters
// look at Broughton et al., Modeling of MF wave mode conversion
double se_percent=0.001;
double kappa=1.584;
double T_se=18.2*eV;
double theta_se=sqrt((kappa-1.5)/kappa*2*T_se/me);
double se_cte= 0.7397/pow(theta_se,3);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////    Simulation parameters

const double endTime=100.0e-3;    // simulation ends (seconds)
const double Tstep=0.5e-7;     // simulation time steps
double TT=endTime/Tstep;              //floor(endTime/Tstep)+2;
const int res=20;
const int TT_res=floor(endTime/Tstep/res);			//floor(endTime/Tstep/res);
double L=70.0;           // simulation box length (meter)j
const int N=2046;          // number of samples in L; should be devidable by 6
double Xstep=L/N;
const int QW=1;     //number of realizations
unsigned int SEED=600;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double eta=(Te+3*Ti)/Te;
double ve=sqrt(Kb*Te/me);
double Cs=sqrt(eta*me/mi)*ve;
double omegae=sqrt(n0*pow(electroncharge,2)/me/epsilon0);
double lambdaD=ve/omegae;

std::string odir;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////   main

int max1 (int, int);
int min1 (int, int);
int sign(double);
void Xsection(double&,double&,double k);

int main(int argc, char ** argv)
{


// argparse
    namespace po = boost::program_options; 
    po::options_description desc("Options"); 
    desc.add_options() 
      ("help", "Print help messages")
      ("outdir,o",po::value<std::string>(&odir)->required(), "Output directory");
 
    po::variables_map vm; 
    try 
    { 
      po::store(po::parse_command_line(argc, argv, desc),  
                vm); // can throw 
 
      /** --help option 
       */ 
      if ( vm.count("help")  ) 
      { 
        std::cout << "Basic Command Line Parameter App" << std::endl 
                  << desc << std::endl; 
        return EXIT_SUCCESS; 
      } 
 
      po::notify(vm); // throws on error, so do after help in case 
                      // there are any problems 
    } 
    catch(po::error& e) 
    { 
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
      std::cerr << desc << std::endl; 
      return EXIT_FAILURE; 
    } 

std::string outDir = odir + boost::filesystem::path::preferred_separator;

// create output directory if it doesn't exist (previously, program segfaulted...)
    boost::filesystem::path odir(outDir);

    if(!(boost::filesystem::exists(odir))){
        if (boost::filesystem::create_directory(odir))
            std::cout << "created output directory " << outDir << std::endl;
    }
    else{
        std::cout << "using output directory " << outDir << std::endl;
    }
//-------------------------------------------------------------------------------------
    printf("Nnbeam=%i \n",Nnbeam);
	printf("Nvbeam=%i \n",Nvbeam);
	printf("TT=%0.1f time steps \n",TT);

	time_t now;
	struct tm *current;
	now = time(0);
	current = localtime(&now);
	int StartTime[3]={current->tm_hour, current->tm_min, current->tm_sec};

/////////////////////////////////////////////////////////////////////////////////////////////////   initialization


std::vector<double> vbeam (Nvbeam);
std::vector<double> tetabeam (Nvbeam);

for (int i=0; i<Nvbeam; i++){
 vbeam[i]=sqrt(eV*vbeam_ev.at(i)*2/me);
 tetabeam[i]=0.3*vbeam.at(i);
			printf("vbeam[%i]=%0.6f \n",i,vbeam.at(i));
            printf("tetabeam[%i]=%0.6f \n",i,tetabeam.at(i));
}


for (int beami=0;beami<Nnbeam;beami++){
for (int beamj=0;beamj<Nvbeam;beamj++){


	int p[N];
	double k[N];
	double Xsection_ion=0.0;
	double Xsection_pl=0.0;
	double E_thermal_k_squared;
	double n_thermal_k_squared;
	double Source_factor_E[N];
	double Source_factor_n[N];
	double omegaL[N];
	double gamas;
	double nui[N];
	double gamal1;
	double gamal2;
	double gamal3;
	double gamal;
	double nue[N];
	double parameters[32];
	double output1[N][12];

	parameters[0]=pi;
	parameters[1]=me;
	parameters[2]=electroncharge;
	parameters[3]=mi;
	parameters[4]=Kb;
	parameters[5]=eV;
	parameters[6]=epsilon0;
	parameters[7]=Z;
	parameters[8]=Te;
	parameters[9]=Ti;
	parameters[10]=nuic;
	parameters[11]=nuec;
	parameters[12]=n0;
	parameters[13]=nbeam.at(beami);
	parameters[14]=vbeam_ev.at(beamj);
	parameters[15]=vbeam.at(beamj);
	parameters[16]=tetabeam.at(beamj);
	parameters[17]=endTime;
	parameters[18]=Tstep;
	parameters[19]=TT;
	parameters[20]=res;
	parameters[21]=TT_res;
	parameters[22]=L;
	parameters[23]=N;
	parameters[24]=Xstep;
	parameters[25]=QW;
	parameters[26]=SEED;
	parameters[27]=eta;
	parameters[28]=ve;
	parameters[29]=Cs;
	parameters[30]=omegae;
	parameters[31]=lambdaD;
	//have to include other parameters regarding the Kappa distribution
	
	std::string fn = outDir + "parameters_n" + std::to_string(beami) + "_v" + std::to_string(beamj);
	const char* fnp = fn.c_str();
	FILE* parameters_out;
	parameters_out = fopen(fnp, "wb");
	int bout=fwrite(parameters, 1, sizeof(parameters), parameters_out);
	fclose(parameters_out);
	std::cout << "Wrote " << bout << " bytes to " << fn << std::endl;

	for (int ii=0;ii<N;ii++){
		p[ii]=ii-N/2;
		if (ii==N/2){
			k[ii]=0;
			Xsection_ion=0;
			Xsection_pl=0;
			n_thermal_k_squared=0.0;
			E_thermal_k_squared=0.0;
		}
		else{
		k[ii]=2*pi*p[ii]/N/Xstep;
		Xsection(Xsection_ion,Xsection_pl,k[ii]);
		Xsection_ion=Xsection_ion/N/N;
		Xsection_pl=Xsection_pl/N/N;
		n_thermal_k_squared=Xsection_ion*n0;
		E_thermal_k_squared=Xsection_pl *n0*pow(electroncharge/epsilon0/k[ii],2);
		}
		omegaL[ii]=sqrt(pow(omegae,2)+3*pow(k[ii]*ve,2));
		gamas= (-1)*sqrt(pi/8)*(sqrt(me/mi)+pow(Te/Ti,2)/sqrt(Te/Ti)*exp((-1)*(Te/2.0/Ti)-1.5))*std::fabs(k[ii])*Cs;
		//gamas= (-1)*sqrt(pi/2)*(sqrt(me/mi)+4*pow(Te/2/Ti,2)/sqrt(Te/2/Ti)*exp((-1)*(Te*4/Ti)))*abs(k[ii])*Cs*10;   //based on Robinson 2002
		//gamas= (-1)*sqrt(pi/8)*pow(1/(1+k[ii]*k[ii]*lambdaD*lambdaD)+3*Ti/Te,2)/sqrt(1/(1+k[ii]*k[ii]*lambdaD*lambdaD)+3*Ti/Te)*(sqrt(me/mi)+pow(Te/Ti,2)/sqrt(Te/Ti)*exp((-1)*(Te/2.0/Ti)/(1+k[ii]*k[ii]*lambdaD*lambdaD)-1.5))*abs(k[ii])*Cs;   //Based on some Chinese paper!!
		nui[ii]=(nuic/2-gamas);
		if (ii==N/2){
			gamal=0.0;           // this one is Nan due to division by zero
			gamal1=0.0;
			nue[ii]=nuec/2-gamal1;
			Source_factor_n[ii]=0.0;
			Source_factor_E[ii]=0.0;
		}
		else{
			gamal1=(-1)*sqrt(pi/8)*pow(omegae/k[ii]/ve,2)*sign(k[ii])*pow(omegaL[ii],2)/(k[ii]*ve)*exp((-1)*pow(omegaL[ii]/k[ii]/ve,2)/2);    //Landau damping due to the thermal electrons
			gamal2=(-1)*sqrt(pi/8)*pow(omegae/k[ii]/tetabeam.at(beamj),2)*sign(k[ii])*nbeam.at(beami)/n0*omegaL[ii]*(omegaL[ii]-k[ii]*vbeam.at(beamj))/(k[ii]*tetabeam.at(beamj))*exp((-1)*pow((omegaL[ii]-k[ii]*vbeam.at(beamj))/k[ii]/tetabeam.at(beamj),2)/2);  //Landau damping due to the beam
			gamal3=(-1)*sqrt(pi)*pow(omegae*omegaL[ii],2)/pow(k[ii],3)*sign(k[ii])*se_cte*pow(1+pow(omegaL[ii],2)/kappa/pow(k[ii]*theta_se,2),(-1)*(kappa+1));
			gamal=gamal1*(1-se_percent)+gamal2+se_percent*gamal3;    // here decide to include the beam and Kappa distribution
			nue[ii]=nuec/2-gamal1;
			Source_factor_n[ii]=2*nui[ii]*sqrt(4*nui[ii]*k[ii]*k[ii]/(4*nui[ii]*nui[ii]+k[ii]*k[ii])*n_thermal_k_squared*power_n_cte);
			Source_factor_E[ii]=sqrt(2*nue[ii]*E_thermal_k_squared*power_E_cte);                                                         // source factor is the factor by which we balance the thermal source intensity
			nue[ii]=nuec/2-gamal;
		}


		output1[ii][0]=p[ii];
		output1[ii][1]=k[ii];
		output1[ii][2]=Xsection_ion;
		output1[ii][3]=Xsection_pl;
		output1[ii][4]=E_thermal_k_squared;
		output1[ii][5]=n_thermal_k_squared;
		output1[ii][6]=omegaL[ii];
		output1[ii][7]=gamas;
		output1[ii][8]=nui[ii];
		output1[ii][9]=nue[ii];
		output1[ii][10]=Source_factor_E[ii];
		output1[ii][11]=Source_factor_n[ii];
	}

	fn = outDir + "output1_n" + std::to_string(beami) + "_v" + std::to_string(beamj);
	const char* fno = fn.c_str();
	FILE* output1_out;
	output1_out = fopen(fno, "wb");
	bout = fwrite(output1, 1, sizeof(output1), output1_out);
	fclose(output1_out);
	std::cout << "Wrote " << bout << " bytes to " << fn << std::endl;


	static double EE [3][N][2];
	static double nn [3][N][2];
	static double vv [3][N][2];
	int LL,UU,pp;
	double CC[2];
	double SSE [N][2];
	double SSn [2];
	double cte1;
	double cte2=omegae/2.0/n0/1;
	double k1[N][2],k2[N][2],k3[N][2],k4[N][2];
	double kn1[2], kn2[2], kn3[2], kn4[2], kv1[2], kv2[2], kv3[2], kv4[2];
	static double total_EE[20000*N*2];
	static double total_nn[20000*N*2];


for (int realization=0;realization<QW;realization++){

	unsigned long int aabb=SEED+realization;
	//unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator (aabb);
	std::normal_distribution<long double> distribution (0.0,1.0);

	fn = outDir + "EE" + std::to_string(SEED+realization) + "_n" + std::to_string(beami) + "_v" + std::to_string(beamj);
    const char* nameE = fn.c_str();
	FILE* EE_out;
	EE_out = fopen(nameE, "wb");
	
	fn = outDir + "nn" + std::to_string(SEED+realization) + "_n" + std::to_string(beami) + "_v" + std::to_string(beamj);
    const char* namen = fn.c_str();
	FILE* nn_out;
	nn_out = fopen(namen, "wb");

	/////////////////////////////////////////////////////////////////////////////////////////////////   main loops

	for (int iij1=0;iij1<3;iij1++){
	    for (int iij2=0;iij2<N;iij2++){
	        EE [iij1][iij2][0]=sqrt(output1[iij2][4]/2.0)*distribution(generator);
	        EE [iij1][iij2][1]=sqrt(output1[iij2][4]/2.0)*distribution(generator);
	        nn [iij1][iij2][0]=sqrt(output1[iij2][5]/2.0)*distribution(generator);
	        nn [iij1][iij2][1]=sqrt(output1[iij2][5]/2.0)*distribution(generator);
	        vv [iij1][iij2][0]=0;
	        vv [iij1][iij2][1]=0;

        }

        for (int iij2=1;iij2<N/2;iij2++){
            nn [iij1][N-iij2][0]=nn [iij1][iij2][0];
            nn [iij1][N-iij2][1]=-nn [iij1][iij2][1];
        }
        nn [iij1][N/2][0]=0;
        nn [iij1][N/2][1]=0;


    }

	int counter1=0;
	for (int tt1=1;tt1<=TT;tt1++){

//		int c0=(tt1-1) % 3;
		int c1=(tt1) % 3;
		int c2=(tt1+1) % 3;
//		long double omega_off=omegae+2*pi*300000;

					//update display every 50th iteration
		if (tt1 % 50 == 0){
		printf("Realization: %i ,  %0.2f%% complete_n%d_v%d \n",realization+1, tt1*100.0/TT,  beami, beamj);
                }

			for (pp=0;pp<N;pp++){


				LL= max1(p[pp]-N/3,-N/3);
				UU= min1(N/3,p[pp]+N/3);
				CC[0]=0.0;
				CC[1]=0.0;

				for (int q=LL;q<=UU;q++){
					CC[0]=CC[0]+EE[c1][(q+N/2)][0]*nn[c1][(p[pp]-q+N/2)][0]-EE[c1][(q+N/2)][1]*nn[c1][(p[pp]-q+N/2)][1];
					CC[1]=CC[1]+EE[c1][(q+N/2)][0]*nn[c1][(p[pp]-q+N/2)][1]+EE[c1][(q+N/2)][1]*nn[c1][(p[pp]-q+N/2)][0];
				}


				SSE[pp][0]=distribution(generator)*Source_factor_E[pp]/sqrt(Tstep);
				SSE[pp][1]=distribution(generator)*Source_factor_E[pp]/sqrt(Tstep);

				cte1=1.5*omegae*(lambdaD*k[pp])*(lambdaD*k[pp]);
				//cte1=1.5*Kb*Te/me/omega_off*k[pp]*k[pp]-(pow(omega_off,2)-pow(omegae,2))/2.0/omega_off;
				k1[pp][0]=Tstep*(cte1*EE[c1][pp][1]-nue[pp]*EE[c1][pp][0]+cte2*CC[1]);
				k1[pp][1]=Tstep*((-1)*cte1*EE[c1][pp][0]-nue[pp]*EE[c1][pp][1]-cte2*CC[0]);
			}

			for (pp=0;pp<N;pp++){


				LL= max1(p[pp]-N/3,-N/3);
				UU= min1(N/3,p[pp]+N/3);
				CC[0]=0.0;
				CC[1]=0.0;

				for (int q=LL;q<=UU;q++){
					CC[0]=CC[0]+(EE[c1][(q+N/2)][0]+k1[(q+N/2)][0]/2)*nn[c1][(p[pp]-q+N/2)][0]-(EE[c1][(q+N/2)][1]+k1[(q+N/2)][1]/2)*nn[c1][(p[pp]-q+N/2)][1];
					CC[1]=CC[1]+(EE[c1][(q+N/2)][0]+k1[(q+N/2)][0]/2)*nn[c1][(p[pp]-q+N/2)][1]+(EE[c1][(q+N/2)][1]+k1[(q+N/2)][1]/2)*nn[c1][(p[pp]-q+N/2)][0];
				}



				cte1=1.5*omegae*(lambdaD*k[pp])*(lambdaD*k[pp]);
				//cte1=1.5*Kb*Te/me/omega_off*k[pp]*k[pp]-(pow(omega_off,2)-pow(omegae,2))/2.0/omega_off;
				k2[pp][0]=Tstep*(cte1*(EE[c1][pp][1]+k1[pp][1]/2.0-SSE[pp][0]/2.0*Tstep)-nue[pp]*(EE[c1][pp][0]+k1[pp][0]/2.0+SSE[pp][1]/2.0*Tstep)+cte2*CC[1]);
                k2[pp][1]=Tstep*((-1)*cte1*(EE[c1][pp][0]+k1[pp][0]/2.0+SSE[pp][1]/2.0*Tstep)-nue[pp]*(EE[c1][pp][1]+k1[pp][1]/2.0-SSE[pp][0]/2.0*Tstep)-cte2*CC[0]);
			}


			for (pp=0;pp<N;pp++){


				LL= max1(p[pp]-N/3,-N/3);
				UU= min1(N/3,p[pp]+N/3);
				CC[0]=0.0;
				CC[1]=0.0;

				for (int q=LL;q<=UU;q++){
					CC[0]=CC[0]+(EE[c1][(q+N/2)][0]+k2[(q+N/2)][0]/2)*nn[c1][(p[pp]-q+N/2)][0]-(EE[c1][(q+N/2)][1]+k2[(q+N/2)][1]/2)*nn[c1][(p[pp]-q+N/2)][1];
					CC[1]=CC[1]+(EE[c1][(q+N/2)][0]+k2[(q+N/2)][0]/2)*nn[c1][(p[pp]-q+N/2)][1]+(EE[c1][(q+N/2)][1]+k2[(q+N/2)][1]/2)*nn[c1][(p[pp]-q+N/2)][0];
				}



				cte1=1.5*omegae*(lambdaD*k[pp])*(lambdaD*k[pp]);
				//cte1=1.5*Kb*Te/me/omega_off*k[pp]*k[pp]-(pow(omega_off,2)-pow(omegae,2))/2.0/omega_off;
				k3[pp][0]=Tstep*(cte1*(EE[c1][pp][1]+k2[pp][1]/2.0-SSE[pp][0]/2.0*Tstep)-nue[pp]*(EE[c1][pp][0]+k2[pp][0]/2.0+SSE[pp][1]/2.0*Tstep)+cte2*CC[1]);
                k3[pp][1]=Tstep*((-1)*cte1*(EE[c1][pp][0]+k2[pp][0]/2.0+SSE[pp][1]/2.0*Tstep)-nue[pp]*(EE[c1][pp][1]+k2[pp][1]/2.0-SSE[pp][0]/2.0*Tstep)-cte2*CC[0]);
			}

			for (pp=0;pp<N;pp++){


				LL= max1(p[pp]-N/3,-N/3);
				UU= min1(N/3,p[pp]+N/3);
				CC[0]=0.0;
				CC[1]=0.0;

				for (int q=LL;q<=UU;q++){
					CC[0]=CC[0]+(EE[c1][(q+N/2)][0]+k3[(q+N/2)][0])*nn[c1][(p[pp]-q+N/2)][0]-(EE[c1][(q+N/2)][1]+k3[(q+N/2)][1])*nn[c1][(p[pp]-q+N/2)][1];
					CC[1]=CC[1]+(EE[c1][(q+N/2)][0]+k3[(q+N/2)][0])*nn[c1][(p[pp]-q+N/2)][1]+(EE[c1][(q+N/2)][1]+k3[(q+N/2)][1])*nn[c1][(p[pp]-q+N/2)][0];
				}


				cte1=1.5*omegae*(lambdaD*k[pp])*(lambdaD*k[pp]);
				//cte1=1.5*Kb*Te/me/omega_off*k[pp]*k[pp]-(pow(omega_off,2)-pow(omegae,2))/2.0/omega_off;
				k4[pp][0]=Tstep*(cte1*(EE[c1][pp][1]+k3[pp][1]-SSE[pp][0]*Tstep)-nue[pp]*(EE[c1][pp][0]+k3[pp][0]+SSE[pp][1]*Tstep)+cte2*CC[1]);
                k4[pp][1]=Tstep*((-1)*cte1*(EE[c1][pp][0]+k3[pp][0]+SSE[pp][1]*Tstep)-nue[pp]*(EE[c1][pp][1]+k3[pp][1]-SSE[pp][0]*Tstep)-cte2*CC[0]);

                EE[c2][pp][0]=EE[c1][pp][0]+(k1[pp][0]+2.0*k2[pp][0]+2.0*k3[pp][0]+k4[pp][0])/6.0+SSE[pp][1]*Tstep;
                EE[c2][pp][1]=EE[c1][pp][1]+(k1[pp][1]+2.0*k2[pp][1]+2.0*k3[pp][1]+k4[pp][1])/6.0-SSE[pp][0]*Tstep;
				EE[c2][N/2][0]=0.0;
				EE[c2][N/2][1]=0.0;


			}

            for (pp=0;pp<=N/2;pp++){


				SSn[0]=distribution(generator)*Source_factor_n[pp]/sqrt(Tstep);
				SSn[1]=distribution(generator)*Source_factor_n[pp]/sqrt(Tstep);


				LL= max1(p[pp]-N/3,-N/3);
				UU= min1(N/3,p[pp]+N/3);
				CC[0]=0.0;
				CC[1]=0.0;
				for (int q=LL;q<=UU;q++){
					CC[0]=CC[0]+EE[c1][(q+N/2)][0]*EE[c1][(q-p[pp]+N/2)][0]+EE[c1][(q+N/2)][1]*EE[c1][(q-p[pp]+N/2)][1];
					CC[1]=CC[1]+EE[c1][(q+N/2)][1]*EE[c1][(q-p[pp]+N/2)][0]-EE[c1][(q+N/2)][0]*EE[c1][(q-p[pp]+N/2)][1];
				}



				kn1[0]=Tstep*(vv[c1][pp][0]);
				kn1[1]=Tstep*(vv[c1][pp][1]);
				kv1[0]=Tstep*((-2.0)*nui[pp]*vv[c1][pp][0]-pow(Cs*k[pp],2)*nn[c1][pp][0]-k[pp]*k[pp]*epsilon0/4/mi*CC[0]);
				kv1[1]=Tstep*((-2.0)*nui[pp]*vv[c1][pp][1]-pow(Cs*k[pp],2)*nn[c1][pp][1]-k[pp]*k[pp]*epsilon0/4/mi*CC[1]);

				CC[0]=0.0;
				CC[1]=0.0;
				for (int q=LL;q<=UU;q++){
					CC[0]=CC[0]+(EE[c1][(q+N/2)][0]+k1[(q+N/2)][0]/2)*(EE[c1][(q-p[pp]+N/2)][0]+k1[(q-p[pp]+N/2)][0]/2)+(EE[c1][(q+N/2)][1]+k1[(q+N/2)][1]/2)*(EE[c1][(q-p[pp]+N/2)][1]+k1[(q-p[pp]+N/2)][1]/2);
					CC[1]=CC[1]+(EE[c1][(q+N/2)][1]+k1[(q+N/2)][1]/2)*(EE[c1][(q-p[pp]+N/2)][0]+k1[(q-p[pp]+N/2)][0]/2)-(EE[c1][(q+N/2)][0]+k1[(q+N/2)][0]/2)*(EE[c1][(q-p[pp]+N/2)][1]+k1[(q-p[pp]+N/2)][1]/2);
				}

				kn2[0]=Tstep*(vv[c1][pp][0]+kv1[0]/2+SSn[0]/2*Tstep);
				kn2[1]=Tstep*(vv[c1][pp][1]+kv1[1]/2+SSn[1]/2*Tstep);
				kv2[0]=Tstep*((-2.0)*nui[pp]*(vv[c1][pp][0]+kv1[0]/2+SSn[0]/2*Tstep)-pow(Cs*k[pp],2)*(nn[c1][pp][0]+kn1[0]/2)-k[pp]*k[pp]*epsilon0/4/mi*CC[0]);
				kv2[1]=Tstep*((-2.0)*nui[pp]*(vv[c1][pp][1]+kv1[1]/2+SSn[1]/2*Tstep)-pow(Cs*k[pp],2)*(nn[c1][pp][1]+kn1[1]/2)-k[pp]*k[pp]*epsilon0/4/mi*CC[1]);

				CC[0]=0.0;
				CC[1]=0.0;
				for (int q=LL;q<=UU;q++){
					CC[0]=CC[0]+(EE[c1][(q+N/2)][0]+k2[(q+N/2)][0]/2)*(EE[c1][(q-p[pp]+N/2)][0]+k2[(q-p[pp]+N/2)][0]/2)+(EE[c1][(q+N/2)][1]+k2[(q+N/2)][1]/2)*(EE[c1][(q-p[pp]+N/2)][1]+k2[(q-p[pp]+N/2)][1]/2);
					CC[1]=CC[1]+(EE[c1][(q+N/2)][1]+k2[(q+N/2)][1]/2)*(EE[c1][(q-p[pp]+N/2)][0]+k2[(q-p[pp]+N/2)][0]/2)-(EE[c1][(q+N/2)][0]+k2[(q+N/2)][0]/2)*(EE[c1][(q-p[pp]+N/2)][1]+k2[(q-p[pp]+N/2)][1]/2);
				}

				kn3[0]=Tstep*(vv[c1][pp][0]+kv2[0]/2+SSn[0]/2*Tstep);
				kn3[1]=Tstep*(vv[c1][pp][1]+kv2[1]/2+SSn[1]/2*Tstep);
				kv3[0]=Tstep*((-2.0)*nui[pp]*(vv[c1][pp][0]+kv2[0]/2+SSn[0]/2*Tstep)-pow(Cs*k[pp],2)*(nn[c1][pp][0]+kn2[0]/2)-k[pp]*k[pp]*epsilon0/4/mi*CC[0]);
				kv3[1]=Tstep*((-2.0)*nui[pp]*(vv[c1][pp][1]+kv2[1]/2+SSn[1]/2*Tstep)-pow(Cs*k[pp],2)*(nn[c1][pp][1]+kn2[1]/2)-k[pp]*k[pp]*epsilon0/4/mi*CC[1]);


				CC[0]=0.0;
				CC[1]=0.0;
				for (int q=LL;q<=UU;q++){
					CC[0]=CC[0]+(EE[c1][(q+N/2)][0]+k3[(q+N/2)][0]/2)*(EE[c1][(q-p[pp]+N/2)][0]+k3[(q-p[pp]+N/2)][0]/2)+(EE[c1][(q+N/2)][1]+k3[(q+N/2)][1]/2)*(EE[c1][(q-p[pp]+N/2)][1]+k3[(q-p[pp]+N/2)][1]/2);
					CC[1]=CC[1]+(EE[c1][(q+N/2)][1]+k3[(q+N/2)][1]/2)*(EE[c1][(q-p[pp]+N/2)][0]+k3[(q-p[pp]+N/2)][0]/2)-(EE[c1][(q+N/2)][0]+k3[(q+N/2)][0]/2)*(EE[c1][(q-p[pp]+N/2)][1]+k3[(q-p[pp]+N/2)][1]/2);
				}

				kn4[0]=Tstep*(vv[c1][pp][0]+kv3[0]+SSn[0]*Tstep);
				kn4[1]=Tstep*(vv[c1][pp][1]+kv3[1]+SSn[1]*Tstep);
				kv4[0]=Tstep*((-2.0)*nui[pp]*(vv[c1][pp][0]+kv3[0]+SSn[0]*Tstep)-pow(Cs*k[pp],2)*(nn[c1][pp][0]+kn3[0])-k[pp]*k[pp]*epsilon0/4/mi*CC[0]);
				kv4[1]=Tstep*((-2.0)*nui[pp]*(vv[c1][pp][1]+kv3[1]+SSn[1]*Tstep)-pow(Cs*k[pp],2)*(nn[c1][pp][1]+kn3[1])-k[pp]*k[pp]*epsilon0/4/mi*CC[1]);


				vv[c2][pp][0]=vv[c1][pp][0]+(kv1[0]+2*kv2[0]+2*kv3[0]+kv4[0])/6+SSn[0]*Tstep;
				vv[c2][pp][1]=vv[c1][pp][1]+(kv1[1]+2*kv2[1]+2*kv3[1]+kv4[1])/6+SSn[1]*Tstep;
				nn[c2][pp][0]=nn[c1][pp][0]+(kn1[0]+2*kn2[0]+2*kn3[0]+kn4[0])/6;
				nn[c2][pp][1]=nn[c1][pp][1]+(kn1[1]+2*kn2[1]+2*kn3[1]+kn4[1])/6;
				nn[c2][N/2][0]=0.0;
				nn[c2][N/2][1]=0.0;

				if (pp>=1){
					nn[c2][N-pp][0]=nn[c2][pp][0];
					nn[c2][N-pp][1]=-nn[c2][pp][1];
				}

			}

			if ( (tt1)%res == 1){
				for (pp=0;pp<N;pp++){
					total_EE[counter1*N*2+pp*2+0]=EE[c2][pp][0];
					total_EE[counter1*N*2+pp*2+1]=EE[c2][pp][1];
					total_nn[counter1*N*2+pp*2+0]=nn[c2][pp][0];
					total_nn[counter1*N*2+pp*2+1]=nn[c2][pp][1];
				}
				counter1++;
			}

			if (counter1==20000){
				fwrite(total_EE, 1, sizeof(total_EE), EE_out);
				fwrite(total_nn, 1, sizeof(total_nn), nn_out);
				counter1=0;
			}

		}

	if (counter1>0){
		fwrite(total_EE, 1, sizeof(double)*counter1*N*2, EE_out);
		fwrite(total_nn, 1, sizeof(double)*counter1*N*2, nn_out);
	}


	fclose(EE_out);
	fclose(nn_out);
}//relizations
}//Nvbeam
}//Nnbeam


	now = time(0);
	current = localtime(&now);
	int EndTime[3]={current->tm_hour, current->tm_min, current->tm_sec};
	printf("Elapsed Time: %i:%i:%i\n", EndTime[0]-StartTime[0], EndTime[1]-StartTime[1], EndTime[2]-StartTime[2]);


	return EXIT_SUCCESS;
}


int sign (double a) {
	return (a>0)-(a<0);
}


int max1(int a, int b){
	if (a>b)
		return a;
	else
		return b;
}

int min1(int a, int b){
	if (a<b)
		return a;
	else
		return b;
}


void Xsection(double& Xsec_ion, double& Xsec_pl,double k){
	double alpha=1.0/(k*lambdaD);
	Xsec_ion=2*pi/(1+pow(alpha,2))*(Z*pow(alpha,4)/(1+pow(alpha,2)+pow(alpha,2)*(Z*Te/Ti)));
	double XX=2*pi*(1+pow(alpha,2)*Z*Te/Ti)/(1+pow(alpha,2)+pow(alpha,2)*(Z*Te/Ti));
	Xsec_pl=XX-Te/Ti*Xsec_ion;
}

