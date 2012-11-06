#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <queue>
#include <sstream>
#include <vector>

#include "SampleTranslateCheck.h"
#include "AllObjects.h"
#include "EnergyInCluster.h"
#include "ParticlesInCluster.h"

#define PI 3.1415926536

class particle
{
public:
	string type;
	int index;
	double px;
	double py;
	double pz;
	double p0;
	
	void clear()
	{
		type.clear();
		index = -1;
		px = 0.;
		py = 0.;
		pz = 0.;
		p0 = 0.;
	}
};

//      common /pgsrec/ numobj,dumobj,          ! number of reconstructed objects
//     .                indobj(nobjmx),         ! index to HEPEVT particle (where relevant)
//     .                typobj(nobjmx),         ! reconstructed type:  0 = photon
//                                              !                      1 = electron
//                                              !                      2 = muon
//                                              !                      3 = tau (hadronic)
//                                              !                      4 = jet
//                                              !                      5 = heavy charged
//     .                pobj(4,nobjmx),         ! four vector of reconstructed object
//     .                qobj(nobjmx),           ! charge of reconstructed object
//     .                vecobj(10,nobjmx),      ! interesting object quantities (see below)
//     .                unique(nobjmx)          ! true for object if it is uniquely identified
//                                              ! and passes cuts in pgs_object_cuts
//c -----------------------------------------------------------------------------------------------------------
//c  type            1           2          3        4        5        6        7        8        9       10
//c -----------------------------------------------------------------------------------------------------------
//c 0  photon     EM energy  HAD energy  track E   N(trk)     -       ET      etiso    ptiso    hadem     ep
//c 1  electron    "   "      "     "       "        "        -       ET      etiso ptiso-pttrk hadem     ep
//c 2  muon        "   "      "     "       "        "     trkiso E  ptiso      -        -       -        -
//c 3  tau         "   "      "     "       "        "      width    mtau     ptmax   trk iso    npi0   pi0 iso
//c 4  jet         "   "      "     "       "        "        "      flavor   c,b tags   
//c -----------------------------------------------------------------------------------------------------------
//c
//c   taus:        vecobj( 8,numobj) =  sum pt of tracks not in cone
//c                vecobj( 9,numobj) =  number of pi0 in cone
//c                vecobj(10,numobj) =  sum pt of pi0's not in cone 
//c
//c  b, c tagging: vecobj(6,iobj) = 21 (g);  2,1,3 (u,d,s);  4 (c) 5 (b)
//c                vecobj(7,iobj) non-zero if loose b tag
//c                vecobj(8,iobj) non-zero if tight b tag
//c
//c    --> tight and loose tags include fake rates for gluon, uds, c and b jets
class recoJet
{
public:
	string type;
	int index;
	double px;
	double py;
	double pz;
	double p0;
	double theVec1;
	double theVec2;
	double theVec3;
	double theVec4;
	double theVec5;
	double theVec6;
	double theVec7;
	double theVec8;
	double theVec9;
	double theVec10;

	void clear()
	{
		type.clear();
		index = -1;
		px = 0.;
		py = 0.;
		pz = 0.;
		p0 = 0.;
		theVec1 = 0.;
		theVec2 = 0.;
		theVec3 = 0.;
		theVec4 = 0.;
		theVec5 = 0.;
		theVec6 = 0.;
		theVec7 = 0.;
		theVec8 = 0.;
		theVec9 = 0.;
		theVec10 = 0.;
	}
};

class visTau
{
public:
	string type;
	int index;
	int sign;
	double px;
	double py;
	double pz;
	double p0;
	double pxvis;
	double pyvis;
	double pzvis;
	double p0vis;
	double ptiso;

	void clear()
	{
		type.clear();
		index = -1;
		sign = 0;
		px = 0.;
		py = 0.;
		pz = 0.;
		p0 = 0.;
		pxvis = 0.;
		pyvis = 0.;
		pzvis = 0.;
		p0vis = 0.;
		ptiso = 0.;
	}
};

// squark -> chargino + quark1 -> neutralino1 + W + quark1
//        -> neutralino1 + quark2 + quark3 + quark1
class wjet
{
public:
	particle squark;
	int nquarks;
	particle quarks[10];
	particle chargino;
	particle w;
	particle neutralino1;
	int njetsLeading;
	recoJet jetsLeading[10];
	int njetsFromW;
	recoJet jetsFromW[10];
	
	void clear()
	{
		squark.clear();
		chargino.clear();
		w.clear();
		neutralino1.clear();
		
		nquarks = 0;
		njetsLeading = 0;
		njetsFromW = 0;
		for (int i = 0; i < 10; i++)
		{
			quarks[i].clear();
			jetsLeading[i].clear();
			jetsFromW[i].clear();
		}
	}
};

// stop1 -> top + neutralino1 -> W + b + neutralino1
//       -> quark1 + quark2 + b + neutralino1
class stopChain
{
public:
	particle stop1;
	particle neutralino1;
	particle top;
	particle w;
	particle bottom;
	int nquarksFromW;
	particle quarksFromW[10];
	int njetsFromW;
	recoJet jetsFromW[10];
	int nbJets;
	recoJet bJets[10];

	void clear()
	{
		stop1.clear();
                neutralino1.clear();
		top.clear();
		w.clear();
		bottom.clear();

		nquarksFromW = 0;
		njetsFromW = 0;
                nbJets = 0;
		for (int i = 0; i < 10; i++)
		{
			
			quarksFromW[i].clear();
			jetsFromW[i].clear();
                        bJets[i].clear();
		}
	}
};

class wboson
{
public:
	particle w;
	int nquarks;
	particle quarks[10];
	int nelectrons;
	particle electrons[10];
	int njetsFromW;
	recoJet jetsFromW[10];
	
	void clear()
	{
		w.clear();
		
		njetsFromW = 0;
		nquarks = 0;
		nelectrons = 0;
		for (int i = 0; i < 10; i++) {
			quarks[i].clear();
			electrons[i].clear();
			jetsFromW[i].clear();
		}
	}
};
// squark -> neutralino2 + quark1 -> stau1 + tau1 + quark1 
//        -> neutralino1 + tau2 + tau1 + quark1
class jettautau
{
public:
	particle squark;
	int nquarks;
	particle quarks[10]; //Question7: why more than 1 quark ?
	particle neutralino2;
	visTau tau1;
	visTau tau2;
	particle stau;
	particle neutralino1;
	int njetsLeading;
	recoJet jetsLeading[10];

	void clear()
	{
		squark.clear();
		neutralino2.clear();
		tau1.clear();
		tau2.clear();
		stau.clear();
		neutralino1.clear();
		
		nquarks = 0;
		njetsLeading = 0;
		for (int i = 0; i < 10; i++)
		{
			quarks[i].clear();
			jetsLeading[i].clear();
		}
	}
};

using namespace std;



// Instantiate the class, providing input and output filenames.
SampleTranslateCheck::SampleTranslateCheck(string anInFileName, string anOutFileName)
{
	inFileName = anInFileName;
	outFileName = anOutFileName;
	totalEvents = 0;
//	cout << "SampleTranslateCheck object created." << endl; 
}



int SampleTranslateCheck::GetTotalEvents()
{
	return totalEvents;
}

// returns cos(theta)for two 3-vectors
float SampleTranslateCheck::determineCos_theta(float px1, float py1, float pz1,
                                               float px2, float py2, float pz2 )
{
  float cos_theta;
  cos_theta = (px1*px2+py1*py2+pz1*pz2)
         /sqrt(px1*px1+py1*py1+pz1*pz1)
         /sqrt(px2*px2+py2*py2+pz2*pz2);
  
  return cos_theta;
}


float SampleTranslateCheck::determinePT(float p_x, float p_y)
{
	float pt;
	pt = sqrt(p_x*p_x + p_y*p_y);
	return pt;
}


float SampleTranslateCheck::determineEta(float p_x, float p_y, float p_z)
{
  float theta, eta, p_t;
  
  p_t = sqrt(p_x*p_x + p_y*p_y);

  theta = atan(fabs(p_t / p_z));
  if(p_z < 0.)
    theta = PI - theta;
  eta = -log(tan(theta / 2));
  
  return eta;
}


// This subroutine performs the matching.
int SampleTranslateCheck::TranslateEvents(int numberOfEvents)
{
// Dummy string, int of double in the text data file
	string dummyString;
	int dummyInt;
        int totalEvents;
        int totalEvtPass;
        int numTotalElectrons;
        int numTotalMuons;
        int numTotalTaus;
        int numTotalTrueElectrons;
        int numTotalFalseElectrons;
        int numTotalTrueMuons;
        int numTotalFalseMuons;

        int numElec;
        int numFalseElec6et;
        int numFalseElec06et;
        int numFalseElec7etiso;
        int numFalseElec67etiso;
        int numFalseElec8ptiso;
        int numFalseElec678ptiso;
        int numFalseElec10epmin;
        int numFalseElec67810epmin;
        int numFalseElec10epmax;
        int numFalseElec67810epminmax;
        int numFalseElec6781011cos_theta;
        int numTrueElec;

        int numTotalWs;
        int numTotalElectronsFromW;
        int numTotalMuonsFromW;
        int numTotalTausFromW;
        int numTotalQuarksFromW;
        int numTotalW2Electron;
        int numTotalW2Muon;
        int numTotalW2Tau;
        int numTotalW2Quark;

        int numTotalElecWpt20;
        int numTotalElecWeta3pt20;
        int numTotalElecWpt15;
        int numTotalElecWeta3pt15;
        int numTotalElecWpt10;
        int numTotalElecWeta3pt10;
        int numTotalElecWpt5;
        int numTotalElecWeta3pt5;
        int numTotalElecWpt0;
        int numTotalElecWeta3pt0;
        int numTotalElecWeta3;

        int numTotalPgsElecsFromW;
        int numTotalTruePgsElecsFromW;
        int numTotalFalsePgsElecsFromW;
        int numTotalTruePGSElecsFromW;

        int nevhep;
        int nevhep_old;
	double dummyDouble;

	queue<int> dummyParentsQ;

	wjet wj[10];
        stopChain wbj[10];
	wboson justW[10];
	jettautau jtt[10];
	int nwj, nwbj, njustW, njtt;
	int cur, cur2;
	bool newW;

// Number of generator particles, particles depositing energy, clusters and reco objects
	int nhep, nhit, nclu, nobj;

// Cluster map - each entry is the cluster index for the cell; -1 if not assigned
	int hitsMap[320][200];

// Mapping of index of generator particle or reco object to index in AllObjects
	int object2myio[4200];

// Missing energy 
	double metcal, phimetcal, metcor, phimetcor;

// COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
//      &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
// DOUBLE PRECISION PHEP, VHEP
// Generator particles
	int istat[4000];
	int itype[4000];
	int m1[4000];
	int m2[4000];
	int d1[4000];
	int d2[4000];
	double px[4000];
	double py[4000];
	double pz[4000];
	double p0[4000];
	double m[4000];
	double v1[4000];
	double v2[4000];
	double v3[4000];
	double v4[4000];

// Cell and energy deposited of generator particles
	int hitParticle[4000];
	int ieta[4000];
	int iphi[4000];
	double emcal[4000];
	double hcal[4000];
// Energy deposited by each particle in the detector emcal + hcal
	double energyDeposited[4000];

// Cluster data - only one entry because the info goes into hitsMap
	int mulclu;
	int indclu;
	int ietaclu;
	int iphiclu;

// Variables for MC jets ???
	int firstParton[200];
	int secondParton[200];
	int thirdParton[200];
	int firstPartonFromW[10];
	int secondPartonFromW[10];
	int thirdPartonFromW[10];
	int firstPartonLeading[10];
	int secondPartonLeading[10];
	int thirdPartonLeading[10];
	int firstPartonB[10];
	int secondPartonB[10];
	int thirdPartonB[10];
	double dr, dphi;

// Variables for MC taus
	double p0vis;
	double pxvis;
	double pyvis;
	double pzvis;
	bool isLeptonic;
	bool notStatus2;
	
// To save my reco objects
	int nmyJets;
	recoJet myJets[45];
	int nmyTaus;
	visTau myTaus[45];
	int nmyLeptons;
	visTau myLeptons[45];

//      common /pgsrec/ numobj,dumobj,          ! number of reconstructed objects
//     .                indobj(nobjmx),         ! index to HEPEVT particle (where relevant)
//     .                typobj(nobjmx),         ! reconstructed type:  0 = photon
//                                              !                      1 = electron
//                                              !                      2 = muon
//                                              !                      3 = tau (hadronic)
//                                              !                      4 = jet
//                                              !                      5 = heavy charged
//     .                pobj(4,nobjmx),         ! four vector of reconstructed object
//     .                qobj(nobjmx),           ! charge of reconstructed object
//     .                vecobj(10,nobjmx),      ! interesting object quantities (see below)
//     .                unique(nobjmx)          ! true for object if it is uniquely identified
//                                              ! and passes cuts in pgs_object_cuts
// Reco objects
	int indobj[200];
	int typobj[200];
	double pobjx[200];
	double pobjy[200];
	double pobjz[200];
	double pobj0[200];
	double q[200];
	double vec1[200];
	double vec2[200];
	double vec3[200];
	double vec4[200];
	double vec5[200];
	double vec6[200];
	double vec7[200];
	double vec8[200];
	double vec9[200];
	double vec10[200];
	bool unique[200];
	bool checkUnique[200];
// isTrue will be converted to unique
	string isTrue;
// Which object is for a given cluster
	int objectOfCluster[200];
// Which particles are for a given cluster
	ParticlesInCluster thePIC[200];

	AllObjects theObjects[4200];

	EnergyInCluster theEnergyDeposited[4000];
	
//	ofstream outFile(outFileName.data(), ios::trunc);
//	ifstream inFile(inFileName.data(), ios::in);
        string inFileNamePart;
        inFileNamePart.clear();
        inFileNamePart.assign(inFileName);

        string outFileNamePart;
        outFileNamePart.clear();
        outFileNamePart.assign(outFileName);

	char tempChar;
        totalEvents = 0;
        totalEvtPass = 0;
        numTotalElectrons = 0;
        numTotalMuons = 0;
        numTotalTaus = 0;
        numTotalTrueElectrons = 0;
        numTotalFalseElectrons = 0;
        numTotalTrueMuons = 0;
        numTotalFalseMuons = 0;

        numElec = 0;
        numFalseElec6et = 0;
        numFalseElec06et = 0;
        numFalseElec7etiso = 0;
        numFalseElec67etiso = 0;
        numFalseElec8ptiso = 0;
        numFalseElec678ptiso = 0;
        numFalseElec10epmin = 0;
        numFalseElec67810epmin = 0;
        numFalseElec10epmax = 0;
        numFalseElec67810epminmax = 0;
        numFalseElec6781011cos_theta = 0;
        numTrueElec = 0;

        numTotalWs = 0;
        numTotalElectronsFromW = 0;
        numTotalMuonsFromW = 0;
        numTotalTausFromW = 0;
        numTotalQuarksFromW = 0;
        numTotalW2Electron = 0;
        numTotalW2Muon = 0;
        numTotalW2Tau = 0;
        numTotalW2Quark = 0;

        numTotalElecWpt20 = 0;
        numTotalElecWeta3pt20 = 0;
        numTotalElecWpt15 = 0;
        numTotalElecWeta3pt15 = 0;
        numTotalElecWpt10 = 0;
        numTotalElecWeta3pt10 = 0;
        numTotalElecWpt5 = 0;
        numTotalElecWeta3pt5 = 0;
        numTotalElecWpt0 = 0;
        numTotalElecWeta3pt0 = 0;
        numTotalElecWeta3 = 0;

        numTotalPgsElecsFromW = 0;
        numTotalTruePgsElecsFromW = 0;
        numTotalFalsePgsElecsFromW = 0;
        numTotalTruePGSElecsFromW = 0;

// Loop file by file
        int totalNumFiles = 0;



        for(int fileNumber=1; fileNumber<=300; fileNumber++){

                int maxNumFiles = 20;
                int fileNumbers[100]={4, 9, 13, 20, 21, 36, 45, 50, 58, 63, 66, 72, 78, 82, 85, 99, 
                                     105, 117, 128, 132, 133, 137, 151, 165, 166, 171, 180, 184, 
                                     218, 225, 226, 248, 260, 262, 263, 264, 265, 267, 268, 269, 
                                     271, 272, 273, 275, 277, 278, 279, 282, 285, 291, 294, 295, 
                                     296, 297, 298};
                bool checkFile = false;
                for(int i=0; i<100; i++)
                  if(fileNumber==fileNumbers[i]){
                     checkFile = true;
                     break;
                  }
                if(!checkFile) continue;
                if(totalNumFiles>=maxNumFiles) return 2;
                totalNumFiles++;





                inFileName.clear();
                inFileName.assign("/x1/data/nkolev/project10/output/ttjnew/");
                inFileName.append(inFileNamePart);
//                inFileName.append("-");

                outFileName.clear();
                outFileName.assign("./ttj006f_check_txt/");
                outFileName.append(outFileNamePart);
                outFileName.append("-");

		if(fileNumber >= 100) {
		   tempChar = 48 + fileNumber/100;
		   inFileName.append(1,tempChar);
		   outFileName.append(1,tempChar);
		}
		if(fileNumber >= 10) {
			tempChar = 48 + (fileNumber%100)/10; // ascii 48 -> 0
			inFileName.append(1,tempChar);
			outFileName.append(1,tempChar);
		}
		tempChar = 48 + fileNumber%10;

		inFileName.append(1,tempChar);
		outFileName.append(1,tempChar);

                inFileName.append(".dat");
//                inFileName.append(".txt");	
		cout<<inFileName<<"\n\n";
		ifstream inFile(inFileName.data(), ifstream::in);

                outFileName.append(".txt");	
		cout<<outFileName<<"\n\n";
                ofstream outFile(outFileName.data(), ofstream::trunc);

        nevhep = 0;
        nevhep_old = 0;

	if(!inFile.good() ) { // The file was bad or did not exist. Write an error message.
		cout<<"\n\nThat file does not exist: "<<inFileName<<"\n\n";
		return 0;
	}
	
//	inFile >> dummyInt; // to catch the first number of the event file, if it is the total number of events.
//	if(dummyInt > 1) dummyInt = 0;
//	if(dummyInt == 1) nevhep = dummyInt;

// Loop over event by event: entire file
	while (!inFile.eof() && (numberOfEvents > 0))//while:383
	{	
		nwj = 0;
                nwbj = 0;
		njustW = 0;
		njtt = 0;
		nmyJets = 0;
		nmyTaus = 0;
		nmyLeptons = 0;
		newW = false;
		
// Initialize my jets and taus
		for(int j=0; j<45; j++)
		{
			myJets[j].clear();
			myTaus[j].clear();
			myLeptons[j].clear();
		}
		
// Initialize theObjects for each event
//              int mcorrec;
//		vector<int> parents;
//		vector<int> children;
//		vector<double> energyDeposited;
		for (int j = 0; j < 4200; j++)
		{
			theObjects[j].clear();
		} 

// Initialize hitsMap for each event
		for (int j = 0; j < 320; j++)
		{
			for (int k = 0; k < 200; k++)
			{
				hitsMap[j][k] = -1;
			}
		}

//		double theEnergy;
//		int theCluster;
// Initialize energy depositions
		for (int j = 0; j < 4000; j++)
		{
			energyDeposited[j] = 0.;//Energy deposited by each particle in the detector emcal + hcal
			theEnergyDeposited[j].theEnergy = -1.;
			theEnergyDeposited[j].theCluster = -1;
		}	

//		int theClusterN;
//		vector<int> theParticlesInCluster;
// Initialize object of cluster relationship
		for (int j = 0; j < 200; j++)
		{
			objectOfCluster[j] = -1;//Which object is for a given cluster
			thePIC[j].theClusterN = 0;
			thePIC[j].theParticlesInCluster.clear();
			firstParton[j] = -1;
			secondParton[j] = -1;
			thirdParton[j] = -1;
		}	 

// Initialize decay types
		for (int j = 0; j < 10; j++)
		{
			wj[j].clear();
                        wbj[j].clear();
			justW[j].clear();
			jtt[j].clear();
		}


// Read generator particles
		// catch the files which have total number of events as their first entry.
                if(nevhep !=0) nevhep_old = nevhep;
                inFile>>nevhep;
//		if(dummyInt == 0){
//			inFile >> dummyInt;
//                        nevhep = dummyInt;
//                }
                 if(nevhep!=nevhep_old+1){
                 cout<<"nevhep_old = "<<nevhep_old<<endl;
                 cout<<"nevhep = "<<nevhep<<endl;
                 break;
                 }

                if(nevhep % 100 == 0)
                cout<<"nevhep = "<<nevhep<<endl;

//                if(totalEvents==99999)
                if(totalEvents==numberOfEvents)
                { cout<<"###############"<<endl;
                  cout<<"totalEvents  = "<<totalEvents<<endl;
                  cout<<"totalEvtPass = "<<totalEvtPass<<endl;
                  cout<<"###############"<<endl;
                  cout<<"numTotalElectrons      = "<<numTotalElectrons<<endl;
                  cout<<"numTotalMuons          = "<<numTotalMuons<<endl;
                  cout<<"numTotalTaus           = "<<numTotalTaus<<endl<<endl;
                  cout<<"numTotalTrueElectrons  = "<<numTotalTrueElectrons<<endl;
                  cout<<"numTotalFalseElectrons = "<<numTotalFalseElectrons<<endl;
                  cout<<"numTotalTrueMuons      = "<<numTotalTrueMuons<<endl;
                  cout<<"numTotalFalseMuons     = "<<numTotalFalseMuons<<endl<<endl;

                  cout<<"###############"<<endl;
                  cout<<"numElec     = "<<numElec<<endl;
                  cout<<"numTrueElec = "<<numTrueElec<<endl;
                  cout<<"-------"<<endl;
                  cout<<"numFalseElec06et             = "<<numFalseElec06et<<endl;
                  cout<<"numFalseElec67etiso          = "<<numFalseElec67etiso<<endl;
                  cout<<"numFalseElec678ptiso         = "<<numFalseElec678ptiso<<endl;
                  cout<<"numFalseElec67810epmin       = "<<numFalseElec67810epmin<<endl;
                  cout<<"numFalseElec67810epminmax    = "<<numFalseElec67810epminmax<<endl;
                  cout<<"numFalseElec6781011cos_theta = "<<numFalseElec6781011cos_theta<<endl;
                  cout<<"-------"<<endl;
                  cout<<"numFalseElec6et     = "<<numFalseElec6et<<endl;
                  cout<<"numFalseElec7etiso  = "<<numFalseElec7etiso<<endl;
                  cout<<"numFalseElec8ptiso  = "<<numFalseElec8ptiso<<endl;
                  cout<<"numFalseElec10epmin = "<<numFalseElec10epmin<<endl;
                  cout<<"numFalseElec10epmax = "<<numFalseElec10epmax<<endl<<endl;

                  cout<<"###############"<<endl;
                  cout<<"numTotalWs             = "<<numTotalWs<<endl;
                  cout<<"numTotalElectronsFromW = "<<numTotalElectronsFromW<<endl;
                  cout<<"numTotalMuonsFromW     = "<<numTotalMuonsFromW<<endl;
                  cout<<"numTotalTausFromW      = "<<numTotalTausFromW<<endl;
                  cout<<"numTotalQuarksFromW    = "<<numTotalQuarksFromW<<endl;
                  cout<<"numTotalW2Electron     = "<<numTotalW2Electron<<endl;
                  cout<<"numTotalW2Muon         = "<<numTotalW2Muon<<endl;
                  cout<<"numTotalW2Tau          = "<<numTotalW2Tau<<endl;
                  cout<<"numTotalW2Quark        = "<<numTotalW2Quark<<endl<<endl;

                  cout<<"numTotalElecWeta3      = "<<numTotalElecWeta3<<endl;
                  cout<<"numTotalElecWeta3pt0   = "<<numTotalElecWeta3pt0<<endl;
                  cout<<"numTotalElecWeta3pt5   = "<<numTotalElecWeta3pt5<<endl;
                  cout<<"numTotalElecWeta3pt10  = "<<numTotalElecWeta3pt10<<endl;
                  cout<<"numTotalElecWeta3pt15  = "<<numTotalElecWeta3pt15<<endl;
                  cout<<"numTotalElecWeta3pt20  = "<<numTotalElecWeta3pt20<<endl;
                  cout<<"numTotalElecWpt0       = "<<numTotalElecWpt0<<endl;
                  cout<<"numTotalElecWpt5       = "<<numTotalElecWpt5<<endl;
                  cout<<"numTotalElecWpt10      = "<<numTotalElecWpt10<<endl;
                  cout<<"numTotalElecWpt15      = "<<numTotalElecWpt15<<endl;
                  cout<<"numTotalElecWpt20      = "<<numTotalElecWpt20<<endl<<endl;

                  cout<<"###############"<<endl;
                  cout<<"numTotalPgsElecsFromW      = "<<numTotalPgsElecsFromW<<endl;
                  cout<<"numTotalTruePgsElecsFromW  = "<<numTotalTruePgsElecsFromW<<endl;
                  cout<<"numTotalFalsePgsElecsFromW = "<<numTotalFalsePgsElecsFromW<<endl;
                  cout<<"numTotalTruePGSElecsFromW  = "<<numTotalTruePGSElecsFromW<<endl<<endl;

                  return 1; 
                }

		inFile >> dummyString >> dummyString;
		inFile >> nhep;
		for (int j = 0; j < nhep; j++)
		{
			inFile >> dummyInt >> istat[j] >> itype[j] >> m1[j] 
                               >> m2[j] >> d1[j] >> d2[j] >> px[j] >> py[j] 
                               >> pz[j] >> p0[j] >> m[j] >> v1[j] >> v2[j] >> v3[j] >> v4[j]; 
			m1[j]--; //because the first entry is saved as index 0
			m2[j]--; //because the first entry is saved as index 0
			d1[j]--; //because the first entry is saved as index 0
			d2[j]--; //because the first entry is saved as index 0
                        //set these values to the class members
                        //void AllObjects::Initialize(int anMcorrec, int anIndex, int aStatus, int aType,
                        //    double aPx, double aPy, double aPz, double aP0, double aMass,
                        //    double aCharge, double aVec6, double aVec7, double aVec8)
                        // if not -1 is either 1 - particle or 31 - reco object
			theObjects[j].Initialize(1, j, istat[j], itype[j], px[j], py[j], pz[j], p0[j], m[j], -99,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, true);
			
			if (m1[j] != -1 && m2[j] != -1) //the m1 and m2 numbers read from file are not 0 
			{
				for (int k = 0; k <= m2[j] - m1[j]; k++)
				{
			           theObjects[j].parents.push_back(m1[j] + k);//??? fill in numbers m1-1,m1,...,m2-1
				}                                             //fill in indices of parents for each entry
			}
			else if (m1[j] != -1)//the m1 read from file is not 0 while m2 is 0
			{
				theObjects[j].parents.push_back(m1[j]);//??? fill in number m1-1
			}
			
		}
//                cout<<"!!!debug605 "<<endl;

		for (int j = 0; j < nhep; j++)
		{       //Returns the number of elements in the vector container.
			for (int k = 0; k < theObjects[j].parents.size(); k++)
			{
				theObjects[theObjects[j].parents[k]].children.push_back(j);//fill in the children index for each entry
			}
		}

//Question3: "cluster"?? Read clusters and fill in the cluster map hitsMap
		inFile >> dummyString >> dummyString;//read "EndStdhep" and "BeginCluster"
		inFile >> nclu;
		for (int j = 0; j < nclu; j++)
		{
                        // Cluster data - only one entry because the info goes into hitsMap
	                // int mulclu;
	                // int indclu;
	                // int ietaclu;
	                // int iphiclu;
			inFile >> mulclu;
			thePIC[j].theClusterN = mulclu;
			for (int k = 0; k < mulclu; k++)
			{
				inFile >> indclu >> ietaclu >> iphiclu;
				ietaclu--;//save as the index position, so -1
				iphiclu--;//save as the index position, so -1
				indclu--; //save as the index position, so -1
				hitsMap[ietaclu][iphiclu] = indclu;
			} 
		}
//                cout<<"!!!debug636 "<<endl;

// Read the cells and energy deposited by each particle
		inFile >> dummyString >> dummyString;//read "EndCluster" and "BeginCalhit"
		inFile >> nhit;
		for (int j = 0; j < nhit; j++)
		{
			inFile >> hitParticle[j] >> ieta[j] >> iphi[j] >> emcal[j] >> hcal[j];
			hitParticle[j] = hitParticle[j] - 1;//save as the index position, so -1
			ieta[j] = ieta[j] - 1;//save as the index position, so -1
			iphi[j] = iphi[j] - 1;//save as the index position, so -1
			if (emcal[j] + hcal[j] > 0.)
			{ 
				theEnergyDeposited[hitParticle[j]].theEnergy = emcal[j] + hcal[j];
				theEnergyDeposited[hitParticle[j]].theCluster = hitsMap[ieta[j]][iphi[j]];
		                // int theClusterN;
		                // vector<int> theParticlesInCluster;
				// hitsMap[ietaclu][iphiclu] = indclu;
				thePIC[hitsMap[ieta[j]][iphi[j]]].theParticlesInCluster.push_back(hitParticle[j]);
			}
		}
//                cout<<"!!!debug657 "<<endl;

// Read reco objects
		inFile >> dummyString >> dummyString;//read " EndCalhit" and "BeginReco"
		inFile >> nobj; 
		for (int j = 0; j < nobj; j++)
		{
             		inFile >> dummyInt >> indobj[j] >> typobj[j] >> pobjx[j] >> pobjy[j] >> pobjz[j] >> pobj0[j]//Question2: "indobj[j]"??
                               >> q[j] >> vec1[j] >> vec2[j] >> vec3[j] >> vec4[j] >> vec5[j]
                               >> vec6[j] >> vec7[j] >> vec8[j] >> vec9[j] >> vec10[j] >> isTrue;
			indobj[j] = indobj[j] - 1;//save as the index position, so -1
			if (isTrue == "T") 
			{
				unique[j] = true; //bool unique[200];
			}
			else
			{
				unique[j] = false; //bool unique[200];
			}

                        //!!! int objectOfCluster[200]; Which object is for a given cluster
                        // 0 = photon; 1 = electron; 2 = muon; 3 = tau (hadronic); 4 = jet; 5 = heavy charged
			if ((typobj[j] == 3 || typobj[j] == 4) && unique[j]) objectOfCluster[indobj[j]] = j;//???

			if (typobj[j]==1) //if electron, does require it must be true.
			{       //set these values to the class members
                                //void AllObjects::Initialize(int anMcorrec, int anIndex, int aStatus, int aType,
                                //    double aPx, double aPy, double aPz, double aP0, double aMass,
                                //    double aCharge, double aVec6, double aVec7, double aVec8)
                                // if not -1 is either 1 - particle or 31 - reco object
				theObjects[4000 + j].Initialize(31, 4000 + j, indobj[j], typobj[j], pobjx[j], pobjy[j], pobjz[j],
					pobj0[j], -1., q[j], vec1[j], vec2[j], vec3[j], vec4[j], vec5[j], vec6[j],
                                        vec7[j], vec8[j], vec9[j], vec10[j], unique[j]);
			}
			else if (unique[j]) //if not electron, requires must be true.
			{       //set these values to the class members
                                //void AllObjects::Initialize(int anMcorrec, int anIndex, int aStatus, int aType,
                                //    double aPx, double aPy, double aPz, double aP0, double aMass,
                                //    double aCharge, double aVec6, double aVec7, double aVec8)
                                // if not -1 is either 1 - particle or 31 - reco object
				theObjects[4000 + j].Initialize(31, 4000 + j, indobj[j], typobj[j], pobjx[j], pobjy[j], pobjz[j],
					pobj0[j], -1., q[j], vec1[j], vec2[j], vec3[j], vec4[j], vec5[j], vec6[j],
                                        vec7[j], vec8[j], vec9[j], vec10[j], unique[j]);
			}

                        //!!!check electrons efficiency
                        if(typobj[j]==1 ){//if:740
                           numTotalElectrons++;
                           if(unique[j]) numTotalTrueElectrons++;
                           else          numTotalFalseElectrons++;


                        }//endif:740

                        if(typobj[j]==2 ){//if:747
                           numTotalMuons++;
                           if(unique[j]) numTotalTrueMuons++;
                           else          numTotalFalseMuons++;

                        }//endif:747


		}
//                cout<<"!!!debug689 "<<endl;

// Read missing energy 
		inFile >> metcal >> phimetcal >> metcor >> phimetcor;

// Read final strings for event
		inFile >> dummyString >> dummyString;//read "EndReco" and "EndEvent"

		totalEvents++;
// ================================================================================
// identify generator taus
		for (int j = 0; j < nhep; j++) //for:925
		{       //tau
			isLeptonic = false;
			// Choose only hadronic taus.
			if(abs(theObjects[j].type) == 15)//tau
                            for(int k=0; k<theObjects[j].children.size(); k++)//check each child of this tau
				if (abs(theObjects[theObjects[j].children[k]].type) == 11 ||//elctron
				    abs(theObjects[theObjects[j].children[k]].type) == 13)//muon
						isLeptonic = true;			
			// While i'm at it, let's now fix the tau momentum to be always visible.
			if (abs(theObjects[j].type) == 15 && theObjects[j].status == 2 && !isLeptonic)//if:936
			{
				p0vis = 0.;
				pxvis = 0.;
				pyvis = 0.;
				pzvis = 0.;
				//consider the children of this tau: should be no eletron or muon
				for (int k = 0; k < theObjects[j].children.size(); k++)
				{
					if (theObjects[j].children[k] != -1 &&//not initial value
					    abs(theObjects[theObjects[j].children[k]].type) != 16)//not neutrino_tau
					{
						p0vis += theObjects[theObjects[j].children[k]].p0;//subtract neutrino_tau momentum
						pxvis += theObjects[theObjects[j].children[k]].px;
						pyvis += theObjects[theObjects[j].children[k]].py;
						pzvis += theObjects[theObjects[j].children[k]].pz;
					}
				}

				if(nmyTaus < 45)
				{
					myTaus[nmyTaus].index = j;//nhep
					myTaus[nmyTaus].type.assign("tau");
                                        //+1 for tau^+ or -1 for tau^-
					myTaus[nmyTaus].sign = (theObjects[j].type)/(abs(theObjects[j].type));
					myTaus[nmyTaus].p0 = theObjects[j].p0;
					myTaus[nmyTaus].px = theObjects[j].px;
					myTaus[nmyTaus].py = theObjects[j].py;
					myTaus[nmyTaus].pz = theObjects[j].pz;
					myTaus[nmyTaus].p0vis = p0vis;
					myTaus[nmyTaus].pxvis = pxvis;
					myTaus[nmyTaus].pyvis = pyvis;
					myTaus[nmyTaus].pzvis = pzvis;
					nmyTaus++;
                                        numTotalTaus++;
				}
				else
				{
					outFile<<"CATASTROPHE!!! out of space for taus!!\n\n";
				}
			}//endif:936

                     }//endfor:925

// ================================================================================
// set cuts
                int numEventJets = 0;
                int numEventElectrons = 0;
                int numEventMuons = 0;
                int numEventTaus = 0;

		for(int i = 4000; i < 4000 + nobj; i++){

                  if(theObjects[i].type == 4 && theObjects[i].theUnique && theObjects[i].mcorrec != -1){//jet
                    float ptJet  = determinePT( theObjects[i].px, theObjects[i].py);
                    float etaJet = determineEta(theObjects[i].px, theObjects[i].py, theObjects[i].pz);
                    if(ptJet>30. && fabs(etaJet)<=2.5) numEventJets++;
//                    if(ptJet>20.) numEventJets++;
                  }

                  if(theObjects[i].type == 1 && theObjects[i].theUnique && theObjects[i].mcorrec != -1){
                    float ptElectron = determinePT( theObjects[i].px, theObjects[i].py);
                    float etaElectron = determineEta(theObjects[i].px, theObjects[i].py, theObjects[i].pz);
                    float ptisoElectron = theObjects[i].theVec8;
                    if(ptElectron>10. && ptisoElectron<5. && fabs(etaElectron)<2.5)
                       numEventElectrons++;

                  }

                  if(theObjects[i].type == 2 && theObjects[i].theUnique && theObjects[i].mcorrec != -1){
                    float ptMuon = determinePT( theObjects[i].px, theObjects[i].py);
                    float etaMuon = determineEta(theObjects[i].px, theObjects[i].py, theObjects[i].pz);
                    float ptisoMuon = theObjects[i].theVec6;
                    if(ptMuon>10. && ptisoMuon<5. && fabs(etaMuon)<2.5)
                       numEventMuons++;

                  }

                }

              for(int i=0; i<nmyTaus; i++ ){
                   float ptvisTau = determinePT( myTaus[i].pxvis, myTaus[i].pyvis);
                   float etavisTau = determineEta( myTaus[i].pxvis, myTaus[i].pyvis, myTaus[i].pzvis);

                   if(ptvisTau>=10. && etavisTau<=2.1) numEventTaus++;

             }
//

//                if(metcal < 200. ) continue;
//                if(numEventJets<5) continue;

//                if(numEventMuons>0) continue;
//                if(numEventTaus>0) continue;
//                if(numEventElectrons>0) continue;

                totalEvtPass++;



// ================================================================================
// loop over Reco objects in object list...
// 0 = photon; 1 = electron; 2 = muon; 3 = tau (hadronic); 4 = jet; 5 = heavy charged
		for(int i = 4000; i < 4000 + nobj; i++){//for:808

                   checkUnique[i] = true;

                   if(theObjects[i].type==1){//if:889 type1
                     numElec++;

                     if(theObjects[i].theVec6 < 5.0){//if:6et
                        checkUnique[i] = false;
                        numFalseElec6et++;
                     }
                     if(!checkUnique[i]) numFalseElec06et++;

                     if(theObjects[i].theVec7 / theObjects[i].theVec6 > 0.10){//if:7etiso
                        checkUnique[i] = false;
                        numFalseElec7etiso++;
                     }
                     if(!checkUnique[i]) numFalseElec67etiso++;

                     if(theObjects[i].theVec8 > 5.0){//if:8ptiso
                        checkUnique[i] = false;
                        numFalseElec8ptiso++;
                     }
                     if(!checkUnique[i]) numFalseElec678ptiso++;

                     if(theObjects[i].theVec10 <= 0.50){//if:10epmin
                        checkUnique[i] = false;
                        numFalseElec10epmin++;
                     }
                     if(!checkUnique[i]) numFalseElec67810epmin++;

                     if(theObjects[i].theVec10 > 1.50){//if:10epmax
                        checkUnique[i] = false;
                        numFalseElec10epmax++;
                     }
                     if(!checkUnique[i]) numFalseElec67810epminmax++;

                     if(checkUnique[i]){

                        for(int j=4000; j < i; j++){//loop over previous objects
                            if( theObjects[j].type!=theObjects[i].type
                             && theObjects[i].type!=2
                             && theObjects[j].type!=2){
                                if(determineCos_theta(theObjects[i].px, theObjects[i].py, theObjects[i].pz,
                                                      theObjects[j].px, theObjects[j].py, theObjects[j].pz  ) > 0.9848){
                                   if(checkUnique[j]) checkUnique[i]=false;
                                }

                            }
                        }
                     }
                     if(!checkUnique[i]) numFalseElec6781011cos_theta++;

                     if(checkUnique[i]) numTrueElec++;

                   }//endif:889 type1

                }//endfor:808


// ================================================================================
// Loop over the entire generated particles: hepevt
		outFile << "Event  " <<nevhep<<endl;
		for (int j = 0; j < nhep; j++)//for:578
		{
			// This looks for any W boson
			if (abs(theObjects[j].type) == 24 && theObjects[j].status == 3) {//if:581
				justW[njustW].w.index = j;
                                numTotalWs++;

				for (int k = 0; k < theObjects[j].children.size(); k++)
				{       // 1=d; 2=u; 3=s; 4=c; 5=b; 6=t;
                                        // only record W which decays into quarks(hardronly decaying )
					if (abs(theObjects[theObjects[j].children[k]].type) <= 4)
					{
						justW[njustW].quarks[justW[njustW].nquarks].index = theObjects[j].children[k];
						justW[njustW].nquarks++;
                                                numTotalQuarksFromW++;
					}
					else if (abs(theObjects[theObjects[j].children[k]].type) == 11) //status==3
					{//elseif:1004
                                           numTotalElectronsFromW++;
				           justW[njustW].electrons[justW[njustW].nelectrons].index = theObjects[j].children[k];
					   justW[njustW].nelectrons++;

                                           int indElc = theObjects[j].children[k];
                                           //if not Final state, find final state
                                           if(theObjects[indElc].status!=1){
                                              for(int iChildElc = 0; iChildElc < theObjects[indElc].children.size(); iChildElc++){
                                                  if(theObjects[theObjects[indElc].children[iChildElc]].type == theObjects[indElc].type
                                                  && theObjects[theObjects[indElc].children[iChildElc]].status == 1)
                                                     indElc =  theObjects[indElc].children[iChildElc];
                                              }

                                           }
                                           if(abs(theObjects[indElc].status)!=1) cout<<"!!! event #"<<nevhep<<": status != 1"<<endl;

                                           float ptElectronsFromW  = determinePT(theObjects[indElc].px, theObjects[indElc].py);
                                           float etaElectronsFromW = determineEta(theObjects[indElc].px, theObjects[indElc].py, theObjects[indElc].pz);

                                           if(ptElectronsFromW > 0)  numTotalElecWpt0++;
                                           if(ptElectronsFromW > 5)  numTotalElecWpt5++;
                                           if(ptElectronsFromW > 10) numTotalElecWpt10++;
                                           if(ptElectronsFromW > 15) numTotalElecWpt15++;
                                           if(ptElectronsFromW > 20) numTotalElecWpt20++;

                                           if(fabs(etaElectronsFromW) < 3){//if:1009
                                             numTotalElecWeta3++;

                                             if(ptElectronsFromW > 10){//if:1013
                                                numTotalElecWeta3pt10++;

               		                        outFile << "JustW"<<endl;
               		                        outFile << "elctronFromW"<<endl;
			                        outFile  << " " 
                                                                << theObjects[indElc].p0 << "  "
					                        << theObjects[indElc].px << "  "
					                        << theObjects[indElc].py << "  "
					                        << theObjects[indElc].pz << "  "
					                        << theObjects[indElc].theVec8 << "  "
					                        << theObjects[indElc].type
                                                             / abs(theObjects[indElc].type) << endl;


                                                for(int nn = 4000; nn < 4000 + nobj; nn++){//for:1032

                                                    if(theObjects[nn].type == 1 && theObjects[nn].mcorrec != -1){//if:1034

                                                       bool findElcFromW = false;
                                                       if(theObjects[nn].status == indElc){
                                                          findElcFromW = true;
                                                       }
                                                       else for(int mm = 0; mm < theObjects[indElc].children.size(); mm++){
                                                           if(theObjects[nn].status == theObjects[indElc].children[mm]){
                                                              findElcFromW = true;
                                                           }
                                                       }

//				theObjects[4000 + j].Initialize(31, 4000 + j, indobj[j], typobj[j], pobjx[j], pobjy[j], pobjz[j],
//					pobj0[j], -1., q[j], vec1[j], vec2[j], vec3[j], vec4[j], vec5[j], vec6[j],
//                                        vec7[j], vec8[j], vec9[j], vec10[j], unique[j]);
                                                       if(findElcFromW){//if:1049
               		                                  outFile << "PGSelectronFromW"<<endl;
			                                  outFile << " " 
                                                                << theObjects[nn].status+1 << "  "//debug only
                                                                << theObjects[nn].p0 << "  "
					                        << theObjects[nn].px << "  "
					                        << theObjects[nn].py << "  "
					                        << theObjects[nn].pz << "  "
					                        << theObjects[nn].charge << "  "
					                        << theObjects[nn].theVec1 << "  "
					                        << theObjects[nn].theVec2 << "  "
					                        << theObjects[nn].theVec3 << "  "
					                        << theObjects[nn].theVec4 << "  "
					                        << theObjects[nn].theVec5 << "  "
					                        << theObjects[nn].theVec6 << "  "
					                        << theObjects[nn].theVec7 << "  "
					                        << theObjects[nn].theVec8 << "  "
					                        << theObjects[nn].theVec9 << "  "
					                        << theObjects[nn].theVec10 << "  "
					                        << theObjects[nn].theUnique << "  "
					                        << checkUnique[nn]<< endl;

                                                           numTotalPgsElecsFromW++;
                                                           if(checkUnique[nn]) numTotalTruePgsElecsFromW++;
                                                           else numTotalFalsePgsElecsFromW++;
                                                           if(theObjects[nn].theUnique) numTotalTruePGSElecsFromW++;

                                                       }//endif:1049

                                                    }//endif:1034

                                                }//endfor:1032

               		                        outFile << "EndJustW"<<endl;

                                            }//if:1013
                                          }//endif:1009
					}//endelseif:1004
					else if (abs(theObjects[theObjects[j].children[k]].type) == 13)
					{
                                                numTotalMuonsFromW++;
					}
					else if (abs(theObjects[theObjects[j].children[k]].type) == 15)
					{
                                                numTotalTausFromW++;
					}
				}


				for (int k = 0; k < theObjects[j].children.size(); k++)
				{       // 1=d; 2=u; 3=s; 4=c; 5=b; 6=t;
                                        // only record W which decays into quarks(hardronly decaying )
					if (abs(theObjects[theObjects[j].children[k]].type) <= 4)
					{
                                            numTotalW2Quark++;
                                            break;
					}
					else if (abs(theObjects[theObjects[j].children[k]].type) == 11)
					{
                                            numTotalW2Electron++;
                                            break;
					}
					else if (abs(theObjects[theObjects[j].children[k]].type) == 13)
					{
                                            numTotalW2Muon++;
                                            break;
					}
					else if (abs(theObjects[theObjects[j].children[k]].type) == 15)
					{
                                            numTotalW2Tau++;
                                            break;
					}

				}


			}//endif:581

//========================
// squark -> chargino + quark1 -> neutralino1 + W + quark1
//        -> neutralino1 + quark2 + quark3 + quark1
// squark -> neutralino2 + quark1 -> stau1 + tau1 + quark1 
//        -> neutralino1 + tau2 + tau1 + quark1
// This looks at the current particle and down the cascade to see if this event is W+jet or Jet+2tau
			if (abs(theObjects[j].type) >= 1000001 && abs(theObjects[j].type) <= 1000004 ||//if:599 regular left squarks
			    abs(theObjects[j].type) >= 2000001 && abs(theObjects[j].type) <= 2000004) //regular right squarks
			{
				wj[nwj].squark.index = j;
				jtt[njtt].squark.index = j;
				for (int k = 0; k < theObjects[j].children.size(); k++)//for:604 check children
				{
                                        //quarks from squarks, may be more than 1
					if (abs(theObjects[theObjects[j].children[k]].type) <= 4)
					{
						jtt[njtt].quarks[jtt[njtt].nquarks].index = theObjects[j].children[k];
						jtt[njtt].nquarks++;
						wj[nwj].quarks[wj[nwj].nquarks].index = theObjects[j].children[k];
						wj[nwj].nquarks++;
					}

                                        // squark -> neutralino2 + quark1 -> stau1 + tau1 + quark1 
                                        //        -> neutralino1 + tau2 + tau1 + quark1
                                        //neutralino2
					else if (abs(theObjects[theObjects[j].children[k]].type) == 1000023) 
					{
						jtt[njtt].neutralino2.index = theObjects[j].children[k];
						cur = theObjects[j].children[k];
                                                //check children of neutralino2
						for (int l = 0; l < theObjects[cur].children.size(); l++)
						{
							if (abs(theObjects[theObjects[cur].children[l]].type) == 15) //tau1
							{
								jtt[njtt].tau1.index = theObjects[cur].children[l];
							}
							else if (abs(theObjects[theObjects[cur].children[l]].type) == 1000015)//stau1
							{
								jtt[njtt].stau.index = theObjects[cur].children[l];
								cur2 = theObjects[cur].children[l];
                                                                //check children of stau1
								for (int kk = 0; kk < theObjects[cur2].children.size(); kk++)
								{
									if (abs(theObjects[theObjects[cur2].children[kk]].type) == 15)//tau2
									{
										jtt[njtt].tau2.index = theObjects[cur2].children[kk];
									}
                                                                        //neutralino1
									else if (abs(theObjects[theObjects[cur2].children[kk]].type) == 1000022)
									{
										jtt[njtt].neutralino1.index = theObjects[cur2].children[kk];
									}
								}
							}
						}
					}

                                        // squark -> chargino1 + quark1 -> neutralino1 + W + quark1
                                        //        -> neutralino1 + quark2 + quark3 + quark1
                                        //chargino1
					else if (abs(theObjects[theObjects[j].children[k]].type) == 1000024)
					{
						wj[nwj].chargino.index = theObjects[j].children[k];
						cur = theObjects[j].children[k];
                                                //check children of chargino1
						for (int l = 0; l < theObjects[cur].children.size(); l++)
						{
							if (abs(theObjects[theObjects[cur].children[l]].type) == 24)// W
							{
								wj[nwj].w.index = theObjects[cur].children[l];
							}
                                                        //neutralino1
							else if (abs(theObjects[theObjects[cur].children[l]].type) == 1000022)
							{
								wj[nwj].neutralino1.index = theObjects[cur].children[l];
							}
						}
					}
				}//endfor:604
				
			}//endif:599
//                cout<<"!!!debug796 "<<endl;

//========================
                        // stop1 -> top + neutralino1 -> W + b + neutralino1
                        //       -> quark1 + quark2 + b + neutralino1
                        // identify this decay chain
			if (abs(theObjects[j].type) == 1000006 )//if:716 stop1
			{
				wbj[nwbj].stop1.index = j;
				for (int k = 0; k < theObjects[j].children.size(); k++)//for:722 check children
				{
                                        //neutralino1 from stop1, so far consider only 1
					if (abs(theObjects[theObjects[j].children[k]].type) == 1000022)//neutralino1
					{
						wbj[nwbj].neutralino1.index = theObjects[j].children[k];
					}

                                        // stop1 -> top + neutralino1 -> W + b + neutralino1
                                        //       -> quark1 + quark2 + b + neutralino1
                                        //top from stop1
					else if (abs(theObjects[theObjects[j].children[k]].type) == 6)//top
					{
						wbj[nwbj].top.index = theObjects[j].children[k];
						cur = theObjects[j].children[k];
                                                //check children of top
						for (int l = 0; l < theObjects[cur].children.size(); l++)
						{       // W
							if (abs(theObjects[theObjects[cur].children[l]].type) == 24)
							{
								wbj[nwbj].w.index = theObjects[cur].children[l];
							}
                                                        //bottom
							else if (abs(theObjects[theObjects[cur].children[l]].type) == 5)
							{
								wbj[nwbj].bottom.index = theObjects[cur].children[l];
							}
						}
					}
				}//endfor:722
				
			}//endif:716
//                cout<<"!!!debug837 "<<endl;


//========================
                        // top -> W + b -> quark1 + quark2 + b
                        // identify this decay chain
			if (theObjects[j].type == 6 && theObjects[j].status == 3)//if:1430 top
			{
//			theObjects[j].Initialize(1, j, istat[j], itype[j], px[j], py[j], pz[j], p0[j], m[j], -99,
//                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, true);
               		    outFile << "top"<<endl;
			    outFile  <<" "<< theObjects[j].index + 1 << "  "
                                          << theObjects[j].status << "  "
                                          << theObjects[j].type << "  "
                                          << theObjects[j].p0 << "  "
					  << theObjects[j].px << "  "
					  << theObjects[j].py << "  "
					  << theObjects[j].pz << "  "
					  << theObjects[j].mass << endl;
				
			}//endif:1430

			if (theObjects[j].type == -6 && theObjects[j].status == 3)//if:1446 antitop
			{
//			theObjects[j].Initialize(1, j, istat[j], itype[j], px[j], py[j], pz[j], p0[j], m[j], -99,
//                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, true);
               		    outFile << "antitop"<<endl;
			    outFile  <<" "<< theObjects[j].index + 1 << "  "
                                          << theObjects[j].status << "  "
                                          << theObjects[j].type << "  "
                                          << theObjects[j].p0 << "  "
					  << theObjects[j].px << "  "
					  << theObjects[j].py << "  "
					  << theObjects[j].pz << "  "
					  << theObjects[j].mass << endl;
				
			}//endif:1446
//                cout<<"!!!debug1461 "<<endl;

//========================
			// If the above process found a W+jet event:
			if (wj[nwj].squark.index != -1 
			 && wj[nwj].chargino.index != -1 
	       		 && wj[nwj].w.index != -1 
			 && wj[nwj].neutralino1.index != -1 
			 && wj[nwj].nquarks > 0)
			{//if:676
				cur = wj[nwj].squark.index;
				wj[nwj].squark.type.assign("squark");
				wj[nwj].squark.px = theObjects[cur].px;
				wj[nwj].squark.py = theObjects[cur].py;
				wj[nwj].squark.pz = theObjects[cur].pz;
				wj[nwj].squark.p0 = theObjects[cur].p0;
				
				cur = wj[nwj].chargino.index;
				wj[nwj].chargino.type.assign("chargino");
				wj[nwj].chargino.px = theObjects[cur].px;
				wj[nwj].chargino.py = theObjects[cur].py;
				wj[nwj].chargino.pz = theObjects[cur].pz;
				wj[nwj].chargino.p0 = theObjects[cur].p0;
				
				cur = wj[nwj].w.index;
				wj[nwj].w.type.assign("w");
				wj[nwj].w.px = theObjects[cur].px;
				wj[nwj].w.py = theObjects[cur].py;
				wj[nwj].w.pz = theObjects[cur].pz;
				wj[nwj].w.p0 = theObjects[cur].p0;
				
				cur = wj[nwj].neutralino1.index;
				wj[nwj].neutralino1.type.assign("neutralino1");
				wj[nwj].neutralino1.px = theObjects[cur].px;
				wj[nwj].neutralino1.py = theObjects[cur].py;
				wj[nwj].neutralino1.pz = theObjects[cur].pz;
				wj[nwj].neutralino1.p0 = theObjects[cur].p0;
				
				for(int k = 0; k < wj[nwj].nquarks; k++) {
					cur = wj[nwj].quarks[k].index;
					wj[nwj].quarks[k].type.assign("quark");
					wj[nwj].quarks[k].px = theObjects[cur].px;
					wj[nwj].quarks[k].py = theObjects[cur].py;
					wj[nwj].quarks[k].pz = theObjects[cur].pz;
					wj[nwj].quarks[k].p0 = theObjects[cur].p0;
				}
				nwj++;
			}//endif:676
			else 
				wj[nwj].clear();//endif:676
//                cout<<"!!!debug887 "<<endl;

//========================
			// If the above process found a W+bJet event:
			if (wbj[nwbj].stop1.index != -1 
			 && wbj[nwbj].neutralino1.index != -1
			 && wbj[nwbj].top.index != -1 
	       		 && wbj[nwbj].w.index != -1 
			 && wbj[nwbj].bottom.index != -1 )
			{//if:806
				cur = wbj[nwbj].stop1.index;
				wbj[nwbj].stop1.type.assign("stop1");
				wbj[nwbj].stop1.px = theObjects[cur].px;
				wbj[nwbj].stop1.py = theObjects[cur].py;
				wbj[nwbj].stop1.pz = theObjects[cur].pz;
				wbj[nwbj].stop1.p0 = theObjects[cur].p0;

				cur = wbj[nwbj].neutralino1.index;
				wbj[nwbj].neutralino1.type.assign("neutralino1");
				wbj[nwbj].neutralino1.px = theObjects[cur].px;
				wbj[nwbj].neutralino1.py = theObjects[cur].py;
				wbj[nwbj].neutralino1.pz = theObjects[cur].pz;
				wbj[nwbj].neutralino1.p0 = theObjects[cur].p0;

				cur = wbj[nwbj].top.index;
				wbj[nwbj].top.type.assign("top");
				wbj[nwbj].top.px = theObjects[cur].px;
				wbj[nwbj].top.py = theObjects[cur].py;
				wbj[nwbj].top.pz = theObjects[cur].pz;
				wbj[nwbj].top.p0 = theObjects[cur].p0;

				cur = wbj[nwbj].w.index;//chidren of W ??
				wbj[nwbj].w.type.assign("w");
				wbj[nwbj].w.px = theObjects[cur].px;
				wbj[nwbj].w.py = theObjects[cur].py;
				wbj[nwbj].w.pz = theObjects[cur].pz;
				wbj[nwbj].w.p0 = theObjects[cur].p0;
				
				cur = wbj[nwbj].bottom.index;
				wbj[nwbj].bottom.type.assign("bottom");
				wbj[nwbj].bottom.px = theObjects[cur].px;
				wbj[nwbj].bottom.py = theObjects[cur].py;
				wbj[nwbj].bottom.pz = theObjects[cur].pz;
				wbj[nwbj].bottom.p0 = theObjects[cur].p0;

				nwbj++;
			}//endif:806
			else 
				wbj[nwbj].clear();//endif:806
//                cout<<"!!!debug936 "<<endl;

//========================
			// If the above process found a W boson:
			if (justW[njustW].w.index != -1
			 && justW[njustW].nquarks > 0) {//if:724
				// Check that the W is NOT in any W+Jet event.
				newW = true;
				for(int k = 0; k < nwj; k++) if(wj[k].w.index == justW[njustW].w.index)
					newW = false;
				// Check that the W is NOT in any W+bJet event.
				for(int k = 0; k < nwbj; k++) if(wbj[k].w.index == justW[njustW].w.index)
					newW = false;
				
				if(newW) {//W is NOT from any W+Jet event or W+bJet event.
					cur = justW[njustW].w.index;
					justW[njustW].w.type.assign("w");
					justW[njustW].w.px = theObjects[cur].px;
					justW[njustW].w.py = theObjects[cur].py;
					justW[njustW].w.pz = theObjects[cur].pz;
					justW[njustW].w.p0 = theObjects[cur].p0;
					
					for(int k = 0; k < justW[njustW].nquarks; k++) {
						cur = justW[njustW].quarks[k].index;
						justW[njustW].quarks[k].type.assign("quark");
						justW[njustW].quarks[k].px = theObjects[cur].px;
						justW[njustW].quarks[k].py = theObjects[cur].py;
						justW[njustW].quarks[k].pz = theObjects[cur].pz;
						justW[njustW].quarks[k].p0 = theObjects[cur].p0;
					}

					for(int k = 0; k < justW[njustW].nelectrons; k++) {
						cur = justW[njustW].electrons[k].index;
						justW[njustW].electrons[k].type.assign("electron");
						justW[njustW].electrons[k].px = theObjects[cur].px;
						justW[njustW].electrons[k].py = theObjects[cur].py;
						justW[njustW].electrons[k].pz = theObjects[cur].pz;
						justW[njustW].electrons[k].p0 = theObjects[cur].p0;
//                                           numTotalElecFromW++;
					}

					njustW++;
				}
				else
					justW[njustW].clear();
			}//endif:724 
			else
				justW[njustW].clear();//endif:724
//                cout<<"!!!debug937 "<<endl;

//========================	
			// If the above process found a Jet+2tau event:
			isLeptonic = false;  // Still need to make sure taus are hadronic decaying.
			if (jtt[njtt].squark.index != -1 &&
			    jtt[njtt].nquarks > 0 &&
			    jtt[njtt].neutralino2.index != -1 &&
		 	    jtt[njtt].tau1.index != -1 &&
			    jtt[njtt].stau.index != -1 &&
		            jtt[njtt].tau2.index != -1 &&
			    jtt[njtt].neutralino1.index != -1)//if:758
			{
				
				// check tau for final state...  if not, choose the tau that IS final state.
                                //= 1 : an existing entry, which has not decayed or fragmented. This is the main
                                //      class of entries, which represents the ‘ﬁnal state’ given by the generator.
                                //= 2 : an entry which has decayed or fragmented and is therefore not appearing in
                                //      the ﬁnal state, but is retained for event history information.
                                //= 3 : a documentation line, deﬁned separately from the event history. This could
                                //      include the two incoming reacting particles, etc.
                                // require the status number for the tau must be 2
				if(theObjects[jtt[njtt].tau1.index].status != 2)
					notStatus2 = true;
				while(notStatus2) {
                                        //check children of tau1: whether it is a tau
					for(int k=0; k<theObjects[jtt[njtt].tau1.index].children.size(); k++)
						if( abs(theObjects[ theObjects[jtt[njtt].tau1.index].children[k] ].type) == 15 )//tau
						{
							jtt[njtt].tau1.index = theObjects[jtt[njtt].tau1.index].children[k];
							if(theObjects[jtt[njtt].tau1.index].status == 2)
								notStatus2 = false;
						}
				}
				if(theObjects[jtt[njtt].tau2.index].status != 2)
					notStatus2 = true;
				while(notStatus2) {
					for(int k=0; k<theObjects[jtt[njtt].tau2.index].children.size(); k++)
						if( abs(theObjects[ theObjects[jtt[njtt].tau2.index].children[k] ].type) == 15 )//tau
						{
							jtt[njtt].tau2.index = theObjects[jtt[njtt].tau2.index].children[k];
							if(theObjects[jtt[njtt].tau2.index].status == 2)
								notStatus2 = false;
						}
				}
				
				for (int k = 0; k < theObjects[jtt[njtt].tau1.index].children.size(); k++)
				{
					if (abs(theObjects[theObjects[jtt[njtt].tau1.index].children[k]].type) == 11 ||//electron
							abs(theObjects[theObjects[jtt[njtt].tau1.index].children[k]].type) == 13)//muon
					{
						isLeptonic = true;
					}
				} 
				for (int k = 0; k < theObjects[jtt[njtt].tau2.index].children.size(); k++)
				{
					if (abs(theObjects[theObjects[jtt[njtt].tau2.index].children[k]].type) == 11 ||//electron
							abs(theObjects[theObjects[jtt[njtt].tau2.index].children[k]].type) == 13)//muon
					{
						isLeptonic = true;
					}
				}
				
				if (!isLeptonic)//requir both tau1 and tau2 are hadronic decaying tau's
				{
					cur = jtt[njtt].squark.index;
					jtt[njtt].squark.type.assign("squark");
					jtt[njtt].squark.px = theObjects[cur].px;
					jtt[njtt].squark.py = theObjects[cur].py;
					jtt[njtt].squark.pz = theObjects[cur].pz;
					jtt[njtt].squark.p0 = theObjects[cur].p0;
					
					cur = jtt[njtt].neutralino2.index;
					jtt[njtt].neutralino2.type.assign("neutralino2");
					jtt[njtt].neutralino2.px = theObjects[cur].px;
					jtt[njtt].neutralino2.py = theObjects[cur].py;
					jtt[njtt].neutralino2.pz = theObjects[cur].pz;
					jtt[njtt].neutralino2.p0 = theObjects[cur].p0;
					
					cur = jtt[njtt].tau1.index;
					p0vis = 0.;
					pxvis = 0.;
					pyvis = 0.;
					pzvis = 0.;
                                        //check children of tau1
					for (int k = 0; k < theObjects[cur].children.size(); k++)
					{
						if (theObjects[cur].children[k] != -1 &&//not initial value
								abs(theObjects[theObjects[cur].children[k]].type) != 16)//not neutrino_tau
						{
							p0vis += theObjects[theObjects[cur].children[k]].p0;//visible momentum
							pxvis += theObjects[theObjects[cur].children[k]].px;//subtract the neutrino momentum
							pyvis += theObjects[theObjects[cur].children[k]].py;
							pzvis += theObjects[theObjects[cur].children[k]].pz;
						}
					}
					jtt[njtt].tau1.type.assign("tau1");//Question6: matching tau ??
                                        //+1 for tau^+ or -1 for tau^-
					jtt[njtt].tau1.sign = (theObjects[cur].type)/(abs(theObjects[cur].type));
					jtt[njtt].tau1.px = theObjects[cur].px;
					jtt[njtt].tau1.py = theObjects[cur].py;
					jtt[njtt].tau1.pz = theObjects[cur].pz;
					jtt[njtt].tau1.p0 = theObjects[cur].p0;
					jtt[njtt].tau1.pxvis = pxvis;
					jtt[njtt].tau1.pyvis = pyvis;
					jtt[njtt].tau1.pzvis = pzvis;
					jtt[njtt].tau1.p0vis = p0vis;
					
					cur = jtt[njtt].stau.index;
					jtt[njtt].stau.type.assign("stau");
					jtt[njtt].stau.px = theObjects[cur].px;
					jtt[njtt].stau.py = theObjects[cur].py;
					jtt[njtt].stau.pz = theObjects[cur].pz;
					jtt[njtt].stau.p0 = theObjects[cur].p0;
					
					cur = jtt[njtt].tau2.index;
					p0vis = 0.;
					pxvis = 0.;
					pyvis = 0.;
					pzvis = 0.;
                                        //check children of tau2
					for (int k = 0; k < theObjects[cur].children.size(); k++)
					{
						if (theObjects[cur].children[k] != -1 &&//not initial value
						    abs(theObjects[theObjects[cur].children[k]].type) != 16)//not neutrino_tau
						{
							p0vis += theObjects[theObjects[cur].children[k]].p0;//visible momentum
							pxvis += theObjects[theObjects[cur].children[k]].px;//subtract neutrino momentum
							pyvis += theObjects[theObjects[cur].children[k]].py;
							pzvis += theObjects[theObjects[cur].children[k]].pz;
						}
					}
					jtt[njtt].tau2.type.assign("tau2");
                                        //+1 for tau^+ or -1 for tau^-
					jtt[njtt].tau2.sign = (theObjects[cur].type)/(abs(theObjects[cur].type));
					jtt[njtt].tau2.px = theObjects[cur].px;
					jtt[njtt].tau2.py = theObjects[cur].py;
					jtt[njtt].tau2.pz = theObjects[cur].pz;
					jtt[njtt].tau2.p0 = theObjects[cur].p0;
					jtt[njtt].tau2.pxvis = pxvis;
					jtt[njtt].tau2.pyvis = pyvis;
					jtt[njtt].tau2.pzvis = pzvis;
					jtt[njtt].tau2.p0vis = p0vis;
					
					cur = jtt[njtt].neutralino1.index;
					jtt[njtt].neutralino1.type.assign("neutralino1");
					jtt[njtt].neutralino1.px = theObjects[cur].px;
					jtt[njtt].neutralino1.py = theObjects[cur].py;
					jtt[njtt].neutralino1.pz = theObjects[cur].pz;
					jtt[njtt].neutralino1.p0 = theObjects[cur].p0;

                                        //quarks from squark; may be more than 1
					for(int k = 0; k < jtt[njtt].nquarks; k++) {
						cur = jtt[njtt].quarks[k].index;
						jtt[njtt].quarks[k].type.assign("quark");
						jtt[njtt].quarks[k].px = theObjects[cur].px;
						jtt[njtt].quarks[k].py = theObjects[cur].py;
						jtt[njtt].quarks[k].pz = theObjects[cur].pz;
						jtt[njtt].quarks[k].p0 = theObjects[cur].p0;
					}
					njtt++;
				}
				else
					jtt[njtt].clear();
			}//endif:758
			else
				jtt[njtt].clear();//endif:758
		}//endfor:578 end loop over the entire generated particles: hepevt
//                cout<<"!!!debug1141 "<<endl;


//++++++++++++++++++++++++++++++++
//		outFile << "Event  " <<nevhep<<endl;
//		outFile << "MET" << endl;	
//		outFile <<" "<< metcal <<"  "<< phimetcal <<"  "<< metcor <<"  "<< phimetcor <<endl;

//		for (int k = 0; k < nmyLeptons; k++)
//		{
//			outFile  << myLeptons[k].type << endl;//type=="electron" or "muon"
//			outFile  << " " << myLeptons[k].p0 << "  "
//					<< myLeptons[k].px << "  "
//					<< myLeptons[k].py << "  "
//					<< myLeptons[k].pz << "  "
//					<< myLeptons[k].ptiso << "  "
//					<< myLeptons[k].sign << endl;
//		}

//		for (int k = 0; k < nmyTaus; k++)
//		{
//			outFile  << myTaus[k].type << endl;//type=="tau"
//			outFile  << " " << myTaus[k].p0 << "  "
//					<< myTaus[k].px << "  "
//					<< myTaus[k].py << "  "
//					<< myTaus[k].pz << "  "
//					<< myTaus[k].p0vis << "  "
//					<< myTaus[k].pxvis << "  "
//					<< myTaus[k].pyvis << "  "
//					<< myTaus[k].pzvis << "  "
//					<< myTaus[k].sign << endl;
//		}

//		for (int k = 0; k < nmyJets; k++)
//		{       //vecobj(6,iobj) = 21 (g);  2,1,3 (u,d,s);  4 (c) 5 (b)
                        //vecobj(7,iobj) non-zero if loose b tag
                        //vecobj(8,iobj) non-zero if tight b tag
//			outFile  << myJets[k].type << endl;//type=="jet"
//			outFile  <<" "<< myJets[k].p0 << "  "
//				      << myJets[k].px << "  "
//				      << myJets[k].py << "  "
//				      << myJets[k].pz << "  "
//				      << myJets[k].theVec6 << "  "
//				      << myJets[k].theVec7 << "  "
//				      << myJets[k].theVec8 << endl;
//		}	
	 
		outFile << "EndEvent" << endl;
		
		dummyInt = 0;
//		if (i % 100 == 0) outFile << i << endl;
	}

	inFile.close();
	outFile.close();

                  cout<<"###############"<<endl;
                  cout<<"totalEvents  = "<<totalEvents<<endl;
                  cout<<"totalEvtPass = "<<totalEvtPass<<endl;
                  cout<<"###############"<<endl;
                  cout<<"numTotalElectrons      = "<<numTotalElectrons<<endl;
                  cout<<"numTotalMuons          = "<<numTotalMuons<<endl;
                  cout<<"numTotalTaus           = "<<numTotalTaus<<endl<<endl;
                  cout<<"numTotalTrueElectrons  = "<<numTotalTrueElectrons<<endl;
                  cout<<"numTotalFalseElectrons = "<<numTotalFalseElectrons<<endl;
                  cout<<"numTotalTrueMuons      = "<<numTotalTrueMuons<<endl;
                  cout<<"numTotalFalseMuons     = "<<numTotalFalseMuons<<endl<<endl;

                  cout<<"###############"<<endl;
                  cout<<"numElec     = "<<numElec<<endl;
                  cout<<"numTrueElec = "<<numTrueElec<<endl;
                  cout<<"-------"<<endl;
                  cout<<"numFalseElec06et             = "<<numFalseElec06et<<endl;
                  cout<<"numFalseElec67etiso          = "<<numFalseElec67etiso<<endl;
                  cout<<"numFalseElec678ptiso         = "<<numFalseElec678ptiso<<endl;
                  cout<<"numFalseElec67810epmin       = "<<numFalseElec67810epmin<<endl;
                  cout<<"numFalseElec67810epminmax    = "<<numFalseElec67810epminmax<<endl;
                  cout<<"numFalseElec6781011cos_theta = "<<numFalseElec6781011cos_theta<<endl;
                  cout<<"-------"<<endl;
                  cout<<"numFalseElec6et     = "<<numFalseElec6et<<endl;
                  cout<<"numFalseElec7etiso  = "<<numFalseElec7etiso<<endl;
                  cout<<"numFalseElec8ptiso  = "<<numFalseElec8ptiso<<endl;
                  cout<<"numFalseElec10epmin = "<<numFalseElec10epmin<<endl;
                  cout<<"numFalseElec10epmax = "<<numFalseElec10epmax<<endl<<endl;

                  cout<<"###############"<<endl;
                  cout<<"numTotalWs             = "<<numTotalWs<<endl;
                  cout<<"numTotalElectronsFromW = "<<numTotalElectronsFromW<<endl;
                  cout<<"numTotalMuonsFromW     = "<<numTotalMuonsFromW<<endl;
                  cout<<"numTotalTausFromW      = "<<numTotalTausFromW<<endl;
                  cout<<"numTotalQuarksFromW    = "<<numTotalQuarksFromW<<endl;
                  cout<<"numTotalW2Electron     = "<<numTotalW2Electron<<endl;
                  cout<<"numTotalW2Muon         = "<<numTotalW2Muon<<endl;
                  cout<<"numTotalW2Tau          = "<<numTotalW2Tau<<endl;
                  cout<<"numTotalW2Quark        = "<<numTotalW2Quark<<endl<<endl;

                  cout<<"numTotalElecWeta3      = "<<numTotalElecWeta3<<endl;
                  cout<<"numTotalElecWeta3pt0   = "<<numTotalElecWeta3pt0<<endl;
                  cout<<"numTotalElecWeta3pt5   = "<<numTotalElecWeta3pt5<<endl;
                  cout<<"numTotalElecWeta3pt10  = "<<numTotalElecWeta3pt10<<endl;
                  cout<<"numTotalElecWeta3pt15  = "<<numTotalElecWeta3pt15<<endl;
                  cout<<"numTotalElecWeta3pt20  = "<<numTotalElecWeta3pt20<<endl;
                  cout<<"numTotalElecWpt0       = "<<numTotalElecWpt0<<endl;
                  cout<<"numTotalElecWpt5       = "<<numTotalElecWpt5<<endl;
                  cout<<"numTotalElecWpt10      = "<<numTotalElecWpt10<<endl;
                  cout<<"numTotalElecWpt15      = "<<numTotalElecWpt15<<endl;
                  cout<<"numTotalElecWpt20      = "<<numTotalElecWpt20<<endl<<endl;

                  cout<<"###############"<<endl;
                  cout<<"numTotalPgsElecsFromW      = "<<numTotalPgsElecsFromW<<endl;
                  cout<<"numTotalTruePgsElecsFromW  = "<<numTotalTruePgsElecsFromW<<endl;
                  cout<<"numTotalFalsePgsElecsFromW = "<<numTotalFalsePgsElecsFromW<<endl;
                  cout<<"numTotalTruePGSElecsFromW  = "<<numTotalTruePGSElecsFromW<<endl<<endl;

        }//endfor:418 File has ended.  Next File.

	return 0;
}


