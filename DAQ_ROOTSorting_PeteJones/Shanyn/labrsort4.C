/* Advanced sorting code for TDR format data */
/* Produces ROOT histograms */

/* Tested up to ROOT 6.18/00 */

/* Changes to implement 8 detectors and generate event data */
/* based on time window */

/* (c) P.M. Jones */
/*   (13/07/18)   */

/*    Ver. 4.2    */
/*   (02/11/20)   */

#include <stdio.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TH2F.h"
#include <stdint.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <inttypes.h>


#define SIZE  16384

int32_t buffer[SIZE];

uint64_t TStop, TSbot, data;
uint64_t TS, TSlast, TSfirst, counter,  SYNC, SYNClast=0;
uint64_t TimeStamp1[100]={0}, TimeStamp2[100]={0}, TimeStamp3[100]={0} ,TimeStamp4[100]={0}, TimeStamp5[100]={0}, TimeStamp6[100]={0}, TimeStamp7[100]={0}, TimeStamp8[100]={0};
uint16_t Energy[100]={0}, energy;
uint16_t QDC[100]={0}, qdc;
uint64_t TS2;

double T1, T2;

int64_t TSdiff, SYNCdiff;
int64_t count=0;

unsigned long int pos = 0;

int verbose = 255;

int main (int argc, char** argv)
{
  FILE *f;
  int i, j, k;
  int blocks_in=0;
  int twidset=0;
  char file_in[50];
  char file_out[50];
  int ident, card=0, adcdata, detectorID;
  int tick, tock;
  double sick;
  int RunEnd=999999, Run=0;
  
  int m1,m2;
  int sum1, sum2;
  
  // Data arrays for events

  int energyF[9] = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
  int energyS[9] = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
  double Time[9] = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
  double Td;
  

  //Energy cablibrations
  double a[9], b[9], c[9], d[9];

  a[1] =  2.2882810120290898e-11;
  b[1] =  -2.8226334147027305e-07;
  c[1] =  0.13217769618776284;
  d[1] =  -2.665958231711416; 
  
  a[2] =  -2.6206302270812858e-11;
  b[2] =  5.20508451009155e-07;
  c[2] =  0.12454191180950766;
  d[2] =  -2.8712415583274207;
  
  a[3] =  -5.8179826116966215e-12;
  b[3] =  2.882971831716912e-07;
  c[3] =  0.1319119871835276;
  d[3] =  1.7809131369665836;
  
  a[4] =  -1.9907752086200708e-12;
  b[4] =  1.138740322817261e-07;
  c[4] =  0.10055457304852329;
  d[4] =  0.09031549663462841;
  
  a[5] =  -5.61935578847695e-12;
  b[5] =  2.1772648634078458e-07;
  c[5] =  0.12846705019381832;
  d[5] =  -0.8463216578107895; 
  
  a[6] =  -3.9525192061812874e-12;
  b[6] =  1.992195829282599e-07;
  c[6] =  0.11108313431790826;
  d[6] =  0.6410335417880054; 
  
  a[7] =  -2.0089658591995725e-12;
  b[7] =  1.5796558109971302e-07;
  c[7] =  0.1158528534483307;
  d[7] =  1.0762093304010865;
  
  a[8] =  -1.193427933259776e-11;
  b[8] =  3.481567315783377e-07;
  c[8] =  0.13492630347249404;
  d[8] =  0.8299451107974863; 


  //Ttme calibrations
  
  double Ta[9];

  Ta[1]=   0.00;
  Ta[2]= -71.60;
  Ta[3]=  32.04;
  Ta[4]= -41.40;
  Ta[5]=  48.42;
  Ta[6]=-114.35;
  Ta[7]= 105.52;
  Ta[8]= -84.83;

  
  int no_read;
  std::vector<int> time_stamp;
  
  
  if ((argc < 2) || (argc > 3))
    {
	fprintf(stderr, "Usage: %s data [last_subrun]\n", argv[0]);
	exit(1);
    }

  /* ROOT STUFF */

  sprintf(file_out, "%s.root", argv[1]);
  if (argc > 2) sscanf(argv[2], "%d", &RunEnd);

  TFile *g = new TFile(file_out,"recreate");

  printf("ROOT file %s opened...\n", file_out);

  //   TTree *T = new TTree("T","Labr tree");
  //   T->Branch("TS",&TS,"TS/l");
  //   T->Branch("labre1",&labre1,"labre1/I");
  //   T->Branch("labre2",&labre2,"labre2/I");
  //   T->Branch("labrt1",&labrt1,"labrt1/I");
  //   T->Branch("labrt2",&labrt2,"labrt2/I");

  TH1D** h1=new TH1D*[9];
  TH1D** h2=new TH1D*[9];
  TH1D** c1=new TH1D*[9];
  TH1D** c2=new TH1D*[9];

  TH1D** k1=new TH1D*[9];
  TH1D** k2=new TH1D*[9];
  TH1D** k3=new TH1D*[9];
  TH1D** k4=new TH1D*[9];
  TH1D** k5=new TH1D*[9];
  TH1D** k6=new TH1D*[9];
  TH1D** k7=new TH1D*[9];
  TH1D** k8=new TH1D*[9];

  TH1D** n1=new TH1D*[9];


  for (j = 1; j <=8; j++)
    {
      h1[j] = new TH1D(TString::Format("Fast_%02d", j),"Spectrum",2048,0,2047);
      h2[j] = new TH1D(TString::Format("Slow_%02d", j),"Spectrum",2048,0,2047);

      c1[j] = new TH1D(TString::Format("Coinc1_%02d", j),"Spectrum",2048,0,2047);
      c2[j] = new TH1D(TString::Format("Coinc2_%02d", j),"Spectrum",2048,0,2047);

      k1[j] = new TH1D(TString::Format("Time_01_%02d", j),"Spectrum",16384,0,16383);
      k2[j] = new TH1D(TString::Format("Time_02_%02d", j),"Spectrum",16384,0,16383);
      k3[j] = new TH1D(TString::Format("Time_03_%02d", j),"Spectrum",16384,0,16383);
      k4[j] = new TH1D(TString::Format("Time_04_%02d", j),"Spectrum",16384,0,16383);
      k5[j] = new TH1D(TString::Format("Time_05_%02d", j),"Spectrum",16384,0,16383);
      k6[j] = new TH1D(TString::Format("Time_06_%02d", j),"Spectrum",16384,0,16383);
      k7[j] = new TH1D(TString::Format("Time_07_%02d", j),"Spectrum",16384,0,16383);
      k8[j] = new TH1D(TString::Format("Time_08_%02d", j),"Spectrum",16384,0,16383);

      n1[j] = new TH1D(TString::Format("Fast_Multiplicity_%02d", j),"Spectrum",2048,0,2047);

    }
  
   TH1D *mul1 = new TH1D("Multiplicity 1","spectrum",16,0,15);
   TH1D *mul2 = new TH1D("Multiplicity 2","spectrum",16,0,15);
   TH1D *mul3 = new TH1D("Multiplicity 3","spectrum",16,0,15);

   TH1D *hit1 = new TH1D("Hitpattern 1","spectrum",16,0,15);
   TH1D *hit2 = new TH1D("Hitpattern 2","spectrum",16,0,15);
   TH1D *hit3 = new TH1D("Hitpattern 3","spectrum",16,0,15);

   TH1D *sumspec1 = new TH1D("Sum 1","Spectrum",65536,0,65535);
   TH1D *sumspec11 = new TH1D("Sum 11","Spectrum",65536,0,65535);
   TH1D *sumspec12 = new TH1D("Sum 12","Spectrum",65536,0,65535);
   TH1D *sumspec13 = new TH1D("Sum 13","Spectrum",65536,0,65535);
   TH1D *sumspec14 = new TH1D("Sum 14","Spectrum",65536,0,65535);
   TH1D *sumspec2 = new TH1D("Sum 2","Spectrum",65536,0,65535);


   
   //File stuff

   while(Run < RunEnd+1)

     {
       sprintf(file_in, "./%s_%d", argv[1], Run);
       
       f = fopen(file_in, "rb");
       if (!f)
	 {
	   fprintf(stderr, "Can't open file '%s'\n", file_in);
	   goto finish;
	 }
       
       printf("File %s opened...\n", file_in);

       
       
       while (!feof(f))
	 
    {

  no_read = fread(buffer, sizeof(buffer[0]), SIZE, f);
  blocks_in++;
  
  if ( (feof(f)) && no_read <= 0 ) goto end;

  pos+=24;

      for (i=6; i< SIZE; i+=2)
	{
	  data = buffer[i+1];
	  TSbot = buffer[i];

	  if ((TSbot == 0) && (data == 0)) goto loop;
	  if ((TSbot == 0x5e5e5e5e) && (data == 0x5e5e5e5e)) goto loop;
	  if ((TSbot == 0x5e5e5e5e) && (data == 0xffffffff)) goto loop;

	  /* DATA */

	  if ((data & 0xc0000000) == 0xc0000000)
	    {
	      ident = (data & 0x0fff0000) >> 16;
	      adcdata = (data & 0x0000ffff);
	      card = (ident / 32);

	      if (card == 99) goto skip;

	      TS = (TS & 0x0000fffff0000000ULL);
	      TS = ((TS | TSbot));
		  if (counter == 0)
		  {
			  TSfirst = TS;
			  counter = 1;
		  }
	      TSdiff = TS-TSlast;

	      if ( (TSbot >= 0x0) && (TSbot <= 0x1a) )
		{
		  if (twidset == 0) 
		    {
		      TS=TS+0x10000000;
		      twidset=1;
		    }

		  TSdiff = TS-TSlast;
		  std::cout << " TSdiff: " << TSdiff << " TS: " << TS << " TSlast: " << TSlast << std::endl;
		}



	      /* MORE ROOT STUFF */
	      energy = (adcdata);

	      //	      printf("yyy %llx %lli %lli\n", TS, TSdiff, TSlast);

	      //Define window here in 10's of NS

	      if (TSdiff > 5)
		{
		  //		  printf("\n");
		  T1=0;
		  T2=0;

		  //		  printf("Event:\n");

		  for (j = 1; j <= 8; j++)
		   { 
		     //		     if (Time[j] != -1) printf("Detector %i --> %08d %08d %g\n", j, energyS[j], energyF[j], Time[j]);
		     

		     //Singles spectra

		     if (Time[j] > 0)
		       {
			 h1[j]->Fill(energyF[j]);
			 h2[j]->Fill(energyS[j]);

			 hit1->Fill(j);
		       }


		     //Time difference spectra
		     for (k = 1; k <= 8; k++)
		       { 
			 if ( (Time[j] > 0) && (Time[k] > 0) )
			   {
			     if (j != k )
			       {
				 Td=Time[j]-Time[k]+5000.0;
				 //std::cout << Time[1] << std::endl;
				 Td=Td+(Ta[j]-Ta[k]);
				 if (j == 1) k1[k]->Fill(Td);
				 if (j == 2) k2[k]->Fill(Td);
				 if (j == 3) k3[k]->Fill(Td);
				 if (j == 4) k4[k]->Fill(Td);
				 if (j == 5) k5[k]->Fill(Td);
				 if (j == 6) k6[k]->Fill(Td);
				 if (j == 7) k7[k]->Fill(Td);
				 if (j == 8) k8[k]->Fill(Td);
			       }
			   }
		       }

		   }


		     // Data analysis here

		     m1=0;
		     m2=0;
		     sum1=0;
		     sum2=0;
		     
		     for (j = 1; j <=8; j++)
		       { 
			 for (k = j+1; k <=8; k++)
			   { 
			     Td=Time[j]-Time[k]+8192.0;
			     Td=Td+(Ta[j]-Ta[k]);

				 			     
			     if ( (Td > 8182) && (Td < 8202) )
			       {
				 c1[j]->Fill(energyF[j]);
				 c2[j]->Fill(energyF[k]);

				 hit2->Fill(j);
				 hit2->Fill(k);

				 m1=m1+1;

				 sum1=sum1+energyF[j]+energyF[k];
				 
				 if ((energyF[j] > 471) && (energyF[j] < 551))
				   {
				     if ((energyF[k] > 471) && (energyF[k] < 551))
				       {
					 hit3->Fill(j);
					 hit3->Fill(k);

					 m2=m2+1;
					 sum2=sum2+energyF[j]+energyF[k];
				       }
				   }
				 

			       }
			   }
		       }

		     mul1->Fill(m1);
		     mul2->Fill(m2);

		     if ((sum1 > 942) && (sum1 < 1102))
		       {
			 mul3->Fill(m1);
		       }

		     
		     if (sum1 > 0) sumspec1->Fill(sum1);
		     if (sum2 > 0) sumspec2->Fill(sum2);

		     if (m1 == 1)  sumspec11->Fill(sum1);
		     if (m1 == 2)  sumspec12->Fill(sum1);
		     if (m1 == 3)  sumspec13->Fill(sum1);
		     if (m1 == 4)  sumspec14->Fill(sum1);

		     
		     for (j = 1; j <=8; j++)
		       { 
			 if ((m1>0) && (m1<9)) n1[m1]->Fill(energyF[j]);
		       }
		     

		  //Reinitialise arrays

		  for (j = 1; j <=8; j++)
		   { 
		     energyS[j] = -1;
		     energyF[j] = -1;
		     Time[j] = -1;
		   }
		}



	      //Ident 64,65,66,67,68,69,70,71 are slow signals 
	      if ((ident > 63) && (ident < 68))
		{
		  detectorID = ident-63;
		  energy = (a[detectorID]*adcdata*adcdata*adcdata)+(b[detectorID]*adcdata*adcdata)+(c[detectorID]*adcdata)+d[detectorID];
		  if (energyS[detectorID] == -1) energyS[detectorID] = energy;
		}
	      //Ident 72,73,74,75,76,77,78,79 are fast signals
	      if ((ident > 71) && (ident < 76))
		{
		  detectorID = ident-71;
		  energy = (a[detectorID]*adcdata*adcdata*adcdata)+(b[detectorID]*adcdata*adcdata)+(c[detectorID]*adcdata)+d[detectorID];
		  if (energyF[detectorID] == -1) energyF[detectorID] = energy;
		}
	      //Ident 80,81,82,83,84,85,86,87 are timing signals of slow signals (not used)

      	      //Ident 88,89,90,91,92,93,94,95 are timing signals of fast signals

	      if ((ident > 87) && (ident < 92))
		{
		  detectorID = ident-87;
		  
		  tock = ((adcdata & 0x0000e000) >> 13);
		  tick = (adcdata & 0x00001fff);

		  sick = (tick/8192.0) * 20.0;

		  if (tock != 7)
		    {
		      if (Time[detectorID] == -1) Time[detectorID] = ((TS * 100.0) + (tock * 20.0) + sick);
		    }

		}
	      TSlast = TS;
	    }

	  /* SYNC */
	  
	  if ((data & 0xc0f00000) == 0x80400000)
	    {
	      card = (data & 0x3f000000) >> 24;
	      TStop = (data & 0x000fffff);
	      TS = (TStop);
	      TS = ((TS << 28));
	      TS = ((TS | TSbot));

	      TSdiff = TS-TSlast;

	      SYNC = TS;
              SYNCdiff = SYNC-SYNClast;

	      twidset=0;

	      SYNClast = SYNC;

	    }

	skip:

	  count++;

	  //	  TSlast = TS;

	loop:
	  pos+=8;
	  
	}
      
    }

 end:  

       fclose(f);
       
       Run++;
       
     }


 
  
 finish:

 
  printf("The first time stamp is (e-8 s): %" PRIu64 "\n", TSfirst);
  printf("The last time stamp is (e-8 s): %" PRIu64 "\n", TS);
  printf("The difference is (e-8 s): %" PRIu64 "\n", TS-TSfirst);
  

  for (j = 1; j <=8; j++)
    {
      h1[j]->Write();
      h2[j]->Write();

      c1[j]->Write();
      c2[j]->Write();

      k1[j]->Write();
      k2[j]->Write();
      k3[j]->Write();
      k4[j]->Write();
      k5[j]->Write();
      k6[j]->Write();
      k7[j]->Write();
      k8[j]->Write();

      n1[j]->Write();

    }
  
  mul1->Write();
  mul2->Write();
  mul3->Write();
  hit1->Write();
  hit2->Write();
 

  sumspec1->Write();
  sumspec11->Write();
  sumspec12->Write();
  sumspec13->Write();
  sumspec14->Write();
  sumspec2->Write();

  
delete g;

}
