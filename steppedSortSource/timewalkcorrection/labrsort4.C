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
#include <cmath>
#include <complex>

#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH1D.h"
#include <stdint.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <inttypes.h>
#include "TRandom3.h"
#include "Math/WrappedTF1.h"
#include "Math/Integrator.h"
#include "TF1.h"
#include "TCanvas.h"

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

double rndUnif;
TRandom3* randy = new TRandom3();

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

  int energyF[2] = {-1,-1};
  int energyS[2] = {-1,-1};
  double Time[2] = {-1,-1};
  double TimeS[2] = {-1,-1};
  double Td;
  
    //Energy cablibrations
  double as[2], bs[2], cs[2], ds[2], af[2], bf[2], cf[2], df[2];
    // LaBr slow 0,1 (Labelled 8) R2
  /* as[0] =  0.000002190000;
  bs[0] =  0.433381500000;
  cs[0] =  0.000000000000;

    // LaBr slow ch2,3 L2  (Labelled 3) R2
  as[1] =  0.000001726752;
  bs[1] =  0.430148300000;
  cs[1] =  0.000000000000;

    // LaBr fast ch8 L0 (Labelled 8) R2
    af[0] =  0.000000105753;
    bf[0] =  0.180526000000;
    cf[0] =  0.000000000000;

    // LaBr fast ch9 L2 (Labelled 9) R2
    af[1] =  0.000000224606;
    bf[1] =  0.186934400000;
    cf[1] =  0.000000000000; */

	// LaBr slow ch0 L0 (Labelled 8) R12
	as[0] = 0;
	bs[0] = 0.2900635;
	cs[0] = 0;

	// LaBr slow ch2 L2  (Labelled 3) R12
	as[1] = 0;
	bs[1] = 0.284559400000;
	cs[1] = 0;

	// LaBr fast ch8 L0 (Labelled 8) R12
	af[0] = 0;
	bf[0] = 0.184797800000;
	cf[0] = 0;

	// LaBr fast ch9 L2 (Labelled 3) R12
	af[1] = 0;
	bf[1] = 0.190535700000;
	cf[1] = 0;
	


  //Ttme calibrations
  
  double Ta[9];

  Ta[0]=   0.00;
  Ta[1]=  -11.39;


  
  int no_read;
  std::vector<int> time_stamp;
  
  
  if ((argc < 2) || (argc > 3))
    {
	fprintf(stderr, "Usage: %s data [last_subrun]\n", argv[0]);
	exit(1);
    }

  /* ROOT STUFF */

  sprintf(file_out, "%s_labrsort4.root", argv[1]);
  if (argc > 2) sscanf(argv[2], "%d", &RunEnd);

  TFile *g = new TFile(file_out,"recreate");

  printf("ROOT file %s opened...\n", file_out);

  //   TTree *T = new TTree("T","Labr tree");
  //   T->Branch("TS",&TS,"TS/l");
  //   T->Branch("labre1",&labre1,"labre1/I");
  //   T->Branch("labre2",&labre2,"labre2/I");
  //   T->Branch("labrt1",&labrt1,"labrt1/I");
  //   T->Branch("labrt2",&labrt2,"labrt2/I");

  TH1D** h1=new TH1D*[2];
  TH1D** h2=new TH1D*[2];
  TH1D** c1=new TH1D*[2];
  TH1D** c2=new TH1D*[2];

  TH1D** kf=new TH1D*[2];
  TH1D** ks=new TH1D*[2];


  TH1D** n1=new TH1D*[2];

  TH1D** t1=new TH1D*[2];

  TH2D** t1e = new TH2D*[2];


  for (j = 0; j <2; j++)
    {
      h1[j] = new TH1D(TString::Format("Fast_%02d", j),"Spectrum",2000,0,2000);
      h2[j] = new TH1D(TString::Format("Slow_%02d", j),"Spectrum",2000,0,2000);

	  t1[j] = new TH1D(TString::Format("Time_Slow-Fast_%02d", j),"Spectrum",400 ,-200,200);
	  t1e[j] = new TH2D(TString::Format("Time_Slow-Fast_vs_Energy_%02d", j),"Spectrum",400 ,-200,200, 2000,0,2000);

      c1[j] = new TH1D(TString::Format("Coinc1_%02d", j),"Spectrum",2048,0,2047);
      c2[j] = new TH1D(TString::Format("Coinc2_%02d", j),"Spectrum",2048,0,2047);

      kf[j] = new TH1D(TString::Format("Time_Fast_00_%02d", j),"Spectrum",2000,-100,100);
	  ks[j] = new TH1D(TString::Format("Time_Slow_00_%02d", j),"Spectrum",20000,-1000,1000);

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

   TH1D *timediffcoinc = new TH1D("timediffcoinc","Spectrum",2000,-100,100);

   TH2D * energyL0L1 = new TH2D("Energy L0 vs L1","Spectrum",2000,0,2000,2000,0,2000);

   TH2D*  energygatedbackgroundL0L1 = new TH2D("Energy background gated L0 vs L1","Spectrum",2000,0,2000,2000,0,2000);
   TH2D * energygatedpeakL0L1 = new TH2D("Energy Coincidence peak L0 vs L1","Spectrum",2000,0,2000,2000,0,2000);

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

	      if (TSdiff > 50)
		{

			// Print
		  //		  printf("\n");
		  T1=0;
		  T2=0;
		  energyL0L1->Fill(energyF[0],energyF[1]);

		  //		  printf("Event:\n");

		  for (j = 0; j < 2; j++)
		   { 
		     //if (Time[j] != -1) printf("Detector %i --> %08d %08d %g\n", j, energyS[j], energyF[j], Time[j]); 
		     //Singles spectra
		     if ((Time[j]> 0) && (TimeS[j] > 0))
		       {
					h1[j]->Fill(energyF[j]);
					h2[j]->Fill(energyS[j]);

					hit1->Fill(j);
				
					t1[j]->Fill((TimeS[j]-Time[j])/10);
					t1e[j]->Fill((TimeS[j]-Time[j])/10,energyS[j]);
			  }

			   //Time difference spectra
		     for (k = 1; k < 2; k++)
		    { 
				if ( (Time[j] > 0) && (Time[k] > 0) )
				{
					if (j != k )
					{
						Td = 0;
						Td=(Time[j]-Time[k])/10;
						Td=Td+(Ta[j]-Ta[k]);
						if (j == 0) kf[k]->Fill(Td);
						//printf("Time difference %i %i %g\n", j, k, Td);
					}
				}
				if ( (TimeS[j] > 0) && (TimeS[k] > 0) )
				{
					if (j != k )
					{
						Td = 0;
						Td=(TimeS[j]-TimeS[k])/10;
						//Td=Td+(Ta[j]-Ta[k]);
						if (j == 0) ks[k]->Fill(Td);

						//printf("Time difference %i %i %g\n", j, k, Td);
					}
				}
		    }

			  

		   }

			// Data analysis here
			m1=0;
			m2=0;
			sum1=0;
			sum2=0;
		     
			for (j = 0; j <2; j++)
			{ 
				for (k = j+1; k <2; k++)
				{ 
					if ((Time[j] > 0) && (Time[k] > 0)) 
					{
						Td=(Time[j]-Time[k])/10;
						Td=Td+(Ta[j]-Ta[k]);	

						//printf("Time difference %i %i %g\n", j, k, Td);								
						//if ( (Td > 8142) && (Td < 8242) ) // 20 bins so 2 ns on either side
						if ((Td< -2) || (Td > 2))
						{
							energygatedbackgroundL0L1->Fill(energyF[0],energyF[1]);
						}
						if ((Td >= -2) && (Td <= 2))
						{
							energygatedpeakL0L1->Fill(energyF[0],energyF[1]);
							c1[j]->Fill(energyF[j]);
							c2[j]->Fill(energyF[k]);

							hit2->Fill(j);
							hit2->Fill(k);

							m1=m1+1;

							sum1=sum1+energyF[j]+energyF[k];
							
							//if ((energyF[j] > 1133) && (energyF[j] < 1213))// 40 bins on either side of 1173
							if ((energyF[j] > 471) && (energyF[j] < 551))
							{
								//if ((energyF[k] > 1292) && (energyF[k] < 1372)) // 40 bins on either side of 1332
								if ((energyF[k] > 1074) && (energyF[k] < 1474))
								{
									hit3->Fill(j);
									hit3->Fill(k);

									timediffcoinc->Fill(Td);
									//printf("Time difference %i %i %g\n", j, k, Td);

									m2=m2+1;
									sum2=sum2+energyF[j]+energyF[k];
								}
							}
							//if ((energyF[j] > 1292) && (energyF[j] < 1372))// 40 bins on either side of 1173
							if ((energyF[j] > 1074) && (energyF[j] < 1474))
							{
								//if ((energyF[k] > 1133) && (energyF[k] < 1213)) // 40 bins on either side of 1332
								if ((energyF[k] > 471) && (energyF[k] < 551))
								{
									hit3->Fill(j);
									hit3->Fill(k);

									timediffcoinc->Fill(Td);
									//printf("Time difference %i %i %g\n", j, k, Td);


									m2=m2+1;
									sum2=sum2+energyF[j]+energyF[k];
								}
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

		     
		    for (j = 0; j <2; j++)
			{ 
				if ((m1>=0) && (m1<2)) n1[m1]->Fill(energyF[j]);
			}
		     

		  //Reinitialise arrays

		  for (j = 0; j <2; j++)
		   { 
				energyS[j] = -1;
				energyF[j] = -1;
				Time[j] = -1;
				TimeS[j] = -1;
		   }
		}
		// LaBr3 Detector 8 (L0) and LaBr3 Detector 3 (L2) but for convenience labelled 0 and 1
		//if (ident < 100) printf("Ident: %i, Energy: %i\n", ident, energy);
	    if ((ident == 64) || (ident == 65)) //slow energy
		{
		  detectorID = ident-64;
		  energy = 0;
		  rndUnif = 0;
		  rndUnif = randy->Uniform(-2.8, 2.8);
		  energy = (as[detectorID]*adcdata*adcdata)+(bs[detectorID]*adcdata)+(cs[detectorID]);
		  if (energyS[detectorID] == -1) energyS[detectorID] = energy + rndUnif;
		}
		else if ((ident == 80) || (ident == 81)) //slow time
        { 
			detectorID = ident-80;
			if (TimeS[detectorID] == -1) TimeS[detectorID] = TS*100;
		}
	    else if ((ident == 72) || (ident == 73)) //fast energy
		{
		  detectorID = ident-72;
		  energy = 0;
		  rndUnif = 0;
		  rndUnif = randy->Uniform(-5.4, 5.4);
		  energy = (af[detectorID]*adcdata*adcdata)+(bf[detectorID]*adcdata)+(cf[detectorID]);
		  if (energyF[detectorID] == -1) energyF[detectorID] = energy + rndUnif;
		}
	    else if ((ident == 88) || (ident == 89)) //fast time
		{
		  detectorID = ident-88;
		  
		  tock = ((adcdata & 0x0000e000) >> 13); // TS raw clock time 
		  tick = (adcdata & 0x00001fff); // 500 MHz sampling speed. Time difference between 2ns samples bits 0-12. 0 is start of tick, 8191 is start of next tick (2,4,6,8,10)
		  sick = (tick/8192.0) * 20.0; // between the 2 ns samples

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

	double duration = (TS-TSfirst)*1.0e-8/60.0;
 
  printf("The first time stamp is (e-8 s): %" PRIu64 "\n", TSfirst);
  printf("The last time stamp is (e-8 s): %" PRIu64 "\n", TS);
  printf("The difference is (min): %f\n", duration);
  
	timediffcoinc->Write();
  for (j = 0; j < 2; j++)
    {
      h1[j]->Write();
      h2[j]->Write();

      c1[j]->Write();
      c2[j]->Write();

      kf[j]->Write();
	  ks[j]->Write();

      n1[j]->Write();

	  t1[j]->Write();
	  t1e[j]->Write();

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

  

  energyL0L1->Write();
  energygatedbackgroundL0L1->Write();
  energygatedpeakL0L1->Write();

g->Close();

  
delete g;

}
