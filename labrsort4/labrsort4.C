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
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
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
  int i, j, k, w;
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

  int energyF[5] = {-1,-1,-1,-1,-1};
  int energyS[5] = {-1,-1,-1,-1,-1};
  double Time[5] = {-1,-1,-1,-1,-1};
  double Td, par0[4], par1[4];

  
  double TD15 = {};
  double calibenergy1 = {};

  //Energy cablibrations
  double a[4], b[4], c[4], d[4];
  a[1] = 0.75263;
  b[1] = 0.71741;
  a[2] = 0.81401;
  b[2] = -0.54088;
  a[3] = 0.51735;
  b[3] = 0.41531;
  a[4] = 0.71452;
  b[4] = 1.14457;

  double bw = 1, energycalib, rndchannel, rndUnif;
  TRandom* randy = new TRandom();

  //Ttme calibrations

  double Ta[5];

  Ta[1]= 0.00;
  Ta[2]= 120.00;
  Ta[3]= 150.00;
  Ta[4]= 80.00;
  Ta[5]= 0.00;
  
  
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


  TH1D** h1=new TH1D*[5];
  TH1D** h2=new TH1D*[5];
  TH1D** h4 = new TH1D*[5];
  TH1D** t = new TH1D*[5];
  TH1D** c1=new TH1D*[5];
  TH1D** c2=new TH1D*[5];
  TH1D** k1=new TH1D*[5];
  TH1D** k2=new TH1D*[5];
  TH1D** k3=new TH1D*[5];
  TH1D** k4=new TH1D*[5];
  TH1D** k5=new TH1D*[5];
  TH1D** n1=new TH1D*[5];
  TH2D** th1l = new TH2D*[5];
  TH1D** InSync = new TH1D*[5];
  TH1D** OutSync = new TH1D*[5];
  TH1D** FullSpec = new TH1D*[5];
  TH2D** OutSync2D = new TH2D*[5];
  TH2D** FullSpec2D = new TH2D*[5];
  TH2D** InSync2D = new TH2D*[5];
  TH1D** OutSyncTd = new TH1D*[5];
  TH1D** FullSpecTd = new TH1D*[5];
  TH1D** InSyncTd = new TH1D*[5];

   long int mintime = 5321429603046705;
   long int maxtime = 5321429603052705;

  	for (j = 1; j <=5; j++)
    {
      h1[j] = new TH1D(TString::Format("Fast_Energy_%02d", j),"Fast Signal Energy Spectrum",8000,0,8000);
      h2[j] = new TH1D(TString::Format("Slow_Energy_%02d", j),"Slow Signal Energy Spectrum",8000,0,8000);

	  t[j] = new TH1D(TString::Format("Fast_Time_%02d", j),"Fast Signal Time Spectrum",1000,5333864072647,5333864062647);

      h4[j] = new TH1D(TString::Format("Calibrated_Slow_Energy_%02d", j),"1^{st} Order Calibrated Slow Signal Energy Spectrum",8000,0,8000);

      c1[j] = new TH1D(TString::Format("Coinc1_%02d", j),"Coincidence Spectrum",2048,0,2047);
      c2[j] = new TH1D(TString::Format("Coinc2_%02d", j),"Coincidence Spectrum",2048,0,2047);

      k1[j] = new TH1D(TString::Format("Time_01_%02d", j),"Time Difference Spectrum",8000,0,8000); //16383 
      k2[j] = new TH1D(TString::Format("Time_02_%02d", j),"Time Difference Spectrum",8000,0,8000); 
      k3[j] = new TH1D(TString::Format("Time_03_%02d", j),"Time Difference Spectrum",8000,0,8000);
      k4[j] = new TH1D(TString::Format("Time_04_%02d", j),"Time Difference Spectrum",8000,0,8000);
	  k5[j] = new TH1D(TString::Format("Time_RF_%02d", j),"Time Difference Spectrum",8000,0,8000);

      n1[j] = new TH1D(TString::Format("Fast_Multiplicity_%02d", j),"Spectrum",2048,0,2047);

   	  th1l[j] = new TH2D(TString::Format("Time_L%2d_Vs_Energy", j), TString::Format("Time L%2d - Vs Energy Spectrum", j),1000,5333864072647,5333864062647,8000,0,8000);
	
	  OutSync2D[j] = new TH2D(TString::Format("OutSync2D%2d", j),TString::Format("Out Sync Spectrum L%2d", j), 1000,0,1000, 8000,0,8000);
	  FullSpec2D[j] = new TH2D(TString::Format("FullSpec2D%2d", j),TString::Format("Full Spectrum L%2d", j), 1000,0,1000, 8000,0,8000);
	  InSync2D[j] = new TH2D(TString::Format("InSync2D%2d", j),TString::Format("In Sync Spectrum L%2d", j), 1000,0,1000, 8000,0,8000);

	  OutSyncTd[j] = new TH1D(TString::Format("OutSyncTd%2d", j),TString::Format("Out Sync Spectrum L%2d", j), 1000,0,1000);
	  FullSpecTd[j] = new TH1D(TString::Format("FullSpecTd%2d", j),TString::Format("Full Spectrum L%2d", j), 1000,0,1000);
	  InSyncTd[j] = new TH1D(TString::Format("InSyncTd%2d", j),TString::Format("In Sync Spectrum L%2d", j), 1000,0,1000);

	  InSync[j] = new TH1D(TString::Format("InSync%2d", j),TString::Format("In Sync Spectrum L%2d", j), 8000,0,8000);
	  OutSync[j] = new TH1D(TString::Format("OutSync%2d", j),TString::Format("Out Sync Spectrum L%2d", j), 8000,0,8000);
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

	      			TS = (TS & 0x0000fffff0000000ULL); // sets TS as 20 bits long (0x0000fffff0000000ULL)
	      			TS = ((TS | TSbot)); // sets TS|TSbot meaning TS is now 40 bits long
		  			
					if (counter == 0)
		  			{
			  			TSfirst = TS;
			  			counter = 1;
		  			}
	      			TSdiff = TS-TSlast; // TSlast is the last TS value from the previous loop

	      			if ( (TSbot >= 0x0) && (TSbot <= 0x1a) ) // TSbot needs to be greater/equal to the first TS value and less than/equal to 0x1a which is the last TS value
					{
		  				if (twidset == 0) // if time window is set to 0
		    			{
		      				TS=TS+0x10000000; // TS is defined as 40 bits long, so TS+0x10000000 is 0x10000000 added to the 40 bit TS value
		      				twidset=1; // time window is set to 1
		    			}

		  				TSdiff = TS-TSlast;

						//printf("TSDifference = %lu\n", TSdiff);
					}



	      			/* MORE ROOT STUFF */
	      			energy = (adcdata);

	      			// printf("yyy %llx %lli %lli\n", TS, TSdiff, TSlast);

	      			//Define window here in 10's of NS
					

	      			if (TSdiff >= 9) 
					{
						//std::cout  << "TSdiff: " << TSdiff<< "| TS: " << TS << "| TSlast: " << TSlast << std::endl;
		  				
		  				// printf("\n");
		 				T1=0;
		  				T2=0;

		  				// printf("Event:\n");

		  				for (j = 1; j <= 5; j++)
		   				{ 
		     				// if (Time[j] != -1) printf("Detector %i --> %08d %08d %g\n", j, energyS[j], energyF[j], Time[j]);

		    				//Singles spectra

		     				if (Time[j] > 0)
		       				{
			 					h1[j]->Fill(energyF[j]);
			 					h2[j]->Fill(energyS[j]);

			 					hit1->Fill(j);
		       				}


		     				//Time difference spectra
		     				for (k = 1; k <= 5; k++)
		       				{ 
			 					if ( (Time[j] > 0) && (Time[k] > 0) )
								{
									if (j != k )
									{
										Td=Time[j]-Time[k]+5000.0+(Ta[k]-Ta[j]);
										if (j == 1) k1[k]->Fill(Td/10);
										if (j == 2) k2[k]->Fill(Td/10);
										if (j == 3) k3[k]->Fill(Td/10);
										if (j == 4) k4[k]->Fill(Td/10);
										if (j == 5) 
										{
											k5[k]->Fill(Td/10);
											rndUnif = randy->Rndm()*bw-(bw/2.);
											rndchannel = energy+rndUnif;
											th1l[k]->Fill(Time[k], rndchannel*a[k]+b[k]);
											FullSpec2D[k]->Fill(Td/10,rndchannel*a[k]+b[k]);
											FullSpecTd[k]->Fill(Td/10);

											FullSpecTd[k]->Fit("pol1","R", "", 434, 455);
											par0[k] = FullSpecTd[k]->GetFunction("pol1")->GetParameter(0);
											par1[k] = FullSpecTd[k]->GetFunction("pol1")->GetParameter(1); //This is the slope of the linear fit

											for (i = 1; i <= FullSpecTd[k]->GetNbinsX(); i++)
											{
												if (FullSpecTd[k]->GetBinContent(i) > (par0[k]+par1[k]*FullSpecTd[k]->GetBinCenter(i)))
												{
													OutSyncTd[k]->SetBinContent(i,0);
												}
											}
											InSyncTd[k]->Add(FullSpecTd[k],OutSyncTd[k],1,-1);
										}
			  						}
		       					}
							}
		   				}
		    		 	// Data analysis here

						m1=0;
						m2=0;
						sum1=0;
						sum2=0;
		     
						for (j = 1; j <=5; j++)
						{ 
							for (k = j+1; k <=5; k++)
							{ 
								Td=Time[j]-Time[k]+5000.0+(Ta[k]-Ta[j]);
								
								if ( (Td > 4340) && (Td < 4550) ) // 200 ns window
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

						
						for (j = 1; j <=5; j++)
						{ 
							if ((m1>0) && (m1<5)) n1[m1]->Fill(energyF[j]);
						}
		     

						//Reinitialise arrays

						for (j = 1; j <=5; j++)
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
						rndUnif = randy->Rndm()*bw-(bw/2.);
						rndchannel = energy+rndUnif; 
						//std::cout << "energy " << energy << " rndchannel " << rndchannel << std::endl;
						energycalib = (a[detectorID]*rndchannel)+(b[detectorID]);
						//if (energyS[detectorID] == -1) 
						energyS[detectorID] = energycalib;
						h4[detectorID]->Fill(energyS[detectorID]);
					}
					
					//Ident 72,73,74,75,76,77,78,79 are fast signals
					if ((ident > 71) && (ident < 77))
					{
						detectorID = ident-71;
						rndUnif = randy->Rndm()*bw-(bw/2.);
						rndchannel = energy+rndUnif; 
						energycalib = (a[detectorID]*rndchannel)+(b[detectorID]);
						//if (energyS[detectorID] == -1) 
						energyF[detectorID] = energycalib;
					}

      	      		//Ident 88,89,90,91,92,93,94,95 are timing signals of fast signals
					if ((ident > 87) && (ident < 93))
					{
						detectorID = ident-87;
						
						tock = ((adcdata & 0x0000e000) >> 13);
						tick = (adcdata & 0x00001fff);

						sick = (tick/5000.0) * 20.0;

						if (tock != 7)
						{
							Time[detectorID] = ((TS * 100.0) + (tock * 20.0) + sick);
							t[detectorID]->Fill(Time[detectorID]/10);
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

				// TSlast = TS;

				loop:
				pos+=8;      
			}
		}
		end:  

		fclose(f);
	
		Run++;
	}
       
	finish:

	double TimeDiff = (double) (TS - TSfirst)*(1e-8)/60.0;
	
	printf("The first time stamp is (e-8 s): %" PRIu64 "\n", TSfirst);
	printf("The last time stamp is (e-8 s): %" PRIu64 "\n", TS);
	printf("The time difference is (min): %f\n", TimeDiff);


	for (j = 1; j <=5; j++)
	{
		h1[j]->Write();
		h2[j]->Write();
		h4[j]->Write();

		c1[j]->Write();
		c2[j]->Write();

		k1[j]->Write();
		k2[j]->Write();
		k3[j]->Write();
		k4[j]->Write();
		k5[j]->Write();

		n1[j]->Write();

		t[j]->Write();

		th1l[j]->Write();

		FullSpec2D[j]->Write();
		FullSpecTd[j]->Write();
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
	return 0;

}
