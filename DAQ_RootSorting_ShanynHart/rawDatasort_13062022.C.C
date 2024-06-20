/* Advanced sorting code for TDR format data */
/* Produces ROOT TTrees */

/* Adapted from P.M Jones code labrsort4 */

/* (c) S. Hart */
/* (10/06/22)  */


#include <stdio.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

// ROOT libraries
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TTree.h"
#include "TObject.h"
#include "TString.h"
#include "TRandom.h"
#include "TH2F.h"

#include <stdint.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <inttypes.h>
#include <initializer_list>

#define SIZE  16384

int32_t buffer[SIZE];
uint64_t TimeStamp1[100]={0}, TimeStamp2[100]={0}, TimeStamp3[100]={0} ,TimeStamp4[100]={0}, TimeStamp5[100]={0}, TimeStamp6[100]={0}, TimeStamp7[100]={0}, TimeStamp8[100]={0};
uint16_t QDC[100]={0}, qdc;
uint64_t TS2;

uint64_t TStop, TSbot, data;
double T1, T2;
int64_t count=0;

unsigned long int pos = 0;
int verbose = 255;

int main (int argc, char** argv)
{
  FILE *f;
  int i;
  char file_in[50];
  char file_out[50];
  int  identity, card=0, adcdata;
  int  fastLaBr1, slowLaBr1, fasttimeLaBr1, slowtimeLaBr1;
  int  fastLaBr2, slowLaBr2, fasttimeLaBr2, slowtimeLaBr2;
  int polaristime;
  int tick, tock;
  double sick;
  uint64_t TS;
  int RunEnd=999999, Run=0;

  int no_read;
  std::vector<int> time_stamp;
  
  
  if ((argc < 2) || (argc > 3))
    {
	fprintf(stderr, "Usage: %s data [last_subrun]\n", argv[0]);
	exit(1);
    }

  /* ROOT STUFF */

    sprintf(file_out, "%s_rawData.root", argv[1]);
    if (argc > 2) sscanf(argv[2], "%d", &RunEnd);

    TFile *g = new TFile(file_out,"recreate");

    printf("ROOT file %s opened...\n", file_out);

    TTree *LaBrData = new TTree("LaBrData", "LaBrRawData"); 
    LaBrData->Branch("Identity", &identity , "identity/I");
    LaBrData->Branch("TimeStamp", &TS , "TS");
    LaBrData->Branch("ADCEnergyData", &adcdata , "adcdata/I");
    LaBrData->Branch("Card", &card , "card/I");

    TTree *LaBr1Data = new TTree("LaBr1Data", "LaBr1RawData"); 
    LaBr1Data->Branch("FastEnergyLaBr1", &fastLaBr1 , "fastLaBr1/I");
    LaBr1Data->Branch("SlowEnergyLaBr1", &slowLaBr1 , "slowLaBr1/I");
    LaBr1Data->Branch("FastTIMINGLaBr1", &fasttimeLaBr1 , "fasttimeLaBr1/I");
    LaBr1Data->Branch("SlowTIMINGLaBr1", &slowtimeLaBr1 , "slowtimeLaBr1/I");
    LaBr1Data->Branch("POLARISCh15Time", &polaristime , "polaristime/I");

    TTree *LaBr2Data = new TTree("LaBr2Data", "LaBr2RawData"); 
    LaBr2Data->Branch("FastEnergyLaBr2", &fastLaBr2 , "fastLaBr2/I");
    LaBr2Data->Branch("SlowEnergyLaBr2", &slowLaBr2 , "slowLaBr2/I");
    LaBr2Data->Branch("FastTIMINGLaBr2", &fasttimeLaBr2 , "fasttimeLaBr2/I");
    LaBr2Data->Branch("POLARISCh15Time", &polaristime , "polaristime/I");
    LaBr2Data->Branch("SlowTIMINGLaBr2", &slowtimeLaBr2 , "slowtimeLaBr2/I");


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
						    int ident = (data & 0x0fff0000) >> 16;
							//printf("  ident: %d /  ", ident);  // resets at 8161
						
						if ((ident == 64) || (ident == 65) || (ident == 72) || (ident == 73) ||(ident == 79) ||(ident == 80) ||(ident == 81) ||(ident == 88) ||(ident == 89) || (ident == 95))   
						{
              identity = ident;
							adcdata = (data & 0x0000ffff);
							card = (ident / 32);
							TS = (TS & 0x0000fffff0000000ULL);

							//printf(" === i = %d / ident: %d / energy: %d / TS: %ld --- \n", i, ident, adcdata, TS);  // resets at 8161

							LaBrData->Fill();

              if (identity == 64) { slowLaBr1 =  (data & 0x0000ffff); LaBr1Data->Fill(); }
              if (identity == 80) { slowtimeLaBr1 =  (data & 0x0000ffff); LaBr1Data->Fill(); }
              if (identity == 72) { fastLaBr1 =  (data & 0x0000ffff); LaBr1Data->Fill(); }
              if (identity == 88) { fasttimeLaBr1 =  (data & 0x0000ffff); LaBr1Data->Fill(); }

              if (identity == 65) { slowLaBr2 =  (data & 0x0000ffff); LaBr2Data->Fill(); }
              if (identity == 81) { slowtimeLaBr2 =  (data & 0x0000ffff); LaBr2Data->Fill(); }
              if (identity == 73) { fastLaBr2 =  (data & 0x0000ffff); LaBr2Data->Fill(); }
              if (identity == 89) { fasttimeLaBr2 =  (data & 0x0000ffff); LaBr2Data->Fill(); }

              if (identity == 95) { polaristime =  (data & 0x0000ffff);  LaBr1Data->Fill(); LaBr2Data->Fill(); }

						}

                                       
                        
                    }

                    skip:
		            count++;

		            loop:
		            pos+=8;
                }
        }
    
     	end:  

       fclose(f);
       
       Run++;
       
    }

	finish:

	LaBrData->Write();
	LaBr1Data->Write();
	LaBr2Data->Write();

	delete g;
}


