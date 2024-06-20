      //int N0first = s0first->GetXaxis()->GetNbins();
      //double bw0first = s0first->GetBinWidth(1);
      N0second = s0second->GetXaxis()->GetNbins();
      bw0second = s0second->GetBinWidth(1);
      //cout << "bw0second: " << bw0second << endl;
      //int N0third = s0third->GetXaxis()->GetNbins();
      //double bw0third = s0third->GetBinWidth(1);

      //int N1first = s1first->GetXaxis()->GetNbins();
      //double bw1first = s1first->GetBinWidth(1);
      N1second = s1second->GetXaxis()->GetNbins();
      bw1second = s1second->GetBinWidth(1);
      //int N1third = s1third->GetXaxis()->GetNbins();
      //double bw1third = s1third->GetBinWidth(1);

      //N2first = s2first->GetXaxis()->GetNbins();
      // bw2first = s2first->GetBinWidth(1);

      N2second = s2second->GetXaxis()->GetNbins();
      bw2second = s2second->GetBinWidth(1);

      //N2third = s2third->GetXaxis()->GetNbins();
      //bw2third = s2third->GetBinWidth(1);

      //N3first = s3first->GetXaxis()->GetNbins();
      //bw3first = s3first->GetBinWidth(1);

      N3second = s3second->GetXaxis()->GetNbins();
      bw3second = s3second->GetBinWidth(1);

      //N3third = s3third->GetXaxis()->GetNbins();
      //bw3third = s3third->GetBinWidth(1);

      // LaBr Det 0
      for(int i=1; i<=N; i++)

      {
        n = s0none->GetBinContent(i);
        ch = s0none->GetBinCenter(i);

        for(int j=1; j<=n; j++)
        {
          //cout << "n: " << n << endl;
          rndUnif = randy->Rndm()*bw-(bw/2.); //generate uniform rnd number inside binwidth
          //rndUnif = 0.158079*bw-(bw/2.); //generate uniform rnd number inside binwidth
          //cout << "random number: " << rndUnif << endl;
          rndchannel = ch+rndUnif; //sum to bincenter
          energy0first = (a10*rndchannel)+(b10);
          //s0first->Fill(energy0first);
          //energy0second = (a20*rndchannel*rndchannel)+(b20*rndchannel)+(c20);
          //s0second->Fill(energy0second);
          //energy0third = (a30*rndchannel*rndchannel*rndchannel)+(b30*rndchannel*rndchannel)+(c30*rndchannel)+d30;
          //s0third->Fill(energy0third);
        }
      }

      // LaBr Det 1
      /*for(int i=1; i<=N; i++)
      {
        n = s1none->GetBinContent(i);
        ch = s1none -> GetBinCenter(i);;

        for(int j=1; j<=n; j++)
        {
          rndUnif = randy->Rndm()*bw-(bw/2.); //generate uniform rnd number iside binwidth
          rndchannel = ch+rndUnif; //sum to bincenter
          energy1first = (a11*rndchannel)+(b11);
          s1first->Fill(energy1first);
          energy1second = (a21*rndchannel*rndchannel)+(b21*rndchannel)+(c21);
          s1second->Fill(energy1second);
          energy1third = (a31*rndchannel*rndchannel*rndchannel)+(b31*rndchannel*rndchannel)+(c31*rndchannel)+d31;
          s1third->Fill(energy1third);
        }
      }*/
      
      // LaBr Det 2
      /*for(int i=1; i<=N; i++)
      {
        n = s2none->GetBinContent(i);
        ch = s2none -> GetBinCenter(i);;

        for(int j=1; j<=n; j++)
        {
          rndUnif = randy->Rndm()*bw-(bw/2.); //generate uniform rnd number iside binwidth
          rndchannel = ch+rndUnif; //sum to bincenter
          energy2first = (a12*rndchannel)+(b12);
          s2first->Fill(energy2first);
          energy2second = (a22*rndchannel*rndchannel)+(b22*rndchannel)+(c22);
          s2second->Fill(energy2second);
          energy2third = (a32*rndchannel*rndchannel*rndchannel)+(b32*rndchannel*rndchannel)+(c32*rndchannel)+d32;
          s2third->Fill(energy2third);
        }
      }*/
      
      // LaBr Det3
      /*for(int i=1; i<=N; i++)
      {
        n = s3none->GetBinContent(i);
        ch = s3none -> GetBinCenter(i);;

        for(int j=1; j<=n; j++)
        {
          rndUnif = randy->Rndm()*bw-(bw/2.); //generate uniform rnd number iside binwidth
          rndchannel = ch+rndUnif; //sum to bincenter
          energy3first = (a13*rndchannel)+(b13);
          s3first->Fill(energy3first);
          energy3second = (a23*rndchannel*rndchannel)+(b23*rndchannel)+(c23);
          s3second->Fill(energy3second);
          energy3third = (a33*rndchannel*rndchannel*rndchannel)+(b33*rndchannel*rndchannel)+(c33*rndchannel)+d33;
          s3third->Fill(energy3third);
        }
      }*/

      // RF
      /*for (int i = 0; i < RFTiming->GetNbinsX(); i++) 
      {
        double RFsignaltime[i] = {RFTiming->GetBinContent(i+1)}; // for each entry in the histo, fill an array (type double) with the values of the histo
      }
      
      for (int i = 0; i < LaBr0SlowTime->GetNbinsX(); i++) // for each entry in the histo, fill an array (type double) with the values of the histo
      {
        double slowtimeLaBr0[i] = {LaBr0SlowTime->GetBinContent(i+1)};
      }

      for (int i = 0; i < LaBr1SlowTime->GetNbinsX(); i++) // for each entry in the histo, fill an array (type double) with the values of the histo
      {
        double slowtimeLaBr1[i] = {LaBr1SlowTime->GetBinContent(i+1)};
      }

      for (int i = 0; i < LaBr2SlowTime->GetNbinsX(); i++) // for each entry in the histo, fill an array (type double) with the values of the histo
      {
        double slowtimeLaBr2[i] = {LaBr2SlowTime->GetBinContent(i+1)};
      }

      for (int i = 0; i < LaBr3SlowTime->GetNbinsX(); i++) // for each entry in the histo, fill an array (type double) with the values of the histo
      {
        double slowtimeLaBr3[i] = {LaBr3SlowTime->GetBinContent(i+1)};
      }*/
  