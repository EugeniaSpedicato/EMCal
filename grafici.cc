#include <iostream>
#include <TCanvas.h>
#include <TGraphErrors.h>

using namespace std;

void grafici()
{
    double x[]{10,40,75,100,130,150};
    
   /*//Sigma gaussiane
    double y1_0[]{83.397,83.164,82.834,82.602,82.403,82.280};
    double s1_0[]{1.932,1.053,0.8463,0.7724,0.7208,0.690};
    
    double y1_sp[]{81.70,81.93,81.92,81.79,81.71,81.760};
    double s1_sp[]{2.17,1.393,0.928,0.886,0.78,0.642};
    
    double y9_0[]{95.974,95.939,95.780,95.617,95.554,95.488};
    double s9_0[]{1.283,0.640,0.528,0.502,0.467,0.461};
    
    double y9_sp[]{95.443,95.586,95.443,95.251,95.263,95.159};
    double s9_sp[]{1.366,0.686,0.572,0.556,0.494,0.48};
    
    //RMS
    double y1_0RMS[]{83.397,83.170,82.812,82.615,82.419,82.280};
    double s1_0RMS[]{1.932,1.065,0.8574,0.7753,0.7124,0.6987};
    
    double y1_spRMS[]{70.31,70.23,69.75,69.53,69.26,69.26};
    double s1_spRMS[]{12.23,11.98,11.99,12,12.26,11.98};*/
    
    double y9_0RMS[]{0.95837,0.95871,0.95658,0.95507,0.95372,0.95266};
    double s9_0RMS[]{0.01417,0.00753,0.006788,0.006849,0.007031,0.00711};
    
    double y9_spRMS[]{0.94641,0.94801,0.9459,0.94395,0.94234,0.94217};
    double s9_spRMS[]{0.02752,0.01942,0.01885,0.01982,0.02047,0.01846};
    
    
    /*TGraphErrors* perc_1_0 = new TGraphErrors(6,x,y1_0,0,s1_0);
    TGraphErrors* perc_1_spread = new TGraphErrors(6,x,y1_sp,0,s1_sp);
    TGraphErrors* perc_3x3_0 = new TGraphErrors(6,x,y9_0,0,s9_0);
    TGraphErrors* perc_3x3_spread = new TGraphErrors(6,x,y9_sp,0,s9_sp);
    
    TGraphErrors* perc_1_0RMS = new TGraphErrors(6,x,y1_0RMS,0,s1_0RMS);
    TGraphErrors* perc_1_spreadRMS = new TGraphErrors(6,x,y1_spRMS,0,s1_spRMS);*/
    TGraphErrors* perc_3x3_0RMS = new TGraphErrors(6,x,y9_0RMS,0,s9_0RMS);
    TGraphErrors* perc_3x3_spreadRMS = new TGraphErrors(6,x,y9_spRMS,0,s9_spRMS);
    
    
    /*RMS
    double sigma_sqrt[]{0.01401*sqrt(10),0.00753*sqrt(40),0.006788*sqrt(75),0.006849*sqrt(100),0.007031*sqrt(130),0.00711*sqrt(150)};*/
    
    //Gauss
     double sigma_sqrt[]{0.01283*sqrt(10)/0.95945,0.00640*sqrt(40)/0.95955,0.00528*sqrt(75)/0.95776,0.00502*sqrt(100)/0.95652,0.00467*sqrt(130)/0.95568,0.00461*sqrt(150)/0.95485};
    
    TGraph* sigma = new TGraph(6,x,sigma_sqrt);
    
    //double sigmaNORM[]{1.283,0.640,0.528,0.502,0.467,0.461};
    double sigmaNORM[]{0.0133,0.00667,0.00551,0.00523,0.00487,0.00483};
    TGraph* sigmainv = new TGraph(6,x,sigmaNORM);
    TF1 *f1 = new TF1("f1","sqrt(([0]/sqrt(x))*([0]/sqrt(x))+([1]/x)*([1]/x)+([2])*([2]))",0,150);

    /*TCanvas * perc= new TCanvas("a","a",1000,100,2500,2000); 
    perc->Divide(2,2);
    perc->cd(1);
    perc_1_0->SetLineWidth(2);
    perc_1_0->SetMarkerStyle(8);
    perc_1_0->SetMarkerSize(1.3);
    perc_1_0->SetMarkerColor(46);
    perc_1_0->SetMaximum(87);
    perc_1_0->SetMinimum(78);
    perc_1_0->SetTitle("Percentage energy 1 cell (0,0)");
    perc_1_0->GetXaxis()->SetTitle("E [GeV]");
    perc_1_0->GetYaxis()->SetTitle("%E_reco/E_in");
    perc_1_0->Draw("ACP");
    
    perc->cd(2);
    perc_1_0RMS->SetLineWidth(2);
    perc_1_0RMS->SetMarkerStyle(8);
    perc_1_0RMS->SetMarkerSize(1.3);
    perc_1_0RMS->SetMarkerColor(42);
    perc_1_0RMS->SetMaximum(87);
    perc_1_0RMS->SetMinimum(78);
    perc_1_0RMS->SetTitle("Percentage energy 1 cell (0,0) RMS");
    perc_1_0RMS->GetXaxis()->SetTitle("E [GeV]");
    perc_1_0RMS->GetYaxis()->SetTitle("%E_reco/E_in");
    perc_1_0RMS->Draw("ACP");
    
    
    perc->cd(3);
    perc_1_spread->SetLineWidth(2);
    perc_1_spread->SetMarkerStyle(8);
    perc_1_spread->SetMarkerSize(1.3);
    perc_1_spread->SetMarkerColor(46);
    perc_1_spread->SetMaximum(87);
    perc_1_spread->SetMinimum(78);
    perc_1_spread->SetTitle("Percentage energy 1 cell (2.6,2.7)");
    perc_1_spread->GetXaxis()->SetTitle("E [GeV]");
    perc_1_spread->GetYaxis()->SetTitle("%E_reco/E_in");
    perc_1_spread->Draw("ACP");
    perc_1_spread->Draw("ACP SAME");
    
    
    perc->cd(4);
    perc_1_spreadRMS->SetLineWidth(2);
    perc_1_spreadRMS->SetMarkerStyle(8);
    perc_1_spreadRMS->SetMarkerSize(1.3);
    perc_1_spreadRMS->SetMarkerColor(42);
    perc_1_spreadRMS->SetMaximum(87);
    perc_1_spreadRMS->SetMinimum(50);
    perc_1_spreadRMS->SetTitle("Percentage energy 1 cell (2.6,2.7) RMS");
    perc_1_spreadRMS->GetXaxis()->SetTitle("E [GeV]");
    perc_1_spreadRMS->GetYaxis()->SetTitle("%E_reco/E_in");
    perc_1_spreadRMS->Draw("ACP");
    
    perc->SaveAs("/Users/eugenia/desktop/EMCal2/grafici/percE1cell.png");  */
    
    
    TCanvas * perc3= new TCanvas("a","a",5000,1000,3500,2500); 
    perc3->Divide(2,1);
    /*perc3->cd(1);
    perc_3x3_0->SetLineWidth(2);
    perc_3x3_0->SetMarkerStyle(8);
    perc_3x3_0->SetMarkerSize(1.3);
    perc_3x3_0->SetMarkerColor(46);
    perc_3x3_0->SetMaximum(100);
    perc_3x3_0->SetMinimum(90);
    perc_3x3_0->SetTitle("Percentage energy 3x3 cell (0,0)");
    perc_3x3_0->GetXaxis()->SetTitle("E [GeV]");
    perc_3x3_0->GetYaxis()->SetTitle("%E_reco/E_in");
    perc_3x3_0->Draw("ACP");*/
    
    perc3->cd(1);
    perc_3x3_0RMS->SetLineWidth(2);
    perc_3x3_0RMS->SetMarkerStyle(8);
    perc_3x3_0RMS->SetMarkerSize(1.3);
    perc_3x3_0RMS->SetMarkerColor(46);
    perc_3x3_0RMS->SetMaximum(1);
    perc_3x3_0RMS->SetMinimum(0.90);
    perc_3x3_0RMS->SetTitle("Fraction of energy 3x3 cell (0,0) RMS");
    perc_3x3_0RMS->GetXaxis()->SetTitle("E_true [GeV]");
    perc_3x3_0RMS->GetYaxis()->SetTitle("E_reco/E_true");
    perc_3x3_0RMS->Draw("ACP");
    
    
    /*perc3->cd(3);
    perc_3x3_spread->SetLineWidth(2);
    perc_3x3_spread->SetMarkerStyle(8);
    perc_3x3_spread->SetMarkerSize(1.3);
    perc_3x3_spread->SetMarkerColor(46);
    perc_3x3_spread->SetMaximum(100);
    perc_3x3_spread->SetMinimum(90);
    perc_3x3_spread->SetTitle("Percentage energy 3x3 cell (2.6,2.7)");
    perc_3x3_spread->GetXaxis()->SetTitle("E [GeV]");
    perc_3x3_spread->GetYaxis()->SetTitle("%E_reco/E_in");
    perc_3x3_spread->Draw("ACP");*/
    
    perc3->cd(2);
    perc_3x3_spreadRMS->SetLineWidth(2);
    perc_3x3_spreadRMS->SetMarkerStyle(8);
    perc_3x3_spreadRMS->SetMarkerSize(1.3);
    perc_3x3_spreadRMS->SetMarkerColor(kBlue);
    perc_3x3_spreadRMS->SetMaximum(1);
    perc_3x3_spreadRMS->SetMinimum(0.90);
    perc_3x3_spreadRMS->SetTitle("Fraction of energy 3x3 cell (2.6,2.7) RMS");
    perc_3x3_spreadRMS->GetXaxis()->SetTitle("E_true [GeV]");
    perc_3x3_spreadRMS->GetYaxis()->SetTitle("E_reco/E_true");
    perc_3x3_spreadRMS->Draw("ACP");
    
    
    
    perc3->SaveAs("/Users/eugenia/desktop/EMCal2/grafici/percE3x3.png");  

    TCanvas * sig= new TCanvas("ab","ab",1000,100,2500,2000); 
    sigma->SetLineWidth(2);
    sigma->SetMarkerStyle(8);
    sigma->SetMarkerSize(1.3);
    sigma->SetMarkerColor(36);
    sigma->SetMaximum(0.5);
    sigma->SetMinimum(0);
    sigma->SetTitle("k = sigma(E)/E * sqrt(E) 3x3 cells (0,0)");
    sigma->GetXaxis()->SetTitle("E [GeV]");
    sigma->GetYaxis()->SetTitle("k");
    sigma->Draw("ACP");
    sig->SaveAs("/Users/eugenia/desktop/EMCal2/grafici/k.png");

    TCanvas * ris= new TCanvas("b","b",1000,100,2500,2000); 
    sigmainv->SetLineWidth(2);
    sigmainv->SetMarkerStyle(8);
    sigmainv->SetMarkerSize(1.3);
    sigmainv->SetMarkerColor(36);
    sigmainv->SetMaximum(0.02);
    sigmainv->SetMinimum(0);
    sigmainv->SetTitle("Sigma(E)/(E)");
    sigmainv->GetXaxis()->SetTitle("E [GeV]");
    sigmainv->GetYaxis()->SetTitle("Sigma(E)/(E)");
    sigmainv->Fit("f1", "R");
    gStyle->SetOptFit();
    sigmainv->Draw("AP");

    ris->SaveAs("/Users/eugenia/desktop/EMCal2/grafici/res.png");
    
    



    
}