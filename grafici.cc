#include <iostream>
#include <TCanvas.h>
#include <TGraphErrors.h>

using namespace std;

void grafici()
{
    double x[]{10,40,75,100,130,150};
    double y1_0[]{83.397,83.164,82.834,82.602,82.403,82.280};
    double s1_0[]{1.932,1.053,0.8463,0.7724,0.7208,0.690};
    
    double y1_sp[]{81.70,81.93,81.92,81.79,81.71,81.760};
    double s1_sp[]{2.17,1.393,0.928,0.886,0.78,0.642};
    
    double y9_0[]{95.974,95.939,95.780,95.617,95.554,95.488};
    double s9_0[]{1.276,0.650,0.526,0.510,0.462,0.4506};
    
    double y9_sp[]{95.443,95.586,95.443,95.251,95.263,95.159};
    double s9_sp[]{1.366,0.686,0.572,0.556,0.494,0.48};
    
    
    TGraphErrors* perc_1_0 = new TGraphErrors(6,x,y1_0,0,s1_0);
    TGraphErrors* perc_1_spread = new TGraphErrors(6,x,y1_sp,0,s1_sp);
    TGraphErrors* perc_3x3_0 = new TGraphErrors(6,x,y9_0,0,s9_0);
    TGraphErrors* perc_3x3_spread = new TGraphErrors(6,x,y9_sp,0,s9_sp);
    
    
    
    double sigma_sqrt[]{1.276*sqrt(10)/100,0.650*sqrt(40)/100,0.526*sqrt(75)/100,0.510*sqrt(100)/100,0.462*sqrt(130)/100,0.4506*sqrt(150)/100};
    TGraph* sigma = new TGraph(6,x,sigma_sqrt);
    
    double sigmaE[]{(0.4506/100)*(0.4506/100),(0.462/100)*(0.462/100),(0.510/100)*(0.510/100),(0.526/100)*(0.526/100),(0.650/100)*(0.650/100),(1.276/100)*(1.276/100)};
    double invx[]{0,1,2,3,4,5,6};
    TGraph* paper = new TGraph(6,invx,sigmaE);
    
    

    TCanvas * perc= new TCanvas("a","a",1000,100,2500,2000); 
    perc->Divide(2,2);
    perc->cd(1);
    perc_1_0->SetLineWidth(2);
    perc_1_0->SetMarkerStyle(8);
    perc_1_0->SetMarkerSize(1.3);
    perc_1_0->SetMarkerColor(46);
    perc_1_0->SetMaximum(87);
    perc_1_0->SetMinimum(79);
    perc_1_0->SetTitle("Percentage energy 1 cell (0,0)");
    perc_1_0->GetXaxis()->SetTitle("E [GeV]");
    perc_1_0->GetYaxis()->SetTitle("%E");
    perc_1_0->Draw("ACP");
    
    perc->cd(2);
    perc_1_spread->SetLineWidth(2);
    perc_1_spread->SetMarkerStyle(8);
    perc_1_spread->SetMarkerSize(1.3);
    perc_1_spread->SetMarkerColor(46);
    perc_1_spread->SetMaximum(87);
    perc_1_spread->SetMinimum(79);
    perc_1_spread->SetTitle("Percentage energy 1 cell (2.6,2.7)");
    perc_1_spread->GetXaxis()->SetTitle("E [GeV]");
    perc_1_spread->GetYaxis()->SetTitle("%E");
    perc_1_spread->Draw("ACP");
    
    perc->cd(3);
    perc_3x3_0->SetLineWidth(2);
    perc_3x3_0->SetMarkerStyle(8);
    perc_3x3_0->SetMarkerSize(1.3);
    perc_3x3_0->SetMarkerColor(42);
    perc_3x3_0->SetMaximum(100);
    perc_3x3_0->SetMinimum(93);
    perc_3x3_0->SetTitle("Percentage energy 3x3 cell (0,0)");
    perc_3x3_0->GetXaxis()->SetTitle("E [GeV]");
    perc_3x3_0->GetYaxis()->SetTitle("%E");
    perc_3x3_0->Draw("ACP");
    
    perc->cd(4);
    perc_3x3_spread->SetLineWidth(2);
    perc_3x3_spread->SetMarkerStyle(8);
    perc_3x3_spread->SetMarkerSize(1.3);
    perc_3x3_spread->SetMarkerColor(42);
    perc_3x3_spread->SetMaximum(100);
    perc_3x3_spread->SetMinimum(93);
    perc_3x3_spread->SetTitle("Percentage energy 3x3 cell (2.6,2.7)");
    perc_3x3_spread->GetXaxis()->SetTitle("E [GeV]");
    perc_3x3_spread->GetYaxis()->SetTitle("%E");
    perc_3x3_spread->Draw("ACP");
    
    
    perc->SaveAs("/Users/eugenia/desktop/EMCal2/grafici/percE.png");  

    TCanvas * sig= new TCanvas("b","b",1000,100,2500,2000); 
    sigma->SetLineWidth(2);
    sigma->SetMarkerStyle(8);
    sigma->SetMarkerSize(1.3);
    sigma->SetMarkerColor(36);
    sigma->SetMaximum(5);
    sigma->SetMinimum(0);
    sigma->SetTitle("k = sigma(E)/E * sqrt(E) 3x3 cells (0,0)");
    sigma->GetXaxis()->SetTitle("E [GeV]");
    sigma->GetYaxis()->SetTitle("k");
    sigma->Draw("ACP");
    
    sig->SaveAs("/Users/eugenia/desktop/EMCal2/grafici/k.png");
    
    TCanvas * pap= new TCanvas("c","c",1000,100,2500,2000); 
    paper->SetLineWidth(2);
    paper->SetMarkerStyle(8);
    paper->SetMarkerSize(1.3);
    paper->SetMarkerColor(36);
    paper->SetTitle("(sigma(E)/E)^2");
    paper->GetXaxis()->SetTitle("E [GeV]");
    paper->GetYaxis()->SetTitle("(sigma(E)/E)^2");
    paper->Draw("ACP");
    
    pap->SaveAs("/Users/eugenia/desktop/EMCal2/grafici/paper.png");


    
}