#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <vector>
#include <iostream>

void makePurityPlot(std::vector<TString> dirs, TString suffix) {

    // IMPORTANTISSIMO
    TH1::AddDirectory(false);
    gStyle->SetOptStat(0);

    TCanvas* c = new TCanvas("c", "", 1600, 1200);
    gPad->SetLeftMargin(0.13);   // More space for Y axis label
    gPad->SetRightMargin(0.05);  // Less space on right
    gPad->SetBottomMargin(0.13); // More space for X axis label
    gPad->SetTopMargin(0.05);    // Less space on top
//    c->SetGrid();

    TLegend* leg = new TLegend(0.40, 0.15, 0.88, 0.65);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    int colors[] = {kRed+1, kBlue+1, kGreen+2, kMagenta+1,
                    kOrange+7, kCyan+2, kViolet+1};

    int ic = 0;
    bool first = true;

    for (auto& dirname : dirs) {

        TString filepath = dirname + "/PurityMuonSpec.root";

        std::cout << "Processing file: " << filepath << std::endl;
        TFile f(filepath, "READ");

        if (f.IsZombie()) {
            std::cout << "Cannot open " << filepath << std::endl;
            continue;
        }

        TH1* h = dynamic_cast<TH1*>(f.Get("hPurityVsP"));
        if (!h) {
            std::cout << "Histogram not found in " << dirname << std::endl;
            continue;
        }

        h->SetName(Form("h_%d", ic));
        h->SetLineColor(colors[ic % 7]);
        h->SetLineWidth(3);

        if (first) {
            h->SetTitle("Purity vs Momentum");
            h->GetYaxis()->SetRangeUser(0.0, 1.0);
            h->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
            h->GetYaxis()->SetTitle("Purity");
            h->Draw("HISTO");
            first = false;
        } else {
            h->Draw("HISTO SAME");
        }

        leg->AddEntry(h, dirname, "l");

        ic++;
    }

    leg->Draw();
    c->SaveAs(suffix + "_PurityVsP.pdf");
    c->SaveAs(suffix + "_PurityVsP.png");
}


void overlayPurity(){
    
    std::vector<TString> dirs = {
        "msDoOutwardPropagationfalsemsDoInwardRefitfalse/field_scaled0.2",
        "msDoOutwardPropagationfalsemsDoInwardRefitfalse/field_scaled0.4",
        "msDoOutwardPropagationfalsemsDoInwardRefitfalse/field_scaled0.6",
        "msDoOutwardPropagationfalsemsDoInwardRefitfalse/field_scaled0.8",
        "msDoOutwardPropagationfalsemsDoInwardRefitfalse/field_scaled1.0",
        "msDoOutwardPropagationfalsemsDoInwardRefitfalse/field_scaled1.2",
        "msDoOutwardPropagationfalsemsDoInwardRefitfalse/field_scaled1.4"
    };

    TString suffix = "field_scaled";
    makePurityPlot(dirs, suffix);
    
    std::vector<TString> dirs2 = {
        "msDoOutwardPropagationfalsemsDoInwardRefitfalse/field_scaled1.0_resX0.05_resY0.25",
        "msDoOutwardPropagationfalsemsDoInwardRefitfalse/field_scaled1.0_resX0.02_resY0.10",
        "msDoOutwardPropagationfalsemsDoInwardRefitfalse/field_scaled1.0_resX0.01_resY0.05",
        "msDoOutwardPropagationfalsemsDoInwardRefitfalse/field_scaled1.0_resX0.005_resY0.025",
        "msDoOutwardPropagationfalsemsDoInwardRefitfalse/field_scaled1.0_resX0.002_resY0.010",
        "msDoOutwardPropagationfalsemsDoInwardRefitfalse/field_scaled1.0_resX0.0005_resY0.0005"
    
    };

    TString suffix2 = "resolution_scaled";
    makePurityPlot(dirs2, suffix2);
    
    std::vector<TString> dirs3 = {
        "msDoOutwardPropagationfalsemsDoInwardRefitfalse/field_scaled1.0",
        "msDoOutwardPropagationtruemsDoInwardRefitfalse/field_scaled1.0",
        "msDoOutwardPropagationtruemsDoInwardRefittrue/field_scaled1.0"
    
    };

    TString suffix3 = "ambiguity4planes";
    makePurityPlot(dirs3, suffix3);

    std::vector<TString> dirs4 = {
        "std_rec_6planes/field_scaled1.0",
        "outer_6planes/field_scaled1.0",
        "refit_6planes/field_scaled1.0"
    
    };

    TString suffix4 = "ambiguity6planes";
    makePurityPlot(dirs4, suffix4);
    

}