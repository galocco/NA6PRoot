// NA6PCCopyright

#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMemFile.h>
#include <TLine.h>
#include <TString.h>
#include <TStyle.h>

#include <algorithm>
#include <filesystem>
#include <array>
#include <span>
#include <stdexcept>
#include <string>

namespace
{
template <typename T>
T& getPlotInput(TFile& input, const char* name)
{
  auto* object = input.Get<T>(name);
  if (!object) {
    throw std::runtime_error(std::string("Missing plot input ") + name + " in " + input.GetName());
  }
  return *object;
}

void configureCanvas(TCanvas& canvas)
{
  canvas.SetMargin(.13, .05, .14, .17);
  canvas.SetTicks(1, 1);
  canvas.SetGridy();
}

void drawHeading(const char* heading, const std::string& sample)
{
  TLatex label;
  label.SetNDC();
  label.SetTextFont(62);
  label.SetTextSize(.043);
  label.DrawLatex(.13, .95, heading);
  label.SetTextFont(42);
  label.SetTextSize(.028);
  label.DrawLatex(.13, .89, sample.c_str());
}

void saveCanvas(TCanvas& canvas, TFile& output, const std::string& dir, const std::string& pdf, bool savePNG = true)
{
  canvas.RedrawAxis();
  canvas.SaveAs((dir + "/" + canvas.GetName() + ".png").c_str());
  output.cd();
  canvas.Write();
}

void drawRate(TFile& input, TFile& output, const std::string& dir, const std::string& pdf,
              const std::string& sample, const std::string& name, const char* heading, const char* axisTitle,
              bool complement = false, const std::string& variable = "")
{
  struct Series {
    int cut, hits, color, marker;
    double offset;
    const char* label;
  };
  constexpr std::array<Series, 3> inclusive{{
    {3, 0, kAzure + 2, 20, -.11, "VT threshold #geq3 hit"},
    {5, 0, kOrange + 7, 25, .15, "VT threshold #geq5 hit"},
    {4, 0, kRed + 1, 24, 0., "VT threshold #geq4 hit"},
  }};
  constexpr std::array<Series, 6> byHits{{
    {3, 3, kAzure + 2, 20, -.21, "Threshold #geq3, 3 VT hits"},
    {3, 4, kTeal + 3, 21, -.07, "Threshold #geq3, 4 VT hits"},
    {3, 5, kViolet + 1, 22, .07, "Threshold #geq3, 5 VT hits"},
    {5, 5, kOrange + 7, 27, .24, "Threshold #geq5, 5 VT hits"},
    {4, 4, kRed + 1, 24, -.28, "Threshold #geq4, 4 VT hits"},
    {4, 5, kBlack, 26, .32, "Threshold #geq4, 5 VT hits"},
  }};
  const bool differential = !variable.empty();
  const std::span<const Series> series = differential ? std::span<const Series>(byHits) : std::span<const Series>(inclusive);
  const std::string suffix = differential ? "_" + variable : "";
  const std::string plotName = (complement ? "fakeRate" : name) + suffix;
  const auto& total = *getPlotInput<TEfficiency>(input, (name + suffix + (differential ? "_hits3" : "") + "_minVT3").c_str()).GetTotalHistogram();
  std::array<TGraphAsymmErrors, byHits.size()> graphs;
  TCanvas canvas(plotName.c_str(), heading, 850, 650);
  configureCanvas(canvas);
  canvas.SetTopMargin(.29);
  TH1D frame((plotName + "Frame").c_str(), (std::string(";VT track hits;") + axisTitle + " (%)").c_str(),
             total.GetNbinsX(), differential ? total.GetXaxis()->GetXmin() : .5,
             differential ? total.GetXaxis()->GetXmax() : total.GetNbinsX() + .5);
  frame.SetDirectory(nullptr);
  frame.SetMinimum(0);
  frame.SetMaximum(100);
  if (differential) {
    const std::string xTitle = variable == "pt" ? "p_{T}^{MS} (GeV/c)" : (variable == "y" ? "y_{lab}^{MS}" : "p^{MS} (GeV/c)");
    frame.GetXaxis()->SetTitle(xTitle.c_str());
    frame.GetXaxis()->SetNdivisions(-total.GetNbinsX());
  } else {
    for (int bin = 1; bin <= total.GetNbinsX(); ++bin) {
      frame.GetXaxis()->SetBinLabel(bin, TString::Format("%d", bin + 2).Data());
    }
  }
  frame.GetYaxis()->SetNdivisions(-504);
  frame.Draw("AXIS");

  TLegend legend(.15, .72, .88, differential ? .865 : .82);
  legend.SetNColumns(2);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(.024);
  for (size_t i = 0; i < series.size(); ++i) {
    const auto& style = series[i];
    const auto selection = (differential ? "_hits" + std::to_string(style.hits) : "") + "_minVT" + std::to_string(style.cut);
    if (!input.GetKey((name + suffix + selection).c_str())) continue;
    auto& rate = getPlotInput<TEfficiency>(input, (name + suffix + selection).c_str());
    auto& graph = graphs[i];
    graph.SetName((plotName + selection).c_str());
    graph.SetMarkerColor(style.color);
    graph.SetLineColor(style.color);
    graph.SetLineWidth(2);
    graph.SetMarkerStyle(style.marker);
    graph.SetMarkerSize(1.4);
    for (int bin = 1; bin <= total.GetNbinsX(); ++bin) {
      if (rate.GetTotalHistogram()->GetBinContent(bin) == 0) {
        continue;
      }
      const int point = graph.GetN();
      const double value = rate.GetEfficiency(bin);
      const double errorLow = rate.GetEfficiencyErrorLow(bin);
      const double errorUp = rate.GetEfficiencyErrorUp(bin);
      const double x = frame.GetBinCenter(bin) + style.offset * frame.GetBinWidth(bin);
      graph.SetPoint(point, x, 100. * (complement ? 1. - value : value));
      graph.SetPointError(point, 0., 0., 100. * (complement ? errorUp : errorLow), 100. * (complement ? errorLow : errorUp));
    }
    if (graph.GetN() > 0) {
      graph.Draw("PZ SAME");
    }
    std::string legendLabel = style.label;
    if (differential && graph.GetN() == 0) legendLabel += " (N/A)";
    legend.AddEntry(&graph, legendLabel.c_str(), "lep");
    if (!differential && name == "purity" && style.cut == 5) {
      TLatex label;
      label.SetTextSize(.026);
      label.SetTextAlign(22);
      label.SetTextColor(style.color);
      for (int bin = 1; bin <= 2; ++bin) {
        if (rate.GetTotalHistogram()->GetBinContent(bin) == 0) {
          label.DrawLatex(bin + style.offset, 5., "N/A");
        }
      }
    }
  }
  legend.Draw();
  drawHeading(heading, sample + "  |  Exact binomial intervals, 68.3%");
  saveCanvas(canvas, output, dir, pdf, !differential);
}

void drawPurityRatio(TFile& input, TFile& output, const std::string& dir, const std::string& pdf,
                     const std::string& sample, const std::string& variable, const std::string& metric = "purity")
{
  const std::string name = metric + "Ratio_" + variable;
  const std::string label = metric == "purity" ? "Purity" : (metric == "fakeRate" ? "Fake rate" : "Efficiency");
  const std::array<std::string, 2> names{metric + "Ratio_" + variable, metric + "Ratio4over3_" + variable};
  std::array<TGraphAsymmErrors*, 2> ratios{input.Get<TGraphAsymmErrors>(names[0].c_str()), input.Get<TGraphAsymmErrors>(names[1].c_str())};
  const auto& axes = *getPlotInput<TEfficiency>(input, ((metric == "efficiency" ? "efficiency_" : (metric == "fakeRate" ? "fakeRate_" : "purity_")) + variable + "_hits5_minVT3").c_str()).GetTotalHistogram();
  TCanvas canvas(name.c_str(), "Ratio, 5-hit VT tracks", 850, 650);
  configureCanvas(canvas);
  canvas.SetTopMargin(.25);
  TH1D frame((name + "Frame").c_str(), "", axes.GetNbinsX(), axes.GetXaxis()->GetXmin(), axes.GetXaxis()->GetXmax());
  frame.SetDirectory(nullptr);
  const std::string xTitle = variable == "pt" ? "p_{T}^{MS} (GeV/c)" : (variable == "y" ? "y_{lab}^{MS}" : "p^{MS} (GeV/c)");
  frame.GetXaxis()->SetTitle(xTitle.c_str());
  frame.GetXaxis()->SetNdivisions(-axes.GetNbinsX());
  frame.GetYaxis()->SetTitle((label + " (threshold) / " + label + " (threshold #geq3)").c_str());
  frame.GetYaxis()->SetTitleOffset(1.25);
  double minimum = .7, maximum = 1.3;
  for (auto* ratio : ratios) {
    if (!ratio) continue;
    for (int i = 0; i < ratio->GetN(); ++i) {
      const double y = ratio->GetPointY(i);
      minimum = std::min(minimum, y - ratio->GetErrorYlow(i) - .05);
      maximum = std::max(maximum, y + ratio->GetErrorYhigh(i) + .05);
    }
  }
  frame.SetMinimum(std::max(0., minimum));
  frame.SetMaximum(maximum);
  frame.Draw("AXIS");
  TLine unity(axes.GetXaxis()->GetXmin(), 1., axes.GetXaxis()->GetXmax(), 1.);
  unity.SetLineStyle(2);
  unity.SetLineColor(kGray + 2);
  unity.Draw();
  TLegend legend(.15, .765, .83, .823);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetNColumns(2);
  legend.SetTextSize(.029);
  const std::array<int, 2> colors{kAzure + 2, kOrange + 7};
  const std::array<int, 2> markers{20, 25};
  std::array<TGraphAsymmErrors, 2> displayed;
  for (size_t j = 0; j < ratios.size(); ++j) {
    if (!ratios[j]) continue;
    displayed[j] = *ratios[j];
    displayed[j].SetName((names[j] + "_display").c_str());
    displayed[j].SetMarkerStyle(markers[j]);
    displayed[j].SetMarkerSize(1.3);
    displayed[j].SetMarkerColor(colors[j]);
    displayed[j].SetLineColor(colors[j]);
    displayed[j].SetLineWidth(2);
    for (int i = 0; i < displayed[j].GetN(); ++i) {
      const double x = displayed[j].GetPointX(i);
      const int bin = axes.GetXaxis()->FindFixBin(x);
      displayed[j].SetPoint(i, x + (j == 0 ? -.08 : .08) * axes.GetBinWidth(bin), displayed[j].GetPointY(i));
    }
    if (displayed[j].GetN()) displayed[j].Draw("P SAME");
    legend.AddEntry(&displayed[j], j == 0 ? "Threshold #geq5 / #geq3" : "Threshold #geq4 / #geq3", "lep");
  }
  legend.Draw();
  drawHeading(("Ratio of " + label + ", 5-hit VT tracks").c_str(), sample + "  |  Paired event bootstrap intervals, 68.3%");
  saveCanvas(canvas, output, dir, pdf);
}

void drawTransitions(TFile& input, TFile& output, const std::string& dir, const std::string& pdf, const std::string& sample)
{
  const std::array<const char*, 3> names{"shortUnmatched", "shortCorrect", "shortFalse"};
  const std::array<const char*, 3> titles{"Unmatched", "Correct match", "Fake match"};
  const std::array<int, 3> colors{kAzure - 9, kTeal + 2, TColor::GetColor("#a53543")};
  std::array<TBox, 3> keys;
  TCanvas canvas("shortTransitions", "Esito dei match corretti da 3 o 4 hit", 850, 650);
  configureCanvas(canvas);
  canvas.SetTopMargin(.25);
  TH1D frame("transitionFrame", ";VT hits for correct match at threshold #geq3;MS tracks (%)", 2, 2.5, 4.5);
  frame.SetDirectory(nullptr);
  frame.GetXaxis()->SetBinLabel(1, "3");
  frame.GetXaxis()->SetBinLabel(2, "4");
  frame.SetMaximum(105.);
  frame.GetYaxis()->SetNdivisions(505);
  frame.Draw("AXIS");
  TLegend legend(.15, .765, .91, .823);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetNColumns(3);
  legend.SetTextSize(.028);
  std::array<TH1D*, 3> counts{};
  for (size_t i = 0; i < counts.size(); ++i) {
    counts[i] = &getPlotInput<TH1D>(input, names[i]);
    keys[i].SetFillColor(colors[i]);
    if (counts[i]->Integral() > 0) {
      legend.AddEntry(&keys[i], titles[i], "f");
    }
  }
  TLatex label;
  label.SetTextAlign(22);
  label.SetTextSize(.035);
  for (int bin = 1; bin <= 2; ++bin) {
    double total = 0.;
    for (auto* hist : counts) total += hist->GetBinContent(bin);
    double bottom = 0.;
    for (size_t i = 0; i < counts.size(); ++i) {
      const double count = counts[i]->GetBinContent(bin);
      if (count == 0.) continue;
      const double height = 100. * count / total;
      TBox bar;
      bar.SetFillColor(colors[i]);
      bar.DrawBox(bin + 1.72, bottom, bin + 2.28, bottom + height);
      label.SetTextColor(i == 0 ? kBlack : kWhite);
      label.DrawLatex(bin + 2, bottom + .5 * height, std::to_string(static_cast<int>(count)).c_str());
      bottom += height;
    }
  }
  legend.Draw();
  drawHeading("Outcome for VT threshold #geq5", sample + "  |  Correct 3- or 4-hit matches at threshold #geq3");
  saveCanvas(canvas, output, dir, pdf);
}

void drawOrigins(TFile& input, TFile& output, const std::string& dir, const std::string& pdf, const std::string& sample)
{
  auto& origins = getPlotInput<TH1D>(input, "falseOrigins");
  TCanvas canvas("falseOrigins", "Provenienza dei falsi a 5 hit", 850, 650);
  configureCanvas(canvas);
  canvas.SetBottomMargin(.18);
  TH1D frame("originFrame", ";Outcome at threshold #geq3;False matches", 3, .5, 3.5);
  frame.SetDirectory(nullptr);
  const std::array<int, 5> order{5, 2, 4, 1, 3};
  const std::array<const char*, 5> labels{"Fake, #geq5 hits", "Correct, 3#minus4 hits", "Fake, 3#minus4 hits", "Unmatched", "Correct, #geq5 hits"};
  std::array<int, 5> shown{};
  int nShown = 0;
  for (size_t i = 0; i < order.size(); ++i) {
    if (origins.GetBinContent(order[i]) > 0.) shown[nShown++] = i;
  }
  frame.SetBins(std::max(1, nShown), .5, std::max(1, nShown) + .5);
  for (int i = 0; i < nShown; ++i) frame.GetXaxis()->SetBinLabel(i + 1, labels[shown[i]]);
  frame.SetMaximum(std::max(1., 1.19 * origins.GetMaximum()));
  frame.GetXaxis()->SetLabelSize(.033);
  frame.Draw("AXIS");
  TLatex label;
  label.SetTextAlign(22);
  label.SetTextSize(.034);
  for (int i = 0; i < nShown; ++i) {
    const double count = origins.GetBinContent(order[shown[i]]);
    TBox bar;
    bar.SetFillColor(shown[i] == 1 ? TColor::GetColor("#a53543") : kAzure + 2);
    bar.DrawBox(i + .72, 0., i + 1.28, count);
    label.DrawLatex(i + 1, count + .039 * frame.GetMaximum(), std::to_string(static_cast<int>(count)).c_str());
  }
  drawHeading("Origin of false matches at threshold #geq5", sample + "  |  " + std::to_string(static_cast<int>(origins.Integral())) + " false matches");
  saveCanvas(canvas, output, dir, pdf);
}
} // namespace

void plotMatchingPNG(const char* runDir = ".", int nEvents = 1000)
{
  const std::string dir = runDir;
  const std::string sample = "Matching VT#minusMS  |  " + std::to_string(nEvents) + " simulated events";
  const std::string pdf = dir + "/matchingComparisonROOT.pdf";
  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(1.);
  TFile input((dir + "/matchingQuality.root").c_str(), "READ");
  if (input.IsZombie()) throw std::runtime_error("Cannot open matchingQuality.root");
  TMemFile output("matchingPlots.root", "RECREATE");
  for (const auto* folder : {"purity", "efficiency", "fake"}) std::filesystem::create_directories(dir + "/" + folder);
  drawRate(input, output, dir + "/efficiency", pdf, sample, "efficiency", "Efficiency of thermal muons", "Efficiency");
  drawRate(input, output, dir + "/purity", pdf, sample, "purity", "Purity for thermal muons", "Purity");
  drawRate(input, output, dir + "/fake", pdf, sample, "fakeRate", "Fake rate (not a correct thermal match)", "Fake rate");
  drawTransitions(input, output, dir + "/fake", pdf, sample);
  drawOrigins(input, output, dir + "/fake", pdf, sample);

  if (input.GetKey("efficiency_pt_hits3_minVT3")) {
    for (const std::string variable : {"pt", "y", "p"}) {
      const auto differentialPDF = dir + "/matchingDifferential_" + variable + ".pdf";
      drawRate(input, output, dir + "/efficiency", differentialPDF, sample, "efficiency", "Efficiency of thermal muons", "Efficiency", false, variable);
      drawRate(input, output, dir + "/purity", differentialPDF, sample, "purity", "Purity for thermal muons", "Purity", false, variable);
      drawRate(input, output, dir + "/fake", differentialPDF, sample, "fakeRate", "Fake rate (not a correct thermal match)", "Fake rate", false, variable);
      if (input.GetKey(("purityRatio_" + variable).c_str())) {
        drawPurityRatio(input, output, dir + "/purity", differentialPDF, sample, variable);
        drawPurityRatio(input, output, dir + "/fake", differentialPDF, sample, variable, "fakeRate");
        drawPurityRatio(input, output, dir + "/efficiency", differentialPDF, sample, variable, "efficiency");

      }
    }
  }
}
