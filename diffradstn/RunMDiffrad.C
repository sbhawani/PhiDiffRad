// ============================================================
//  RunMDiffrad.C  –  ROOT macro to run the MDIFFRAD FORTRAN code
//                   and produce a CSV summary:
//                   xB, Q, W, t, rad_corr, rad_corr_err
//
//  Uses:
//    - Input  : inmdi.dat  (copied to INMDI.DAT for Fortran)
//    - Output : ALLmc.dat (per-sample). etamc.dat optional.
//
//  NOTE:
//   - Assumes output "main ..." lines appear in order:
//       point1 sample1..sampleNev, point2 sample1..sampleNev, ...
//   - rad_corr     = mean of de over samples for that point
//   - rad_corr_err = sample standard deviation (sigma) of de over samples
// ============================================================

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"

static constexpr double Mp = 0.9382720813;

// -------------------------------
// helper: split numeric tokens
// -------------------------------
static std::vector<double> ParseNumbersFromLine(const std::string& line) {
  std::vector<double> out;
  std::istringstream iss(line);
  std::string tok;
  while (iss >> tok) {
    char* endptr = nullptr;
    double v = std::strtod(tok.c_str(), &endptr);
    if (endptr != tok.c_str() && *endptr == '\0') out.push_back(v);
  }
  return out;
}

// -------------------------------------------------------
//  Helper: write inmdi.dat for MDIFFRAD (MC version)
// -------------------------------------------------------
void WriteInmdiDat(const char* fname,
                   double bmom,
                   double tmom,
                   int lepton,
                   int ivec,
                   double ann1,
                   double ann2,
                   double ann3,
                   double vcut,
                   int nev,
                   long seed,
                   const std::vector<double>& W2min,
                   const std::vector<double>& W2max,
                   const std::vector<double>& Q2min,
                   const std::vector<double>& Q2max,
                   const std::vector<double>& tmin,
                   const std::vector<double>& tmax)
{
  int npoi = (int)W2min.size();
  auto check = [&](const std::vector<double>& v, const char* name){
    if ((int)v.size() != npoi) {
      Error("WriteInmdiDat", "Vector %s has size %d but expected %d", name, (int)v.size(), npoi);
      return false;
    }
    return true;
  };
  if (!check(W2max,"W2max") || !check(Q2min,"Q2min") || !check(Q2max,"Q2max") ||
      !check(tmin,"tmin")   || !check(tmax,"tmax")) return;

  std::ofstream f(fname);
  if (!f.is_open()) {
    Error("WriteInmdiDat", "Cannot open %s for writing", fname);
    return;
  }

  f << bmom << "     !  bmom - lepton momentum\n";
  f << tmom << "     !  tmom - proton/nucleon momentum\n";
  f << lepton << "         !  lepton - 1 electron, 2 muon\n";
  f << ivec   << "         !  ivec -(1) rho, (2) omega, (3) phi, (4) J/Psi\n";
  f << ann1 << " " << ann2 << " " << ann3 << "   ! ann1 ann2 ann3 -- numbers of events\n";
  f << vcut << "        !  vcut - cut on v (0 = none)\n";
  f << nev  << "         !  nev - number of samples for each point\n";
  f << seed << "    !  random seed\n";
  f << npoi << "         !  npoi - number of kinematical point\n";

  auto writeRow = [&](const std::vector<double>& v){
    for (int i=0;i<npoi;i++) f << v[i] << (i+1<npoi ? " " : "\n");
  };

  // README convention: lines are "-W2_min", "-W2_max", "-Q2_min", "-Q2_max", "-t_min", "-t_max"
  writeRow(W2min);
  writeRow(W2max);
  writeRow(Q2min);
  writeRow(Q2max);
  writeRow(tmin);
  writeRow(tmax);

  f.close();
  Printf("  [WriteInmdiDat] Written %s  (npoi=%d)", fname, npoi);
}

// -------------------------------------------------------
// Parse ALLmc.dat (per sample).
// We treat nums as:
//   [born, dev, der, de, mean, err]  (common)
// -------------------------------------------------------
struct AllMcRow {
  double born = NAN;
  double dev  = NAN;
  double der  = NAN;
  double de   = NAN;  // total RC factor for that sample
  double mean = NAN;
  double err  = NAN;
};

static std::vector<AllMcRow> ParseAllMc(const char* fname)
{
  std::vector<AllMcRow> rows;
  std::ifstream f(fname);
  if (!f.is_open()) {
    Warning("ParseAllMc", "Cannot open %s", fname);
    return rows;
  }
  std::string line;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    auto nums = ParseNumbersFromLine(line);
    if ((int)nums.size() < 4) continue;

    AllMcRow r;
    r.born = nums[0];
    r.dev  = nums.size() > 1 ? nums[1] : NAN;
    r.der  = nums.size() > 2 ? nums[2] : NAN;
    r.de   = nums.size() > 3 ? nums[3] : NAN;

    if (nums.size() == 6) {
      r.mean = nums[4];
      r.err  = nums[5];
    } else if (nums.size() >= 7) {
      r.mean = nums[nums.size()-2];
      r.err  = nums[nums.size()-1];
    }
    rows.push_back(r);
  }
  Printf("  [ParseAllMc] Read %zu row(s) from %s", rows.size(), fname);
  return rows;
}

// -------------------------------------------------------
// Utilities: mean/std (ignoring NaN)
// -------------------------------------------------------
static bool IsFinite(double x){ return std::isfinite(x); }

static double MeanFinite(const std::vector<double>& v){
  double s=0; int n=0;
  for (auto x: v){ if (IsFinite(x)){ s+=x; n++; } }
  return (n>0)? s/n : NAN;
}

static double StdFiniteSample(const std::vector<double>& v){
  // sample std (N-1 in denom), ignoring NaN
  std::vector<double> u;
  u.reserve(v.size());
  for (auto x: v) if (IsFinite(x)) u.push_back(x);
  int n = (int)u.size();
  if (n <= 1) return NAN;
  double m = 0;
  for (auto x: u) m += x;
  m /= n;
  double ss = 0;
  for (auto x: u) ss += (x-m)*(x-m);
  return std::sqrt(ss / (n-1));
}

// -------------------------------------------------------
// Convert (W2,Q2) -> xB using: W2 = Mp^2 + Q2*(1/xB - 1)
// => xB = Q2 / (W2 - Mp^2 + Q2)
// -------------------------------------------------------
static double ComputeXB(double W2, double Q2){
  double denom = (W2 - Mp*Mp + Q2);
  if (!IsFinite(denom) || std::fabs(denom) < 1e-12) return NAN;
  return Q2 / denom;
}

// -------------------------------------------------------
// Main entry
// -------------------------------------------------------
void RunMDiffrad(const char* inFile     = "inmdi.dat",
                 bool        regenInput = true,
                 bool        recompile  = true,
                 const char* fortranSrc = "mdiffrad.f",
                 const char* executable = "./mdiffrad_exec",
                 const char* outCsv     = "diffrad_rc_summary.csv")
{
  Printf("=======================================================");
  Printf("  RunMDiffrad  –  MDIFFRAD (Monte Carlo) wrapper        ");
  Printf("=======================================================");

  // -----------------------
  // 1) Write inmdi.dat
  // -----------------------
  // Keep these vectors accessible for CSV (we need to know kinematics).
  // We'll define them here and only fill if regenInput=true.
  double bmom = 10.60;
  double tmom = 0.0;
  int lepton  = 1;
  int ivec    = 3;

  double ann1 = 1e6;
  double ann2 = 1e6;
  double ann3 = 5e5;
  double vcut = 0.0;
  int    nev  = 3;
  long   seed = 333522;

  std::vector<double> W2min, W2max, Q2min, Q2max, tmin, tmax;

  if (regenInput) {
    Printf("\n[Step 1] Writing input file: %s", inFile);

    // Example: 10 bins
    W2min = {-12.25, -12.25, -12.25, -12.25, -12.25, -12.25, -12.25, -12.25, -12.25, -12.25};
    W2max = { -6.25,  -6.25,  -6.25,  -6.25,  -6.25,  -6.25,  -6.25,  -6.25,  -6.25,  -6.25};
    Q2min = { -0.25,  -0.45,  -0.65,  -0.85,  -1.05,  -1.25,  -1.45,  -1.65,  -1.85,  -2.05};
    Q2max = { -0.15,  -0.35,  -0.55,  -0.75,  -0.95,  -1.15,  -1.35,  -1.55,  -1.75,  -1.95};
    tmin  = { -0.25,  -0.25,  -0.25,  -0.25,  -0.25,  -0.25,  -0.25,  -0.25,  -0.25,  -0.25};
    tmax  = { -0.15,  -0.15,  -0.15,  -0.15,  -0.15,  -0.15,  -0.15,  -0.15,  -0.15,  -0.15};

    WriteInmdiDat(inFile, bmom, tmom, lepton, ivec,
                  ann1, ann2, ann3, vcut, nev, seed,
                  W2min, W2max, Q2min, Q2max, tmin, tmax);
  } else {
    Printf("\n[Step 1] Using existing input file: %s", inFile);
    // If you set regenInput=false, and still want CSV with correct kinematics,
    // then you must also load the kinematics from inFile here.
    // For simplicity, we keep regenInput=true in your workflow.
  }

  // -----------------------
  // 2) Compile Fortran
  // -----------------------
  if (recompile) {
    Printf("\n[Step 2] Compiling %s ...", fortranSrc);
    TString cmd = TString::Format(
      "gfortran -O0 -g -std=legacy -ffixed-line-length-none "
      "-fbacktrace -ffpe-trap=invalid,overflow,zero -o %s %s -lm 2>&1",
      executable, fortranSrc);
    Printf("  cmd: %s", cmd.Data());
    int ret = gSystem->Exec(cmd.Data());
    if (ret != 0) {
      Error("RunMDiffrad", "Compilation failed (exit %d). Check %s", ret, fortranSrc);
      return;
    }
    Printf("  Compilation successful.");
  }

  // -----------------------
  // 3) Run executable
  // -----------------------
  Printf("\n[Step 3] Running MDIFFRAD ...");
  gSystem->Exec(TString::Format("cp -f %s INMDI.DAT", inFile));
  int retRun = gSystem->Exec(TString::Format("%s 2>&1", executable));
  if (retRun != 0) {
    Warning("RunMDiffrad", "MDIFFRAD returned exit code %d. Will still try to parse outputs.", retRun);
  } else {
    Printf("  MDIFFRAD finished successfully.");
  }

  // -----------------------
  // 4) Parse ALLmc.dat
  // -----------------------
  Printf("\n[Step 4] Parsing outputs ALLmc.dat ...");
  auto allmc = ParseAllMc("ALLmc.dat");
  if (allmc.empty()) {
    Error("RunMDiffrad", "No rows found in ALLmc.dat. Aborting.");
    return;
  }

  // -----------------------
  // 5) Group by point and compute RC mean/sigma
  // -----------------------
  int npoi = (int)W2min.size();
  if (npoi <= 0) {
    Error("RunMDiffrad", "npoi is 0 (no kinematics vectors). Keep regenInput=true or implement inFile reading.");
    return;
  }

  // Expect allmc.size() ~ npoi*nev (often exactly)
  // We'll group sequentially in chunks of nev.
  if ((int)allmc.size() < npoi * nev) {
    Warning("RunMDiffrad",
            "ALLmc has %zu rows but expected at least npoi*nev=%d. Will group as much as possible.",
            allmc.size(), npoi*nev);
  }

  struct RcSummary {
    double xB= NAN, Q=NAN, W=NAN, t=NAN;
    double rc_mean=NAN;
    double rc_sigma=NAN;     // "single-point multiple-sample sigma"
    double rc_sem=NAN;       // sigma/sqrt(N)
    int    nsamp=0;
  };

  std::vector<RcSummary> summary;
  summary.reserve(npoi);

  for (int i=0; i<npoi; i++) {
    // pull this point's samples
    std::vector<double> de_samples;
    de_samples.reserve(nev);

    for (int s=0; s<nev; s++) {
      int idx = i*nev + s;
      if (idx >= (int)allmc.size()) break;
      de_samples.push_back(allmc[idx].de);
    }

    RcSummary r;
    r.nsamp = (int)de_samples.size();

    // Compute kinematic centers from the bin ranges.
    // W2/Q2 lines are stored as NEGATIVE of the limits in README convention.
    // So physical W2 range is [-W2min, -W2max] (then take min/max).
    auto w2a = -W2min[i];
    auto w2b = -W2max[i];
    double W2lo = std::min(w2a, w2b);
    double W2hi = std::max(w2a, w2b);
    double W2c  = 0.5*(W2lo + W2hi);

    auto q2a = -Q2min[i];
    auto q2b = -Q2max[i];
    double Q2lo = std::min(q2a, q2b);
    double Q2hi = std::max(q2a, q2b);
    double Q2c  = 0.5*(Q2lo + Q2hi);

    // t is negative already; treat (tmin,tmax) as the actual t-range
    double tc = 0.5*(tmin[i] + tmax[i]);

    r.Q = std::sqrt(std::max(0.0, Q2c));
    r.W = std::sqrt(std::max(0.0, W2c));
    r.t = tc;
    r.xB = ComputeXB(W2c, Q2c);

    r.rc_mean  = MeanFinite(de_samples);
    r.rc_sigma = StdFiniteSample(de_samples);                 // <-- your requested "sigma"
    if (IsFinite(r.rc_sigma) && r.nsamp>0) r.rc_sem = r.rc_sigma/std::sqrt((double)r.nsamp);

    summary.push_back(r);
  }

  // -----------------------
  // 6) Write CSV
  // -----------------------
  Printf("\n[Step 5] Writing CSV: %s", outCsv);
  std::ofstream csv(outCsv);
  if (!csv.is_open()) {
    Error("RunMDiffrad", "Cannot open %s for writing", outCsv);
    return;
  }

  csv << "xB,Q,W,t,rad_corr,rad_corr_err\n";
  for (auto& r : summary) {
    // rad_corr_err per your requirement: sigma over samples (not SEM)
    csv << r.xB << ","
        << r.Q  << ","
        << r.W  << ","
        << r.t  << ","
        << r.rc_mean  << ","
        << r.rc_sigma << "\n";
  }
  csv.close();
  Printf("  CSV written: %s", outCsv);

  // --------------------------
  // 6) Plots with ERROR BARS from our computed per-point sigma
  // ----------------------------------------------------------
  Printf("\n[Step 6] Drawing summary plots (with error bars) ...");
  gStyle->SetOptStat(0);

  int M = (int)summary.size();
  std::vector<double> xIdx(M), yRC(M), ex0(M,0.0), eySig(M);
  std::vector<double> xQ2(M), yRC2(M), exQ2(M,0.0), eySig2(M);

  for (int i=0;i<M;i++){
    xIdx[i]  = i+1;
    yRC[i]   = summary[i].rc_mean;
    eySig[i] = summary[i].rc_sigma;   // <-- your rad_corr_err (sigma)

    double Q2c = summary[i].Q * summary[i].Q;
    xQ2[i]   = Q2c;
    yRC2[i]  = summary[i].rc_mean;
    eySig2[i]= summary[i].rc_sigma;
  }

  TCanvas* c = new TCanvas("cMDiffrad", "MDIFFRAD summary (error bars)", 1200, 600);
  c->Divide(2,1);

  // Left: RC vs point index with sigma error bars
  c->cd(1);
  {
    auto gr = new TGraphErrors(M, xIdx.data(), yRC.data(), ex0.data(), eySig.data());
    gr->SetTitle("MDIFFRAD: RC mean per point (error = sigma over samples);point index;RC (mean de)");
    gr->SetMarkerStyle(20);
    gr->Draw("AP");
  }

  // Right: RC vs Q2 with sigma error bars
  c->cd(2);
  {
    auto gr = new TGraphErrors(M, xQ2.data(), yRC2.data(), exQ2.data(), eySig2.data());
    gr->SetTitle("MDIFFRAD: RC mean vs Q^{2} center (error = sigma);Q^{2} [GeV^{2}];RC (mean de)");
    gr->SetMarkerStyle(20);
    gr->Draw("AP");
  }
  c->Update();
  c->SaveAs("mdiffrad_summary.png");

  // -----------------------
  // 7) (Optional) Save ROOT file too
  // -----------------------
  Printf("\n[Step 6] Writing ROOT output (optional) ...");
  TFile* fout = new TFile("mdiffrad_output.root","RECREATE");

  // Store per-sample tree
  TTree* tAll = new TTree("allmc", "MDIFFRAD per-sample output");
  Double_t born, dev, der, de, mean, err;
  tAll->Branch("born", &born, "born/D");
  tAll->Branch("dev",  &dev,  "dev/D");
  tAll->Branch("der",  &der,  "der/D");
  tAll->Branch("de",   &de,   "de/D");
  tAll->Branch("mean", &mean, "mean/D");
  tAll->Branch("err",  &err,  "err/D");

  for (auto& rr : allmc) {
    born = rr.born; dev = rr.dev; der = rr.der; de = rr.de; mean = rr.mean; err = rr.err;
    tAll->Fill();
  }

  // Store summary tree
  TTree* tSum = new TTree("rcsum", "RC summary per kinematic point");
  Double_t xB, Q, W, tt, rc, rc_sigma, rc_sem;
  Int_t nsamp;
  tSum->Branch("xB", &xB, "xB/D");
  tSum->Branch("Q",  &Q,  "Q/D");
  tSum->Branch("W",  &W,  "W/D");
  tSum->Branch("t",  &tt, "t/D");
  tSum->Branch("rad_corr",      &rc,       "rad_corr/D");
  tSum->Branch("rad_corr_sigma",&rc_sigma, "rad_corr_sigma/D");
  tSum->Branch("rad_corr_sem",  &rc_sem,   "rad_corr_sem/D");
  tSum->Branch("nsamp", &nsamp, "nsamp/I");

  for (auto& r : summary) {
    xB = r.xB; Q=r.Q; W=r.W; tt=r.t;
    rc=r.rc_mean; rc_sigma=r.rc_sigma; rc_sem=r.rc_sem;
    nsamp=r.nsamp;
    tSum->Fill();
  }

  tAll->Write();
  tSum->Write();
  fout->Close();

  Printf("[RunMDiffrad] Done.\n");
}
