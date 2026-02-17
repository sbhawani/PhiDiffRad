// ============================================================
//  RunDiffrad.C  –  ROOT macro to run the DIFFRAD FORTRAN code
//                  and visualise the radiative-correction output
//
//  Program:  DIFFRAD (Akushevich et al., Eur.Phys.J.C8:457, 1999)
//  Purpose:  QED radiative corrections in diffractive vector-meson
//            electroproduction.
//
//  Usage (from the ROOT prompt or command line):
//    root -l RunDiffrad.C                    <- uses defaults
//    root -l 'RunDiffrad.C("myinput.dat")'  <- custom input file
//    root -b -q RunDiffrad.C                 <- batch / no GUI
//
//  The macro will:
//    1. Write (or overwrite) input.dat from configurable parameters.
//    2. Compile idiffrad.f with gfortran and run the executable.
//    3. Parse allu.dat (the main output file) into a TTree.
//    4. Draw a summary canvas with four physics plots.
// ============================================================

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdio>
#include <cstdlib>

// -------------------------------------------------------
//  Helper: write an input.dat file from parameters
// -------------------------------------------------------
void WriteInputDat(const char* fname,
                   double bmom,   // lepton beam momentum  [GeV]
                   double tmom,   // nucleon momentum      [GeV]
                   int    lepton, // 1=electron  2=muon
                   int    ivec,   // 1=rho 2=omega 3=phi 4=J/Psi
                   int    intphi, // 1=integrate over phi_h
                   double vcut,   // inelasticity cut (0 = none)
                   std::vector<double> xv,   // x  or -W2 per point
                   std::vector<double> yv,   // y  or -Q2 per point
                   std::vector<double> tv,   // t  (must be negative)
                   std::vector<double> phiv) // phi_h (if intphi=0)
{
    int npoi = (int)xv.size();
    std::ofstream f(fname);
    if (!f.is_open()) {
        Error("WriteInputDat", "Cannot open %s for writing", fname);
        return;
    }
    f << bmom   << "    !  bmom  - lepton momentum\n";
    f << tmom   << "    !  tmom  - momentum per nucleon\n";
    f << lepton << "    !  lepton - 1 electron, 2 muon\n";
    f << ivec   << "    !  ivec  - (1) rho, (2) omega, (3) phi, (4) J/Psi\n";
    f << intphi << "    !  intphi - integrate over phi_h (1) or not (0)\n";
    f << vcut   << "    !  vcut  - cut on v (0. = no cut)\n";
    f << npoi   << "    !  npoi  - number of kinematical points\n";

    // x (or -W2) row
    for (int i = 0; i < npoi; ++i)  f << xv[i] << (i+1<npoi ? " " : "\n");
    // y (or -Q2) row
    for (int i = 0; i < npoi; ++i)  f << yv[i] << (i+1<npoi ? " " : "\n");
    // t row
    for (int i = 0; i < npoi; ++i)  f << tv[i] << (i+1<npoi ? " " : "\n");
    // phi_h row (only if intphi == 0)
    if (intphi == 0 && !phiv.empty()) {
        for (int i = 0; i < npoi; ++i) f << phiv[i] << (i+1<npoi ? " " : "\n");
    }
    f.close();
    Printf("  [WriteInputDat] Written %s  (npoi=%d)", fname, npoi);
}

// -------------------------------------------------------
//  Helper: parse allu.dat into parallel vectors
//  Format (from WRITE statement in idiffrad.f):
//    write(21,'(2f6.3,f8.3,f8.0,f8.3,3f9.3,2g12.4)')
//      xs, ys, q2, w2, tdif, sig/sib, wei, |wei-sig/sib|, sig, sib
// -------------------------------------------------------
struct DiffradPoint {
    double xs, ys, Q2, W2, t;
    double rc_num;   // sig/sib  (full numerical RC)
    double rc_apx;   // wei      (approximate RC)
    double rc_diff;  // |numerical - approximate|
    double sig;      // observed cross section [nb]
    double sib;      // Born cross section     [nb]
};

std::vector<DiffradPoint> ParseAlluDat(const char* fname)
{
    std::vector<DiffradPoint> pts;
    std::ifstream f(fname);
    if (!f.is_open()) {
        Warning("ParseAlluDat", "Cannot open %s – did DIFFRAD run successfully?", fname);
        return pts;
    }
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == ' ' && line.size() < 5) continue;
        DiffradPoint p;
        // The FORTRAN format packs columns tightly; use sscanf for robustness
        int n = sscanf(line.c_str(),
                       "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                       &p.xs, &p.ys, &p.Q2, &p.W2, &p.t,
                       &p.rc_num, &p.rc_apx, &p.rc_diff,
                       &p.sig, &p.sib);
        if (n >= 8) pts.push_back(p);
    }
    Printf("  [ParseAlluDat] Read %zu point(s) from %s", pts.size(), fname);
    return pts;
}

// -------------------------------------------------------
//  Main macro entry point
// -------------------------------------------------------
void RunDiffrad(const char* inputFile = "input.dat",
                bool        regenerateInput = true,
                bool        recompile       = true,
                const char* fortranSrc      = "idiffrad.f",
                const char* executable      = "./diffrad_exec")
{
    Printf("=======================================================");
    Printf("  RunDiffrad  –  DIFFRAD radiative-correction wrapper  ");
    Printf("=======================================================");

    // ----------------------------------------------------------
    // 1.  Write input.dat
    //     Edit the parameters below to change the kinematic grid.
    // ----------------------------------------------------------
    if (regenerateInput) {
        Printf("\n[Step 1] Writing input file: %s", inputFile);

        double bmom   = 10.6;   // lepton beam momentum  [GeV]
        double tmom   = 0.0;   // nucleon momentum      [GeV/c]
        int    lepton = 1;       // 1 = electron
        int    ivec   = 3;       // 1 = rho meson
        int    intphi = 1;       // integrate over phi_h
        double vcut   = 0.0;     // no inelasticity cut

        // Kinematic points: negative values mean -W2 or -Q2
        std::vector<double> xv  = {-9.,  -9.,  -9.,  -9.,  -9.,  -9.,  -9.,  -9.,  -9.,  -9.}; //-W^2
        std::vector<double> yv  = {-0.2, -0.4, -0.6, -0.8, -1.0, -1.2, -1.4, -1.6, -1.8, -2.0};//-Q^2
        std::vector<double> tv  = {-0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2};//t (must be negative)
        std::vector<double> phv;  // not needed when intphi=1

        WriteInputDat(inputFile, bmom, tmom, lepton, ivec,
                      intphi, vcut, xv, yv, tv, phv);
    } else {
        Printf("\n[Step 1] Using existing input file: %s", inputFile);
    }

    // ----------------------------------------------------------
    // 2.  Compile FORTRAN source with gfortran
    // ----------------------------------------------------------
    if (recompile) {
        Printf("\n[Step 2] Compiling %s ...", fortranSrc);
        TString compileCmd = TString::Format(
            "gfortran -O2 -o %s %s -lm 2>&1", executable, fortranSrc);
        Printf("  cmd: %s", compileCmd.Data());
        int ret = gSystem->Exec(compileCmd.Data());
        if (ret != 0) {
            Error("RunDiffrad",
                  "Compilation failed (exit %d). "
                  "Make sure gfortran is installed and %s is present.",
                  ret, fortranSrc);
            return;
        }
        Printf("  Compilation successful.");
    } else {
        Printf("\n[Step 2] Skipping recompilation – using %s", executable);
    }

    // ----------------------------------------------------------
    // 3.  Run the DIFFRAD executable
    // ----------------------------------------------------------
    Printf("\n[Step 3] Running DIFFRAD ...");
    TString runCmd = TString::Format("%s 2>&1", executable);
    Printf("  cmd: %s", runCmd.Data());
    int ret = gSystem->Exec(runCmd.Data());
    if (ret != 0) {
        Warning("RunDiffrad",
                "DIFFRAD returned non-zero exit code %d. "
                "Output may still be available in allu.dat.", ret);
    } else {
        Printf("  DIFFRAD finished successfully.");
    }

    // ----------------------------------------------------------
    // 4.  Parse output file allu.dat
    // ----------------------------------------------------------
    Printf("\n[Step 4] Parsing output file allu.dat ...");
    std::vector<DiffradPoint> pts = ParseAlluDat("allu.dat");
    if (pts.empty()) {
        Error("RunDiffrad", "No data points found in allu.dat. Aborting.");
        return;
    }

    int N = (int)pts.size();

    // ----------------------------------------------------------
    // 5.  Fill a TTree (makes data available for further analysis)
    // ----------------------------------------------------------
    Printf("\n[Step 5] Filling TTree ...");
    Double_t t_xs, t_ys, t_Q2, t_W2, t_t;
    Double_t t_rc_num, t_rc_apx, t_rc_diff, t_sig, t_sib;

    TTree* tree = new TTree("diffrad", "DIFFRAD output");
    tree->Branch("xs",      &t_xs,      "xs/D");
    tree->Branch("ys",      &t_ys,      "ys/D");
    tree->Branch("Q2",      &t_Q2,      "Q2/D");
    tree->Branch("W2",      &t_W2,      "W2/D");
    tree->Branch("t",       &t_t,       "t/D");
    tree->Branch("rc_num",  &t_rc_num,  "rc_num/D");   // full RC  sig/sib
    tree->Branch("rc_apx",  &t_rc_apx,  "rc_apx/D");   // approx  RC
    tree->Branch("rc_diff", &t_rc_diff, "rc_diff/D");   // |num - apx|
    tree->Branch("sig",     &t_sig,     "sig/D");       // observed [nb]
    tree->Branch("sib",     &t_sib,     "sib/D");       // Born     [nb]

    for (auto& p : pts) {
        t_xs = p.xs;  t_ys = p.ys;  t_Q2 = p.Q2;
        t_W2 = p.W2;  t_t  = p.t;
        t_rc_num  = p.rc_num;
        t_rc_apx  = p.rc_apx;
        t_rc_diff = p.rc_diff;
        t_sig = p.sig;  t_sib = p.sib;
        tree->Fill();
    }
    Printf("  TTree filled with %d entries.", N);

    // ----------------------------------------------------------
    // 6.  Draw summary canvas
    // ----------------------------------------------------------
    Printf("\n[Step 6] Drawing summary plots ...");

    // Arrays for manual graph filling
    std::vector<Double_t> vW   (N), vQ2   (N);
    std::vector<Double_t> vRCn (N), vRCa  (N);
    std::vector<Double_t> vSig (N), vSib  (N);
    std::vector<Double_t> vDiff(N), vidx  (N);

    for (int i = 0; i < N; ++i) {
        auto& p = pts[i];
        vW[i]    = TMath::Sqrt(p.W2);
        vQ2[i]   = p.Q2;
        vRCn[i]  = p.rc_num;
        vRCa[i]  = p.rc_apx;
        vSig[i]  = p.sig;
        vSib[i]  = p.sib;
        vDiff[i] = p.rc_diff;
        vidx[i]  = (Double_t)(i + 1);
    }

    // ------ canvas ------
    TCanvas* c = new TCanvas("cDiffrad",
                             "DIFFRAD – Radiative corrections summary",
                             1200, 900);
    c->Divide(2, 2);
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.06);

    // ---- Pad 1: RC factor (numerical vs approximate) vs Q2 ----
    c->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.13);

    TGraph* grNum = new TGraph(N, vQ2.data(), vRCn.data());
    grNum->SetTitle("RC factor vs Q^{2};Q^{2} [GeV^{2}];RC factor  #sigma/#sigma_{Born}");
    grNum->SetMarkerStyle(kFullCircle);
    grNum->SetMarkerColor(kBlue+1);
    grNum->SetLineColor(kBlue+1);
    grNum->SetLineWidth(2);
    grNum->SetMarkerSize(1.2);

    TGraph* grApx = new TGraph(N, vQ2.data(), vRCa.data());
    grApx->SetMarkerStyle(kOpenSquare);
    grApx->SetMarkerColor(kRed+1);
    grApx->SetLineColor(kRed+1);
    grApx->SetLineWidth(2);
    grApx->SetLineStyle(2);
    grApx->SetMarkerSize(1.2);

    grNum->Draw("APL");
    grApx->Draw("PL same");

    TLegend* leg1 = new TLegend(0.55, 0.75, 0.88, 0.90);
    leg1->SetBorderSize(0);
    leg1->AddEntry(grNum, "Full numerical", "pl");
    leg1->AddEntry(grApx, "Approximate",    "pl");
    leg1->Draw();

    // ---- Pad 2: Cross sections vs W ----
    c->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.13);
    gPad->SetLogy();

    TGraph* grSig = new TGraph(N, vW.data(), vSig.data());
    grSig->SetTitle("Cross sections vs W;W [GeV];#sigma [nb]");
    grSig->SetMarkerStyle(kFullCircle);
    grSig->SetMarkerColor(kBlue+1);
    grSig->SetLineColor(kBlue+1);
    grSig->SetLineWidth(2);
    grSig->SetMarkerSize(1.2);

    TGraph* grSib = new TGraph(N, vW.data(), vSib.data());
    grSib->SetMarkerStyle(kOpenSquare);
    grSib->SetMarkerColor(kRed+1);
    grSib->SetLineColor(kRed+1);
    grSib->SetLineWidth(2);
    grSib->SetLineStyle(2);
    grSib->SetMarkerSize(1.2);

    grSig->Draw("APL");
    grSib->Draw("PL same");

    TLegend* leg2 = new TLegend(0.55, 0.75, 0.88, 0.90);
    leg2->SetBorderSize(0);
    leg2->AddEntry(grSig, "#sigma (observed)", "pl");
    leg2->AddEntry(grSib, "#sigma_{Born}",     "pl");
    leg2->Draw();

    // ---- Pad 3: |RC_num - RC_apx| vs point index ----
    c->cd(3);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.13);

    TGraph* grDiff = new TGraph(N, vidx.data(), vDiff.data());
    grDiff->SetTitle("|#DeltaRC| per kinematic point;"
                     "Point index;|#sigma/#sigma_{Born} - wei|");
    grDiff->SetMarkerStyle(kFullTriangleUp);
    grDiff->SetMarkerColor(kGreen+2);
    grDiff->SetLineColor(kGreen+2);
    grDiff->SetLineWidth(2);
    grDiff->SetMarkerSize(1.3);
    grDiff->Draw("APL");

    // draw zero line
    Double_t xlo = 0.5, xhi = N + 0.5;
    TLine* zeroline = new TLine(xlo, 0., xhi, 0.);
    zeroline->SetLineStyle(3);
    zeroline->SetLineColor(kGray+1);
    zeroline->Draw();

    // ---- Pad 4: RC factor vs W (colour-coded by Q2) ----
    c->cd(4);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.13);

    // Use a multi-graph so each point can have its own colour
    TMultiGraph* mg = new TMultiGraph();
    mg->SetTitle("RC factor vs W;W [GeV];RC factor #sigma/#sigma_{Born}");

    // group by unique Q2 value (within 0.1 GeV^2 tolerance)
    std::vector<int>    colIdx;
    std::vector<double> uniqueQ2;
    int palette[] = {kBlue+1, kRed+1, kGreen+2, kOrange+1,
                     kMagenta+1, kCyan+2, kViolet+2, kPink+1};
    int nPal = 8;

    for (int i = 0; i < N; ++i) {
        int found = -1;
        for (int j = 0; j < (int)uniqueQ2.size(); ++j) {
            if (TMath::Abs(vQ2[i] - uniqueQ2[j]) < 0.1) { found = j; break; }
        }
        if (found < 0) {
            found = (int)uniqueQ2.size();
            uniqueQ2.push_back(vQ2[i]);
        }
        colIdx.push_back(found);
    }

    TLegend* leg4 = new TLegend(0.58, 0.60, 0.92, 0.90);
    leg4->SetBorderSize(0);
    leg4->SetHeader("Q^{2} [GeV^{2}]", "C");

    for (int j = 0; j < (int)uniqueQ2.size(); ++j) {
        std::vector<Double_t> wx, wy;
        for (int i = 0; i < N; ++i) {
            if (colIdx[i] == j) { wx.push_back(vW[i]); wy.push_back(vRCn[i]); }
        }
        if (wx.empty()) continue;
        TGraph* g = new TGraph((int)wx.size(), wx.data(), wy.data());
        int col = palette[j % nPal];
        g->SetMarkerStyle(kFullCircle);
        g->SetMarkerColor(col);
        g->SetLineColor(col);
        g->SetLineWidth(2);
        g->SetMarkerSize(1.2);
        mg->Add(g, "PL");
        leg4->AddEntry(g, TString::Format("Q^{2}=%.1f", uniqueQ2[j]), "pl");
    }
    mg->Draw("A");
    leg4->Draw();

    c->Update();

    // ----------------------------------------------------------
    // 7.  Save canvas and TTree
    // ----------------------------------------------------------
    Printf("\n[Step 7] Saving outputs ...");
    c->SaveAs("diffrad_summary.png");
    Printf("  Canvas saved to diffrad_summary.png");

    TFile* fout = new TFile("diffrad_output.root", "RECREATE");
    tree->Write();
    c->Write();
    fout->Close();
    Printf("  TTree and canvas written to diffrad_output.root");

    // ----------------------------------------------------------
    // 8.  Print a text summary table to the terminal
    // ----------------------------------------------------------
    Printf("\n-------------------------------------------------------");
    Printf("  DIFFRAD Results Summary");
    Printf("-------------------------------------------------------");
    Printf("  %-6s %-6s %-8s %-8s %-6s | %-8s %-8s %-8s",
           "xs", "ys", "Q2", "W2", "t",
           "RC_num", "RC_apx", "|diff|");
    Printf("  %s", TString('-', 72).Data());
    for (auto& p : pts) {
        Printf("  %-6.3f %-6.3f %-8.3f %-8.1f %-6.3f | %-8.4f %-8.4f %-8.4f",
               p.xs, p.ys, p.Q2, p.W2, p.t,
               p.rc_num, p.rc_apx, p.rc_diff);
    }
    Printf("-------------------------------------------------------\n");
    Printf("[RunDiffrad] Done.\n");
}
