// ============================================================
//  RunMDiffrad.C
//  Drop-in replacement for RunMDiffrad_backup.C
//
//  What this macro does (end-to-end):
//    1. Expands your Q² × t' bin edges into flat kinematic points
//    2. Converts each t'-bin into a physical |t| range using |tmin(Q²)|
//    3. Writes MDIFFRAD's inmdi.dat (FORTRAN sign convention)
//    4. Compiles + runs the FORTRAN executable
//    5. Parses ALLmc.dat → RC mean ± σ per bin
//    6. Writes two CSVs:
//         a) diffrad_rc_summary.csv        – RC results per bin
//         b) diffrad_bins_metadata.csv     – full bin geometry for RDF
//    7. Produces ROOT publication-quality plots:
//         – δ_RC vs t'  (one curve per Q² bin, coloured)
//         – δ_RC vs |t| (one curve per Q² bin, coloured)
//         – 2D δ_RC(Q², t') heatmap
//         – δ_RC vs point index
//    8. Writes mdiffrad_output.root (per-sample tree + summary tree)
//
//  Column conventions (from DISANAMMUtils.h / DISANAMath):
//    t        → NEGATIVE  (GetT)
//    tmin     → NEGATIVE  (GetTmin)
//    mtprime  → POSITIVE  = |t| − |tmin| = |t + tmin|
//    tprime   → NEGATIVE  = −mtprime
//    |t|      → POSITIVE  = mtprime + |tmin|
//
//  FORTRAN inmdi.dat sign convention:
//    all kinematic limits stored as their NEGATIVE physical values
//    (i.e. −W², −Q², −|t|)
//
//  Usage:
//    1. Edit the "── Edit only this section ──" block at the top of this file
//       to set kFortranSrc, kExecutable, kInmdiDat, kBinningCsv, etc.
//    2. Run:
//         root -b -q './RunMDiffradNew.C++'
//
//    Use kRecompile = false after the first successful compile to save time.
//    Use kRegenInput = false to reuse an existing inmdi.dat unchanged.
// ============================================================

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <future>
#include <iomanip>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TPad.h"
#include "TAxis.h"
#include "TExec.h"

// ============================================================
//  ██████╗ ██╗███╗   ██╗███╗   ██╗██╗███╗   ██╗ ██████╗
//  ██╔══██╗██║████╗  ██║████╗  ██║██║████╗  ██║██╔════╝
//  ██████╔╝██║██╔██╗ ██║██╔██╗ ██║██║██╔██╗ ██║██║  ███╗
//  ██╔══██╗██║██║╚██╗██║██║╚██╗██║██║██║╚██╗██║██║   ██║
//  ██████╔╝██║██║ ╚████║██║ ╚████║██║██║ ╚████║╚██████╔╝
//  ╚═════╝ ╚═╝╚═╝  ╚═══╝╚═╝  ╚═══╝╚═╝╚═╝  ╚═══╝ ╚═════╝
//
//  ── Edit only this section ───────────────────────────────
// ============================================================

// ── Beam / meson ─────────────────────────────────────────────────────────────
static constexpr double kBmom   = 10.60;   // beam momentum (GeV)
static constexpr double kTmom   = 0.0;     // target momentum (0 = fixed target)
static constexpr int    kLepton = 1;       // 1=electron, 2=muon
static constexpr int    kIvec   = 3;       // 1=rho, 2=omega, 3=phi, 4=J/psi
static constexpr double kMphi   = 1.0195;  // phi meson mass (GeV) — used for |tmin|

// ── MDIFFRAD Monte Carlo settings ────────────────────────────────────────────
static constexpr double kAnn1 = 1e6;
static constexpr double kAnn2 = 1e6;
static constexpr double kAnn3 = 5e5;
static constexpr double kVcut = 0.02;
static constexpr int    kNev  = 3;        // MC samples per kinematic point
static constexpr long   kSeed = 333522;

// ── Parallelism ───────────────────────────────────────────────────────────────
// Max number of MDIFFRAD subprocesses to run simultaneously.
// Set to the number of physical cores available on your node.
// Each job uses ~1 CPU core and negligible memory.
static constexpr int    kMaxJobs = 36;    // e.g. 36 for a 36-core JLab node

// ── W² range (fixed for all bins) ────────────────────────────────────────────
static constexpr double kW2lo = 4.00;    // W > 2.0 GeV
static constexpr double kW2hi = 12.25;   // W < 3.5 GeV

// ── Q² bin edges (5 edges → 4 bins) ─────────────────────────────────────────
static const std::vector<double> kQ2Edges = {
    0.1, 1.0667, 1.4667, 2.0000, 8.0000
};

// ── t' (mtprime = |t|−|tmin|) bin edges (9 edges → 8 bins) ──────────────────
static const std::vector<double> kTpEdges = {
    0.001, 0.3000, 0.3750, 0.4500, 0.5250,
    0.6000, 0.6750, 0.8250, 4.5000
};

// ── Paths ─────────────────────────────────────────────────────────────────────
// Set these to match your directory layout before running.
static const char* kFortranSrc    = "./mdiffrad.f";       // FORTRAN source
static const char* kExecutable    = "./mdiffrad_exec";    // compiled binary
static const char* kBinningCsv    = "";                   // CSV from DISANA_PhiAnalysisPlotter
                                                          // (leave "" to use hard-coded edges above)

// ── Output directory ──────────────────────────────────────────────────────────
// All outputs (CSV, ROOT, PNGs, work/ subdirs, inmdi.dat) go here.
// The directory is created automatically if it does not exist.
// Relative paths are resolved from the directory where you run root.
static const char* kOutDir = "./outputs";

// ── Run control ───────────────────────────────────────────────────────────────
static constexpr bool kRegenInput = true;   // write a fresh inmdi.dat each run
static constexpr bool kRecompile  = true;   // recompile FORTRAN each run
                                            // (set false after first successful compile)

// ============================================================
//  END OF USER CONFIG — do not edit below this line
// ============================================================

// ── Physics constants ─────────────────────────────────────────────────────────
static constexpr double kMp = 0.9382720813;

// ── Derived: s = Mp² + 2·Mp·Ebeam ────────────────────────────────────────────
// ── |tmin|(Q²) at LO: |tmin| = (Q²+Mv²)² / (4s) ────────────────────────────
static double ComputeS(double beamMom = kBmom)  { return kMp * kMp + 2.0 * kMp * beamMom; }
static double TminAbs(double Q2c, double beamMom = kBmom) {
    double num = Q2c + kMphi * kMphi;
    return (num * num) / (4.0 * ComputeS(beamMom));
}

// ── xB from (W², Q²) ─────────────────────────────────────────────────────────
static double ComputeXB(double W2, double Q2) {
    double d = W2 - kMp * kMp + Q2;
    return (std::fabs(d) > 1e-12) ? Q2 / d : NAN;
}

static bool IsFinite(double x) { return std::isfinite(x); }

// ── Robust mean / sample-std ignoring NaN ────────────────────────────────────
static double MeanFinite(const std::vector<double>& v) {
    double s = 0; int n = 0;
    for (auto x : v) if (IsFinite(x)) { s += x; ++n; }
    return n > 0 ? s / n : NAN;
}
static double StdFiniteSample(const std::vector<double>& v) {
    std::vector<double> u;
    for (auto x : v) if (IsFinite(x)) u.push_back(x);
    int n = (int)u.size();
    if (n <= 1) return NAN;
    double m = MeanFinite(u);
    double ss = 0;
    for (auto x : u) ss += (x - m) * (x - m);
    return std::sqrt(ss / (n - 1));
}

// ─────────────────────────────────────────────────────────────────────────────
//  Bin descriptor — one entry per (iQ2, itp) cell
// ─────────────────────────────────────────────────────────────────────────────
struct BinDesc {
    int    iQ2, itp;
    // Q² bin
    double Q2lo, Q2hi, Q2c;
    // t' (mtprime) bin
    double tplo, tphi, tpc;
    // |tmin| for this Q² centre
    double tmin_abs;
    // |t| = t' + |tmin|
    double t_lo_abs, t_hi_abs, t_c_abs;
    // kinematics
    double W2c, xBc, Wc;
    // FORTRAN inmdi.dat values (all negative)
    double fort_W2min, fort_W2max;
    double fort_Q2min, fort_Q2max;
    double fort_tmin,  fort_tmax;
};

// ─────────────────────────────────────────────────────────────────────────────
//  CSV binning reader  (reads CSV produced by DISANA_PhiAnalysisPlotter)
//
//  Format — pure CSV, key is first field:
//    beam_momentum,<value>
//    W2_range,<W2lo>,<W2hi>
//    Q2_edges,e0,e1,...,eN
//    tprime_edges,e0,e1,...,eM
//  Comment lines starting with # are skipped.
//  tmin / |t| are NOT stored — BuildBins derives them internally.
// ─────────────────────────────────────────────────────────────────────────────
struct CsvBinningResult {
    std::vector<double> q2Edges;
    std::vector<double> tpEdges;
    double W2lo    = kW2lo;
    double W2hi    = kW2hi;
    double beamMom = kBmom;
    bool   valid   = false;
};

// Split a CSV line into fields, stripping whitespace from each token.
static std::vector<std::string> SplitCSV(const std::string& line) {
    std::vector<std::string> fields;
    std::istringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        // strip leading/trailing whitespace and CR
        auto b = tok.find_first_not_of(" \t\r");
        auto e = tok.find_last_not_of(" \t\r");
        fields.push_back((b == std::string::npos) ? "" : tok.substr(b, e - b + 1));
    }
    return fields;
}

static CsvBinningResult ReadBinningCSV(const char* path) {
    CsvBinningResult res;
    std::ifstream f(path);
    if (!f.is_open()) {
        Warning("ReadBinningCSV", "Cannot open %s — falling back to hard-coded edges.", path);
        return res;
    }

    std::string line;
    while (std::getline(f, line)) {
        // strip trailing CR/LF, skip blank lines and comments
        while (!line.empty() && (line.back()=='\r'||line.back()=='\n')) line.pop_back();
        if (line.empty() || line[0]=='#') continue;

        auto fields = SplitCSV(line);
        if (fields.empty()) continue;
        const std::string& key = fields[0];

        // Helper: parse fields[1..] as doubles
        auto parseValues = [&](int startIdx) {
            std::vector<double> vals;
            for (int i = startIdx; i < (int)fields.size(); ++i) {
                if (fields[i].empty()) continue;
                char* ep = nullptr;
                double v = std::strtod(fields[i].c_str(), &ep);
                if (ep != fields[i].c_str()) vals.push_back(v);
            }
            return vals;
        };

        if (key == "beam_momentum") {
            auto v = parseValues(1);
            if (!v.empty()) res.beamMom = v[0];
        } else if (key == "W2_range") {
            auto v = parseValues(1);
            if (v.size() >= 2) { res.W2lo = v[0]; res.W2hi = v[1]; }
        } else if (key == "Q2_edges") {
            res.q2Edges = parseValues(1);
        } else if (key == "tprime_edges") {
            res.tpEdges = parseValues(1);
        }
        // unknown keys silently ignored — forward compatible
    }

    if (res.q2Edges.size() < 2 || res.tpEdges.size() < 2) {
        Warning("ReadBinningCSV", "Missing or too-short edge lists in %s — falling back.", path);
        return res;
    }

    res.valid = true;
    Printf("  [ReadBinningCSV] %s", path);
    Printf("    beam_mom=%.4f GeV  W2=[%.4f, %.4f] GeV²", res.beamMom, res.W2lo, res.W2hi);
    Printf("    Q² bins=%d  edges:", (int)res.q2Edges.size()-1);
    for (auto e : res.q2Edges) Printf("      %.6f", e);
    Printf("    t' bins=%d  edges:", (int)res.tpEdges.size()-1);
    for (auto e : res.tpEdges) Printf("      %.6f", e);
    return res;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Build the flat list of kinematic points from the bin-edge vectors
// ─────────────────────────────────────────────────────────────────────────────
static std::vector<BinDesc> BuildBins(
    const std::vector<double>& q2e,
    const std::vector<double>& tpe,
    double W2lo    = kW2lo,
    double W2hi    = kW2hi,
    double beamMom = kBmom)
{
    std::vector<BinDesc> bins;
    int nQ2 = (int)q2e.size() - 1;
    int nTp  = (int)tpe.size() - 1;

    for (int iQ2 = 0; iQ2 < nQ2; ++iQ2) {
        double Q2lo = q2e[iQ2], Q2hi = q2e[iQ2+1], Q2c = 0.5*(Q2lo+Q2hi);
        double tmin = TminAbs(Q2c, beamMom);

        for (int itp = 0; itp < nTp; ++itp) {
            BinDesc b;
            b.iQ2 = iQ2;  b.itp = itp;
            b.Q2lo = Q2lo; b.Q2hi = Q2hi; b.Q2c = Q2c;
            b.tplo = tpe[itp]; b.tphi = tpe[itp+1];
            b.tpc  = 0.5*(b.tplo + b.tphi);
            b.tmin_abs  = tmin;
            b.t_lo_abs  = b.tplo + tmin;
            b.t_hi_abs  = b.tphi + tmin;
            b.t_c_abs   = b.tpc  + tmin;
            b.W2c  = 0.5*(W2lo + W2hi);
            b.xBc  = ComputeXB(b.W2c, Q2c);
            b.Wc   = std::sqrt(std::max(0.0, b.W2c));
            // FORTRAN sign convention: −(physical value)
            b.fort_W2min = -W2hi;        // more negative = higher W²
            b.fort_W2max = -W2lo;
            b.fort_Q2min = -Q2hi;
            b.fort_Q2max = -Q2lo;
            b.fort_tmin  = -b.t_hi_abs;   // most negative t edge (largest |t|)
            b.fort_tmax  = -b.t_lo_abs;   // least negative t edge
            bins.push_back(b);
        }
    }
    return bins;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Write inmdi.dat  (FORTRAN sign convention throughout)
// ─────────────────────────────────────────────────────────────────────────────
static void WriteInmdiDat(const char* fname,
                          const std::vector<BinDesc>& bins,
                          double beamMom = kBmom,
                          long   seedOverride = 0)   // 0 = use kSeed
{
    int npoi = (int)bins.size();
    std::ofstream f(fname);
    if (!f.is_open()) {
        Error("WriteInmdiDat", "Cannot open %s", fname);
        return;
    }
    long seed = (seedOverride != 0) ? seedOverride : kSeed;
    f << std::fixed << std::setprecision(6);
    f << beamMom << "     !  bmom  - lepton momentum\n";
    f << kTmom   << "     !  tmom  - proton momentum (0 = fixed target)\n";
    f << kLepton << "         !  lepton - 1=electron 2=muon\n";
    f << kIvec   << "         !  ivec  - 1=rho 2=omega 3=phi 4=J/psi\n";
    f << kAnn1   << " " << kAnn2 << " " << kAnn3
      << "   ! ann1 ann2 ann3 - event counts\n";
    f << kVcut   << "        !  vcut  - v cut (0=none)\n";
    f << kNev    << "         !  nev   - MC samples per point\n";
    f << seed    << "    !  seed\n";
    f << npoi    << "         !  npoi  - number of kinematic points\n";

    // Each of the 6 rows: one value per point, space-separated
    auto row = [&](auto getter) {
        for (int i = 0; i < npoi; ++i)
            f << getter(bins[i]) << (i+1<npoi ? " " : "\n");
    };
    row([](const BinDesc& b){ return b.fort_W2min; });
    row([](const BinDesc& b){ return b.fort_W2max; });
    row([](const BinDesc& b){ return b.fort_Q2min; });
    row([](const BinDesc& b){ return b.fort_Q2max; });
    row([](const BinDesc& b){ return b.fort_tmin;  });
    row([](const BinDesc& b){ return b.fort_tmax;  });
    f.close();
    Printf("  [WriteInmdiDat] %s  (%d points)", fname, npoi);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Parse ALLmc.dat  — one row per (point × sample)
// ─────────────────────────────────────────────────────────────────────────────
struct AllMcRow { double born=NAN, dev=NAN, der=NAN, de=NAN, mean=NAN, err=NAN; };

static std::vector<double> ParseTokensFromLine(const std::string& line) {
    std::vector<double> out;
    std::istringstream iss(line);
    std::string tok;
    while (iss >> tok) {
        char* ep = nullptr;
        double v = std::strtod(tok.c_str(), &ep);
        if (ep != tok.c_str() && *ep == '\0') out.push_back(v);
    }
    return out;
}

static std::vector<AllMcRow> ParseAllMc(const char* fname) {
    std::vector<AllMcRow> rows;
    std::ifstream f(fname);
    if (!f.is_open()) { Warning("ParseAllMc","Cannot open %s", fname); return rows; }
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        auto n = ParseTokensFromLine(line);
        if ((int)n.size() < 4) continue;
        AllMcRow r;
        r.born = n[0]; r.dev = n.size()>1?n[1]:NAN;
        r.der  = n.size()>2?n[2]:NAN; r.de = n.size()>3?n[3]:NAN;
        if      (n.size()==6)  { r.mean=n[4]; r.err=n[5]; }
        else if (n.size()>=7)  { r.mean=n[n.size()-2]; r.err=n[n.size()-1]; }
        rows.push_back(r);
    }
    Printf("  [ParseAllMc] %zu rows from %s", rows.size(), fname);
    return rows;
}

// ─────────────────────────────────────────────────────────────────────────────
//  RC summary per bin (filled after grouping ALLmc rows)
// ─────────────────────────────────────────────────────────────────────────────
struct RcBin {
    BinDesc bin;
    double rc_mean=NAN, rc_sigma=NAN, rc_sem=NAN;
    int nsamp=0;
};

// ─────────────────────────────────────────────────────────────────────────────
//  Style helpers
// ─────────────────────────────────────────────────────────────────────────────
static void ApplyPhysicsStyle() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTopMargin(0.07);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetPadRightMargin(0.04);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetLineWidth(2);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42,"xyz");
    gStyle->SetTitleFont(42,"xyz");
    gStyle->SetLabelSize(0.046,"xyz");
    gStyle->SetTitleSize(0.050,"xyz");
    gStyle->SetTitleOffset(1.05,"x");
    gStyle->SetTitleOffset(1.15,"y");
    gStyle->SetMarkerSize(1.2);
    gStyle->SetEndErrorSize(5);
    gStyle->SetTickLength(0.025,"xyz");
}

// Colour palette — one per Q² bin (max 8 supported)
static Int_t BinColor(int iQ2) {
    static const Int_t pal[] = {
        TColor::GetColor("#2166ac"),   // deep blue
        TColor::GetColor("#d6604d"),   // red-orange
        TColor::GetColor("#33a02c"),   // green
        TColor::GetColor("#8856a7"),   // purple
        TColor::GetColor("#f4a582"),   // salmon
        TColor::GetColor("#4dac26"),   // lime
        TColor::GetColor("#e08214"),   // amber
        TColor::GetColor("#01665e"),   // teal
    };
    return pal[iQ2 % 8];
}
static Int_t BinMarker(int iQ2) {
    static const Int_t mk[] = {20, 21, 22, 23, 29, 33, 34, 47};
    return mk[iQ2 % 8];
}

// ─────────────────────────────────────────────────────────────────────────────
//  Draw one TGraphErrors band (filled) + markers on top
// ─────────────────────────────────────────────────────────────────────────────
static TGraphErrors* MakeGraph(const std::vector<RcBin>& summary,
                                int iQ2, bool useT) {
    std::vector<double> xs, ys, exs, eys;
    for (auto& rb : summary) {
        if (rb.bin.iQ2 != iQ2) continue;
        xs.push_back(useT ? rb.bin.t_c_abs : rb.bin.tpc);
        ys.push_back(rb.rc_mean);
        double ex = 0.3 * 0.5 * (useT
            ? (rb.bin.t_hi_abs - rb.bin.t_lo_abs)
            : (rb.bin.tphi     - rb.bin.tplo));
        exs.push_back(ex);
        eys.push_back(IsFinite(rb.rc_sigma) ? rb.rc_sigma : 0.0);
    }
    if (xs.empty()) return nullptr;
    auto* gr = new TGraphErrors((int)xs.size(),
                                xs.data(), ys.data(),
                                exs.data(), eys.data());
    gr->SetMarkerColor(BinColor(iQ2));
    gr->SetLineColor  (BinColor(iQ2));
    gr->SetFillColorAlpha(BinColor(iQ2), 0.20);
    gr->SetMarkerStyle(BinMarker(iQ2));
    gr->SetMarkerSize(1.3);
    gr->SetLineWidth(2);
    return gr;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Format a Q² bin label
// ─────────────────────────────────────────────────────────────────────────────
static TString Q2Label(int iQ2, const std::vector<double>& q2e) {
    return TString::Format("%.4f #leq Q^{2} < %.4f GeV^{2}",
                           q2e[iQ2], q2e[iQ2+1]);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Colour / marker / label for t' bins (warm palette — distinct from Q² blues)
// ─────────────────────────────────────────────────────────────────────────────
static Int_t TpBinColor(int itp) {
    static const Int_t pal[] = {
        TColor::GetColor("#b2182b"),  // deep red
        TColor::GetColor("#e34a33"),  // orange-red
        TColor::GetColor("#fc8d59"),  // salmon
        TColor::GetColor("#fdbb84"),  // peach
        TColor::GetColor("#9ecae1"),  // light blue
        TColor::GetColor("#4292c6"),  // mid blue
        TColor::GetColor("#08519c"),  // dark blue
        TColor::GetColor("#2c7bb6"),  // steel blue
        TColor::GetColor("#1a9641"),  // green (9th bin)
    };
    return pal[itp % 9];
}
static Int_t TpBinMarker(int itp) {
    static const Int_t mk[] = {20, 21, 22, 23, 24, 25, 26, 27, 28};
    return mk[itp % 9];
}
static TString TpLabel(int itp, const std::vector<double>& tpe) {
    return TString::Format("%.4f #leq t' < %.4f GeV^{2}",
                           tpe[itp], tpe[itp+1]);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Build TGraphErrors: δ_RC vs Q² centre, for a fixed t' bin (itp)
//    x  = Q²_centre,    xerr = 0.3 × half Q²-bin-width
//    y  = rc_mean,      yerr = rc_sigma
// ─────────────────────────────────────────────────────────────────────────────
static TGraphErrors* MakeGraph_vsQ2(const std::vector<RcBin>& summary, int itp) {
    std::vector<double> xs, ys, exs, eys;
    for (auto& rb : summary) {
        if (rb.bin.itp != itp) continue;
        if (!IsFinite(rb.rc_mean)) continue;
        xs.push_back(rb.bin.Q2c);
        ys.push_back(rb.rc_mean);
        exs.push_back(0.3 * 0.5 * (rb.bin.Q2hi - rb.bin.Q2lo));
        eys.push_back(IsFinite(rb.rc_sigma) ? rb.rc_sigma : 0.0);
    }
    if (xs.empty()) return nullptr;
    auto* gr = new TGraphErrors((int)xs.size(),
                                xs.data(), ys.data(),
                                exs.data(), eys.data());
    gr->SetMarkerColor(TpBinColor(itp));
    gr->SetLineColor  (TpBinColor(itp));
    gr->SetFillColorAlpha(TpBinColor(itp), 0.18);
    gr->SetMarkerStyle(TpBinMarker(itp));
    gr->SetMarkerSize(1.3);
    gr->SetLineWidth(2);
    return gr;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Draw bin-edge tick lines on a pad
// ─────────────────────────────────────────────────────────────────────────────
static void DrawEdgeLines(const std::vector<double>& edges,
                          double ylo, double yhi, bool vertical=true) {
    TLine ln;
    ln.SetLineStyle(3);
    ln.SetLineColor(kGray+1);
    ln.SetLineWidth(1);
    for (double e : edges) {
        if (vertical) ln.DrawLine(e, ylo, e, yhi);
        else          ln.DrawLine(ylo, e, yhi, e);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  MAIN
// ─────────────────────────────────────────────────────────────────────────────
void RunMDiffradNew()
{
    // ── All settings are at the top of this file ──────────────────────────────
    const char* fortranSrc     = kFortranSrc;
    const char* executable     = kExecutable;
    const char* binningCsvPath = kBinningCsv;
    bool        regenInput     = kRegenInput;
    bool        recompile      = kRecompile;

    // ── Create output directory and build output-path helper ──────────────────
    std::string outDir = kOutDir;
    // Trim trailing slash for consistency
    while (outDir.size() > 1 && outDir.back() == '/') outDir.pop_back();
    ::system(("mkdir -p " + outDir).c_str());
    ::system(("mkdir -p " + outDir + "/work").c_str());

    // O(path) — prepend outDir to a filename
    auto O = [&](const std::string& name) -> std::string {
        return outDir + "/" + name;
    };

    // All output paths
    const std::string kInmdiDatPath  = O("inmdi.dat");
    const std::string kOutCsvRC      = O("diffrad_rc_results.csv");
    const std::string kOutRoot       = O("mdiffrad_output.root");
    const std::string kOutPng        = O("mdiffrad_RC_summary.png");
    const std::string kOutPngQ2      = O("mdiffrad_RC_vsQ2.png");
    const std::string kWorkDir       = outDir + "/work";

    const char* inFile = kInmdiDatPath.c_str();

    Printf("\n══════════════════════════════════════════════════════════");
    Printf("  RunMDiffrad  –  φ electroproduction  (Q²×t' binning)");
    Printf("  Output directory: %s", outDir.c_str());

    // ── Step 0: Optionally read bin edges from DISANA_PhiAnalysisPlotter CSV ──
    // If binningCsvPath is non-empty and readable, it overrides the hard-coded
    // kQ2Edges / kTpEdges.  The tmin values are taken directly from the CSV so
    // the t ↔ t' mapping is identical in both macros.
    std::vector<double> activeQ2Edges;
    std::vector<double> activeTpEdges;
    double activeW2lo    = kW2lo;
    double activeW2hi    = kW2hi;
    double activeBeamMom = kBmom;         // overridden by CSV when provided

    bool useCsv = (binningCsvPath && binningCsvPath[0] != '\0');
    if (useCsv) {
        Printf("[Step 0] Reading bin edges from CSV: %s", binningCsvPath);
        CsvBinningResult csvResult = ReadBinningCSV(binningCsvPath);
        if (csvResult.valid) {
            activeQ2Edges   = csvResult.q2Edges;
            activeTpEdges   = csvResult.tpEdges;
            activeW2lo      = csvResult.W2lo;
            activeW2hi      = csvResult.W2hi;
            activeBeamMom   = csvResult.beamMom;
            Printf("  → Using CSV binning: %d Q² bins × %d t' bins",
                   (int)activeQ2Edges.size()-1, (int)activeTpEdges.size()-1);
            Printf("  → W² range from CSV: [%.4f, %.4f] GeV²", activeW2lo, activeW2hi);
            Printf("  → Beam momentum from CSV: %.4f GeV  (overrides kBmom=%.4f)",
                   activeBeamMom, kBmom);
        } else {
            Warning("RunMDiffradNew","CSV read failed — falling back to hard-coded edges.");
            useCsv = false;
        }
    }
    if (!useCsv) {
        Printf("[Step 0] Using hard-coded bin edges (no CSV provided or CSV invalid)");
        activeQ2Edges = kQ2Edges;
        activeTpEdges = kTpEdges;
        // activeBeamMom stays kBmom; BuildBins will compute tmin via TminAbs(Q2c, kBmom)
    }

    Printf("  %d Q² bins × %d t' bins = %d kinematic points",
           (int)activeQ2Edges.size()-1, (int)activeTpEdges.size()-1,
           (int)((activeQ2Edges.size()-1)*(activeTpEdges.size()-1)));
    Printf("══════════════════════════════════════════════════════════\n");

    // ── Print active edges ───────────────────────────────────────────────────
    Printf("=== Q^{2} bin edges (%d bins) ===", (int)activeQ2Edges.size()-1);
    for (int i = 0; i < (int)activeQ2Edges.size(); ++i)
        Printf("  edge[%d] = %.4f", i, activeQ2Edges[i]);
    Printf("\n=== t' bin edges (%d bins) ===", (int)activeTpEdges.size()-1);
    for (int i = 0; i < (int)activeTpEdges.size(); ++i)
        Printf("  edge[%d] = %.4f", i, activeTpEdges[i]);
    Printf("");

    // ── Step 1: Build bins ───────────────────────────────────────────────────
    auto bins = BuildBins(activeQ2Edges, activeTpEdges,
                          activeW2lo, activeW2hi, activeBeamMom);
    int npoi  = (int)bins.size();
    Printf("[Step 1] %d kinematic points built", npoi);

    // ── Step 2: Write inmdi.dat ──────────────────────────────────────────────
    if (regenInput) {
        Printf("[Step 2] Writing %s  (beam momentum = %.4f GeV)...", inFile, activeBeamMom);
        WriteInmdiDat(inFile, bins, activeBeamMom);
    } else {
        Printf("[Step 2] Skipping inmdi.dat write (regenInput=false)");
    }

    // ── Step 3: Compile FORTRAN ───────────────────────────────────────────────
    if (recompile) {
        Printf("[Step 3] Compiling %s ...", fortranSrc);
        TString cmd = TString::Format(
            "gfortran -O0 -g -std=legacy -ffixed-line-length-none "
            "-fbacktrace -ffpe-trap=invalid,overflow,zero "
            "-o %s %s -lm 2>&1",
            executable, fortranSrc);
        int ret = gSystem->Exec(cmd.Data());
        if (ret != 0) {
            Error("RunMDiffrad","Compilation failed (exit %d). Aborting.", ret);
            return;
        }
        Printf("  Compilation OK.");
    } else {
        Printf("[Step 3] Skipping compilation (recompile=false)");
    }

    // ── Step 4: Run MDIFFRAD in parallel — one subprocess per kinematic point ──
    //
    //  Strategy:
    //    • Create work/pt_<iQ2>_<itp>/ for every bin
    //    • Write a single-point INMDI.DAT (npoi=1) in each directory
    //    • Launch the FORTRAN executable concurrently with std::async
    //      (each process runs in its own directory → no file conflicts)
    //    • Throttle to kMaxJobs concurrent processes so we don't overload
    //      the node; set kMaxJobs = number of available cores (or less)
    //    • Collect ALLmc.dat from each directory in original bin order
    //
    Printf("[Step 4] Launching %d MDIFFRAD jobs in parallel (max %d concurrent)...",
           npoi, kMaxJobs);

    // ── Per-point work directory helper ──────────────────────────────────────
    auto PointDir = [&](const BinDesc& b) -> std::string {
        return kWorkDir + "/pt_" + std::to_string(b.iQ2) + "_" + std::to_string(b.itp);
    };

    // ── Resolve executable to absolute path BEFORE any cd into subdirs ─────────
    // The user may pass a relative path like "./../mdiffrad_exec".
    // Once we cd into work/pt_X_Y/ that relative path breaks.
    std::string absExec;
    {
        char* rp = realpath(executable, nullptr);
        if (rp) { absExec = rp; free(rp); }
        else     { absExec = std::string(gSystem->WorkingDirectory()) + "/" + executable; }
        Printf("  Executable (absolute): %s", absExec.c_str());
        if (gSystem->AccessPathName(absExec.c_str())) {
            Error("RunMDiffrad",
                  "Executable not found: %s  (from arg '%s', cwd '%s')",
                  absExec.c_str(), executable, gSystem->WorkingDirectory());
            return;
        }
    }

    // ── Write one INMDI.DAT per point ─────────────────────────────────────────
    // All filesystem ops here are single-threaded — safe to use system().
    // kWorkDir already created at startup; create per-point subdirs
    for (int i = 0; i < npoi; ++i) {
        const auto& b = bins[i];
        std::string dir = PointDir(b);
        ::system(("mkdir -p " + dir).c_str());
        std::string inmdi = dir + "/inmdi.dat";   // FORTRAN opens lowercase "inmdi.dat"
        // Unique seed per point so MC streams are independent
        long ptSeed = kSeed + static_cast<long>(i) * 1000L;
        WriteInmdiDat(inmdi.c_str(), {b}, activeBeamMom, ptSeed);
    }

    // ── Parallel launch ───────────────────────────────────────────────────────
    // RunPoint is the only function called from worker threads.
    // It captures absExec by value and touches nothing ROOT-related —
    // only std::string ops and ::system(), both of which are thread-safe.
    auto RunPoint = [absExec, kWorkDir](const BinDesc& b) -> int {
        std::string dir = kWorkDir + "/pt_" + std::to_string(b.iQ2)
                                     + "_" + std::to_string(b.itp);
        std::string cmd = "cd " + dir + " && " + absExec + " > mdiffrad.log 2>&1";
        return ::system(cmd.c_str());
    };

    // ── Throttled async pool ──────────────────────────────────────────────────
    // Invariant: every future in `running` is valid (not yet get()-ed).
    // Flush() collects completed futures safely without invalidating iterators:
    //   - builds a list of ready indices first, then harvests them in reverse
    //     order so erase() indices stay valid.
    // All Printf/logging happens on the main thread only.
    {
        using FutType = std::future<int>;
        std::vector<FutType> running;
        running.reserve(kMaxJobs);

        int submitted = 0;
        int failed    = 0;

        // Harvest all futures that are already done (or all if waitAll=true).
        auto Flush = [&](bool waitAll) {
            std::vector<int> doneIdx;
            for (int i = 0; i < (int)running.size(); ++i) {
                if (!running[i].valid()) { doneIdx.push_back(i); continue; }
                bool ready = waitAll ||
                    running[i].wait_for(std::chrono::seconds(0))
                        == std::future_status::ready;
                if (ready) {
                    int rc = running[i].get();   // blocks if waitAll and not yet done
                    if (rc != 0) ++failed;
                    doneIdx.push_back(i);
                }
            }
            // Erase in reverse order so indices stay valid
            for (int i = (int)doneIdx.size()-1; i >= 0; --i)
                running.erase(running.begin() + doneIdx[i]);
        };

        for (auto& b : bins) {
            // Block until a slot is free
            while ((int)running.size() >= kMaxJobs) {
                std::this_thread::sleep_for(std::chrono::milliseconds(200));
                Flush(false);
            }
            running.push_back(std::async(std::launch::async, RunPoint, b));
            ++submitted;
            if (submitted % 4 == 0 || submitted == npoi)
                Printf("  [Step 4] Submitted %d / %d jobs  (%d running)",
                       submitted, npoi, (int)running.size());
        }
        // Drain remaining jobs
        while (!running.empty()) {
            Flush(true);
            if (!running.empty())
                std::this_thread::sleep_for(std::chrono::milliseconds(200));
        }

        if (failed > 0) {
            Printf("Warning in <RunMDiffrad>: %d / %d jobs returned non-zero exit code", failed, npoi);
            Printf("  Check logs in: %s/pt_*/mdiffrad.log", kWorkDir.c_str());
            // Print the first failure log to help diagnose the problem
            std::string firstLog = kWorkDir + "/pt_0_0/mdiffrad.log";
            Printf("  Printing first failure log (%s):", firstLog.c_str());
            gSystem->Exec(("head -30 " + firstLog + " 2>/dev/null || echo '(log not found)'").c_str());
        } else
            Printf("  [Step 4] All %d jobs finished OK.", npoi);
    }

    // ── Step 5: Collect ALLmc.dat from each work directory ────────────────────
    Printf("[Step 5] Collecting ALLmc.dat results ...");
    std::vector<AllMcRow> allmc;
    allmc.reserve(npoi * kNev);

    for (auto& b : bins) {
        std::string almcPath = PointDir(b) + "/ALLmc.dat";
        auto rows = ParseAllMc(almcPath.c_str());
        if (rows.empty())
            Warning("RunMDiffrad",
                    "Empty ALLmc.dat for iQ2=%d itp=%d — bin RC will be NaN",
                    b.iQ2, b.itp);
        // Append exactly kNev rows (pad with NaN if short)
        for (int s = 0; s < kNev; ++s) {
            if (s < (int)rows.size())
                allmc.push_back(rows[s]);
            else
                allmc.push_back(AllMcRow{});  // NaN-filled
        }
    }
    Printf("  Collected %zu rows total (%d points × %d samples)",
           allmc.size(), npoi, kNev);

    // ── Step 6: Group by point → RC mean/σ ───────────────────────────────────
    Printf("[Step 6] Computing RC per bin ...");
    std::vector<RcBin> summary;
    summary.reserve(npoi);

    for (int i = 0; i < npoi; ++i) {
        std::vector<double> de;
        de.reserve(kNev);
        for (int s = 0; s < kNev; ++s) {
            int idx = i * kNev + s;
            if (idx < (int)allmc.size()) de.push_back(allmc[idx].de);
        }
        RcBin rb;
        rb.bin     = bins[i];
        rb.nsamp   = (int)de.size();
        rb.rc_mean = MeanFinite(de);
        rb.rc_sigma= StdFiniteSample(de);
        if (IsFinite(rb.rc_sigma) && rb.nsamp > 0)
            rb.rc_sem = rb.rc_sigma / std::sqrt((double)rb.nsamp);
        summary.push_back(rb);
    }

    // ── Step 7: Write Python-friendly RC results CSV ──────────────────────────
    //
    //  Single output file kOutCsvRC containing everything needed for plotting
    //  and RDF postprocessing:
    //    Bin indices, Q² edges & centre, t' edges & centre,
    //    |t| edges & centre, tmin, signed-t edges & centre,
    //    W², xB, W at bin centre,
    //    rad_corr (δ_RC), its σ and SEM, number of MC samples,
    //    and a ready-made RDF filter string.
    //
    //  Conventions (matching DISANAMMUtils):
    //    tprime_*   positive  (mtprime = |t| − |tmin|)
    //    t_abs_*    positive  (= tprime + tmin_abs)
    //    t_phys_*   negative  (= −t_abs,  matches the "t" column in the RDF)
    // ─────────────────────────────────────────────────────────────────────────
    Printf("[Step 7] Writing %s ...", kOutCsvRC.c_str());
    {
        std::ofstream csv(kOutCsvRC.c_str());
        csv << std::fixed << std::setprecision(8);

        // ── Header comments ──────────────────────────────────────────────────
        csv << "# MDIFFRAD radiative correction results\n"
            << "# Generated by RunMDiffradNew\n"
            << "#\n"
            << "# Conventions:\n"
            << "#   tprime_*  = positive  (mtprime = |t| - |tmin|)\n"
            << "#   t_abs_*   = positive  (= tprime + tmin_abs)\n"
            << "#   t_phys_*  = negative  (= -t_abs, matches RDF 't' column)\n"
            << "#   W2c, xBc, Wc are evaluated at (Q2c, W2_centre)\n"
            << "#   rad_corr = delta_RC  (Born / observed cross-section ratio)\n"
            << "#   rdf_filter can be pasted directly into RDataFrame::Filter()\n"
            << "#\n"
            << "# beam_momentum_GeV = " << activeBeamMom << "\n"
            << "# W2_range_GeV2     = " << activeW2lo << " " << activeW2hi << "\n"
            << "#\n";

        // ── Column header ────────────────────────────────────────────────────
        csv << "iQ2,itp,"
               "Q2_lo,Q2_hi,Q2_c,"
               "tprime_lo,tprime_hi,tprime_c,"
               "t_abs_lo,t_abs_hi,t_abs_c,"
               "t_phys_lo,t_phys_hi,t_phys_c,"
               "tmin_abs,"
               "W2_lo,W2_hi,W2_c,"
               "xBc,Wc,"
               "rad_corr,rad_corr_sigma,rad_corr_sem,"
               "nsamp,"
               "rdf_filter\n";

        // ── Data rows ────────────────────────────────────────────────────────
        for (auto& rb : summary) {
            const auto& b = rb.bin;
            TString filt = TString::Format(
                "Q2 >= %.6f && Q2 < %.6f && mtprime >= %.6f && mtprime < %.6f",
                b.Q2lo, b.Q2hi, b.tplo, b.tphi);
            csv << b.iQ2        << "," << b.itp        << ","
                << b.Q2lo       << "," << b.Q2hi       << "," << b.Q2c       << ","
                << b.tplo       << "," << b.tphi       << "," << b.tpc       << ","
                << b.t_lo_abs   << "," << b.t_hi_abs   << "," << b.t_c_abs   << ","
                << -b.t_hi_abs  << "," << -b.t_lo_abs  << "," << -b.t_c_abs  << ","
                << b.tmin_abs   << ","
                << activeW2lo   << "," << activeW2hi   << "," << b.W2c       << ","
                << b.xBc        << "," << b.Wc         << ","
                << rb.rc_mean   << "," << rb.rc_sigma  << "," << rb.rc_sem   << ","
                << rb.nsamp     << ","
                << "\"" << filt.Data() << "\"\n";
        }
        Printf("  Written: %s  (%d rows, %d Q2 bins x %d t' bins)",
               kOutCsvRC.c_str(), npoi,
               (int)activeQ2Edges.size()-1, (int)activeTpEdges.size()-1);
    }

    // ── Step 8: Publication-quality plots ────────────────────────────────────
    Printf("[Step 8] Drawing plots ...");
    ApplyPhysicsStyle();

    // Aliases so the plot code works identically whether edges came from CSV
    // or from the hard-coded constants.
    const std::vector<double>& plotQ2Edges = activeQ2Edges;
    const std::vector<double>& plotTpEdges = activeTpEdges;
    const double plotW2lo = activeW2lo;
    const double plotW2hi = activeW2hi;

    int nQ2 = (int)plotQ2Edges.size() - 1;
    int nTp = (int)plotTpEdges.size()  - 1;

    // Compute y-range from data
    double ymin=1e9, ymax=-1e9;
    for (auto& rb : summary) {
        if (!IsFinite(rb.rc_mean)) continue;
        double lo = rb.rc_mean - (IsFinite(rb.rc_sigma) ? rb.rc_sigma : 0.0);
        double hi = rb.rc_mean + (IsFinite(rb.rc_sigma) ? rb.rc_sigma : 0.0);
        ymin = std::min(ymin, lo);
        ymax = std::max(ymax, hi);
    }
    if (ymin > ymax) { ymin=0.90; ymax=1.15; }
    double ypad  = 0.10 * (ymax - ymin);
    double ylo   = ymin - ypad;
    double yhi   = ymax + ypad;
    double tp_hi = plotTpEdges.back() * 1.03;
    double t_hi  = summary.back().bin.t_hi_abs * 1.03;

    // ─── Canvas: 2×2 layout ──────────────────────────────────────────────────
    //   [0,0] δ_RC vs t'          [0,1] δ_RC vs |t|
    //   [1,0] 2D heatmap          [1,1] δ_RC vs point index
    TCanvas* c = new TCanvas("cMDiffrad","MDIFFRAD RC Summary",1700,1200);
    c->Divide(2, 2, 0.004, 0.004);

    // ──────────────────────────────
    //  Panel (1,1): δ_RC vs t'
    // ──────────────────────────────
    c->cd(1);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);

    TLegend* leg1 = new TLegend(0.50, 0.62, 0.96, 0.92);
    leg1->SetTextFont(42); leg1->SetTextSize(0.036);
    leg1->SetBorderSize(0); leg1->SetFillStyle(0);

    bool first1 = true;
    for (int iQ2 = 0; iQ2 < nQ2; ++iQ2) {
        auto* gr = MakeGraph(summary, iQ2, false);
        if (!gr) continue;
        if (first1) {
            gr->GetXaxis()->SetTitle("t' = |t| #minus |t_{min}| (GeV^{2})");
            gr->GetYaxis()->SetTitle("Radiative Correction  #delta_{RC}");
            gr->GetXaxis()->SetLimits(-0.02, tp_hi);
            gr->SetMinimum(ylo); gr->SetMaximum(yhi);
            gr->Draw("A P E2");
            first1 = false;
        } else {
            gr->Draw("P E2 same");
        }
        leg1->AddEntry(gr, Q2Label(iQ2, plotQ2Edges), "lp");
    }
    // bin edge guide lines
    DrawEdgeLines(plotTpEdges, ylo, yhi, true);
    // unity line
    TLine unity1;
    unity1.SetLineStyle(2); unity1.SetLineColor(kGray+2); unity1.SetLineWidth(1);
    unity1.DrawLine(0, 1.0, tp_hi, 1.0);
    leg1->Draw();
    // labels
    TLatex lat;
    lat.SetNDC(); lat.SetTextFont(62); lat.SetTextSize(0.045);
    lat.DrawLatex(0.17, 0.90, "#delta_{RC} vs t'");
    lat.SetTextFont(42); lat.SetTextSize(0.035);
    lat.DrawLatex(0.17, 0.84,
        TString::Format("E_{beam}=%.2f GeV, iv=%d (#phi), W^{2}#in[%.1f,%.1f]",
                        activeBeamMom, kIvec, plotW2lo, plotW2hi));

    // ──────────────────────────────
    //  Panel (1,2): δ_RC vs |t|
    // ──────────────────────────────
    c->cd(2);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);

    TLegend* leg2 = new TLegend(0.50, 0.62, 0.96, 0.92);
    leg2->SetTextFont(42); leg2->SetTextSize(0.036);
    leg2->SetBorderSize(0); leg2->SetFillStyle(0);

    bool first2 = true;
    for (int iQ2 = 0; iQ2 < nQ2; ++iQ2) {
        auto* gr = MakeGraph(summary, iQ2, true);
        if (!gr) continue;
        if (first2) {
            gr->GetXaxis()->SetTitle("|t| (GeV^{2})");
            gr->GetYaxis()->SetTitle("Radiative Correction  #delta_{RC}");
            gr->GetXaxis()->SetLimits(0, t_hi);
            gr->SetMinimum(ylo); gr->SetMaximum(yhi);
            gr->Draw("A P E2");
            first2 = false;
        } else {
            gr->Draw("P E2 same");
        }
        leg2->AddEntry(gr, Q2Label(iQ2, plotQ2Edges), "lp");
    }
    TLine unity2;
    unity2.SetLineStyle(2); unity2.SetLineColor(kGray+2); unity2.SetLineWidth(1);
    unity2.DrawLine(0, 1.0, t_hi, 1.0);
    leg2->Draw();
    TLatex lat2;
    lat2.SetNDC(); lat2.SetTextFont(62); lat2.SetTextSize(0.045);
    lat2.DrawLatex(0.17, 0.90, "#delta_{RC} vs |t|");

    // ──────────────────────────────
    //  Panel (2,1): 2D heatmap δ_RC(Q², t')
    // ──────────────────────────────
    c->cd(3);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
    gPad->SetRightMargin(0.14);

    // Convert edge vectors to arrays for TH2D
    std::vector<Double_t> q2arr(plotQ2Edges.begin(), plotQ2Edges.end());
    std::vector<Double_t> tparr(plotTpEdges.begin(),  plotTpEdges.end());

    TH2D* h2 = new TH2D("h2rc",
        ";Q^{2} (GeV^{2});t' (GeV^{2})",
        nQ2, q2arr.data(), nTp, tparr.data());

    for (auto& rb : summary) {
        if (!IsFinite(rb.rc_mean)) continue;
        h2->Fill(rb.bin.Q2c, rb.bin.tpc, rb.rc_mean);
    }

    // Use a diverging palette centred on 1.0
    gStyle->SetPalette(kTemperatureMap);
    h2->GetZaxis()->SetTitle("#delta_{RC}");
    h2->GetZaxis()->SetTitleOffset(1.4);
    h2->Draw("COLZ");

    // Annotate each cell with its value
    TLatex cellLat;
    cellLat.SetTextAlign(22);
    cellLat.SetTextFont(42);
    cellLat.SetTextSize(0.028);
    for (auto& rb : summary) {
        if (!IsFinite(rb.rc_mean)) continue;
        cellLat.DrawLatex(rb.bin.Q2c, rb.bin.tpc,
            TString::Format("%.3f", rb.rc_mean));
    }

    TLatex lat3;
    lat3.SetNDC(); lat3.SetTextFont(62); lat3.SetTextSize(0.045);
    lat3.DrawLatex(0.17, 0.93, "#delta_{RC} (Q^{2}, t')  heatmap");

    // ──────────────────────────────
    //  Panel (2,2): δ_RC vs flat point index
    // ──────────────────────────────
    c->cd(4);
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.14);

    int M = (int)summary.size();
    std::vector<double> xi(M), yi(M), exi(M, 0.0), eyi(M);
    for (int i = 0; i < M; ++i) {
        xi[i]  = i + 1;
        yi[i]  = summary[i].rc_mean;
        eyi[i] = IsFinite(summary[i].rc_sigma) ? summary[i].rc_sigma : 0.0;
    }
    auto* grIdx = new TGraphErrors(M, xi.data(), yi.data(), exi.data(), eyi.data());
    grIdx->SetMarkerStyle(20);
    grIdx->SetMarkerColor(kBlue+1);
    grIdx->SetLineColor(kBlue+1);
    grIdx->SetMarkerSize(0.9);
    grIdx->SetTitle(";Point index (iQ2 #times nTp + itp);#delta_{RC}");
    grIdx->GetYaxis()->SetRangeUser(ylo, yhi);
    grIdx->Draw("AP");

    // Vertical lines between Q² groups
    TLine vln;
    vln.SetLineStyle(3); vln.SetLineColor(kGray+1);
    for (int iQ2 = 1; iQ2 < nQ2; ++iQ2)
        vln.DrawLine(iQ2 * nTp + 0.5, ylo, iQ2 * nTp + 0.5, yhi);

    // Q² group labels
    TLatex grpLat;
    grpLat.SetTextFont(42); grpLat.SetTextSize(0.030); grpLat.SetTextAlign(22);
    for (int iQ2 = 0; iQ2 < nQ2; ++iQ2) {
        double xmid = iQ2 * nTp + 0.5 * nTp + 1;
        grpLat.DrawLatex(xmid, yhi - 0.4*ypad,
            TString::Format("Q^{2}_{%d}", iQ2));
    }

    TLine unity3;
    unity3.SetLineStyle(2); unity3.SetLineColor(kGray+2); unity3.SetLineWidth(1);
    unity3.DrawLine(0.5, 1.0, M + 0.5, 1.0);

    TLatex lat4;
    lat4.SetNDC(); lat4.SetTextFont(62); lat4.SetTextSize(0.045);
    lat4.DrawLatex(0.15, 0.93, "#delta_{RC} vs point index");

    c->Update();
    c->SaveAs(kOutPng.c_str());
    Printf("  Saved: %s", kOutPng.c_str());

    // ─────────────────────────────────────────────────────────────────────────
    //  Canvas 2: δ_RC vs Q²  (one curve per t' bin)
    //
    //  Layout  2 × 2:
    //    [1,1]  δ_RC vs Q²  — all t' bins overlaid, coloured by t' bin
    //    [1,2]  δ_RC vs Q²  — split into two sub-groups (low / high t')
    //           to avoid legend overcrowding with 8 bins
    //    [2,1]  same as [1,1] but plotted vs Q²_centre on a log-x axis
    //    [2,2]  δ_RC vs t' bin index for each Q² bin
    //           (transposed view of the point-index panel from Canvas 1)
    // ─────────────────────────────────────────────────────────────────────────
    TCanvas* cQ2 = new TCanvas("cMDiffrad_vsQ2",
                               "MDIFFRAD  #delta_{RC} vs Q^{2}", 1700, 1200);
    cQ2->Divide(2, 2, 0.004, 0.004);

    double q2_lo_axis = plotQ2Edges.front();
    double q2_hi_axis = plotQ2Edges.back() * 1.05;

    // ── Panel (1,1): all t' bins, δ_RC vs Q² ─────────────────────────────
    cQ2->cd(1);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);

    // Legend: split into two columns to keep it compact with 8 t' bins
    TLegend* legQ1 = new TLegend(0.14, 0.13, 0.96, 0.42);
    legQ1->SetNColumns(2);
    legQ1->SetTextFont(42); legQ1->SetTextSize(0.030);
    legQ1->SetBorderSize(0); legQ1->SetFillStyle(0);

    bool firstQ1 = true;
    for (int itp = 0; itp < nTp; ++itp) {
        auto* gr = MakeGraph_vsQ2(summary, itp);
        if (!gr) continue;
        if (firstQ1) {
            gr->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
            gr->GetYaxis()->SetTitle("Radiative Correction  #delta_{RC}");
            gr->GetXaxis()->SetLimits(q2_lo_axis, q2_hi_axis);
            gr->SetMinimum(ylo); gr->SetMaximum(yhi);
            gr->Draw("A P E2");
            firstQ1 = false;
        } else {
            gr->Draw("P E2 same");
        }
        legQ1->AddEntry(gr, TpLabel(itp, plotTpEdges), "lp");
    }
    DrawEdgeLines(plotQ2Edges, ylo, yhi, true);
    TLine uQ1;
    uQ1.SetLineStyle(2); uQ1.SetLineColor(kGray+2); uQ1.SetLineWidth(1);
    uQ1.DrawLine(q2_lo_axis, 1.0, q2_hi_axis, 1.0);
    legQ1->Draw();
    {
        TLatex lx; lx.SetNDC(); lx.SetTextFont(62); lx.SetTextSize(0.045);
        lx.DrawLatex(0.17, 0.93, "#delta_{RC} vs Q^{2}  (all t' bins)");
        lx.SetTextFont(42); lx.SetTextSize(0.033);
        lx.DrawLatex(0.17, 0.87,
            TString::Format("E_{beam}=%.2f GeV,  W^{2}#in[%.1f,%.1f] GeV^{2}",
                            activeBeamMom, plotW2lo, plotW2hi));
    }

    // ── Panel (1,2): split low / high t' bins ────────────────────────────
    //  Low  = first  ceil(nTp/2) bins   High = remaining bins
    cQ2->cd(2);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);

    int nLow  = (nTp + 1) / 2;   // ceil(nTp/2)
    int nHigh = nTp - nLow;

    // Two sub-legends side by side
    TLegend* legLow  = new TLegend(0.14, 0.13, 0.55, 0.13 + 0.055*nLow);
    TLegend* legHigh = new TLegend(0.55, 0.13, 0.96, 0.13 + 0.055*nHigh);
    for (auto* lg : {legLow, legHigh}) {
        lg->SetTextFont(42); lg->SetTextSize(0.032);
        lg->SetBorderSize(0); lg->SetFillStyle(0);
    }
    legLow->SetHeader("Low t' bins", "C");
    legHigh->SetHeader("High t' bins", "C");

    bool firstQ2p = true;
    for (int itp = 0; itp < nTp; ++itp) {
        auto* gr = MakeGraph_vsQ2(summary, itp);
        if (!gr) continue;
        if (firstQ2p) {
            gr->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
            gr->GetYaxis()->SetTitle("Radiative Correction  #delta_{RC}");
            gr->GetXaxis()->SetLimits(q2_lo_axis, q2_hi_axis);
            gr->SetMinimum(ylo); gr->SetMaximum(yhi);
            gr->Draw("A P E2");
            firstQ2p = false;
        } else {
            gr->Draw("P E2 same");
        }
        if (itp < nLow)
            legLow->AddEntry(gr, TpLabel(itp, plotTpEdges), "lp");
        else
            legHigh->AddEntry(gr, TpLabel(itp, plotTpEdges), "lp");
    }
    DrawEdgeLines(plotQ2Edges, ylo, yhi, true);
    TLine uQ2p;
    uQ2p.SetLineStyle(2); uQ2p.SetLineColor(kGray+2); uQ2p.SetLineWidth(1);
    uQ2p.DrawLine(q2_lo_axis, 1.0, q2_hi_axis, 1.0);
    legLow->Draw();
    legHigh->Draw();
    {
        TLatex lx; lx.SetNDC(); lx.SetTextFont(62); lx.SetTextSize(0.045);
        lx.DrawLatex(0.17, 0.93, "#delta_{RC} vs Q^{2}  (low / high t')");
    }

    // ── Panel (2,1): log-scale Q² axis ───────────────────────────────────
    cQ2->cd(3);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
    gPad->SetLogx(1);

    double q2_lo_log = std::max(plotQ2Edges[1] * 0.85, 0.05); // skip the 0 edge for log
    bool firstQ3 = true;
    for (int itp = 0; itp < nTp; ++itp) {
        auto* gr = MakeGraph_vsQ2(summary, itp);
        if (!gr) continue;
        if (firstQ3) {
            gr->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
            gr->GetYaxis()->SetTitle("Radiative Correction  #delta_{RC}");
            gr->GetXaxis()->SetLimits(q2_lo_log, q2_hi_axis);
            gr->SetMinimum(ylo); gr->SetMaximum(yhi);
            gr->Draw("A P E2");
            firstQ3 = false;
        } else {
            gr->Draw("P E2 same");
        }
    }
    // Q² edge lines on log scale
    {
        TLine ln; ln.SetLineStyle(3); ln.SetLineColor(kGray+1); ln.SetLineWidth(1);
        for (double e : plotQ2Edges)
            if (e > 0) ln.DrawLine(e, ylo, e, yhi);
    }
    TLine uQ3;
    uQ3.SetLineStyle(2); uQ3.SetLineColor(kGray+2); uQ3.SetLineWidth(1);
    uQ3.DrawLine(q2_lo_log, 1.0, q2_hi_axis, 1.0);
    {
        TLatex lx; lx.SetNDC(); lx.SetTextFont(62); lx.SetTextSize(0.045);
        lx.DrawLatex(0.17, 0.93, "#delta_{RC} vs Q^{2}  (log scale)");
        lx.SetTextFont(42); lx.SetTextSize(0.033);
        lx.DrawLatex(0.17, 0.87, "colours: t' bins (same as panels above)");
    }

    // ── Panel (2,2): δ_RC vs t'-bin index, one curve per Q² bin ─────────
    //  This is the "transposed" view: same data as panel(1,1) of Canvas 1
    //  but now x = t'-bin index and curves = Q² bins, making it easy to
    //  read off Q²-dependence at a glance for any given t' slice.
    cQ2->cd(4);
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.14);
    gPad->SetLogx(0);

    TLegend* legQ4 = new TLegend(0.52, 0.62, 0.96, 0.92);
    legQ4->SetTextFont(42); legQ4->SetTextSize(0.036);
    legQ4->SetBorderSize(0); legQ4->SetFillStyle(0);

    bool firstQ4 = true;
    for (int iQ2 = 0; iQ2 < nQ2; ++iQ2) {
        // Build graph: x = t'-bin centre index (1..nTp), y = rc_mean
        std::vector<double> xs4, ys4, exs4, eys4;
        for (auto& rb : summary) {
            if (rb.bin.iQ2 != iQ2 || !IsFinite(rb.rc_mean)) continue;
            xs4.push_back(rb.bin.itp + 1.0 + iQ2 * 0.15);  // slight x-offset per Q² bin
            ys4.push_back(rb.rc_mean);
            exs4.push_back(0.0);
            eys4.push_back(IsFinite(rb.rc_sigma) ? rb.rc_sigma : 0.0);
        }
        if (xs4.empty()) continue;
        auto* gr4 = new TGraphErrors((int)xs4.size(),
                                     xs4.data(), ys4.data(),
                                     exs4.data(), eys4.data());
        gr4->SetMarkerColor(BinColor(iQ2));
        gr4->SetLineColor  (BinColor(iQ2));
        gr4->SetMarkerStyle(BinMarker(iQ2));
        gr4->SetMarkerSize(1.2);
        gr4->SetLineWidth(2);
        if (firstQ4) {
            gr4->GetXaxis()->SetTitle("t' bin index");
            gr4->GetYaxis()->SetTitle("Radiative Correction  #delta_{RC}");
            gr4->GetXaxis()->SetLimits(0.5, nTp + 1.0);
            gr4->SetMinimum(ylo); gr4->SetMaximum(yhi);
            // label the t'-bin tick marks with t' edges
            for (int itp = 0; itp < nTp; ++itp) {
                gr4->GetXaxis()->SetBinLabel(
                    gr4->GetXaxis()->FindBin(itp + 1.0),
                    TString::Format("%.2f", plotTpEdges[itp]));
            }
            gr4->Draw("A P");
            firstQ4 = false;
        } else {
            gr4->Draw("P same");
        }
        legQ4->AddEntry(gr4, Q2Label(iQ2, plotQ2Edges), "lp");
    }
    // Vertical guide lines at each t'-bin boundary
    {
        TLine ln; ln.SetLineStyle(3); ln.SetLineColor(kGray+1); ln.SetLineWidth(1);
        for (int itp = 1; itp < nTp; ++itp)
            ln.DrawLine(itp + 0.5, ylo, itp + 0.5, yhi);
    }
    TLine uQ4;
    uQ4.SetLineStyle(2); uQ4.SetLineColor(kGray+2); uQ4.SetLineWidth(1);
    uQ4.DrawLine(0.5, 1.0, nTp + 1.0, 1.0);
    legQ4->Draw();
    {
        TLatex lx; lx.SetNDC(); lx.SetTextFont(62); lx.SetTextSize(0.043);
        lx.DrawLatex(0.15, 0.93, "#delta_{RC} vs t' index  (Q^{2} slices)");
        lx.SetTextFont(42); lx.SetTextSize(0.031);
        lx.DrawLatex(0.15, 0.87, "x-axis = t' bin index, offset by Q^{2} bin");
    }

    cQ2->Update();
    cQ2->SaveAs(kOutPngQ2.c_str());
    Printf("  Saved: %s", kOutPngQ2.c_str());

    // ── Step 9: Write ROOT file ────────────────────────────────────────────
    Printf("[Step 9] Writing ROOT file %s ...", kOutRoot.c_str());
    TFile* fout = new TFile(kOutRoot.c_str(), "RECREATE");

    // Per-sample tree
    TTree* tAll = new TTree("allmc","MDIFFRAD per-sample rows");
    Double_t br_born, br_dev, br_der, br_de, br_mean, br_err;
    tAll->Branch("born",&br_born,"born/D"); tAll->Branch("dev",&br_dev,"dev/D");
    tAll->Branch("der",&br_der,"der/D");   tAll->Branch("de",&br_de,"de/D");
    tAll->Branch("mean",&br_mean,"mean/D"); tAll->Branch("err",&br_err,"err/D");
    for (auto& rr : allmc) {
        br_born=rr.born; br_dev=rr.dev; br_der=rr.der;
        br_de=rr.de; br_mean=rr.mean; br_err=rr.err;
        tAll->Fill();
    }

    // Summary tree — complete bin geometry + RC
    TTree* tSum = new TTree("rcsum","RC summary per (Q2,t') bin");
    Int_t    br_iQ2, br_itp, br_nsamp;
    Double_t br_Q2lo, br_Q2hi, br_Q2c;
    Double_t br_tplo, br_tphi, br_tpc;
    Double_t br_tlo,  br_thi,  br_tc, br_tmin;
    Double_t br_xBc,  br_Wc;
    Double_t br_rc,   br_rcsig, br_rcsem;
    tSum->Branch("iQ2",        &br_iQ2,   "iQ2/I");
    tSum->Branch("itp",        &br_itp,   "itp/I");
    tSum->Branch("Q2lo",       &br_Q2lo,  "Q2lo/D");
    tSum->Branch("Q2hi",       &br_Q2hi,  "Q2hi/D");
    tSum->Branch("Q2c",        &br_Q2c,   "Q2c/D");
    tSum->Branch("tprime_lo",  &br_tplo,  "tprime_lo/D");
    tSum->Branch("tprime_hi",  &br_tphi,  "tprime_hi/D");
    tSum->Branch("tprime_c",   &br_tpc,   "tprime_c/D");
    tSum->Branch("t_lo_abs",   &br_tlo,   "t_lo_abs/D");
    tSum->Branch("t_hi_abs",   &br_thi,   "t_hi_abs/D");
    tSum->Branch("t_c_abs",    &br_tc,    "t_c_abs/D");
    tSum->Branch("tmin_abs",   &br_tmin,  "tmin_abs/D");
    tSum->Branch("xBc",        &br_xBc,   "xBc/D");
    tSum->Branch("Wc",         &br_Wc,    "Wc/D");
    tSum->Branch("rad_corr",       &br_rc,    "rad_corr/D");
    tSum->Branch("rad_corr_sigma", &br_rcsig, "rad_corr_sigma/D");
    tSum->Branch("rad_corr_sem",   &br_rcsem, "rad_corr_sem/D");
    tSum->Branch("nsamp",          &br_nsamp, "nsamp/I");
    for (auto& rb : summary) {
        const auto& b = rb.bin;
        br_iQ2=b.iQ2; br_itp=b.itp;
        br_Q2lo=b.Q2lo; br_Q2hi=b.Q2hi; br_Q2c=b.Q2c;
        br_tplo=b.tplo; br_tphi=b.tphi; br_tpc=b.tpc;
        br_tlo=b.t_lo_abs; br_thi=b.t_hi_abs; br_tc=b.t_c_abs;
        br_tmin=b.tmin_abs; br_xBc=b.xBc; br_Wc=b.Wc;
        br_rc=rb.rc_mean; br_rcsig=rb.rc_sigma; br_rcsem=rb.rc_sem;
        br_nsamp=rb.nsamp;
        tSum->Fill();
    }
    h2->Write();
    tAll->Write(); tSum->Write();
    fout->Close();
    Printf("  Written: %s", kOutRoot.c_str());

    // ── Done ──────────────────────────────────────────────────────────────
    Printf("\n══════════════════════════════════════════════════════════");
    Printf("  DONE   %d bins processed", npoi);
    Printf("  Outputs:");
    Printf("    %-42s  RC results + bin geometry (Python-ready)", kOutCsvRC.c_str());
    Printf("    %-42s  RC vs t', |t|, heatmap, index",            kOutPng.c_str());
    Printf("    %-42s  RC vs Q² per t' bin",                      kOutPngQ2.c_str());
    Printf("    %-42s  ROOT trees",                               kOutRoot.c_str());
    Printf("    %-42s  per-point logs",                           (kWorkDir + "/pt_*/mdiffrad.log").c_str());
    Printf("══════════════════════════════════════════════════════════");
    Printf("  Python plotting:");
    Printf("    python3 plot_rad_corrections.py %s", kOutCsvRC.c_str());
    Printf("  RDF postprocessing:");
    Printf("    // read rdf_filter column from %s, then:", kOutCsvRC.c_str());
    Printf("    df.Define(\"sigma_rc\", \"sigma_born / rad_corr\");");
    Printf("══════════════════════════════════════════════════════════\n");
}