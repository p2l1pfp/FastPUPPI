#include <TH1F.h>
#include <TEfficiency.h>
#include <TGraph.h>

unsigned int fillTH1FastGenCut(TH1F *ret, unsigned int n, const float *corrArr, const float *genArr, float genThr) {
    unsigned int nsel = 0;
    for (unsigned int i = 0; i < n; ++i) {
        if (genArr[i] > genThr) { ret->Fill(corrArr[i]); nsel++; }
    }
    unsigned int nbins = ret->GetNbinsX();
    ret->SetBinContent(nbins, ret->GetBinContent(nbins) + ret->GetBinContent(nbins+1));
    ret->GetBinContent(nbins+1, 0);
    return nsel;
}

unsigned int fillTH1Fast(TH1F *ret, unsigned int n, const float *corrArr) {
    for (unsigned int i = 0; i < n; ++i) {
        ret->Fill(corrArr[i]);
    }
    unsigned int nbins = ret->GetNbinsX();
    ret->SetBinContent(nbins, ret->GetBinContent(nbins) + ret->GetBinContent(nbins+1));
    ret->GetBinContent(nbins+1, 0);
    return n;
}

void fillTEffFast(TEfficiency *eff, unsigned int n, const float *refArr, const float *corrArr, float corrThr) {
    for (unsigned int i = 0; i < n; ++i) {
        eff->Fill(corrArr[i] > corrThr, refArr[i]);
    }
}

TGraph *makeROCFast(TH1 *effsig, TH1 *effbkg) {
    TGraph *graph = new TGraph(effsig->GetNbinsX());
    for (unsigned int i = 1, n = graph->GetN()+1; i < n; ++i) {
        graph->SetPoint(i-1, effsig->GetBinContent(i), effbkg->GetBinContent(i));
    }
    return graph;
}
