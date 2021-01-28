#ifndef FUNCTIONS
#define FUNCTIONS

#include "TLorentzVector.h"
#include "TVector2.h"
#include "TF1.h"
#include "TMath.h"
#include "TH1F.h"
#include <iostream>

TLorentzVector VVec(float lep0_pt, float lep0_eta, float lep0_phi, float lep0_m,
                    float lep1_pt, float lep1_eta, float lep1_phi, float lep1_m)
{
    TLorentzVector vlep0;
    vlep0.SetPtEtaPhiM(lep0_pt, lep0_eta, lep0_phi, lep0_m);
    TLorentzVector vlep1;
    vlep1.SetPtEtaPhiM(lep1_pt, lep1_eta, lep1_phi, lep1_m);

    return (vlep0 + vlep1);
}

TLorentzVector VLep(float lep0_pt, float lep0_eta, float lep0_phi, float lep0_m){
    TLorentzVector vlep0;
    vlep0.SetPtEtaPhiM(lep0_pt, lep0_eta, lep0_phi, lep0_m);
    return vlep0;
}

#endif
