#include <SuSyAachen/DiLeptonHistograms/interface/TtFullLepKinSolver.h>
#include "TF2.h"

TtFullLepKinSolver::TtFullLepKinSolver():
  mw(80.4),
  mb(4.8),
  pxmiss_(0),
  pymiss_(0),
  doSmearing(true)
{

}

TtFullLepKinSolver::TtFullLepKinSolver(const std::string ResolutionSmearingFileName, bool doSmearing, const double mW, const double mB):
  mw(mW),
  mb(mB),
  smearingFilename(ResolutionSmearingFileName),
  pxmiss_(0),
  pymiss_(0),
  doSmearing(doSmearing)
{
  
  if(doSmearing){
    TFile smearingFile(smearingFilename.c_str(), "read");
    leptonEnergyHist = (*(TH1F*)smearingFile.Get("leptonGenRecoERatioPlot"));
    leptonEnergyHist.SetDirectory(0);
    leptonEtaHist = *(TH1F*)smearingFile.Get("leptonGenRecoDEtaPlot");
    leptonEtaHist.SetDirectory(0);
    leptonPhiHist = *(TH1F*)smearingFile.Get("leptonGenRecoDPhiPlot");
    leptonPhiHist.SetDirectory(0);
    
    jetEnergyHist = *(TH1F*)smearingFile.Get("jetGenRecoERatioPlot");
    jetEnergyHist.SetDirectory(0);
    jetEtaHist = *(TH1F*)smearingFile.Get("jetGenRecoDEtaPlot");
    jetEtaHist.SetDirectory(0);
    jetPhiHist = *(TH1F*)smearingFile.Get("jetGenRecoDPhiPlot");
    jetPhiHist.SetDirectory(0);
    
    wMassHist = *(TH1F*)smearingFile.Get("genWMass");
    wMassHist.SetDirectory(0);
    lbMassHist = *(TH1F*)smearingFile.Get("genMLB");
    lbMassHist.SetDirectory(0);
    
    smearingFile.Close();
  }
  
}

//
// destructor
//
TtFullLepKinSolver::~TtFullLepKinSolver() 
{
}


void
TtFullLepKinSolver::SetConstraints(const double xx, const double yy)
{
  pxmiss_ = xx;
  pymiss_ = yy;
}

TtFullLepKinSolver::NeutrinoSolution
TtFullLepKinSolver::getNuSolution(const TLorentzVector& LV_l, 
          const TLorentzVector& LV_l_, 
          const TLorentzVector& LV_b, 
          const TLorentzVector& LV_b_,
          const TLorentzVector& met)
{ 
  TLorentzVector minLV_n  = TLorentzVector(0,0,0,0); 
  TLorentzVector minLV_n_ = TLorentzVector(0,0,0,0);   
  TLorentzVector minLV_b  = TLorentzVector(0,0,0,0); 
  TLorentzVector minLV_b_ = TLorentzVector(0,0,0,0);   
  TLorentzVector minLV_l  = TLorentzVector(0,0,0,0); 
  TLorentzVector minLV_l_ = TLorentzVector(0,0,0,0);   
  
  TLorentzVector smearedMet(0,0,0,0);
  
  std::vector<double> sols_weights;
  std::vector<TLorentzVector> sols_LV_n;
  std::vector<TLorentzVector> sols_LV_n_;
  std::vector<TLorentzVector> sols_LV_b;
  std::vector<TLorentzVector> sols_LV_b_;
  std::vector<TLorentzVector> sols_LV_l;
  std::vector<TLorentzVector> sols_LV_l_;
  std::vector<TLorentzVector> sols_LV_t;
  std::vector<TLorentzVector> sols_LV_t_;
  
  math::XYZTLorentzVector LV_avg_n(0,0,0,0);
  math::XYZTLorentzVector LV_avg_n_(0,0,0,0);

  int smearingsWithSols = 0;
  double weightSum = 0;
  
  double mt = 172.5;
  for(int i = 0; i < 100; i++) {
    double energyScale_l = leptonEnergyHist.GetRandom();
    double energyScale_l_ = leptonEnergyHist.GetRandom();
    double etaShift_l = leptonEtaHist.GetRandom();
    double etaShift_l_ = leptonEtaHist.GetRandom();
    double phiShift_l = leptonPhiHist.GetRandom();
    double phiShift_l_ = leptonPhiHist.GetRandom();
    
    double energyScale_b = jetEnergyHist.GetRandom();
    double energyScale_b_ = jetEnergyHist.GetRandom();
    double etaShift_b = jetEtaHist.GetRandom();
    double etaShift_b_ = jetEtaHist.GetRandom();
    double phiShift_b = jetPhiHist.GetRandom();
    double phiShift_b_ = jetPhiHist.GetRandom();
    
    double xconstraint = 0; double yconstraint = 0;
    
    
    mw = wMassHist.GetRandom();
    
    LV_L.SetPtEtaPhiE (LV_l.Pt() *energyScale_l,  LV_l.Eta()- etaShift_l,  LV_l.Phi()- phiShift_l,  LV_l.E() *energyScale_l);
    LV_L_.SetPtEtaPhiE(LV_l_.Pt()*energyScale_l_, LV_l_.Eta()-etaShift_l_, LV_l_.Phi()-phiShift_l_, LV_l_.E()*energyScale_l_);
    
    LV_B.SetPtEtaPhiE (LV_b.Pt() *energyScale_b,  LV_b.Eta() -etaShift_b,  LV_b.Phi() -phiShift_b,  LV_b.E() *energyScale_b);
    LV_B_.SetPtEtaPhiE(LV_b_.Pt()*energyScale_b_, LV_b_.Eta()-etaShift_b_, LV_b_.Phi()-phiShift_b_, LV_b_.E()*energyScale_b_);
    
    // recalculate MET
    smearedMet = met + (LV_l - LV_L) + (LV_l_ - LV_L_) + (LV_b - LV_B) + (LV_b_ - LV_B_);
    
    //xconstraint = smearedMet.Px()+LV_L.Px()+LV_L_.Px()+LV_B.Px()+LV_B_.Px();
    //yconstraint = smearedMet.Py()+LV_L.Py()+LV_L_.Py()+LV_B.Py()+LV_B_.Py();
    xconstraint = smearedMet.Px()+LV_l.Px()+LV_l_.Px()+LV_b.Px()+LV_b_.Px();
    yconstraint = smearedMet.Py()+LV_l.Py()+LV_l_.Py()+LV_b.Py()+LV_b_.Py();
    SetConstraints(xconstraint, yconstraint);
    
    double q_coeff[5], q_sol[4];   
    FindCoeff(LV_L, LV_L_, LV_B, LV_B_, mt, mt, pxmiss_, pymiss_, q_coeff);
    int NSol = quartic(q_coeff, q_sol);
    double ttmassMin = 10000;
    //loop on all solutions
    for (int isol = 0; isol < NSol; isol++) {
      TopRec(LV_L, LV_L_, LV_B, LV_B_, q_sol[isol]);
      double ttmass = WeightSolfromMass();
      if (ttmass < ttmassMin) {
        ttmassMin = ttmass;
        minLV_n.SetPxPyPzE(LV_n.Px(), LV_n.Py(), LV_n.Pz(), LV_n.E());
        minLV_n_.SetPxPyPzE(LV_n_.Px(), LV_n_.Py(), LV_n_.Pz(), LV_n_.E());
        minLV_b.SetPxPyPzE(LV_B.Px(), LV_B.Py(), LV_B.Pz(), LV_B.E());
        minLV_b_.SetPxPyPzE(LV_B_.Px(), LV_B_.Py(), LV_B_.Pz(), LV_B_.E());
        minLV_l.SetPxPyPzE(LV_L.Px(), LV_L.Py(), LV_L.Pz(), LV_L.E());
        minLV_l_.SetPxPyPzE(LV_L_.Px(), LV_L_.Py(), LV_L_.Pz(), LV_L_.E());
      }
    }
    
    if (NSol > 0){
      double massLB1 = (minLV_b+minLV_l_).M();
      double massLB2 = (minLV_b_+minLV_l).M();
      double weight = WeightSolfromMLB(massLB1, massLB2);
      smearingsWithSols++;
      weightSum += weight;
      sols_weights.push_back(weight);
      sols_LV_n.push_back (minLV_n);
      sols_LV_n_.push_back(minLV_n_);
      sols_LV_b.push_back (minLV_b);
      sols_LV_b_.push_back(minLV_b_);
      sols_LV_l.push_back (minLV_l);
      sols_LV_l_.push_back(minLV_l_);
      sols_LV_t.push_back (minLV_l+minLV_b_+minLV_n_);
      sols_LV_t_.push_back(minLV_l_+minLV_b+minLV_n);
    }
    
  }
  
  TLorentzVector LV_avg_t;
  TLorentzVector LV_avg_t_;
  
  TLorentzVector LV_sum_n;
  TLorentzVector LV_sum_n_;
  
  
  
  if (smearingsWithSols > 0){
    for(int i = 0; i < smearingsWithSols; i++){
      LV_avg_t  += (sols_weights[i]*sols_LV_t[i]);
      LV_avg_t_ += (sols_weights[i]*sols_LV_t_[i]);
      LV_sum_n  += (sols_weights[i]*sols_LV_n[i]);
      LV_sum_n_ += (sols_weights[i]*sols_LV_n_[i]);
    }
    
    
    LV_avg_t *= 1.0/weightSum;
    LV_avg_t_ *= 1.0/weightSum;
    LV_sum_n *= 1.0/weightSum;
    LV_sum_n_ *= 1.0/weightSum;
    LV_avg_t.SetE(TMath::Sqrt(LV_avg_t.P()*LV_avg_t.P()+mt*mt));
    LV_avg_t_.SetE(TMath::Sqrt(LV_avg_t_.P()*LV_avg_t_.P()+mt*mt));
    
    LV_avg_n.SetPxPyPzE(LV_sum_n.Px(), LV_sum_n.Py(), LV_sum_n.Pz(), LV_sum_n.E());
    LV_avg_n_.SetPxPyPzE(LV_sum_n_.Px(), LV_sum_n_.Py(), LV_sum_n_.Pz(), LV_sum_n_.E());
  }
  
  TtFullLepKinSolver::NeutrinoSolution nuSol;
  nuSol.neutrino    = reco::LeafCandidate(0, LV_avg_n  );
  nuSol.neutrinoBar = reco::LeafCandidate(0, LV_avg_n_ ); 
  nuSol.LV_t = LV_avg_t;
  nuSol.LV_t_ = LV_avg_t_;
  nuSol.weight = weightSum; 
  nuSol.smearingsWithSols = smearingsWithSols;
  return nuSol;
}

void
TtFullLepKinSolver::FindCoeff(const TLorentzVector& al, 
            const TLorentzVector& l,
            const TLorentzVector& b_al,
            const TLorentzVector& b_l,
            const double mt, 
            const double mat, 
            const double px_miss, 
            const double py_miss,
            double* koeficienty)
{
  double E, apom1, apom2, apom3;
  double k11, k21, k31, k41,  cpom1, cpom2, cpom3, l11, l21, l31, l41, l51, l61, k1, k2, k3, k4, k5,k6;
  double l1, l2, l3, l4, l5, l6, k15, k25, k35, k45;

  C = -al.Px()-b_al.Px()-l.Px()-b_l.Px() + px_miss;
  D = -al.Py()-b_al.Py()-l.Py()-b_l.Py() + py_miss;

  // right side of first two linear equations - missing pT
  
  E = (sqr(mt)-sqr(mw)-sqr(mb))/(2*b_al.E())-sqr(mw)/(2*al.E())-al.E()+al.Px()*b_al.Px()/b_al.E()+al.Py()*b_al.Py()/b_al.E()+al.Pz()*b_al.Pz()/b_al.E();
  F = (sqr(mat)-sqr(mw)-sqr(mb))/(2*b_l.E())-sqr(mw)/(2*l.E())-l.E()+l.Px()*b_l.Px()/b_l.E()+l.Py()*b_l.Py()/b_l.E()+l.Pz()*b_l.Pz()/b_l.E();
  
  m1 = al.Px()/al.E()-b_al.Px()/b_al.E();
  m2 = al.Py()/al.E()-b_al.Py()/b_al.E();
  m3 = al.Pz()/al.E()-b_al.Pz()/b_al.E();
  
  n1 = l.Px()/l.E()-b_l.Px()/b_l.E();
  n2 = l.Py()/l.E()-b_l.Py()/b_l.E();
  n3 = l.Pz()/l.E()-b_l.Pz()/b_l.E();
  
  pom = E-m1*C-m2*D;
  apom1 = sqr(al.Px())-sqr(al.E());
  apom2 = sqr(al.Py())-sqr(al.E());
  apom3 = sqr(al.Pz())-sqr(al.E());
  
  k11 = 1/sqr(al.E())*(pow(mw,4)/4+sqr(C)*apom1+sqr(D)*apom2+apom3*sqr(pom)/sqr(m3)+sqr(mw)*(al.Px()*C+al.Py()*D+al.Pz()*pom/m3)+2*al.Px()*al.Py()*C*D+2*al.Px()*al.Pz()*C*pom/m3+2*al.Py()*al.Pz()*D*pom/m3);
  k21 = 1/sqr(al.E())*(-2*C*m3*n3*apom1+2*apom3*n3*m1*pom/m3-sqr(mw)*m3*n3*al.Px()+sqr(mw)*m1*n3*al.Pz()-2*al.Px()*al.Py()*D*m3*n3+2*al.Px()*al.Pz()*C*m1*n3-2*al.Px()*al.Pz()*n3*pom+2*al.Py()*al.Pz()*D*m1*n3);
  k31 = 1/sqr(al.E())*(-2*D*m3*n3*apom2+2*apom3*n3*m2*pom/m3-sqr(mw)*m3*n3*al.Py()+sqr(mw)*m2*n3*al.Pz()-2*al.Px()*al.Py()*C*m3*n3+2*al.Px()*al.Pz()*C*m2*n3-2*al.Py()*al.Pz()*n3*pom+2*al.Py()*al.Pz()*D*m2*n3);
  k41 = 1/sqr(al.E())*(2*apom3*m1*m2*sqr(n3)+2*al.Px()*al.Py()*sqr(m3)*sqr(n3)-2*al.Px()*al.Pz()*m2*m3*sqr(n3)-2*al.Py()*al.Pz()*m1*m3*sqr(n3));
  k51 = 1/sqr(al.E())*(apom1*sqr(m3)*sqr(n3)+apom3*sqr(m1)*sqr(n3)-2*al.Px()*al.Pz()*m1*m3*sqr(n3));
  k61 = 1/sqr(al.E())*(apom2*sqr(m3)*sqr(n3)+apom3*sqr(m2)*sqr(n3)-2*al.Py()*al.Pz()*m2*m3*sqr(n3));
  
  cpom1 = sqr(l.Px())-sqr(l.E());
  cpom2 = sqr(l.Py())-sqr(l.E());
  cpom3 = sqr(l.Pz())-sqr(l.E());
  
  l11 = 1/sqr(l.E())*(pow(mw,4)/4+cpom3*sqr(F)/sqr(n3)+sqr(mw)*l.Pz()*F/n3);
  l21 = 1/sqr(l.E())*(-2*cpom3*F*m3*n1/n3+sqr(mw)*(l.Px()*m3*n3-l.Pz()*n1*m3)+2*l.Px()*l.Pz()*F*m3);
  l31 = 1/sqr(l.E())*(-2*cpom3*F*m3*n2/n3+sqr(mw)*(l.Py()*m3*n3-l.Pz()*n2*m3)+2*l.Py()*l.Pz()*F*m3);
  l41 = 1/sqr(l.E())*(2*cpom3*n1*n2*sqr(m3)+2*l.Px()*l.Py()*sqr(m3)*sqr(n3)-2*l.Px()*l.Pz()*n2*n3*sqr(m3)-2*l.Py()*l.Pz()*n1*n3*sqr(m3));
  l51 = 1/sqr(l.E())*(cpom1*sqr(m3)*sqr(n3)+cpom3*sqr(n1)*sqr(m3)-2*l.Px()*l.Pz()*n1*n3*sqr(m3));
  l61 = 1/sqr(l.E())*(cpom2*sqr(m3)*sqr(n3)+cpom3*sqr(n2)*sqr(m3)-2*l.Py()*l.Pz()*n2*n3*sqr(m3));
  
  k1 = k11*k61;
  k2 = k61*k21/k51;
  k3 = k31;
  k4 = k41/k51;
  k5 = k61/k51;
  k6 = 1;
  
  l1 = l11*k61;
  l2 = l21*k61/k51;
  l3 = l31;
  l4 = l41/k51;
  l5 = l51*k61/(sqr(k51));
  l6 = l61/k61;
  
  k15 = k1*l5-l1*k5;
  k25 = k2*l5-l2*k5;
  k35 = k3*l5-l3*k5;
  k45 = k4*l5-l4*k5;
  
  k16 = k1*l6-l1*k6;
  k26 = k2*l6-l2*k6;
  k36 = k3*l6-l3*k6;
  k46 = k4*l6-l4*k6;
  k56 = k5*l6-l5*k6;

  koeficienty[0] = k15*sqr(k36)-k35*k36*k16-k56*sqr(k16);
  koeficienty[1] = 2*k15*k36*k46+k25*sqr(k36)+k35*(-k46*k16-k36*k26)-k45*k36*k16-2*k56*k26*k16;
  koeficienty[2] = k15*sqr(k46)+2*k25*k36*k46+k35*(-k46*k26-k36*k56)-k56*(sqr(k26)+2*k56*k16)-k45*(k46*k16+k36*k26);
  koeficienty[3] = k25*sqr(k46)-k35*k46*k56-k45*(k46*k26+k36*k56)-2*sqr(k56)*k26;
  koeficienty[4] = -k45*k46*k56-pow(k56,3);
  
  // normalization of coefficients
  int moc=(int(log10(fabs(koeficienty[0])))+int(log10(fabs(koeficienty[4]))))/2;
  
  koeficienty[0]=koeficienty[0]/TMath::Power(10,moc);
  koeficienty[1]=koeficienty[1]/TMath::Power(10,moc);
  koeficienty[2]=koeficienty[2]/TMath::Power(10,moc);
  koeficienty[3]=koeficienty[3]/TMath::Power(10,moc);
  koeficienty[4]=koeficienty[4]/TMath::Power(10,moc);
}

void TtFullLepKinSolver::TopRec(const TLorentzVector& al, 
                                const TLorentzVector& l,
                          const TLorentzVector& b_al,
                          const TLorentzVector& b_l, 
        const double sol)
{
  TVector3 t_ttboost;
  TLorentzVector aux;
  double pxp, pyp, pzp, pup, pvp, pwp;
    
  pxp = sol*(m3*n3/k51);   
  pyp = -(m3*n3/k61)*(k56*pow(sol,2) + k26*sol + k16)/(k36 + k46*sol);
  pzp = -1/n3*(n1*pxp + n2*pyp - F);
  pwp = 1/m3*(m1*pxp + m2*pyp + pom);
  pup = C - pxp;
  pvp = D - pyp;
     
  LV_n_.SetXYZM(pxp, pyp, pzp, 0.0);
  LV_n.SetXYZM(pup, pvp, pwp, 0.0);
  
  LV_t_ = b_l + l + LV_n_;
  LV_t = b_al + al + LV_n;  
 
  aux = (LV_t_ + LV_t);
  t_ttboost = -aux.BoostVector();
  LV_tt_t_ = LV_t_;
  LV_tt_t = LV_t;
  LV_tt_t_.Boost(t_ttboost);
  LV_tt_t.Boost(t_ttboost); 
}

double
TtFullLepKinSolver::WeightSolfromMass() const
{
  return (LV_n + LV_n_+ LV_L + LV_L_ + LV_B + LV_B_).M();
}

double
TtFullLepKinSolver::WeightSolfromMLB(double mlb1, double mlb2)
{ 
  double integral = lbMassHist.Integral();
  if(integral > 0){
    return lbMassHist.GetBinContent(lbMassHist.FindBin(mlb1))*lbMassHist.GetBinContent(lbMassHist.FindBin(mlb2))/integral/integral;
  }else{
    return 1;
  }
  
} 
         
int
TtFullLepKinSolver::quartic(double *koeficienty, double* koreny) const
{ 
  double w, b0, b1, b2;
  double c[4];
  double d0, d1, h, t, z;
  double *px;
 
  if (koeficienty[4]==0.0) 
    return cubic(koeficienty, koreny);
  /* quartic problem? */
  w = koeficienty[3]/(4*koeficienty[4]);
  /* offset */
  b2 = -6*sqr(w) + koeficienty[2]/koeficienty[4];
  /* koeficienty. of shifted polynomial */
  b1 = (8*sqr(w) - 2*koeficienty[2]/koeficienty[4])*w + koeficienty[1]/koeficienty[4];
  b0 = ((-3*sqr(w) + koeficienty[2]/koeficienty[4])*w - koeficienty[1]/koeficienty[4])*w + koeficienty[0]/koeficienty[4];

  c[3] = 1.0;
  /* cubic resolvent */
  c[2] = b2;
  c[1] = -4*b0;
  c[0] = sqr(b1) - 4*b0*b2;
  
  cubic(c, koreny);
  z = koreny[0];
  //double z1=1.0,z2=2.0,z3=3.0;
  //TMath::RootsCubic(c,z1,z2,z3);
  //if (z2 !=0) z = z2;
  //if (z1 !=0) z = z1;
  /* only lowermost root needed */

  int nreal = 0;
  px = koreny;
  t = sqrt(0.25*sqr(z) - b0);
  for(int i=-1; i<=1; i+=2) {
    d0 = -0.5*z + i*t;
    /* coeffs. of quadratic factor */
    d1 = (t!=0.0)? -i*0.5*b1/t : i*sqrt(-z - b2);
    h = 0.25*sqr(d1) - d0;
    if (h>=0.0) {
      h = sqrt(h);
      nreal += 2;
      *px++ = -0.5*d1 - h - w;
      *px++ = -0.5*d1 + h - w;
    }
  }

  //  if (nreal==4) {
    /* sort results */
//    if (koreny[2]<koreny[0]) SWAP(koreny[0], koreny[2]);
//    if (koreny[3]<koreny[1]) SWAP(koreny[1], koreny[3]);
//    if (koreny[1]<koreny[0]) SWAP(koreny[0], koreny[1]);
//    if (koreny[3]<koreny[2]) SWAP(koreny[2], koreny[3]);
//    if (koreny[2]<koreny[1]) SWAP(koreny[1], koreny[2]);
//  }
  return nreal;

}

int
TtFullLepKinSolver::cubic(const double *coeffs, double* koreny) const
{
  unsigned nreal;
  double w, p, q, dis, h, phi;
  
  if (coeffs[3]!=0.0) {
    /* cubic problem? */
    w = coeffs[2]/(3*coeffs[3]);
    p = sqr(coeffs[1]/(3*coeffs[3])-sqr(w))*(coeffs[1]/(3*coeffs[3])-sqr(w));
    q = -0.5*(2*sqr(w)*w-(coeffs[1]*w-coeffs[0])/coeffs[3]);
    dis = sqr(q)+p;
    /* discriminant */
    if (dis<0.0) {
      /* 3 real solutions */
      h = q/sqrt(-p);
      if (h>1.0) h = 1.0;
      /* confine the argument of */
      if (h<-1.0) h = -1.0;
      /* acos to [-1;+1] */
      phi = acos(h);
      p = 2*TMath::Power(-p, 1.0/6.0);
      for(unsigned i=0; i<3; i++) 
  koreny[i] = p*cos((phi+2*i*TMath::Pi())/3.0) - w;
      if (koreny[1]<koreny[0]) SWAP(koreny[0], koreny[1]);
      /* sort results */
      if (koreny[2]<koreny[1]) SWAP(koreny[1], koreny[2]);
      if (koreny[1]<koreny[0]) SWAP(koreny[0], koreny[1]);
      nreal = 3;
    }
    else {
      /* only one real solution */
      dis = sqrt(dis);
      h = TMath::Power(fabs(q+dis), 1.0/3.0);
      p = TMath::Power(fabs(q-dis), 1.0/3.0);
      koreny[0] = ((q+dis>0.0)? h : -h) + ((q-dis>0.0)? p : -p) -  w;
      nreal = 1;
    }

    /* Perform one step of a Newton iteration in order to minimize
       round-off errors */
    for(unsigned i=0; i<nreal; i++) {
      h = coeffs[1] + koreny[i] * (2 * coeffs[2] + 3 * koreny[i] * coeffs[3]);
      if (h != 0.0)
  koreny[i] -= (coeffs[0] + koreny[i] * (coeffs[1] + koreny[i] * (coeffs[2] + koreny[i] * coeffs[3])))/h;
    }
  }

  else if (coeffs[2]!=0.0) {
    /* quadratic problem? */
    p = 0.5*coeffs[1]/coeffs[2];
    dis = sqr(p) - coeffs[0]/coeffs[2];
    if (dis>=0.0) {
      /* two real solutions */
      dis = sqrt(dis);
      koreny[0] = -p - dis;
      koreny[1] = -p + dis;
      nreal = 2;
    }
    else
      /* no real solution */
      nreal = 0;
  }

  else if (coeffs[1]!=0.0) {
    /* linear problem? */
    koreny[0] = -coeffs[0]/coeffs[1];
    nreal = 1;
  }

  else
    /* no equation */
    nreal = 0;
  
  return nreal;
}


void
TtFullLepKinSolver::SWAP(double& realone, double& realtwo) const
{
  if (realtwo < realone) {
    double aux = realtwo;
    realtwo = realone;
    realone = aux;
  }
}
