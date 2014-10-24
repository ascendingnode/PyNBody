#ifndef __RKN670_HPP__
#define __RKN670_HPP__

#include <cmath>
#include <vector>
#include <cstdio>

template<class doubtype,class state> class RKN670 { public:
   bool onestep,first,loop,high,last,phase2,inefficient;
   int failed,nsuccess,efficient2,efficient1,nfailed,NATT,neq,nstage,finish_phase2;
   doubtype H,uround,R,RMN1,RMN2,RMN3,tol,exponent,STBRAD,scale,err,direction,rang,
            HUSED,HSTART,TOLD,Tcheck,toosmall, olddir;
   int pointer[17];
   doubtype C[18],A[18][18],BHAT[18],BPHAT[18],B[18],BP[18],BSTAR[7],BPSTAR[7];
   std::vector<doubtype> NTHR,NTHRP,NWT,NWTP,NRT1,NRT2,NRT3,NRT4;
   std::vector<std::vector<doubtype> > JS;

   RKN670() {  setup_hi(1,1e-5); }

   RKN670(int n0, double tol0, bool high0=true) {
       if(high0) setup_hi(n0,tol0);
       else setup_64(n0,tol0);
   }

   void integrate(state &s,doubtype TEND) {
      
      using namespace std;

      if(s.T==TEND) return;

      failed = 0;
      rang = fabs(TEND-s.T);
      if (TEND-s.T<0.) direction = -1.0;
      else             direction =  1.0;
      if(direction!=olddir) {
         first = true;
         H = 0;
      }

      if(first) {
         s.dt(s.T,s.r,s.a);

         if (H==0.) {
            doubtype YPNRM = 0.;
            doubtype YDPNRM = 0.;
            bool DFLT = false;
            bool DFLTP = false;
            phase2 = true;
            for (int K(0);K<neq;K++) {
               doubtype THRES = NTHR[K];
               doubtype WT = max(THRES,fabs(s.r[K]));
               NWT[K] = WT;
               doubtype TRY = fabs(s.v[K]/WT);
               if (TRY>YPNRM) {
                  YPNRM = TRY;
                  DFLT = (WT==THRES);
               }
               doubtype THRESP = NTHRP[K];
               doubtype WTP = max(THRESP,fabs(s.v[K]));
               NWTP[K] = WTP;
               TRY = fabs(s.a[K]/WTP);
               if (TRY>YDPNRM) {
                  YDPNRM = TRY;
                  DFLTP = (WTP==THRESP);
               }
            }
            doubtype ABSH = rang;
            doubtype ABSHP = rang;
            doubtype TAU = pow(tol,exponent);
            if (ABSH*YPNRM>TAU and not DFLT)    ABSH  = TAU/YPNRM;
            if (ABSHP*YDPNRM>TAU and not DFLTP) ABSHP = TAU/YDPNRM;
            H = direction*min(ABSH,ABSHP);
         }
         else
            phase2 = false;
      }

      bool success = false;
      bool sfailed = false;
      bool optimal_reduct = false;
      doubtype alpha = 0.0;
      while (true) {
      
         if (direction*(s.T+H-TEND)>=0.) {
            inefficient = ((TEND-s.T)/H<0.5);
            H = TEND - s.T;
            last = true;
         }
         else {
            last = false;
            if (direction*(s.T+2.*H-TEND)>=0.) {
               H = (TEND-s.T)*0.5;
            }
         }

         if (fabs(H)<toosmall*max(fabs(s.T),fabs(s.T+H))) {
            //printf("# timestep too small! H = %0.3e\n",double(H));
            failed = 3;
            return;
         }

         if (not s.is_sane()) {
             failed = 8;
             return;
         }

         doubtype HSQ = H*H;
         for (int K(0);K<neq;K++)
            NRT2[K] = BPHAT[1]*s.a[K];

         for (int I(2);I<nstage+1;I++) {

            for (int K(0);K<neq;K++) {
               NRT1[K] = A[I][1]*s.a[K];
               for (int J(2);J<I;J++)
                  //if (A[I][J]!=0.0)
                     NRT1[K] += A[I][J]*JS[pointer[J-1]][K];
               NRT3[K] = s.r[K] + C[I]*H*s.v[K] + HSQ*NRT1[K];
            }

            int JSTAGE = pointer[I-1];
            s.dt(s.T+C[I]*H,NRT3,JS[JSTAGE]);

            if (BPHAT[I]!=0.)
               for (int K(0);K<neq;K++)
                  NRT2[K] += BPHAT[I]*JS[JSTAGE][K];

            if (phase2) {
               doubtype YPNRM=0.,YDPNRM=0.,DEL=0.,FDEL=0.,YPDEL=0.,YDPDEL=0.;
               for (int K(0);K<neq;K++) {
                  doubtype YINTK = NRT3[K];
                  doubtype YPINTK = s.v[K] + H*NRT2[K];
                  doubtype WT = max(NWT[K],fabs(YINTK));
                  doubtype WTP = max(NWTP[K],fabs(YPINTK));
                  NWT[K] = WT;
                  NWTP[K] = WTP;
                  YPNRM = max(YPNRM,max(fabs(YINTK)/WT,fabs(s.r[K])/WT));
                  YDPNRM = max(YDPNRM,max(fabs(YPINTK)/WTP,fabs(s.v[K])/WTP));
                  DEL = max(DEL,fabs(YINTK-s.r[K])/WT);
                  FDEL = max(FDEL,fabs(YPINTK-s.v[K])/WT);
                  YPDEL = max(YPDEL,fabs(YPINTK-s.v[K])/WTP);
                  YDPDEL = max(YDPDEL,fabs(JS[JSTAGE][K]-s.a[K])/WTP);
               }
               DEL = max(DEL,fabs(C[I]*H)/rang);
               doubtype SCL = 1.0, SCLP = 1.0;
               if (DEL>10.0*uround*YPNRM)
                  if (fabs(H)*FDEL>STBRAD*DEL)
                     SCL = STBRAD/R*max(DEL/(fabs(H)*FDEL),RMN3);
               if (YPDEL>10.0*uround*YDPNRM)
                  if (fabs(H)*YDPDEL>STBRAD*YPDEL)
                     SCLP = STBRAD/R*max(YPDEL/(fabs(H)*YDPDEL),RMN3);
               doubtype scale1 = min(SCL,SCLP);
               if (scale1<1.0) {
                  NATT = NATT + 1;
                  H = scale1*fabs(H);
                  H = direction*H;
                  last = false;
                  for (int K(0);K<neq;K++) {
                     NWT[K] = max(fabs(s.r[K]),NTHR[K]);
                     NWTP[K] = max(fabs(s.v[K]),NTHRP[K]);
                  }
                  continue;
               }
               phase2 = (I<finish_phase2);
            }
         }

         for (int K(0);K<neq;K++) {
            NRT4[K] = s.v[K] + H*NRT2[K];
            if (high) {
               NRT1[K] = BHAT[1]*s.a[K];
               for (int J(7);J<16;J++) 
                  NRT1[K] += BHAT[J]*JS[pointer[J-1]][K];
               NRT3[K] = s.r[K] + H*s.v[K] + HSQ*NRT1[K];
            }
         }

         err = 0.0;
         doubtype errp = 0.0;
         for (int K(0);K<neq;K++) {
            NRT1[K] -= B[1]*s.a[K];
            NRT2[K] -= BP[1]*s.a[K];
            if (not high) {
               for (int I(2);I<4;I++)
                  NRT1[K] -= B[I] *JS[pointer[I-1]][K];
               for (int I(2);I<nstage+1;I++)
                  NRT2[K] -= BP[I]*JS[pointer[I-1]][K];
            } else 
               for(int I(7);I<17;I++) {
                  NRT1[K] -= B[I] *JS[pointer[I-1]][K];
                  NRT2[K] -= BP[I]*JS[pointer[I-1]][K];
               }
            doubtype SK = H*H*NRT1[K];
            doubtype SPK = H*NRT2[K];
            doubtype WT = 0.5*fabs(s.r[K])+0.25*fabs(NRT3[K])+0.25*fabs(NRT3[K]-SK);
            WT = max(WT,NTHR[K]);
            doubtype WTP = 0.5*fabs(s.v[K])+0.25*fabs(NRT4[K])+0.25*fabs(NRT4[K]-SPK);
            WTP = max(WTP,NTHRP[K]);
            err = max(err,fabs(SK)/WT);
            errp = max(errp,fabs(SPK)/WTP);
         }
         err = max(err,errp);
         
         if (err>tol) {
            last = false;
            if (first) {
               NATT = NATT + 1;
               if (success) {
                  optimal_reduct = true;
                  alpha = pow(scale*(tol/err),exponent);
                  alpha = max(RMN2,alpha);
               }
               else {
                  alpha = RMN1;
                  phase2 = true;
                  for (int K(0);K<neq;K++) {
                     NWT[K] = max(fabs(s.r[K]),NTHR[K]);
                     NWTP[K] = max(fabs(s.v[K]),NTHRP[K]);
                  }
               }
            }
            else {
               nfailed = nfailed + 1;
               if (sfailed)
                  alpha = 0.5;
               else {
                  sfailed = true;
                  alpha = pow(scale*(tol/err),exponent);
                  alpha = max(RMN1,alpha);
               }
            }
            H = alpha*H;
            continue;
         }

         success = true;
         doubtype beta = (pow(err/tol,exponent))/scale;
         doubtype one = 1.0;
         if (first)
            alpha = one/max(RMN3,beta);
         else {
            alpha = one/max(RMN1,beta);
            if (sfailed) alpha = min(one,alpha);
         }
         HUSED = H;
         H = alpha*H;

         if (first and not last and not optimal_reduct and alpha>R) {
            HUSED = 0.;
            NATT = NATT + 1;
            last = false;
            continue;
         }

         if (first) HSTART = HUSED;
         first = false;
         TOLD = s.T;
         s.T += HUSED;
         efficient1 = efficient1 + 1;
         if (last) {
            s.T = TEND;
            if (inefficient)
               efficient2 = efficient2 + 1;
            if (efficient1>100) 
               if (efficient2>0.1*efficient1) {
                  //failed = 4;
                  efficient1 = 0;
                  efficient2 = 0;
               }
         }
         nsuccess = nsuccess + 1;

         NRT1 = s.r; s.r = NRT3;
         NRT2 = s.v; s.v = NRT4;
         NRT3 = s.a;

         if (high) s.dt(s.T,s.r,s.a);
         else      s.a = JS[pointer[5]];

         olddir = direction;

         /*if (s.is_sane()>0) {
            failed = 8;
            success = false;
            sfailed = true;
            //break;
         }*/

         if (last or onestep) break;
         success = false;
         sfailed = false;
      }

      if (failed>0) {
         if (not loop) {
            loop = true;
            Tcheck = s.T;
         }
         else {
            if (Tcheck==s.T) {
               printf(" RKNINT STOPPED DUE TO TWO SUCCESSIVE FAILURES AT");
               printf(" THE SAME VALUE OF T (= %12.5f)\n",double(s.T));
               return;
            }
         }
      }
      else 
         loop = false;
   }

   void evaluate64(const state &s,doubtype TWANT,state &s2) {
      s2 = s;
      s2.T = TWANT;
      if (high) {
         failed = 3;
         return;
      }
      failed = 0;
      if (TWANT==s.T) return;
      if (TWANT==TOLD) {
         s2.r = NRT1; s2.v = NRT2; return; }
      bool INTERP;
      if (HUSED>0.0) INTERP = ((TOLD<TWANT) and (TWANT<s.T));
      else           INTERP = ((TOLD>TWANT) and (TWANT>s.T));
      if (not INTERP) failed = 1;
      doubtype HSIG = TWANT - TOLD;
      doubtype S = HSIG/HUSED;
      
      BSTAR[1] = (((900.0*S-3819.0)*S+6386.0)*S-5244.0)*S + 2106.0;
      BSTAR[1] = BSTAR[1]/4212.0;
      BPSTAR[1] = (((5400.0*S-19095.0)*S+25544.0)*S-15732.0)*S + 4212.0;
      BPSTAR[1] = BPSTAR[1]/4212.0;
      BSTAR[2] = 0.;
      BPSTAR[2] = 0.;

      doubtype U = -.365774278775981018074400638167e-3;
      doubtype V = -.423066852735158392593161679389e0;
      doubtype W = -.357765287229492094385275612595e1;
      doubtype X = .110329844381694622139549133756e0;
      BSTAR[3] = (((1800.0*U*S-6.0*V)*S+W)*S+18.0*X)*S;
      U = -.731548557551962036148801276334e-3;
      V = -.846133705470316785186323358781e0;
      W = -.715530574458984188770551225188e1;
      X = .220659688763389244279098267512e0;
      BPSTAR[3] = (((5400.0*U*S-15.0*V)*S+2.0*W)*S+27.0*X)*S;

      U = .886060093179442480359392334923e-3;
      V = .860688925651348955842548458187e0;
      W = .554758675105232911302563118117e1;
      X = -.103107545982925635135268011828e0;
      BSTAR[4] = (((1800.0*U*S-6.0*V)*S+W)*S+18.0*X)*S;
      U = .177212018635888496071878466985e-2;
      V = .172137785130269791168509691638e1;
      W = .110951735021046582260512623623e2;
      X = -.206215091965851270270536023656e0;
      BPSTAR[4] = (((5400.0*U*S-15.0*V)*S+2.0*W)*S+27.0*X)*S;

      BSTAR[5] = ((225.0*S-651.0)*S+620.0)*S - 195.0;
      BSTAR[5] = -200.0*S*BSTAR[5]/17901.0;
      BPSTAR[5] = ((270.0*S-651.0)*S+496.0)*S - 117.0;
      BPSTAR[5] = -1000.0*S*BPSTAR[5]/17901.0;
      BSTAR[6] = S*(S-1.0)*((300.0*S-523.0)*S+234.0)/220.0;
      BPSTAR[6] = S*(((1800.0*S-4115.0)*S+3028.0)*S-702.0)/220.0;
   
      for(int K(0);K<neq;K++) {
         doubtype SUM = BSTAR[1]*NRT3[K];
         doubtype SUMP = BPSTAR[1]*NRT3[K];
         for(int I(3);I<=6;I++) {
            SUM  = SUM  +  BSTAR[I]*JS[I-2][K];
            SUMP = SUMP + BPSTAR[I]*JS[I-2][K];
         }
         s2.r[K] = NRT1[K] + HSIG*(NRT2[K]+HSIG*SUM);
         s2.v[K] = NRT2[K] + HSIG*SUMP;
      }
   }

   void init(int n,doubtype stol) {
      olddir = 1.0;
      onestep = false;
      uround = 3.0e-17;
      H = 0.;
      doubtype thrdef=50.0*uround;
      neq = n;
      NTHR = std::vector<doubtype>(neq,thrdef);
      NTHRP = std::vector<doubtype>(neq,thrdef);
      NWT = std::vector<doubtype>(neq,0.0);
      NWTP = std::vector<doubtype>(neq,0.0);
      first = true;
      R=10.0;
      RMN1=1.0/R;
      RMN2=1.0/(R*R);
      RMN3=1.0/(R*R*R);
      NATT = 0;
      nfailed = 0;
      efficient1 = 0;
      efficient2 = 0;
      nsuccess = 0;
      tol = stol;
      loop = false;
      failed = 0;
      NRT1 = std::vector<doubtype>(neq,0.0);
      NRT2 = std::vector<doubtype>(neq,0.0);
      NRT3 = std::vector<doubtype>(neq,0.0);
      NRT4 = std::vector<doubtype>(neq,0.0);
      JS.resize(14);
      for(auto it=JS.begin();it!=JS.end();++it)
         *it = std::vector<doubtype>(neq,0.0);
   }

   void setup_64(int n,doubtype stol) {
      this->init(n,stol);

      nstage = 6;
      exponent = 0.2;
      high = false;
      STBRAD = 4.;
      finish_phase2 = 4;
      scale = 0.8;
      toosmall = 20.*uround;

      pointer[1] = 0;
      pointer[2] = 1;
      pointer[3] = 2;
      pointer[4] = 3;
      pointer[5] = 4;

      C[1] = 0.0;
      C[2] = 1.29295903136704415288990053209e-1;
      C[3] = 2.58591806273408830577980106418e-1;
      C[4] = 6.70297082615480058310908782471e-1;
      C[5] = 9.0e-1;
      C[6] = 1.0;

      A[2][1] = 8.35871528396802532822102346743e-3;
      A[3][1] = 1.11449537119573671042946979566e-2;
      A[3][2] = 2.22899074239147342085893959132e-2;
      A[4][1] = 1.45474742801091785895935232316e-1;
      A[4][2] = -2.29860640522647473120261456297e-1;
      A[4][3] = 3.09034987202967536528726080729e-1;
      A[5][1] = -2.07668262950789954335146205252e-1;
      A[5][2] = 6.86366784292514312273571851042e-1;
      A[5][3] = -1.99549277872349252201326358888e-1;
      A[5][4] = 1.25850756530624894262900713098e-1;
      A[6][1] = 7.81101614434947768281101614435e-2;
      A[6][2] = 0.0;
      A[6][3] = 2.882917411897667776841772093e-1;
      A[6][4] = 1.22425537174570410182422422005e-1;
      A[6][5] = 1.1172560192168035305290207251e-2;

      BPHAT[1] = 7.81101614434947768281101614435e-2;
      BPHAT[2] = 0.0;
      BPHAT[3] = 3.88843478705982602715952208523e-1;
      BPHAT[4] = 3.71320757928842267403035557523e-1;
      BPHAT[5] = 1.1172560192168035305290207251e-1;
      BPHAT[6] = 5.0e-2;

      B[1] = 1.05885926037041827822001005886;
      B[2] = -2.40675137192445205319176777632;
      B[3] = 1.84789211155403377497175771746;
      B[4] = 0.0;
      B[5] = 0.0;
      B[6] = 0.0;

      BP[1] = 5.46058879392212725546058879392e-2;
      BP[2] = 0.0;
      BP[3] = 4.6126678590362684429468353488e-1;
      BP[4] = 1.95880859479312656291875875208e-1;
      BP[5] = 3.88246466677839226858834701972e-1;
      BP[6] = -1.0e-1;
   }

   void setup_hi(int n,doubtype tol) {
      this->init(n,tol);

      nstage = 17;
      exponent = 1.0/11.0;
      high = true;
      STBRAD = 8.;
      finish_phase2 = 14;
      scale = 0.75;
      toosmall = 200.*uround;

      pointer[1] = 0;
      pointer[2] = 1;
      pointer[3] = 2;
      pointer[4] = 0;
      pointer[5] = 3;
      pointer[6] = 4;
      pointer[7] = 5;
      pointer[8] = 6;
      pointer[9] = 7;
      pointer[10] = 8;
      pointer[11] = 9;
      pointer[12] = 10;
      pointer[13] = 11;
      pointer[14] = 12;
      pointer[15] = 13;
      pointer[16] = 3;

      C[1] = 0.0;
      C[2] = 2.0e-2;
      C[3] = 4.0e-2;
      C[4] = 1.0e-1;
      C[5] = 1.33333333333333333333333333333e-1;
      C[6] = 1.6e-1;
      C[7] = 5.0e-2;
      C[8] = 2.0e-1;
      C[9] = 2.5e-1;
      C[10] = 3.33333333333333333333333333333e-1;
      C[11] = 5.0e-1;
      C[12] = 5.55555555555555555555555555556e-1;
      C[13] = 7.5e-1;
      C[14] = 8.57142857142857142857142857143e-1;
      C[15] = 9.45216222272014340129957427739e-1;
      C[16] = 1.0;
      C[17] = 1.0;

      A[2][1] = 2.0e-4;

      A[3][1] = 2.66666666666666666666666666667e-4;
      A[3][2] = 5.33333333333333333333333333333e-4;

      A[4][1] = 2.91666666666666666666666666667e-3;
      A[4][2] = -4.16666666666666666666666666667e-3;
      A[4][3] = 6.25e-3;

      A[5][1] = 1.64609053497942386831275720165e-3;
      A[5][2] = 0.0;
      A[5][3] = 5.48696844993141289437585733882e-3;
      A[5][4] = 1.75582990397805212620027434842e-3;

      A[6][1] = 1.9456e-3;
      A[6][2] = 0.0;
      A[6][3] = 7.15174603174603174603174603175e-3;
      A[6][4] = 2.91271111111111111111111111111e-3;
      A[6][5] = 7.89942857142857142857142857143e-4;

      A[7][1] = 5.6640625e-4;
      A[7][2] = 0.0;
      A[7][3] = 8.80973048941798941798941798942e-4;
      A[7][4] = -4.36921296296296296296296296296e-4;
      A[7][5] = 3.39006696428571428571428571429e-4;
      A[7][6] = -9.94646990740740740740740740741e-5;

      A[8][1] = 3.08333333333333333333333333333e-3;
      A[8][2] = 0.0;
      A[8][3] = 0.0;
      A[8][4] = 1.77777777777777777777777777778e-3;
      A[8][5] = 2.7e-3;
      A[8][6] = 1.57828282828282828282828282828e-3;
      A[8][7] = 1.08606060606060606060606060606e-2;

      A[9][1] = 3.65183937480112971375119150338e-3;
      A[9][2] = 0.0;
      A[9][3] = 3.96517171407234306617557289807e-3;
      A[9][4] = 3.19725826293062822350093426091e-3;
      A[9][5] = 8.22146730685543536968701883401e-3;
      A[9][6] = -1.31309269595723798362013884863e-3;
      A[9][7] = 9.77158696806486781562609494147e-3;
      A[9][8] = 3.75576906923283379487932641079e-3;

      A[10][1] = 3.70724106871850081019565530521e-3;
      A[10][2] = 0.0;
      A[10][3] = 5.08204585455528598076108163479e-3;
      A[10][4] = 1.17470800217541204473569104943e-3;
      A[10][5] = -2.11476299151269914996229766362e-2;
      A[10][6] = 6.01046369810788081222573525136e-2;
      A[10][7] = 2.01057347685061881846748708777e-2;
      A[10][8] = -2.83507501229335808430366774368e-2;
      A[10][9] = 1.48795689185819327555905582479e-2;

      A[11][1] = 3.51253765607334415311308293052e-2;
      A[11][2] = 0.0;
      A[11][3] = -8.61574919513847910340576078545e-3;
      A[11][4] = -5.79144805100791652167632252471e-3;
      A[11][5] = 1.94555482378261584239438810411e0;
      A[11][6] = -3.43512386745651359636787167574e0;
      A[11][7] = -1.09307011074752217583892572001e-1;
      A[11][8] = 2.3496383118995166394320161088e0;
      A[11][9] = -7.56009408687022978027190729778e-1;
      A[11][10] = 1.09528972221569264246502018618e-1;

      A[12][1] = 2.05277925374824966509720571672e-2;
      A[12][2] = 0.0;
      A[12][3] = -7.28644676448017991778247943149e-3;
      A[12][4] = -2.11535560796184024069259562549e-3;
      A[12][5] = 9.27580796872352224256768033235e-1;
      A[12][6] = -1.65228248442573667907302673325e0;
      A[12][7] = -2.10795630056865698191914366913e-2;
      A[12][8] = 1.20653643262078715447708832536e0;
      A[12][9] = -4.13714477001066141324662463645e-1;
      A[12][10] = 9.07987398280965375956795739516e-2;
      A[12][11] = 5.35555260053398504916870658215e-3;

      A[13][1] = -1.43240788755455150458921091632e-1;
      A[13][2] = 0.0;
      A[13][3] = 1.25287037730918172778464480231e-2;
      A[13][4] = 6.82601916396982712868112411737e-3;
      A[13][5] = -4.79955539557438726550216254291e0;
      A[13][6] = 5.69862504395194143379169794156e0;
      A[13][7] = 7.55343036952364522249444028716e-1;
      A[13][8] = -1.27554878582810837175400796542e-1;
      A[13][9] = -1.96059260511173843289133255423e0;
      A[13][10] = 9.18560905663526240976234285341e-1;
      A[13][11] = -2.38800855052844310534827013402e-1;
      A[13][12] = 1.59110813572342155138740170963e-1;

      A[14][1] = 8.04501920552048948697230778134e-1;
      A[14][2] = 0.0;
      A[14][3] = -1.66585270670112451778516268261e-2;
      A[14][4] = -2.1415834042629734811731437191e-2;
      A[14][5] = 1.68272359289624658702009353564e1;
      A[14][6] = -1.11728353571760979267882984241e1;
      A[14][7] = -3.37715929722632374148856475521e0;
      A[14][8] = -1.52433266553608456461817682939e1;
      A[14][9] = 1.71798357382154165620247684026e1;
      A[14][10] = -5.43771923982399464535413738556e0;
      A[14][11] = 1.38786716183646557551256778839e0;
      A[14][12] = -5.92582773265281165347677029181e-1;
      A[14][13] = 2.96038731712973527961592794552e-2;

      A[15][1] = -9.13296766697358082096250482648e-1;
      A[15][2] = 0.0;
      A[15][3] = 2.41127257578051783924489946102e-3;
      A[15][4] = 1.76581226938617419820698839226e-2;
      A[15][5] = -1.48516497797203838246128557088e1;
      A[15][6] = 2.15897086700457560030782161561e0;
      A[15][7] = 3.99791558311787990115282754337e0;
      A[15][8] = 2.84341518002322318984542514988e1;
      A[15][9] = -2.52593643549415984378843352235e1;
      A[15][10] = 7.7338785423622373655340014114e0;
      A[15][11] = -1.8913028948478674610382580129e0;
      A[15][12] = 1.00148450702247178036685959248e0;
      A[15][13] = 4.64119959910905190510518247052e-3;
      A[15][14] = 1.12187550221489570339750499063e-2;
      
      A[16][1] = -2.75196297205593938206065227039e-1;
      A[16][2] = 0.0;
      A[16][3] = 3.66118887791549201342293285553e-2;
      A[16][4] = 9.7895196882315626246509967162e-3;
      A[16][5] = -1.2293062345886210304214726509e1;
      A[16][6] = 1.42072264539379026942929665966e1;
      A[16][7] = 1.58664769067895368322481964272e0;
      A[16][8] = 2.45777353275959454390324346975e0;
      A[16][9] = -8.93519369440327190552259086374e0;
      A[16][10] = 4.37367273161340694839327077512e0;
      A[16][11] = -1.83471817654494916304344410264e0;
      A[16][12] = 1.15920852890614912078083198373e0;
      A[16][13] = -1.72902531653839221518003422953e-2;
      A[16][14] = 1.93259779044607666727649875324e-2;
      A[16][15] = 5.20444293755499311184926401526e-3;

      A[17][1] = 1.30763918474040575879994562983e0;
      A[17][2] = 0.0;
      A[17][3] = 1.73641091897458418670879991296e-2;
      A[17][4] = -1.8544456454265795024362115588e-2;
      A[17][5] = 1.48115220328677268968478356223e1;
      A[17][6] = 9.38317630848247090787922177126e0;
      A[17][7] = -5.2284261999445422541474024553e0;
      A[17][8] = -4.89512805258476508040093482743e1;
      A[17][9] = 3.82970960343379225625583875836e1;
      A[17][10] = -1.05873813369759797091619037505e1;
      A[17][11] = 2.43323043762262763585119618787e0;
      A[17][12] = -1.04534060425754442848652456513e0;
      A[17][13] = 7.17732095086725945198184857508e-2;
      A[17][14] = 2.16221097080827826905505320027e-3;
      A[17][15] = 7.00959575960251423699282781988e-3;
      A[17][16] = 0.0;

      BHAT[1] = 1.21278685171854149768890395495e-2;
      BHAT[2] = 0.0;
      BHAT[3] = 0.0;
      BHAT[4] = 0.0;
      BHAT[5] = 0.0;
      BHAT[6] = 0.0;
      BHAT[7] = 8.62974625156887444363792274411e-2;
      BHAT[8] = 2.52546958118714719432343449316e-1;
      BHAT[9] = -1.97418679932682303358307954886e-1;
      BHAT[10] = 2.03186919078972590809261561009e-1;
      BHAT[11] = -2.07758080777149166121933554691e-2;
      BHAT[12] = 1.09678048745020136250111237823e-1;
      BHAT[13] = 3.80651325264665057344878719105e-2;
      BHAT[14] = 1.16340688043242296440927709215e-2;
      BHAT[15] = 4.65802970402487868693615238455e-3;
      BHAT[16] = 0.0;
      BHAT[17] = 0.0;

      BPHAT[1] = 1.21278685171854149768890395495e-2;
      BPHAT[2] = 0.0;
      BPHAT[3] = 0.0;
      BPHAT[4] = 0.0;
      BPHAT[5] = 0.0;
      BPHAT[6] = 0.0;
      BPHAT[7] = 9.08394342270407836172412920433e-2;
      BPHAT[8] = 3.15683697648393399290429311645e-1;
      BPHAT[9] = -2.63224906576909737811077273181e-1;
      BPHAT[10] = 3.04780378618458886213892341513e-1;
      BPHAT[11] = -4.15516161554298332243867109382e-2;
      BPHAT[12] = 2.46775609676295306562750285101e-1;
      BPHAT[13] = 1.52260530105866022937951487642e-1;
      BPHAT[14] = 8.14384816302696075086493964505e-2;
      BPHAT[15] = 8.50257119389081128008018326881e-2;
      BPHAT[16] = -9.15518963007796287314100251351e-3;
      BPHAT[17] = 2.5e-2;

      B[1] = 1.70087019070069917527544646189e-2;
      B[2] = 0.0;
      B[3] = 0.0;
      B[4] = 0.0;
      B[5] = 0.0;
      B[6] = 0.0;
      B[7] = 7.22593359308314069488600038463e-2;
      B[8] = 3.72026177326753045388210502067e-1;
      B[9] = -4.01821145009303521439340233863e-1;
      B[10] = 3.35455068301351666696584034896e-1;
      B[11] = -1.31306501075331808430281840783e-1;
      B[12] = 1.89431906616048652722659836455e-1;
      B[13] = 2.68408020400290479053691655806e-2;
      B[14] = 1.63056656059179238935180933102e-2;
      B[15] = 3.79998835669659456166597387323e-3;
      B[16] = 0.0;
      B[17] = 0.0;

      BP[1] = 1.70087019070069917527544646189e-2;
      BP[2] = 0.0;
      BP[3] = 0.0;
      BP[4] = 0.0;
      BP[5] = 0.0;
      BP[6] = 0.0;
      BP[7] = 7.60624588745593757356421093119e-2;
      BP[8] = 4.65032721658441306735263127583e-1;
      BP[9] = -5.35761526679071361919120311817e-1;
      BP[10] = 5.03182602452027500044876052344e-1;
      BP[11] = -2.62613002150663616860563681567e-1;
      BP[12] = 4.26221789886109468625984632024e-1;
      BP[13] = 1.07363208160116191621476662322e-1;
      BP[14] = 1.14139659241425467254626653171e-1;
      BP[15] = 6.93633866500486770090602920091e-2;
      BP[16] = 2.0e-2;
      BP[17] = 0.0;
   }
};

#endif
