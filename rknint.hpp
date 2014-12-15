#ifndef __RKNINT_HPP__
#define __RKNINT_HPP__

#include <cstdio>
#include <cmath>
#include <vector>

template<class O> class RKN_Integrator { 

    public:

        int iset,ifail,ieval;
        unsigned nsuccessful,nfailed,nattempted,nmaxsteps;
        double uround, Hnext,Hused,Hstart,Told,tolerance;
        bool high,oldhigh,onestep,start,justset,istart;
        std::vector<double> nrt1,nrt2,nrt3,nrt4;
        std::vector<double> threshold,thresholdp,weight,weightp;
        std::vector< std::vector<double> > jwork;

        RKN_Integrator(unsigned NEQ,double TOL=1e-12,bool HIGH=true);
        
        void setup(int NEQ,double H,double TOL,const std::vector<double> &THRES,const std::vector<double> &THRESP,
                int MAXSTP,bool ONESTP,bool HIGH);
        
        void integrate(O &obj,double TEND); 
        
        void evaluate64(const O &obj,unsigned NWANT,double TWANT,std::vector<double> &YWANT,std::vector<double> &YPWANT);

    private:

        double DIR, EXPON, H, RANGE, SCALE, STBRAD, TCHECK, TOL, TOOSML;
        int EFFCT1, EFFCT2, FINPH2, ISTP, NASTEP, NATT, neq, NFLSTP, NSTAGE, NSTMAX;
        bool FAILED, FIRST, HGIVEN, LAST, LOOP, OPTRED, PHASE2, SUCCES;

        double A[18][18], B[18], BHAT[18], BP[18], BPHAT[18], C[18];
        int PTR[17];

        void CFEVAL(double *B,double *BP,double S);
        void coefficients64();
        void coefficientsHi();

};

template<class O> RKN_Integrator<O>::
RKN_Integrator(unsigned NEQ,double TOL,bool HIGH) {
    istart = true;
    std::vector<double> THRES(NEQ),THRESP(NEQ);
    setup(NEQ,0,TOL,THRES,THRESP,0,false,HIGH);
    istart = true;
}

template<class O> void RKN_Integrator<O>::
setup(int NEQ,double H,double TOL,const std::vector<double> &THRES,const std::vector<double> &THRESP,
        int MAXSTP,bool ONESTP,bool HIGH) {
    
    //const int NSTMAX=100000;
    const double UROUND=3.0e-17,THRDEF=50.0*UROUND,U10=10.0*UROUND;
    
    start = istart;
    if(istart) {
        Hstart = fabs(H);
        Hused = 0.;
        Hnext = 0.;
        nsuccessful = 0;
        nfailed = 0;
        nattempted = 0;
        uround = UROUND;
        Told = 0.;
    } else {
        if(H != 0.) Hnext = fabs(H);
    }
   
    onestep = ONESTP; 
    
    /*if(MAXSTP<=0) nmaxsteps = NSTMAX;
    else nmaxsteps = MAXSTP;*/
    if(MAXSTP>0) nmaxsteps = MAXSTP;
    else nmaxsteps = 0;
   
    high = HIGH;
    
    if((TOL>1.) or (TOL<U10)) {
        iset = 3;
        istart = false;
        return;
    } else {
        tolerance = TOL;
    }
   
    neq = NEQ; 
    nrt1.resize(neq);
    nrt2.resize(neq);
    nrt3.resize(neq);
    nrt4.resize(neq);
    threshold.resize(neq);
    thresholdp.resize(neq);
    weight.resize(neq);
    weightp.resize(neq);
    jwork.resize(14);
    for(unsigned i=0;i<jwork.size();i++) jwork[i].resize(neq);
    
    if(THRES[0] <= 0.) {
        for(int K=0;K<neq;K++) {
            threshold[K] = THRDEF;
        }
    } else {
        for(int K=0;K<neq;K++) {
            if(THRES[K] <= 0.) {
                iset = 1;
                istart = false;
                return;
            } else {
                threshold[K] = THRES[K];
            }
        }
    }
   
    if(THRESP[0] <= 0.) {
        for(int K=0;K<neq;K++) {
            thresholdp[K] = THRDEF;
        }
    } else {
        for(int K=0;K<neq;K++) {
            if(THRESP[K] <= 0.) {
                iset = 1;
                istart = false;
                return;
            } else {
                thresholdp[K] = THRESP[K];
            }
        }
    }

    justset = true;
    iset = 0;
    istart = false;
    ifail = 0;
}

template<class O> void RKN_Integrator<O>::
integrate(O &obj,double TEND) {
    
    const double R=10.0,RMN1=1.0/R,RMN2=1.0/(R*R),RMN3=1.0/(R*R*R);
    
    double ABSH, ABSHP, ALPHA, BETA, DEL, ERR, ERRP, FDEL,
           HSQ, HUSED=0, SCALE1, SCL, SCLP, SK, SPK, TAU,
           THRES, THRESP, TOLD=0, TRY, WT, WTP, YDPDEL,
           YDPNRM, YINTK, YPDEL, YPINTK, YPNRM;
    int I, J, JSTAGE, K;
    bool DFLT, DFLTP, INEFFC=false, WASTE, save, checktime, phase2fail=false;

    // BLOCK A.
    // CHECK FOR NONZERO ENTRY VALUE OF IFAIL.

    if(ifail != 0) {
        printf(" RKNINT STOPPED DUE TO NONZERO INPUT VALUE OF IFAIL. ");
        printf("IFAIL = %3d\n",ifail);
        return;
    }
    NASTEP = 0;
    ifail = 0;
    save = true;
    checktime = true;

    if(obj.t == TEND) ifail = 1;
    if(not start) {
        if(unsigned(neq) != obj.neq) ifail = 1;
        if(DIR*(TEND-obj.t) < 0.) ifail = 1;
    }
    if(ifail!=0) {
        save = false;
    } else if(start) {

        // BLOCK B.
        // FIRST CALL IN AN INTEGRATION.
        // EXTRACT DATA FROM SETUP ROUTINE RKNSET.

        TOLD = obj.t;
        H = Hstart;
        HUSED = 0.;
        Hstart = 0.;
        NSTMAX = nmaxsteps;
        justset = false;
        start = false;
        TOL = tolerance;

        // INITIALIZATION.
    
        EFFCT1 = 0;
        EFFCT2 = 0;
        ISTP = 0;
        NFLSTP = 0;
        NATT = 0;
        FIRST = true;
        LOOP = false;
        RANGE = fabs(TEND-obj.t);
        if(TEND>obj.t) DIR = 1.0;
        else DIR = -1.0;

        // THE NEXT THREE LINES RELATE TO THE CODE CHOICE OF H0.

        OPTRED = false;
        HGIVEN = (fabs(H) > uround);
        PHASE2 = (FIRST and (not HGIVEN));

        // PICK UP CONSTANTS DEPENDING ON ORDER OF METHOD SELECTED AND SET
        // POINTERS. THE LOCAL ARRAY PTR CONTAINS POINTERS TO WHERE THE STAGES
        // ARE HELD IN THE ARRAY RWORK. THIS IS TO ALLOW SAVINGS IN STORAGE FOR
        // 12(10) PAIR AND BOTH PAIRS TO BE IMPLEMENTED IN THE SAME BLOCK OF
        // CODE.

        oldhigh = high;
        if(high) coefficientsHi();
        else     coefficients64();

        obj.dt(obj.t,obj.r,obj.a);

        if( not HGIVEN) {

            // PHASE 1 OF CALCULATION OF AN INITIAL STEPSIZE H0 BY THE CODE
            //   IF || YP(T0) || .NE. 0 THEN H0 = (TOL**EXPON) / ||YP(T0)||
            //                          ELSE H0 = TEND - T0
            //   IF || YDP(T0) || .NE. 0 THEN HP0 = (TOL**EXPON) / ||YDP(T0)||
            //                           ELSE HP0 = TEND - T0
            //   H = MIN(H0, HP0)
            // NOTE: 1) EXPON = 1/5  FOR THE 6(4) PAIR,
            //          EXPON = 1/11 FOR THE 12(10) PAIR.
            //       2) THE STRATEGY GIVEN ABOVE HAS BEEN MODIFIED TO TRY TO
            //          AVOID SELECTING TOO SMALL A VALUE FOR H0 WHEN THE
            //          INITIAL VALUES OF Y AND YP ARE SMALL.

            YPNRM = 0.;
            YDPNRM = 0.;
            DFLT = false;
            DFLTP = false;

            // EVALUATE NORMS AS DESCRIBED ABOVE AND INITIALISE THE WEIGHT VECTORS
            // FOR USE IN PHASE2.

            for(K=0;K<neq;K++) {
                THRES = threshold[K];
                WT = std::max(THRES,fabs(obj.r[K]));
                weight[K] = WT;
                TRY = fabs(obj.v[K]/WT);
                if(TRY > YPNRM) {
                    YPNRM = TRY;
                    DFLT = (WT == THRES);
                }
                THRESP = thresholdp[K];
                WTP = std::max(THRESP,fabs(obj.v[K]));
                weightp[K] = WTP;
                TRY = fabs(obj.a[K]/WTP);
                if(TRY > YDPNRM) {
                    YDPNRM = TRY;
                    DFLTP = (WTP == THRESP);
                }
            }
            ABSH = RANGE;
            ABSHP = RANGE;
            TAU = pow(TOL,EXPON);
            if (ABSH*YPNRM > TAU and not DFLT) ABSH = TAU/YPNRM;
            if (ABSHP*YDPNRM > TAU and not DFLTP) ABSHP = TAU/YDPNRM;
            H = std::min(ABSH,ABSHP);

            // END OF PHASE 1.

        }

        // END OF BLOCK B.

    } else {

        // BLOCK C.
        // CONTINUATION CALL.
        // CHECK DATA FOR CONSISTENCY.

        TOLD = Told;
        HUSED = Hused;

        // CHECK FOR OPTIONAL INPUTS.

        if(justset) {
            H = Hnext;
            NSTMAX = nmaxsteps;
            justset = false;
            TOL = tolerance;
            if(oldhigh != high) {
                oldhigh = high;
                if(high) coefficientsHi();
                else     coefficients64();
            }
        }

        // END OF BLOCK C.

    }

    // BLOCK D.
    // THE START OF A NEW STEP.

    if(ifail == 0) {
        H = DIR*fabs(H);
        SUCCES = false;
    }

    while(ifail == 0) {
    
        if(checktime) {
    
            LAST = false;

            // THE CODE IS REQUIRED TO HIT TEND EXACTLY WITH A REASONABLE STEPSIZE.
            // CHECK THAT  T+H .LE. TEND  OR  T+2*H .LE. TEND  IN THE DIRECTION
            // OF INTEGRATION. A CHECK IS INCORPORATED LATER IN THE CODE TO CHECK
            // WHEN H IS REDUCED TOO OFTEN BY A SIGNIFICANT AMOUNT, HENCE IMPAIRING
            // EFFICIENCY.

            if (DIR*(obj.t+H-TEND) >= 0.) {
                INEFFC = ((TEND-obj.t)/H < 0.5);
                H = TEND - obj.t;
                LAST = true;
            } else {
                if (DIR*(obj.t+2.*H-TEND) >= 0.) {
                    H = (TEND-obj.t)*0.5;
                }
            }

            FAILED = false;
        }

        // RESTART HERE AFTER A STEP FAILURE.

        // THE CHECK FOR LIMITING PRECISION, USING  TOOSML,  IS BASED ON
        // THE DISTANCE BETWEEN COEFFICIENTS  C(I)  IN THE RKN FORMULAS,
        // AS DESCRIBED IN
        //     SHAMPINE,L.F.,WATTS,H.A. (1979).  THE ART OF WRITING A
        //     RUNGE-KUTTA CODE.  APPL. MATH. COMPUTATION 5, 93-121.

        /*if(fabs(H)<TOOSML*std::max(fabs(obj.t),fabs(obj.t+H))) {

            // STEP TOO SMALL. CONTROL RETURNS TO USER THROUGH BOTTOM OF CODE.

            ifail = 3;
            break;
        }
        NASTEP = NASTEP + 1;
        if(NSTMAX > 0 and NASTEP > NSTMAX) {
            ifail = 2;
            break;
        }*/

        if (not obj.is_sane()) {
            ifail = 8;
            return;
        }

        // END OF BLOCK D.

        // BLOCK E.
        // FORM THE STAGES AND STORE THEM IN THE WORKING STORAGE ARRAY, RWORK.
        // TEMPORARY STORAGE:
        //      RWORK(NRT1+K) CONTAINS THE SUM OF THE STAGES  SUM A(I,J)*G(J)
        //          USED IN FORMING THE INTERNAL Y VALUES. AT THE END OF
        //          THE BLOCK IT WILL CONTAIN   SUM BHAT(I) * G(I)   AND
        //          WILL BE USED IN ERROR ESTIMATION FOR Y.
        //      RWORK(NRT2+K) WILL CONTAIN  SUM BPHAT(I) * G(I)   AND WILL
        //          BE USED IN ESTIMATING THE ERROR IN YP.
        // AT THE END OF THE BLOCK, THE NEW APPROXIMATION TO Y(K) WILL BE
        // IN RWORK(NRT3+K) AND THE APPROXIMATION TO YP(K) WILL BE IN
        // RWORK(NRT4+K).
        //
        // NOTE.
        // MANY OF THE FOLLOWING LOOPS OVER  K=1,NEQ  HAVE CONSTANT ARRAY VALUES
        // INSIDE. THE CODE IS WRITTEN WITH CLARITY IN MIND. ANY OPTIMIZING
        // COMPILER WILL IDENTIFY THESE OCCURRENCES AND MOVE THESE CONSTANT
        // VALUES OUTSIDE THE LOOPS. TWO LOOPS CONTAIN CHECKS FOR ZERO
        // MULTIPLIERS SO AS TO PREVENT NEEDLESS COMPUTATION.

        HSQ = H*H;
        for(K=0;K<neq;K++) {
            nrt2[K] = BPHAT[1]*obj.a[K];
        }

        for(I=2;I<=NSTAGE;I++) {
    
            for(K=0;K<neq;K++) {
                nrt1[K] = A[I][1]*obj.a[K];
            }

            for(J=2;J<=I-1;J++) {
                if (A[I][J] != 0.) {
                    for(K=0;K<neq;K++) {
                        nrt1[K] += A[I][J]*jwork[PTR[J-1]][K];
                    }
                }
            }

            for(K=0;K<neq;K++) {
                nrt3[K] = obj.r[K] + C[I]*H*obj.v[K] + HSQ*nrt1[K];
            }

            JSTAGE = PTR[I-1];
            obj.dt(obj.t+C[I]*H,nrt3,jwork[JSTAGE]);

            if(BPHAT[I]!=0.) {
                for(K=0;K<neq;K++) {
                    nrt2[K] += BPHAT[I]*jwork[JSTAGE][K];
                }
            }

            // CONTROL USUALLY GOES FROM HERE TO THE START OF THE LOOP (I=2,NSTAGE),
            // UNLESS IT IS THE FIRST STEP AND THE CODE IS COMPUTING H0,
            // IN WHICH CASE THE NEXT BLOCK OF CODE IS EXECUTED.

            phase2fail = false;
            if(PHASE2) {
                YPNRM = 0.;
                YDPNRM = 0.;
                DEL = 0.;
                FDEL = 0.;
                YPDEL = 0.;
                YDPDEL = 0.;
            
                for(K=0;K<neq;K++) {
                    YINTK = nrt3[K];
                    YPINTK = obj.v[K] + H*nrt2[K];
                    WT = std::max(weight[K],fabs(YINTK));
                    WTP = std::max(weightp[K],fabs(YPINTK));
                    weight[K] = WT;
                    weightp[K] = WTP;
                    YPNRM = std::max(YPNRM,std::max(fabs(YINTK)/WT,fabs(obj.r[K])/WT));
                    YDPNRM = std::max(YDPNRM,std::max(fabs(YPINTK)/WTP,fabs(obj.v[K])/WTP));
                    DEL = std::max(DEL,fabs(YINTK-obj.r[K])/WT);
                    FDEL = std::max(FDEL,fabs(YPINTK-obj.v[K])/WT);
                    YPDEL = std::max(YPDEL,fabs(YPINTK-obj.v[K])/WTP);
                    YDPDEL = std::max(YDPDEL,fabs(jwork[JSTAGE][K]-obj.a[K])/WTP);
                }
                DEL = std::max(DEL,fabs(C[I]*H)/RANGE);
                SCL = 1.;
                SCLP = 1.;
                if(DEL > 10.*uround*YPNRM) {
                    if(fabs(H)*FDEL > STBRAD*DEL) {
                        SCL = STBRAD/R*std::max(DEL/(fabs(H)*FDEL),RMN3);
                    }
                }
                if(YPDEL > 10.*uround*YDPNRM) {
                    if(fabs(H)*YDPDEL > STBRAD*YPDEL) {
                        SCLP = STBRAD/R*std::max(YPDEL/(fabs(H)*YDPDEL),RMN3);
                    }
                }
                SCALE1 = std::min(SCL,SCLP);
                if (SCALE1 < 1.) {
                    NATT = NATT + 1;
                    H = SCALE1*fabs(H);
                    H = DIR*H;
                    LAST = false;
            
                    // RESET THE WEIGHT VECTORS

                    for(K=0;K<neq;K++) {
                        weight[K] = std::max(fabs(obj.r[K]),threshold[K]);
                        weightp[K] = std::max(fabs(obj.v[K]),thresholdp[K]);
                    }
                    checktime = false;
                    phase2fail = true;
                    break;
                }

                PHASE2 = (I<FINPH2);
            }
        }
        if(phase2fail) continue;

        for(K=0;K<neq;K++) {
            nrt4[K] = obj.v[K] + H*nrt2[K];
        }

        // THE 12(10) PAIR DOES NOT USE FSAL (FIRST STAGE ON STEP N+1 IS SAME
        // AS LAST STAGE ON STEP N), HENCE MUST ADD STAGES UP TO OBTAIN THE
        // APPROXIMATIONS. UTILIZE FACT THAT
        //     BHAT(2)=BHAT(3)=BHAT(4)=BHAT(5)=BHAT(6)=BHAT(16)=BHAT(17)=0.

        if(oldhigh) {
            for(K=0;K<neq;K++) {
                nrt1[K] = BHAT[1]*obj.a[K];
            }
            for(J=7;J<=15;J++) {
                for(K=0;K<neq;K++) {
                    nrt1[K] += BHAT[J]*jwork[PTR[J-1]][K];
                }
            }
            for(K=0;K<neq;K++) {
                nrt3[K] = obj.r[K] + H*obj.v[K] + HSQ*nrt1[K];
            }
        }

        // END OF BLOCK E.

        // BLOCK F.
        // ERROR ESTIMATION.
        // USE THE TEMPORARY STORAGE VECTORS RWORK(NRT1+K) AND RWORK(NRT2+K)
        // TO COMPUTE  YHAT - Y  AND  YPHAT - YP  RESPECTIVELY.
        // SEPARATE SECTIONS ARE USED FOR THE 12(10) AND 6(4) PAIRS
        // IN ORDER TO TAKE ADVANTAGE OF THE FACT THAT FOR THE 6(4) PAIR
        // CASE, B(4) = B(5) = B(6) = ZERO, WHILST FOR THE 12(10) PAIR
        // B(2) = B(3) = B(4) = B(5) = B(6) = B(17) = ZERO AND SIMILARLY FOR
        // BP.

        for(K=0;K<neq;K++) {
            nrt1[K] -= B[1]*obj.a[K];
            nrt2[K] -= BP[1]*obj.a[K];
        }

        if(not oldhigh) {
        
            // ERROR ESTIMATES FOR THE 6(4) PAIR.

            for(I=2;I<=3;I++) {
                for(K=0;K<neq;K++) {
                    nrt1[K] -= B[I]*jwork[PTR[I-1]][K];
                }
            }
            for(I=2;I<=NSTAGE;I++) {
                for(K=0;K<neq;K++) {
                    nrt2[K] -= BP[I]*jwork[PTR[I-1]][K];
                }
            }
        } else {

            // ERROR ESTIMATES FOR THE 12(10) PAIR.

            for(I=7;I<=16;I++) {
                JSTAGE = PTR[I-1];
                for(K=0;K<neq;K++) {
                    nrt1[K] -= B[I]*jwork[JSTAGE][K];
                    nrt2[K] -= BP[I]*jwork[JSTAGE][K];
                }
            }
        }

        // FORM ERR, THE MAXIMUM OF THE WEIGHTED MAX NORMS OF THE ESTIMATED
        // LOCAL ERRORS IN Y AND YP.

        ERR = 0.;
        ERRP = 0.;
        for(K=0;K<neq;K++) {
            SK = HSQ*nrt1[K];
            SPK = H*nrt2[K];
            WT = 0.5*fabs(obj.r[K]) + 0.25*fabs(nrt3[K]) + 0.25*fabs(nrt3[K]-SK);
            WT = std::max(WT,threshold[K]);
            WTP = 0.5*fabs(obj.v[K]) + 0.25*fabs(nrt4[K]) + 0.25*fabs(nrt4[K]-SPK);
            WTP = std::max(WTP,thresholdp[K]);
            ERR = std::max(ERR,fabs(SK)/WT);
            ERRP = std::max(ERRP,fabs(SPK)/WTP);
        }
        ERR = std::max(ERR,ERRP);

        // END OF BLOCK F.

        if(ERR > TOL) {

            // BLOCK G.
            // FAILED STEP.

            LAST = false;
            if (FIRST and (not HGIVEN)) {

                // PHASE 3 FOR CODE CHOICE OF H0.

                NATT = NATT + 1;
                if (SUCCES) {

                    // THE CODE HAS DISCARDED AN INITIAL STEP IN FAVOUR OF A MUCH LARGER
                    // PREDICTED STEP BUT THIS LARGER STEP HAS FAILED.  THEREFORE CARRY
                    // OUT OPTIMAL REDUCTION.

                    OPTRED = true;
                    ALPHA = SCALE*pow((TOL/ERR),EXPON);
                    ALPHA = std::max(RMN2,ALPHA);

                } else {

                    // NO SUCCESSFUL STEP YET.  REDUCE H0 TO H0/R AND RESET WEIGHT VECTORS.

                    ALPHA = RMN1;
                    PHASE2 = true;
                    for(K=0;K<neq;K++) {
                        weight[K] = std::max(fabs(obj.r[K]),threshold[K]);
                        weightp[K] = std::max(fabs(obj.v[K]),thresholdp[K]);
                    }
                }
            } else {

                // NOT THE FIRST STEP (OR H0 WAS SUPPLIED BY THE USER), SO USE THE
                // NORMAL STEP REDUCTION ALGORITHM.

                NFLSTP = NFLSTP + 1;
                if (FAILED) {
                    ALPHA = 0.5;
                } else {
                    FAILED = true;
                    ALPHA = SCALE*pow((TOL/ERR),EXPON);
                    ALPHA = std::max(RMN1,ALPHA);
                }
            }

            H = ALPHA*H;

            // TRY THE STEP AGAIN.

            checktime = false;
            continue;

            // END OF BLOCK G.

        } else {

            // BLOCK H.
            // SUCCESSFUL STEP.

            SUCCES = true;
            BETA = pow((ERR/TOL),EXPON)/SCALE;
            if (FIRST) {

                // PHASE 3 FOR CODE CHOICE OF H0.

                ALPHA = 1./std::max(RMN3,BETA);
            } else {

                // NORMAL STEPSIZE PREDICTION.

                ALPHA = 1./std::max(RMN1,BETA);
                if (FAILED) ALPHA = std::min(1.,ALPHA);
            }
            HUSED = H;
            H = ALPHA*H;

            if(FIRST and (not LAST) and (not OPTRED) and (ALPHA > R)) {

                // PHASE 3 FOR CODE CHOICE OF H0. THE CODE HAS ATTEMPTED THE FIRST STEP,
                // IS NOT AT THE END OF THE INTEGRATION RANGE, HAS NOT ENCOUNTERED A STEP
                // FAILURE IN PHASE3, AND THE PREDICTED INCREASE IS LARGER THAN THE
                // MAXIMUM PERMITTED ON A NORMAL STEP. RETAKE THE STEP.

                HUSED = 0.;
                NATT = NATT + 1;
                checktime = true;
                continue;
            }

            if(FIRST) Hstart = HUSED;
            FIRST = false;
            TOLD = obj.t;
            obj.t += HUSED;
            EFFCT1 += 1;
            if (LAST) {
                obj.t = TEND;
                if (INEFFC) EFFCT2 += 1;

                // TEST FOR INEFFICIENT USE OF FORMULAS FOR OUTPUT.

                if (EFFCT1 > 100) {
                    WASTE = (double(EFFCT2)/double(EFFCT1) > 0.1);
                    if (WASTE) {
                        ifail = 4;
                        EFFCT1 = 0;
                        EFFCT2 = 0;
                    }
                }
            }
            ISTP = ISTP + 1;

            // ON A SUCCESSFUL STEP, COPY THE NEW APPROXIMATIONS INTO Y AND YP.
            // SAVE THE OLD VALUES OF Y, YP AND YDP IN RWORK,
            // ALONG WITH THE INTERMEDIATE STAGES FOR USE IN CASE DENSE
            // OUTPUT IS REQUIRED.

            nrt1 = obj.r;
            obj.r = nrt3;
            nrt2 = obj.v;
            obj.v = nrt4;
            nrt3 = obj.a;

            // SINCE THE 12(10) PAIR DOES NOT USE FSAL (SEE ABOVE), A FUNCTION
            // EVALUATION IS DONE HERE WHICH WILL BE USED ON THE NEXT STEP. FOR THE
            // 6(4) PAIR SET YDP FROM THE LAST STAGE.

            if(oldhigh) {
                //F(T,Y,YDP);
                obj.dt(obj.t,obj.r,obj.a);
            } else {
                obj.a = jwork[PTR[5]];
            }

            // IF THE CODE HAS NOT REACHED THE END OF THE INTEGRATION RANGE OR IS
            // OPERATING IN INTERVAL MODE, ATTEMPT THE NEXT STEP.

            if(LAST or onestep) break;
        
            SUCCES = false;
            checktime = true;

            // END OF BLOCK H.

        }

    }

    // BLOCK I.
    // SET DIAGNOSTIC INFORMATION AND EXIT.

    if(save) {
        Hused = HUSED;
        Hnext = H;
        nsuccessful = ISTP;
        nfailed = NFLSTP;
        nattempted = NATT;
        Told = TOLD;
    }

    // TEST FOR SUCCESSIVE FAILURES AT THE SAME VALUE OF T.

    if(ifail > 0) {
        if(not LOOP) {
            LOOP = true;
            TCHECK = obj.t;
        } else {
            if(TCHECK==obj.t) {
                printf("\n RKNINT STOPPED DUE TO TWO SUCCESSIVE FAILURES AT THE SAME VALUE OF T (= %12.5E)\n\n",obj.t);
                ifail += 10;
                return;
            } else {
               LOOP = false;
            }
        }
    } else {
        LOOP = false;
    }

    // END OF BLOCK I.
}

template<class O> void RKN_Integrator<O>::
CFEVAL(double *B,double *BP,double S) {
    double U, V, W, X;
    const double R = sqrt(8581.0);
    const double Y = 22529243880.0;
    const double Z = 11264621940.0;

    B[1] = (((900.e0*S-3819.e0)*S+6386.e0)*S-5244.e0)*S + 2106.e0;
    B[1] = B[1]/4212.e0;
    BP[1] = (((5400.e0*S-19095.e0)*S+25544.e0)*S-15732.e0)*S + 4212.e0;
    BP[1] = BP[1]/4212.e0;
    B[2] = 0.0;
    BP[2] = 0.e0;

    U = 5860823.0 - 152228.0*R;
    V = 4929647204.0 - 156109769.0*R;
    W = 22190560391.0 - 1109665151.0*R;
    X = 81356461.0 + 25954829.0*R;
    B[3] = (((1800.*U*S-6.*V)*S+W)*S+18.*X)*S/Y;
    BP[3] = (((5400.*U*S-15.*V)*S+2.*W)*S+27.*X)*S/Z;

    U = 5860823.0 + 152228.0*R;
    V = 4929647204.0 + 156109769.0*R;
    W = 22190560391.0 + 1109665151.0*R;
    X = 81356461.0 - 25954829.0*R;
    B[4] = (((1800.*U*S-6.*V)*S+W)*S+18.*X)*S/Y;
    BP[4] = (((5400.*U*S-15.*V)*S+2.*W)*S+27.*X)*S/Z;

    B[5] = ((225.e0*S-651.e0)*S+620.e0)*S - 195.e0;
    B[5] = -200.e0*S*B[5]/17901.e0;
    BP[5] = ((270.0*S-651.0)*S+496.0)*S - 117.0;
    BP[5] = -1000.0*S*BP[5]/17901.0;
    B[6] = S*(S-1.0)*((300.0*S-523.0)*S+234.0)/220.0;
    BP[6] = S*(((1800.0*S-4115.0)*S+3028.0)*S-702.0)/220.0;
}

template<class O> void RKN_Integrator<O>::
evaluate64(const O &obj,unsigned NWANT,double TWANT,std::vector<double> &YWANT,std::vector<double> &YPWANT) {
    
    double H, HSIG, SIGMA, SUM, SUMP, TOLD;
    unsigned I, JSTAGE, K;
    bool INTERP;
    double BPSTAR[7], BSTAR[7];
    
    if(high) {
        ieval = 3;
        return;
    }
    if (NWANT>unsigned(neq) or NWANT<1) {
        ieval = 2;
        return;
    }

    // RWORK(NRT1+...) CONTAINS THE PREVIOUS VALUE OF Y.
    // RWORK(NRT2+...) CONTAINS THE PREVIOUS VALUE OF YP.
    // RWORK(NRT3+...) CONTAINS THE PREVIOUS VALUE OF YDP.

    // CHECK FOR TRIVIAL CASES

    ieval = 0;
    if (TWANT==obj.t) {
        for(K=0;K<NWANT;K++) {
            YWANT[K] = obj.r[K];
            YPWANT[K] = obj.v[K];
        }
        return;
    }
    H = Hused;
    TOLD = Told;
    if (TWANT==TOLD) {
        for(K=0;K<NWANT;K++) {
            YWANT[K] = nrt1[K];
            YPWANT[K] = nrt2[K];
        }
        return;
    }

    // CHECK FOR EXTRAPOLATION

    if (H>0) {
        INTERP = ((TOLD<TWANT) and (TWANT<obj.t));
    } else {
        INTERP = ((TOLD>TWANT) and (TWANT>obj.t));
    }
    if( not INTERP) ieval = 1;

    HSIG = TWANT - TOLD;
    SIGMA = HSIG/H;
    CFEVAL(BSTAR,BPSTAR,SIGMA);

    // COMPUTE YWANT AND YPWANT.
    // THE FACT THAT BSTAR(2) = 0 = BPSTAR(2) IS USED.

    for(K=0;K<NWANT;K++) {
        SUM = BSTAR[1]*nrt3[K];
        SUMP = BPSTAR[1]*nrt3[K];
        for(I=3;I<=6;I++) {
            JSTAGE = I-2;
            SUM += BSTAR[I]*jwork[JSTAGE][K];
            SUMP += BPSTAR[I]*jwork[JSTAGE][K];
        }
        YWANT[K] = nrt1[K] + HSIG*(nrt2[K]+HSIG*SUM);
        YPWANT[K] = nrt2[K] + HSIG*SUMP;
    }
}

template<class O> void RKN_Integrator<O>::
coefficients64() {

    EXPON = 0.2;
    NSTAGE = 6;
    STBRAD = 4.0;
    FINPH2 = 4;
    TOOSML = 20.*uround;
    SCALE = 0.8;
    
    PTR[1] = 0;
    PTR[2] = 1;
    PTR[3] = 2;
    PTR[4] = 3;
    PTR[5] = 4;

    C[1] = 0.0e0;
    C[2] = 1.29295903136704415288990053209e-1;
    C[3] = 2.58591806273408830577980106418e-1;
    C[4] = 6.70297082615480058310908782471e-1;
    C[5] = 9.0e-1;
    C[6] = 1.0e0;

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
    A[6][2] = 0.0e0;
    A[6][3] = 2.882917411897667776841772093e-1;
    A[6][4] = 1.22425537174570410182422422005e-1;
    A[6][5] = 1.1172560192168035305290207251e-2;

    BPHAT[1] = 7.81101614434947768281101614435e-2;
    BPHAT[2] = 0.0e0;
    BPHAT[3] = 3.88843478705982602715952208523e-1;
    BPHAT[4] = 3.71320757928842267403035557523e-1;
    BPHAT[5] = 1.1172560192168035305290207251e-1;
    BPHAT[6] = 5.0e-2;

    B[1] = 1.05885926037041827822001005886e0;
    B[2] = -2.40675137192445205319176777632e0;
    B[3] = 1.84789211155403377497175771746e0;
    B[4] = 0.0e0;
    B[5] = 0.0e0;
    B[6] = 0.0e0;

    BP[1] = 5.46058879392212725546058879392e-2;
    BP[2] = 0.0e0;
    BP[3] = 4.6126678590362684429468353488e-1;
    BP[4] = 1.95880859479312656291875875208e-1;
    BP[5] = 3.88246466677839226858834701972e-1;
    BP[6] = -1.0e-1;

}

template<class O> void RKN_Integrator<O>::
coefficientsHi() {

    EXPON = 1.0/11.0;
    NSTAGE = 17;
    STBRAD = 8.0;
    FINPH2 = 14;
    TOOSML = 200.*uround;
    SCALE = 0.75;

    PTR[1] = 0;
    PTR[2] = 1;
    PTR[3] = 2;
    PTR[4] = 0;
    PTR[5] = 3;
    PTR[6] = 4;
    PTR[7] = 5;
    PTR[8] = 6;
    PTR[9] = 7;
    PTR[10] = 8;
    PTR[11] = 9;
    PTR[12] = 10;
    PTR[13] = 11;
    PTR[14] = 12;
    PTR[15] = 13;
    PTR[16] = 3;

    C[1] = 0.0e0;
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
    C[16] = 1.0e0;
    C[17] = 1.0e0;

    A[2][1] = 2.0e-4;

    A[3][1] = 2.66666666666666666666666666667e-4;
    A[3][2] = 5.33333333333333333333333333333e-4;

    A[4][1] = 2.91666666666666666666666666667e-3;
    A[4][2] = -4.16666666666666666666666666667e-3;
    A[4][3] = 6.25e-3;

    A[5][1] = 1.64609053497942386831275720165e-3;
    A[5][2] = 0.0e0;
    A[5][3] = 5.48696844993141289437585733882e-3;
    A[5][4] = 1.75582990397805212620027434842e-3;

    A[6][1] = 1.9456e-3;
    A[6][2] = 0.0e0;
    A[6][3] = 7.15174603174603174603174603175e-3;
    A[6][4] = 2.91271111111111111111111111111e-3;
    A[6][5] = 7.89942857142857142857142857143e-4;

    A[7][1] = 5.6640625e-4;
    A[7][2] = 0.0e0;
    A[7][3] = 8.80973048941798941798941798942e-4;
    A[7][4] = -4.36921296296296296296296296296e-4;
    A[7][5] = 3.39006696428571428571428571429e-4;
    A[7][6] = -9.94646990740740740740740740741e-5;

    A[8][1] = 3.08333333333333333333333333333e-3;
    A[8][2] = 0.0e0;
    A[8][3] = 0.0e0;
    A[8][4] = 1.77777777777777777777777777778e-3;
    A[8][5] = 2.7e-3;
    A[8][6] = 1.57828282828282828282828282828e-3;
    A[8][7] = 1.08606060606060606060606060606e-2;

    A[9][1] = 3.65183937480112971375119150338e-3;
    A[9][2] = 0.0e0;
    A[9][3] = 3.96517171407234306617557289807e-3;
    A[9][4] = 3.19725826293062822350093426091e-3;
    A[9][5] = 8.22146730685543536968701883401e-3;
    A[9][6] = -1.31309269595723798362013884863e-3;
    A[9][7] = 9.77158696806486781562609494147e-3;
    A[9][8] = 3.75576906923283379487932641079e-3;

    A[10][1] = 3.70724106871850081019565530521e-3;
    A[10][2] = 0.0e0;
    A[10][3] = 5.08204585455528598076108163479e-3;
    A[10][4] = 1.17470800217541204473569104943e-3;
    A[10][5] = -2.11476299151269914996229766362e-2;
    A[10][6] = 6.01046369810788081222573525136e-2;
    A[10][7] = 2.01057347685061881846748708777e-2;
    A[10][8] = -2.83507501229335808430366774368e-2;
    A[10][9] = 1.48795689185819327555905582479e-2;

    A[11][1] = 3.51253765607334415311308293052e-2;
    A[11][2] = 0.0e0;
    A[11][3] = -8.61574919513847910340576078545e-3;
    A[11][4] = -5.79144805100791652167632252471e-3;
    A[11][5] = 1.94555482378261584239438810411e0;
    A[11][6] = -3.43512386745651359636787167574e0;
    A[11][7] = -1.09307011074752217583892572001e-1;
    A[11][8] = 2.3496383118995166394320161088e0;
    A[11][9] = -7.56009408687022978027190729778e-1;
    A[11][10] = 1.09528972221569264246502018618e-1;

    A[12][1] = 2.05277925374824966509720571672e-2;
    A[12][2] = 0.0e0;
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
    A[13][2] = 0.0e0;
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
    A[14][2] = 0.0e0;
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
    A[15][2] = 0.0e0;
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
    A[16][2] = 0.0e0;
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
    A[17][2] = 0.0e0;
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
    A[17][16] = 0.0e0;


    BHAT[1] = 1.21278685171854149768890395495e-2;
    BHAT[2] = 0.0e0;
    BHAT[3] = 0.0e0;
    BHAT[4] = 0.0e0;
    BHAT[5] = 0.0e0;
    BHAT[6] = 0.0e0;
    BHAT[7] = 8.62974625156887444363792274411e-2;
    BHAT[8] = 2.52546958118714719432343449316e-1;
    BHAT[9] = -1.97418679932682303358307954886e-1;
    BHAT[10] = 2.03186919078972590809261561009e-1;
    BHAT[11] = -2.07758080777149166121933554691e-2;
    BHAT[12] = 1.09678048745020136250111237823e-1;
    BHAT[13] = 3.80651325264665057344878719105e-2;
    BHAT[14] = 1.16340688043242296440927709215e-2;
    BHAT[15] = 4.65802970402487868693615238455e-3;
    BHAT[16] = 0.0e0;
    BHAT[17] = 0.0e0;

    BPHAT[1] = 1.21278685171854149768890395495e-2;
    BPHAT[2] = 0.0e0;
    BPHAT[3] = 0.0e0;
    BPHAT[4] = 0.0e0;
    BPHAT[5] = 0.0e0;
    BPHAT[6] = 0.0e0;
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
    B[2] = 0.0e0;
    B[3] = 0.0e0;
    B[4] = 0.0e0;
    B[5] = 0.0e0;
    B[6] = 0.0e0;
    B[7] = 7.22593359308314069488600038463e-2;
    B[8] = 3.72026177326753045388210502067e-1;
    B[9] = -4.01821145009303521439340233863e-1;
    B[10] = 3.35455068301351666696584034896e-1;
    B[11] = -1.31306501075331808430281840783e-1;
    B[12] = 1.89431906616048652722659836455e-1;
    B[13] = 2.68408020400290479053691655806e-2;
    B[14] = 1.63056656059179238935180933102e-2;
    B[15] = 3.79998835669659456166597387323e-3;
    B[16] = 0.0e0;
    B[17] = 0.0e0;

    BP[1] = 1.70087019070069917527544646189e-2;
    BP[2] = 0.0e0;
    BP[3] = 0.0e0;
    BP[4] = 0.0e0;
    BP[5] = 0.0e0;
    BP[6] = 0.0e0;
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
    BP[17] = 0.0e0;
}

#endif
