#ifndef __CONIC_HPP__
#define __CONIC_HPP__

#include <cmath>
#include <vector>

namespace vec_math {
    
    static std::vector<double> cross(const std::vector<double> &u,const std::vector<double> &v) {
        return { u[1]*v[2] - u[2]*v[1],   u[2]*v[0] - u[0]*v[2],   u[0]*v[1] - u[1]*v[0] };
    }

    static double dot(const std::vector<double> &u,const std::vector<double> &v) {
        return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
    }

    static double norm(const std::vector<double> &v) {
        return sqrt(dot(v,v));
    }

    static std::vector<double> mult(const std::vector<double> &v,double s) {
        return {v[0]*s, v[1]*s, v[2]*s};
    }

    static std::vector<double> mult(double s,const std::vector<double> &v) {
        return mult(v,s);
    }

    static std::vector<double> div(const std::vector<double> &v,double s) {
        return {v[0]/s, v[1]/s, v[2]/s};
    }

    static std::vector<double> add(const std::vector<double> &u,const std::vector<double> &v) {
        return {u[0]+v[0], u[1]+v[1], u[2]+v[2]};
    }

    static std::vector<double> sub(const std::vector<double> &u,const std::vector<double> &v) {
        return {u[0]-v[0], u[1]-v[1], u[2]-v[2]};
    }

    static std::vector<double> normalize(const std::vector<double> &v) {
        return div(v,norm(v));
    }

}

class Conic { public:

    double rp,e,i,O,w,M0,t0,mu;
    double n,p,muh;
    double sO,cO, si,ci, sw,cw;

    Conic() {}

    Conic(const std::vector<double> &elts) {
        setup_elements(elts);
    }

    Conic(double rp0,double e0,double i0,double O0,double w0,double M00,double t00,double mu0) {
        setup_elements(rp0,e0,i0,O0,w0,M00,t00,mu0);
    }

    Conic(const double &t,const double &mu0,const std::vector<double> &sv) {
        setup_state(t,mu0,sv);
    }

    Conic(const double &t,const double &mu0,const std::vector<double> &r,const std::vector<double> &v) {
        setup_rv(t,mu0,r,v);
    }

    void setup_elements(const std::vector<double> &elts) {
        setup_elements(elts[0],elts[1],elts[2],elts[3],elts[4],elts[5],elts[6],elts[7]);
    }

    void setup_elements(double rp0,double e0,double i0,double O0,double w0,double M00,double t00,double mu0) {
        rp=rp0; e=e0; i=i0; O=O0; w=w0; M0=M00; t0=t00; mu=mu0;
        if (e==1.) {
            n = sqrt(mu/(2.*rp*rp*rp));
            p = 2.*rp;
        } else {
            double alpha = (1.-e)/rp;
            n = sqrt(fabs(mu*alpha*alpha*alpha));
            p = (1.-e*e)/alpha;
        }
        double h = sqrt(mu*p);
        muh = mu/h;
        sincos(O,&sO,&cO); sincos(i,&si,&ci); sincos(w,&sw,&cw);
    }

    void setup_state(const double &t,const double &mu0,const std::vector<double> &sv) {
        setup_rv(t,mu0, {sv[0],sv[1],sv[2]}, {sv[3],sv[4],sv[5]} );
    }

    void setup_rv(const double &t,const double &mu0,const std::vector<double> &r,const std::vector<double> &v) {
        using namespace vec_math;
        t0 = t; mu = mu0;
        
        auto h = cross(r,v);
        double hm = norm(h);
        muh = mu/hm;
        p = (hm*hm)/mu;
        double rm = norm(r);
        double alpha = (2.*mu - dot(v,v)*rm)/(mu*rm);
        if (fabs(alpha) < 1e-15) {
            e = 1.0;
            rp = p/2.;
        } else {
            e = sqrt(1.-(p*alpha));
            rp = (1.-e)/alpha;
        }
        
        double ic = ((p/rm)-1)/e;
        if (ic> 1) ic= 1;
        if (ic<-1) ic=-1;
        double f = acos(ic);
        if (dot(r,v)<0) f = 2.*M_PI-f;

        if (e<1) {
            double E = 2.*atan(sqrt((1.-e)/(1.+e))*tan(f/2.));
            M0 = E - e*sin(E);
            n = sqrt(mu*alpha*alpha*alpha);
        } else if (e==1.) {
            double D = tan(f/2.);
            M0 = D + (D*D*D)/3.;
            n = sqrt(mu/(2.*rp*rp*rp));
        } else {
            double F = 2.*atanh(sqrt((e-1.)/(e+1.))*tan(f/2.));
            M0 = e*sinh(F) - F;
            n = sqrt(-mu*alpha*alpha*alpha);
        }
        M0 = fixrad(M0);

        i = acos(h[2]/hm);
        if (i<0.) i = M_PI-i;
        O = fixrad(atan2(h[0],-h[1]));
        double swf;
        if (i==0) swf = (r[1]*cos(O)-r[0]*sin(O))/rm;
        else swf = r[2]/(rm*sin(i));
        double cwf = (r[0]*cos(O)+r[1]*sin(O))/rm;
        w = fixrad(atan2(swf,cwf)-f);

        sincos(O,&sO,&cO); sincos(i,&si,&ci); sincos(w,&sw,&cw);
    }

    static double fixrad(double x) {
        const double tpi = 2.0*M_PI;
        return x-floor(x/tpi)*tpi;
    }

    std::vector<double> state(const double &t) const {
        double f = kepler( n*(t-t0) + M0 );
        double rm = p/(1. + e*cos(f));
        double st,ct; sincos(w+f,&st,&ct);
        std::vector<double> sv(6);
        sv[0] = rm*(cO*ct - sO*st*ci);
        sv[1] = rm*(sO*ct + cO*st*ci);
        sv[2] = rm*st*si;
        sv[3] = -muh*(cO*(st + e*sw) + sO*(ct + e*cw)*ci);
        sv[4] = -muh*(sO*(st + e*sw) - cO*(ct + e*cw)*ci);
        sv[5] =  muh*(ct + e*cw)*si;
        return sv;
    }

    double kepler(const double &M) const {
        if (e==1.) return kepler_parabola(M);
        else if (e<1) return kepler_ellipse(M);
        else return kepler_hyperbola(M);
    }

    double kepler_parabola(const double &M) const {
        double A = 1.5*fixrad(M);
        double B = pow( A + sqrt(A*A + 1.), 1./3.);
        return 2.*atan(B - (1./B));
    }

    double kepler_hyperbola(const double &M) const {
        double F = fixrad(M);
        for (unsigned nc(0);nc<100;nc++) {
            double d = -(e*sinh(F) - F - M)/(e*cosh(F) - F);
            F += d;
            if(fabs(d)<1e-15) break;
        }
        return 2.*atan(sqrt((e+1.)/(e-1.))*tanh(F/2.));
    }

    // Solve Kepler's Equation using an accelerated Newton's method
    // From "Fundamentals of Celestial Mechanics", Danby, 2nd ed., section 6.6
    double kepler_ellipse(const double &M) const {
        double Ms = fixrad(M);
        double x=Ms, es,ec, f,fp,dx;
        if (sin(Ms)>0) x += 0.85*e;
        else           x -= 0.85*e;
        for (unsigned nc(0);nc<10;nc++) {
            sincos(x,&es,&ec);
            es*=e; ec*=e; f = x-es-Ms;
            if (fabs(f)<1e-14) break;
            fp = 1.0-ec; dx = -f/fp;
            dx = -f/(fp + 0.5*dx*es);
            dx = -f/(fp + 0.5*dx*es + dx*dx*ec/6.0);
            x += dx;
        }
        return 2.*atan(sqrt((1.+e)/(1.-e))*tan(x/2.));
    }

    // Recursion FTW
    std::vector<double> elements(double t) const {
        if(t==t0) return {rp,e,i,O,w,M0,t0,mu};
        return Conic(t,mu,this->state(t)).elements(t);
    }

};

#endif
