#ifndef __NBODY_HPP__
#define __NBODY_HPP__

#include "rknint.hpp"
#include <iostream>

class NBody { public:

    unsigned neq, nobj;
    double t, maxdist;
    std::vector<double> r,v,a, u,R, emax_e;
    std::vector<std::string> name;
    std::vector<int> emax_i,emax_j;

    bool oblate,verbose,echeck;
    unsigned iob;
    double J2,J4, ra,dec,et_pole;
    std::vector<double> Cx,Cy,Cz;

    NBody() {
        nobj = neq = 0;
        t = maxdist = 0.0;
        verbose = oblate = echeck = false;
    }

    unsigned add_object(const std::string &name0,double u0,double R0,const std::vector<double> &r0,const std::vector<double> &v0) {
        name.push_back(name0);
        u.push_back(u0); R.push_back(R0);
        for(int i=0;i<3;i++) {
            r.push_back(r0[i]);
            v.push_back(v0[i]);
        }
        nobj++; neq+=3;
        return nobj-1;
    }

    unsigned add_object_state(const std::string &name0,double u0,double R0,const std::vector<double> &s0) {
        name.push_back(name0);
        u.push_back(u0); R.push_back(R0);
        for(int i=0;i<3;i++) {
            r.push_back(s0[i  ]);
            v.push_back(s0[i+3]);
        }
        nobj++; neq+=3;
        return nobj-1;
    }

    void set_oblate(const std::string &name0,double J20,double J40,double ra0,double dec0) {
        oblate = true;
        iob = lookup(name0);
        J2 = J20; J4 = J40; ra = ra0; dec = dec0;
        const double d2r = M_PI/180.;
        double sr,cr; sincos(ra *d2r,&sr,&cr);
        double sd,cd; sincos(dec*d2r,&sd,&cd);
        Cx.resize(3); Cy.resize(3); Cz.resize(3);
        Cx[0]=-sr   ; Cx[1]= cr   ; Cx[2]=0.;
        Cy[0]=-cr*sd; Cy[1]=-sr*sd; Cy[2]=cd;
        Cz[0]= cr*cd; Cz[1]= sr*cd; Cz[2]=sd;
    }

    void add_emax(unsigned i,double e,int j=-1) {
        emax_i.push_back(i); 
        emax_e.push_back(e); 
        emax_j.push_back(j);
        echeck = true;
    }

    int lookup(const std::string &name0) const {
        for(unsigned i=0;i<nobj;i++) {
            if(name[i].compare(name0)==0)
                return i; }
        return -1;
    }

    void dt(double tt,const std::vector<double> &rr,std::vector<double> &aa) const {
        aa = std::vector<double>(neq,0.0);
        particle_grav(tt,rr,aa);
        if(oblate) oblate_grav(tt,rr,aa);
    }

    void particle_grav(double tt,const std::vector<double> &rr,std::vector<double> &aa) const {
        for(unsigned i=0;i<nobj-1;i++) {
            for(unsigned j=i+1;j<nobj;j++) {
                if(u[i]+u[j]==0.) continue;
                const unsigned i2=i*3, j2=j*3;
                const double dx = rr[i2  ]-rr[j2  ];
                const double dy = rr[i2+1]-rr[j2+1];
                const double dz = rr[i2+2]-rr[j2+2];
                const double dr2 = dx*dx + dy*dy + dz*dz;
                const double dr3 = dr2*sqrt(dr2);
                if (u[i]>0.) {
                    aa[j2  ] += (u[i]*dx)/dr3;
                    aa[j2+1] += (u[i]*dy)/dr3;
                    aa[j2+2] += (u[i]*dz)/dr3;
                }
                if (u[j]>0.) {
                    aa[i2  ] -= (u[j]*dx)/dr3;
                    aa[i2+1] -= (u[j]*dy)/dr3;
                    aa[i2+2] -= (u[j]*dz)/dr3;
                }
            }
        }
    }

    void oblate_grav(double tt,const std::vector<double> &rr,std::vector<double> &aa) const {
        const static double k2 = -1.5*J2*R[iob]*R[iob];
        const static double k4 = -0.625*J4*R[iob]*R[iob]*R[iob]*R[iob];
        for(unsigned i=0;i<nobj;i++) {
            if(i==iob) continue;
            const double dx = rr[3*i  ]-rr[3*iob  ];
            const double dy = rr[3*i+1]-rr[3*iob+1];
            const double dz = rr[3*i+2]-rr[3*iob+2];
            const double x = Cx[0]*dx + Cx[1]*dy + Cx[2]*dz;
            const double y = Cy[0]*dx + Cy[1]*dy + Cy[2]*dz;
            const double z = Cz[0]*dx + Cz[1]*dy + Cz[2]*dz;
            const double r2 = x*x + y*y + z*z;
            const double r3 = r2*sqrt(r2), r4 = r2*r2;
            const double z2r2 = (z*z)/r2, z4r4 = z2r2*z2r2;
            const double J2tx = k2*(5.*z2r2 - 1.)/r2;
            const double J2tz = k2*(5.*z2r2 - 3.)/r2;
            const double J4tx = k4*(63.*z4r4 - 42.*z2r2 +  3.)/r4;
            const double J4tz = k4*(63.*z4r4 - 70.*z2r2 + 15.)/r4;
            const double GEx = (J2tx+J4tx)*(x/r3);
            const double GEy = (J2tx+J4tx)*(y/r3);
            const double GEz = (J2tz+J4tz)*(z/r3);
            const double Gx = Cx[0]*GEx + Cy[0]*GEy + Cz[0]*GEz;
            const double Gy = Cx[1]*GEx + Cy[1]*GEy + Cz[1]*GEz;
            const double Gz = Cx[2]*GEx + Cy[2]*GEy + Cz[2]*GEz;
            aa[3*i  ] -= u[iob]*Gx;
            aa[3*i+1] -= u[iob]*Gy;
            aa[3*i+2] -= u[iob]*Gz;
            if(u[i] > 0) {
                aa[3*iob  ] += u[i]*Gx;
                aa[3*iob+1] += u[i]*Gy;
                aa[3*iob+2] += u[i]*Gz;
            }
        }
    }

    std::vector<double> get_position(unsigned i) const {
        std::vector<double> ret(3,0.0);
        if(i<nobj) {
            ret[0]=r[3*i]; ret[1]=r[3*i+1]; ret[2]=r[3*i+2]; }
        return ret;
    }

    std::vector<double> get_velocity(unsigned i) const {
        std::vector<double> ret(3,0.0);
        if(i<nobj) {
            ret[0]=v[3*i]; ret[1]=v[3*i+1]; ret[2]=v[3*i+2]; }
        return ret;
    }

    std::vector<double> get_state(unsigned i) const {
        std::vector<double> ret(6,0.0);
        if(i<nobj) {
            ret[0]=r[3*i  ]; ret[1]=r[3*i+1]; ret[2]=r[3*i+2]; 
            ret[3]=v[3*i  ]; ret[4]=v[3*i+1]; ret[5]=v[3*i+2];
        }
        return ret;
    }

    double get_eccentricity(unsigned i,int j=-1) const {
        double ub; std::vector<double> rb,vb;
        if(j>-1) { ub = u[i]+u[j]; rb = get_position(j); vb = get_velocity(j); } 
        else calc_bary(ub,rb,vb);
        std::vector<double> r = get_position(i);
        std::vector<double> v = get_velocity(i);
        for(int k=0;k<3;k++) { r[k] -= rb[k]; v[k] -= vb[k]; }
        std::vector<double> h(3);
        h[0] = r[1]*v[2] - r[2]*v[1];
        h[1] = r[2]*v[0] - r[0]*v[2];
        h[2] = r[0]*v[1] - r[1]*v[0];
        const double rm = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        const double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
        const double h2 = h[0]*h[0] + h[1]*h[1] + h[2]*h[2];
        const double alpha = (2.*ub - v2*rm)/(ub*rm);
        if(fabs(alpha) > 1e-15) return sqrt(1.-(h2*alpha/ub));
        else return 1.0;
    }

    void calc_bary(double &ub,std::vector<double> &rb,std::vector<double> &vb) const {
        ub = 0; rb = std::vector<double>(3,0.0); vb = std::vector<double>(3,0.0);
        for(unsigned i=0;i<nobj;i++) {
            ub += u[i];
            for(int j=0;j<3;j++) {
                rb[j] += r[3*i+j]*u[i];
                vb[j] += v[3*i+j]*u[i];
            }
        }
        for(int j=0;j<3;j++) {
            rb[j] /= ub;
            vb[j] /= ub;
        }
    }

    void move2bary() {
        double ub; std::vector<double> rb,vb;
        calc_bary(ub,rb,vb);
        for(unsigned i=0;i<nobj;i++) {
            for(int j=0;j<3;j++) {
                r[3*i+j] -= rb[j];
                v[3*i+j] -= vb[j];
            }
        }
    }

    std::vector<double> momentum() const {
        std::vector<double> L(3,0.0);
        for(unsigned i=0;i<nobj;i++)
            for(int j=0;j<3;j++) 
                L[j] += u[i]*v[3*i+j];
        return L;
    }

    std::vector<double> angular_momentum() const {
        std::vector<double> L(3,0.0);
        for(unsigned i=0;i<nobj;i++) {
            L[0] += u[i]*(r[3*i+1]*v[3*i+2] - r[3*i+2]*v[3*i+1]);
            L[1] += u[i]*(r[3*i+2]*v[3*i+0] - r[3*i+0]*v[3*i+2]);
            L[2] += u[i]*(r[3*i+0]*v[3*i+1] - r[3*i+1]*v[3*i+0]);
        }
        return L;
    }

    double hamiltonian() const {
        double T = 0.0, V = 0.0;
        for(unsigned i=0;i<nobj;i++) 
            T += 0.5*u[i]*( v[3*i]*v[3*i] + v[3*i+1]*v[3*i+1] + v[3*i+2]*v[3*i+2] );
        for(unsigned i=0;i<nobj-1;i++) 
            for(unsigned j=i+1;j<nobj;j++) {
                if(u[i]*u[j]==0.0) continue;
                const double dx = r[3*i  ]-r[3*j  ];
                const double dy = r[3*i+1]-r[3*j+1];
                const double dz = r[3*i+2]-r[3*j+2];
                V -= (u[i]*u[j])/sqrt( dx*dx + dy*dy + dz*dz );
            }
        return T+V;
    }

    bool impact_test() const {
        for(unsigned i=0;i<nobj-1;i++) {
            for(unsigned j=i+1;j<nobj;j++) {
                const double r2 = R[i]+R[j];
                if(r2<=0) continue;
                const double dx = r[3*i  ]-r[3*j  ];
                const double dy = r[3*i+1]-r[3*j+1];
                const double dz = r[3*i+2]-r[3*j+2];
                const double dist = dx*dx + dy*dy + dz*dz;
                if(dist<=r2*r2) {
                    if(verbose) std::cout<<"# "<<name[i]<<" and "<<name[j]<<" impacted!\n";
                    return true;
                }
            }
        }
        return false;
    }

    bool distance_test() const {
        for(unsigned i=0;i<nobj;i++) {
            const double dist2 = r[3*i]*r[3*i] + r[3*i+1]*r[3*i+1] + r[3*i+2]*r[3*i+2];
            if(dist2>maxdist*maxdist) {
                if(verbose) std::cout<<"# "<<name[i]<<" reached maximum distance!\n";
                return true;
            }
        }
        return false;
    }

    bool check_e() const {
        for(unsigned i=0;i<emax_i.size();i++) {
            if(get_eccentricity(emax_i[i],emax_j[i])>emax_e[i]) {
                if(verbose) std::cout<<"# "<<name[emax_i[i]]<<" exceeded eccentricity of "<<emax_e[i]<<"!\n";
                return true;
            }
        }
        return false;
    }

    bool is_sane() const {
        if(impact_test()) return false;
        if(maxdist>=0. && distance_test()) return false;
        if(echeck && check_e()) return false;
        return true;
    }

    int evolve_rkn(double tgoal,double tol=1e-12) {
        RKN_Integrator<NBody> rk(neq,tol,true);
        rk.integrate(*this,tgoal);
        return rk.ifail;
    }

};

#endif
