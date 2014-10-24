#ifndef __ORBIT_HPP__
#define __ORBIT_HPP__

#include <iostream>
#include "conic.hpp"
#include "rkn670_4.hpp"

class NBody { public:

    double T, maxdist, mindist;
    unsigned nobj; //, iob;
    int failed;
    std::vector<double> r,v,a, u,R,maxe;
    std::vector<std::string> name;

    bool verbose;
    //Orbit helio_orb;
    //bool oblate, helio;

    NBody() {
        T = 0.;
        nobj = 0;
        failed = 0;
        verbose = false;
        //helio = false;
        maxdist = 0;
        mindist = -1;
    }

    unsigned add_object(const std::string &name0,double u0,double R0,const std::vector<double> &r0,const std::vector<double> &v0) {
        name.push_back(name0);
        u.push_back(u0); R.push_back(R0); maxe.push_back(1);
        for(unsigned i=0;i<3;i++) {
            r.push_back(r0[i]); v.push_back(v0[i]); a.push_back(0.);
        }
        nobj++;
        return nobj-1;
    }

    unsigned add_object_state(const std::string &name0,double u0,double R0,const std::vector<double> &s0) {
        std::vector<double> r0(3),v0(3);
        for(int i=0;i<3;i++) { r0[i]=s0[i]; v0[i]=s0[i+3]; }
        return add_object(name0,u0,R0,r0,v0);
    }

    /*void set_helio(const Orbit &orb) {
        helio = true; helio_orb = orb;
    }*/

    int lookup(const std::string &name0) const {
        for(unsigned i=0;i<nobj;i++) {
            if(name[i].compare(name0)==0)
                return i; }
        return -1;
    }

    std::vector<double> get_position(unsigned i) const {
        return {r[3*i],r[3*i+1],r[3*i+2]};
    }

    std::vector<double> get_velocity(unsigned i) const {
        return {v[3*i],v[3*i+1],v[3*i+2]};
    }

    std::vector<double> get_state(unsigned i) const {
        return {r[3*i],r[3*i+1],r[3*i+2],v[3*i],v[3*i+1],v[3*i+2]};
    }
    
    Conic orbit(unsigned i,int j) const {
        double ut; std::vector<double> rj,vj;
        if(j<0) calc_bary(ut,rj,vj);
        else {
            ut=u[j]+u[i];
            rj = get_position(j);
            vj = get_velocity(j);
        }
        auto ri = vec_math::sub(get_position(i),rj);
        auto vi = vec_math::sub(get_velocity(i),vj);
        return Conic(T,ut,ri,vi);
    }

    void first_gravity(const std::vector<double> &rr, std::vector<double> &aa) const {
        for(unsigned i=0;i<nobj;i++) {
            for(unsigned j=i+1;j<nobj;j++) {
                const unsigned i2=i*3, j2=j*3;
                const double dx = rr[i2  ]-rr[j2  ];
                const double dy = rr[i2+1]-rr[j2+1];
                const double dz = rr[i2+2]-rr[j2+2];
                const double dr = sqrt(dx*dx + dy*dy + dz*dz);
                const double dr3 = dr*dr*dr;
                if (u[i]>0.) {
                    aa[j2  ] += u[i]*dx/dr3;
                    aa[j2+1] += u[i]*dy/dr3;
                    aa[j2+2] += u[i]*dz/dr3;
                }
                if (u[j]>0.) {
                    aa[i2  ] -= u[j]*dx/dr3;
                    aa[i2+1] -= u[j]*dy/dr3;
                    aa[i2+2] -= u[j]*dz/dr3;
                }
            }
        }
    }

    /*void helio_gravity(const double &t, const vector<double> &rr, vector<double> &aa) const {
        Vector3d rh,vh; helio_orb.rv(t,rh,vh);
        double rhm = rh.norm(), hg = -helio_orb.u/(rhm*rhm*rhm);
        for(unsigned j=0;j<nobj;j++) {
            for(unsigned k=0;k<3;k++) 
                aa[3*j+k] -= hg*rh[k];
        }
    }*/

    void dt(const double &t, const std::vector<double> &rr, std::vector<double> &aa) const {
        for(unsigned i=0;i<nobj*3;i++) aa[i] = 0.0;
        first_gravity(rr,aa);
        //if(helio) helio_gravity(t,rr,aa);
    }

    void calc_bary(double &ub,std::vector<double> &rb,std::vector<double> &vb) const {
        using namespace vec_math;
        ub = 0; rb = {0,0,0}; vb = {0,0,0};
        for(unsigned i=0;i<nobj;i++) {
            ub += u[i];
            rb = add(rb, mult(u[i],get_position(i)) );
            vb = add(vb, mult(u[i],get_velocity(i)) );
        }
        rb = div(rb,ub);
        vb = div(vb,ub);
    }

    void move2bary() {
        double ub; std::vector<double> rb,vb;
        calc_bary(ub,rb,vb);
        for(unsigned i=0;i<nobj;i++) {
            for(unsigned j=0;j<3;j++) {
                r[3*i+j] -= rb[j];
                v[3*i+j] -= vb[j];
            }
        }
    }

    void print_state() const {
        std::cout.precision(15);
        std::cout<<"T = "<<std::scientific<<T<<'\n';
        for(unsigned i=0;i<nobj;i++)
            printf("%-10s %.15e %.8e %+.15e %+.15e %+.15e %+.15e %+.15e %+.15e\n",
                    name[i].c_str(),u[i],R[i],
                    r[3*i],r[3*i+1],r[3*i+2],v[3*i],v[3*i+1],v[3*i+2]);
    }

    void calc_dist(double &md,unsigned &id) const {
        double d; md = 0.; id = 0;
        for(unsigned i=0;i<nobj;i++) {
            d = vec_math::norm(get_position(i));
            if(d>md) { md=d; id=i; }
        }
    }

    bool impact_test(unsigned &id,unsigned &jd) {
        id = 0; jd = 0;
        for(unsigned i=0;i<nobj-1;i++) {
            for(unsigned j=i+1;j<nobj;j++) {
                double r2 = R[i]+R[j];
                if(r2<=0) continue;
                double dist = vec_math::norm(vec_math::sub(get_position(i),get_position(j)));
                if(dist<=r2) {
                    id = i; jd = j;
                    if(verbose) std::cout<<"# "<<name[i]<<" and "<<name[j]<<" impacted!\n";
                    return true;
                }
                if(mindist<0 or mindist>dist) mindist = dist;
            }
        }
        return false;
    }

    bool e_test() const {
        for(unsigned i=0;i<nobj;i++) {
            if(maxe[i]>=1.-1e-5) continue;
            auto orb = orbit(i,-1);
            if(orb.e>maxe[i]) {
                if(verbose) std::cout<<"# "<<name[i]<<" reached an eccentricity of "<<orb.e<<"!\n";
                return true;
            }
        }
        return false;
    }

    bool is_sane() {
        double md; unsigned id,jd;
        if(maxdist>0) {
            calc_dist(md,id);
            if(md>maxdist) {
                if(verbose) std::cout<<"# "<<name[id]<<" exceeded maximum distance!\n";
                return false;
            }
        }
        if(impact_test(id,jd)) return false;
        if(e_test()) return false;
        return true;
    }

    void evolve_self(double tgoal) {
        RKN670<double,NBody> rkn(nobj*3,1e-12);
        rkn.integrate(*this,tgoal);
        failed = rkn.failed;
    }

};

void rkn_evolve(NBody &nb,double tgoal) {
    RKN670<double,NBody> rkn(nb.nobj*3,1e-12);
    rkn.integrate(nb,tgoal);
    nb.failed = rkn.failed;
}

#endif
