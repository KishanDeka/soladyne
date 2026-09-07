#include "solar_dynamo/model.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace solar_dynamo {

Model::Model(Parameters p):parameters_(p),
 vp_(halfmax+1),vq_(halfmax+1),vp1_(halfmax+1),vq1_(halfmax+1),
 psi_(halfmax+1),etab2_(halfmax+1),dvp_(halfmax+1),vpb_(halfmax+1),deta_(halfmax+1),
 rho2_(halfmax+1),nut2_(halfmax+1),ft1_(halfmax+1),gt1_(halfmax+1),
 ft2_(halfmax+1),gt2_(halfmax+1),fr1_(halfmax+1),gr1_(halfmax+1),
 fr2_(halfmax+1),gr2_(halfmax+1),v00_(halfmax+1),h1_(halfmax+1),
 radius_(nmax+1),theta_(nmax+1),ra_(nmax+1),smooth_(halfmax+1),
 ss1_p_(lmax+1),sn_(nmax+1),pl_((nmax+1)*(lmax+1),0.0){
  if(parameters_.dt<=0.0||parameters_.tmax<0.0)throw std::invalid_argument("dt must be positive and tmax non-negative");
  pi_=4.0*std::atan(1.0);pm_=6.96;pb_=0.55*pm_;pk_=0.73*pm_;pw_=2.5*pm_;qm_=pi_;
  dp_=(pm_-pb_)/double(n);dq_=-qm_/double(n);fac_=1.0e8/(3600.0*24.0*365.0);
}

double Model::legacy_erf(double x){
  constexpr double p=.3275911,a1=.254829592,a2=-.284496736,a3=1.421413741,a4=-1.453152027,a5=1.061405429;
  const double sig=x<0?-1.0:1.0,xa=std::abs(x);if(xa>20)return sig;
  const double t=1.0/(1.0+p*xa);return sig*(1.0-((((a5*t+a4)*t+a3)*t+a2)*t+a1)*t*std::exp(-xa*xa));
}

double Model::associated_legendre(int l,int m,double x){
  if(m<0||m>l||std::abs(x)>1.0)throw std::invalid_argument("bad Legendre arguments");
  double pmm=1.0;
  if(m>0){double somx2=std::sqrt((1.0-x)*(1.0+x)),fact=1.0;for(int i=1;i<=m;++i){pmm=-pmm*fact*somx2;fact+=2.0;}}
  if(l==m)return pmm;
  double pmmp1=x*(2*m+1)*pmm;
  if(l==m+1)return pmmp1;
  double pll=0;for(int ll=m+2;ll<=l;++ll){pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);pmm=pmmp1;pmmp1=pll;}return pll;
}

void Model::build_profiles(){
  const double pbc=.71*pm_,rhobc=.2e-6,hbc=.58,gmm=5.0/3.0;
  for(int i=1;i<=2*nmax;++i){double r=pb_+double(i-1)/double(2*n)*(pm_-pb_);smooth_[i]=.5*(1+legacy_erf((r-pk_)/(.01*pm_)));}
  for(int i=1;i<=nmax;++i){double r=pb_+double(i-1)/n*(pm_-pb_);radius_[i]=r;ra_[i]=r;
    for(int j=1;j<=nmax;++j)rho_(i,j)=rhobc*std::pow(1+((gmm-1)/gmm)*(pbc/hbc)*((pbc/r)-1),1/(gmm-1));}
  for(int i=1;i<=2*n+1;++i){double r=pb_+double(i-1)*(pm_-pb_)/(2*n);for(int j=1;j<=2*n+1;++j)
    rho2_(i,j)=rhobc*std::pow(1+((gmm-1)/gmm)*(pbc/hbc)*((pbc/r)-1),1/(gmm-1));}
  for(int i=1;i<=nmax;++i){double r=radius_[i];for(int j=1;j<=nmax;++j){double q=qm_-double(j-1)/n*qm_,co=std::cos(q);theta_[j]=q;
    alpha_(i,j)=.25*co*(1+legacy_erf((r-.95*pm_)/(.03*pm_)))*(1-legacy_erf((r-pm_)/(.03*pm_)));
    eta_(i,j)=.000220+(parameters_.et0/2)*(1+legacy_erf((r-.7*pm_)/(.03*pm_)));
    etab_(i,j)=.00022+(parameters_.et1/2)*(1+legacy_erf((r-.725*pm_)/(.03*pm_)))+(parameters_.et0/2)*(1+legacy_erf((r-.975*pm_)/(.03*pm_)));
    nu_d_(i,j)=80*(1+(4.5-1)*smooth_[2*i-1])*etab_(i,j);
    dom_(i,j)=271.9+.5*(1+legacy_erf((r-.7*pm_)/(.03*pm_)))*(289.5-39.4*co*co-42.2*std::pow(co,4)-271.9);
  }}
  for(int i=2;i<=n;++i){double r=radius_[i];for(int j=2;j<=n;++j){double q=theta_[j];
    wr_(i,j)=(std::pow(r+dp_,4)*rho_(i+1,j)*nu_d_(i+1,j)-std::pow(r-dp_,4)*rho_(i-1,j)*nu_d_(i-1,j))/(2*dp_*rho_(i,j)*std::pow(r,4));
    wt_(i,j)=(rho_(i,j+1)*nu_d_(i,j+1)*std::pow(std::sin(q+dq_),3)-rho_(i,j-1)*nu_d_(i,j-1)*std::pow(std::sin(q-dq_),3))/(2*dq_*rho_(i,j)*r*r*std::pow(std::sin(q),3));
  }}
  for(int i=1;i<=2*n+1;++i){double r=pb_+double(i-1)*(pm_-pb_)/(2*n);for(int j=1;j<=2*n+1;++j){
    etab2_(i,j)=.00022+(parameters_.et1/2)*(1+legacy_erf((r-.725*pm_)/(.03*pm_)))+(parameters_.et0/2)*(1+legacy_erf((r-.975*pm_)/(.03*pm_)));
    nut2_(i,j)=(.01+(1-.01)*smooth_[i])*etab2_(i,j);
  }}
  for(int i=2;i<=n;++i)for(int j=2;j<=n;++j){dror_(i,j)=(dom_(i+1,j)-dom_(i-1,j))/(2*dp_);drot_(i,j)=(dom_(i,j+1)-dom_(i,j-1))/(2*dq_);}
  const double beta1=1.5,beta2=1.3,pp=.73*pm_,del=2.0000001,gm=3.47,p0=(pm_-pb_)/3.5;
  for(int i=1;i<=2*n+1;++i){double r=pb_+double(i-1)/double(2*n)*(pm_-pb_);for(int j=1;j<=2*n+1;++j){double q=qm_-double(j-1)/double(2*n)*qm_;
    double exq0,exqm,gau=std::exp(-1.05*std::pow((r-p0)/gm,2));
    if(q<=pi_/2){exq0=std::exp(-beta1*std::pow(q,del));exqm=std::exp(beta2*(q-pi_/2));psi_(i,j)=(r-pp)*std::sin(pi_*(r-pp)/(pm_-pp))*(1-exq0)*(1-exqm)*gau;}
    else{exq0=std::exp(-beta1*std::pow(pi_-q,del));exqm=std::exp(beta2*(pi_/2-q));psi_(i,j)=-(r-pp)*std::sin(pi_*(r-pp)/(pm_-pp))*(1-exq0)*(1-exqm)*gau;}
    if(r<pp)psi_(i,j)=0;
  }}
  for(int i=2;i<=2*n;++i){double r=pb_+double(i-1)*(pm_-pb_)/(2*n);for(int j=2;j<=2*n;++j){double q=qm_-double(j-1)*qm_/(2*n);
    vp1_(i,j)=2.4573*parameters_.v0*(psi_(i,j+1)-psi_(i,j-1))/(r*r*std::sin(q)*dq_*std::pow(pm_/r-.95,1.5));
    vq1_(i,j)=-2.4573*parameters_.v0*(psi_(i+1,j)-psi_(i-1,j))/(r*std::sin(q)*dp_*std::pow(pm_/r-.95,1.5));
    deta_(i,j)=(etab2_(i+1,j)-etab2_(i-1,j))/dp_;
    vpb_(i,j)=(vp1_(i,j)-deta_(i,j))*parameters_.dt/(2*r*dp_);vp_(i,j)=vp1_(i,j)*parameters_.dt/(2*r*dp_);
    vq_(i,j)=vq1_(i,j)*parameters_.dt/(2*r*std::sin(q)*dq_);
  }}
  for(int i=1;i<=2*n+1;++i){double r=pb_+double(i-1)*(pm_-pb_)/(2*n),pr=r/pm_;for(int j=1;j<=2*n+1;++j){double q=qm_-double(j-1)*qm_/(2*n);
    double c0=110.80*std::exp(-602.3235*std::pow(pr,3)+1496.1566*pr*pr-1245.933*pr+344.6711);
    double z0=(9-(2*c0*c0)/(1+c0*c0)-(c0*c0+9)*std::atan(c0)/c0)/(2*std::pow(c0,4));
    double z1=(45+c0*c0-(4*c0*c0)/(1+c0*c0)+(std::pow(c0,4)-12*c0*c0-45)*std::atan(c0)/c0)/(2*std::pow(c0,4));
    double y0=(-19-5/(1+c0*c0)+(3*c0*c0+24)*std::atan(c0)/c0)/(4*std::pow(c0,4));
    double y1=3*(-15+2*c0*c0/(1+c0*c0)+(3*c0*c0+15)*std::atan(c0)/c0)/(4*std::pow(c0,4));
    v00_(i,j)=100*smooth_[i]*(z0+4*y0);h1_(i,j)=900*smooth_[i]*(z1+4*y1);
    double s=std::sin(q);
    if(std::abs(s)>1e-15){ft1_(i,j)=vq1_(i,j)/(r*s*s);ft2_(i,j)=1/(rho2_(i,j)*r*r*std::pow(s,3));}
    gt1_(i,j)=s*s;gt2_(i,j)=rho2_(i,j)*nut2_(i,j)*std::cos(q)*std::pow(s,4)*h1_(i,j);
    fr1_(i,j)=vp1_(i,j)/(r*r);gr1_(i,j)=r*r;fr2_(i,j)=1/(rho2_(i,j)*std::pow(r,4));
    gr2_(i,j)=rho2_(i,j)*nut2_(i,j)*std::pow(r,3)*(v00_(i,j)-h1_(i,j)*std::pow(std::cos(q),2));
  }}
  for(int l=1;l<=lmax;++l){double cl=1/(1+double(l)/(l+1)*std::pow(pm_/pw_,2*l+1));ss1_p_[l]=double(2*l+1)/(l*pm_)*(cl-double(l)/(l+1)*(1-cl));
    for(int j=1;j<=n+1;++j){double q=pi_+double(j-1)*dq_;pl_[j*(lmax+1)+l]=associated_legendre(l,1,std::cos(q));sn_[j]=std::sin(q);}}
}

void Model::load_initial_state(const std::filesystem::path& path){
  if(!parameters_.relaxed_initial_state){for(int i=1;i<=n+1;++i){double r=pb_+double(i-1)/n*(pm_-pb_);for(int j=1;j<=n+1;++j){double q=qm_-double(j-1)/n*qm_;u_(i,j)=0;ub_(i,j)=std::sin(2*q)*std::sin(pi_*(r-pb_)/(pm_-pb_));}}return;}
  std::ifstream in(path);if(!in)throw std::runtime_error("cannot open init.dat: "+path.string());
  double q,r,uin,btor;for(int i=1;i<=n+1;++i)for(int j=1;j<=n+1;++j){if(!(in>>q>>r>>uin>>btor))throw std::runtime_error("init.dat must contain 66049 four-column rows");ub_(i,j)=btor;u_(i,j)=j==n+1?0.0:uin/(r*std::sin(q));}
  double extra;if(in>>extra)throw std::runtime_error("init.dat contains more than 66049 rows");
}

void Model::initialize(const std::filesystem::path& init_file){
  build_profiles();load_initial_state(init_file);for(int i=1;i<=nmax;++i)for(int j=1;j<=nmax;++j)omega_(i,j)=parameters_.omega0;
  build_coefficients();time_=0;step_number_=0;initialized_=true;
}

void Model::build_coefficients(){
  const double dt=parameters_.dt;
  for(int i=2;i<=2*n;++i)for(int j=2;j<=2*n;++j)dvp_(i,j)=(vp1_(i+1,j)-vp1_(i-1,j))/dp_;
  for(int i=2;i<=n;++i){double p=pb_+double(i-1)*(pm_-pb_)/n;for(int j=2;j<=n;++j){double q=qm_-double(j-1)*qm_/n;
    oa_(i,j)=dt*wt_(i,j)/(4*dq_)-dt*nu_d_(i,j)/(2*std::pow(p*dq_,2))
      -(dt*ft1_(2*i-1,2*j-1)*gt1_(2*i-1,2*j-2)/(4*dq_))*(1+dt*ft1_(2*i-1,2*j-2)*gt1_(2*i-1,2*j-3)/(2*dq_))
      -(dt*ft2_(2*i-1,2*j-1)*gt2_(2*i-1,2*j-2)/(4*dq_))*(1+dt*ft2_(2*i-1,2*j-2)*gt2_(2*i-1,2*j-3)/(2*dq_));
    ob_(i,j)=dt*nu_d_(i,j)/std::pow(p*dq_,2)
      +(dt*ft1_(2*i-1,2*j-1)*gt1_(2*i-1,2*j)/(4*dq_))*(1+dt*ft1_(2*i-1,2*j)*gt1_(2*i-1,2*j-1)/(2*dq_))
      -(dt*ft1_(2*i-1,2*j-1)*gt1_(2*i-1,2*j-2)/(4*dq_))*(1-dt*ft1_(2*i-1,2*j-2)*gt1_(2*i-1,2*j-1)/(2*dq_))
      +(dt*ft2_(2*i-1,2*j-1)*gt2_(2*i-1,2*j)/(4*dq_))*(1+dt*ft2_(2*i-1,2*j)*gt2_(2*i-1,2*j-1)/(2*dq_))
      -(dt*ft2_(2*i-1,2*j-1)*gt2_(2*i-1,2*j-2)/(4*dq_))*(1-dt*ft2_(2*i-1,2*j-2)*gt2_(2*i-1,2*j-1)/(2*dq_));
    oc_(i,j)=-dt*wt_(i,j)/(4*dq_)-dt*nu_d_(i,j)/(2*std::pow(p*dq_,2))
      +(dt*ft1_(2*i-1,2*j-1)*gt1_(2*i-1,2*j)/(4*dq_))*(1-dt*ft1_(2*i-1,2*j)*gt1_(2*i-1,2*j+1)/(2*dq_))
      +(dt*ft2_(2*i-1,2*j-1)*gt2_(2*i-1,2*j)/(4*dq_))*(1-dt*ft2_(2*i-1,2*j)*gt2_(2*i-1,2*j+1)/(2*dq_));
    od_(i,j)=dt*wr_(i,j)/(4*dp_)-dt*nu_d_(i,j)/(2*dp_*dp_)
      -(dt*fr1_(2*i-1,2*j-1)*gr1_(2*i-2,2*j-1)/(4*dp_))*(1+dt*fr1_(2*i-2,2*j-1)*gr1_(2*i-3,2*j-1)/(2*dp_))
      -(dt*fr2_(2*i-1,2*j-1)*gr2_(2*i-2,2*j-1)/(4*dp_))*(1+dt*fr2_(2*i-2,2*j-1)*gr2_(2*i-3,2*j-1)/(2*dp_));
    oe_(i,j)=dt*nu_d_(i,j)/(dp_*dp_)
      +(dt*fr1_(2*i-1,2*j-1)*gr1_(2*i,2*j-1)/(4*dp_))*(1+dt*fr1_(2*i,2*j-1)*gr1_(2*i-1,2*j-1)/(2*dp_))
      +(dt*fr2_(2*i-1,2*j-1)*gr2_(2*i,2*j-1)/(4*dp_))*(1+dt*fr2_(2*i,2*j-1)*gr2_(2*i-1,2*j-1)/(2*dp_))
      -(dt*fr1_(2*i-1,2*j-1)*gr1_(2*i-3,2*j-1)/(4*dp_))*(1-dt*fr1_(2*i-3,2*j-1)*gr1_(2*i-1,2*j-1)/(2*dp_))
      -(dt*fr2_(2*i-1,2*j-1)*gr2_(2*i-3,2*j-1)/(4*dp_))*(1-dt*fr2_(2*i-3,2*j-1)*gr2_(2*i-1,2*j-1)/(2*dp_));
    of_(i,j)=-dt*wr_(i,j)/(4*dp_)-dt*nu_d_(i,j)/(2*dp_*dp_)
      +(dt*fr1_(2*i-1,2*j-1)*gr1_(2*i,2*j-1)/(4*dp_))*(1-dt*fr1_(2*i,2*j-1)*gr1_(2*i+1,2*j-1)/(2*dp_))
      +(dt*fr2_(2*i-1,2*j-1)*gr2_(2*i,2*j-1)/(4*dp_))*(1-dt*fr2_(2*i,2*j-1)*gr2_(2*i+1,2*j-1)/(2*dp_));

    ma_(i,j)=-(eta_(i,j)*dt/(2*std::pow(p*dq_,2))-eta_(i,j)*dt/(4*std::tan(q)*p*p*dq_)+vq_(2*i-1,2*j-1)*std::sin(q-dq_/2)*(1+vq_(2*i-1,2*j-2)*std::sin(q-dq_))/2);
    mb_(i,j)=-(-eta_(i,j)*dt/std::pow(p*dq_,2)-eta_(i,j)*dt/(4*std::pow(p*std::sin(q),2))-vq_(2*i-1,2*j-1)*(std::sin(q+dq_/2)*(1+vq_(2*i-1,2*j)*std::sin(q))-std::sin(q-dq_/2)*(1-vq_(2*i-1,2*j-2)*std::sin(q)))/2);
    mc_(i,j)=-(eta_(i,j)*dt/(2*std::pow(p*dq_,2))+eta_(i,j)*dt/(4*std::tan(q)*p*p*dq_)-vq_(2*i-1,2*j-1)*std::sin(q+dq_/2)*(1-vq_(2*i-1,2*j)*std::sin(q+dq_))/2);
    md_(i,j)=-(eta_(i,j)*dt/(2*dp_*dp_)-eta_(i,j)*dt/(2*p*dp_)+vp_(2*i-1,2*j-1)*(p-dp_/2)*(1+vp_(2*i-2,2*j-1)*(p-dp_))/2);
    me_(i,j)=-(-eta_(i,j)*dt/(dp_*dp_)-eta_(i,j)*dt/(4*std::pow(p*std::sin(q),2))-vp_(2*i-1,2*j-1)*(dp_+p*((p+dp_/2)*vp_(2*i,2*j-1)+(p-dp_/2)*vp_(2*i-2,2*j-1)))/2);
    mf_(i,j)=-(eta_(i,j)*dt/(2*dp_*dp_)+eta_(i,j)*dt/(2*p*dp_)-vp_(2*i-1,2*j-1)*(p+dp_/2)*(1-vp_(2*i,2*j-1)*(p+dp_))/2);
    ta_(i,j)=-(etab_(i,j)*dt/(2*std::pow(p*dq_,2))-etab_(i,j)*dt/(4*std::tan(q)*p*p*dq_)+vq_(2*i-1,2*j-2)*std::sin(q-dq_/2)*(1+vq_(2*i-1,2*j-3)*std::sin(q-dq_))/2);
    tb_(i,j)=-(-etab_(i,j)*dt/std::pow(p*dq_,2)-etab_(i,j)*dt/(4*std::pow(p*std::sin(q),2))-(vq_(2*i-1,2*j)*std::sin(q+dq_/2)*(1+vq_(2*i-1,2*j-1)*std::sin(q))-vq_(2*i-1,2*j-2)*std::sin(q-dq_/2)*(1-vq_(2*i-1,2*j-1)*std::sin(q)))/2);
    tc_(i,j)=-(etab_(i,j)*dt/(2*std::pow(p*dq_,2))+etab_(i,j)*dt/(4*std::tan(q)*p*p*dq_)-vq_(2*i-1,2*j)*std::sin(q+dq_/2)*(1-vq_(2*i-1,2*j+1)*std::sin(q+dq_))/2);
    td_(i,j)=-(etab_(i,j)*dt/(2*dp_*dp_)-etab_(i,j)*dt/(2*p*dp_)+vpb_(2*i-1,2*j-1)*(p-dp_/2)*(1+vpb_(2*i-2,2*j-1)*(p-dp_))/2);
    te_(i,j)=-(-etab_(i,j)*dt/(dp_*dp_)-etab_(i,j)*dt/(4*std::pow(p*std::sin(q),2))-dvp_(2*i-1,2*j-1)*dt/2-vpb_(2*i-1,2*j-1)*(dp_+vpb_(2*i,2*j-1)*p*(p+dp_/2)+(p-dp_/2)*p*vpb_(2*i-2,2*j-1))/2);
    tf_(i,j)=-(etab_(i,j)*dt/(2*dp_*dp_)+etab_(i,j)*dt/(2*p*dp_)-vpb_(2*i-1,2*j-1)*(p+dp_/2)*(1-vpb_(2*i,2*j-1)*(p+dp_))/2);
  }}
}

void Model::tridag(const std::vector<double>& a,const std::vector<double>& b,const std::vector<double>& c,const std::vector<double>& r,std::vector<double>& u,int size){
  std::vector<double> gam(size+1);if(b[1]==0)throw std::runtime_error("tridag zero pivot");double bet=b[1];u[1]=r[1]/bet;
  for(int j=2;j<=size;++j){gam[j]=c[j-1]/bet;bet=b[j]-a[j]*gam[j];if(bet==0)throw std::runtime_error("tridag zero pivot");u[j]=(r[j]-a[j]*u[j-1])/bet;}
  for(int j=size-1;j>=1;--j)u[j]-=gam[j+1]*u[j+1];
}

double Model::simpson_ss(int l)const{double t1=0,t2=0;for(int k=2;k<=n;k+=2)t1+=pl_[k*(lmax+1)+l]*u_(n+1,k)*sn_[k];for(int k=3;k<=n-1;k+=2)t2+=pl_[k*(lmax+1)+l]*u_(n+1,k)*sn_[k];
  double y1=pl_[(lmax+1)+l]*u_(n+1,1)*sn_[1],yn=pl_[(n+1)*(lmax+1)+l]*u_(n+1,n+1)*sn_[n+1];return dq_/3*(y1+4*t1+2*t2+yn);}
double Model::upper_boundary_coefficient(const std::vector<double>& ss1,int j)const{double cf=0;for(int l=1;l<=lmax;++l)cf+=ss1[l]*pl_[j*(lmax+1)+l];return cf;}

void Model::advance_one(){
  const double dt=parameters_.dt,mu0=4*pi_;
  std::vector<double>a(nmax+1),b(nmax+1),c(nmax+1),r(nmax+1),x(nmax+1);
  for(int i=2;i<=n;++i){double p=radius_[i];for(int j=2;j<=n;++j){double q=theta_[j];
    lorentz_(i,j)=1/(mu0*rho_(i,j)*std::sin(q)*p*p)*(
      ((std::sin(q+dq_)*u_(i,j+1)-std::sin(q-dq_)*u_(i,j-1))/(2*dq_))*(((p+dp_)*ub_(i+1,j)-(p-dp_)*ub_(i-1,j))/(2*dp_))
      -(((p+dp_)*u_(i+1,j)-(p-dp_)*u_(i-1,j))/(2*dp_))*((std::sin(q+dq_)*ub_(i,j+1)-std::sin(q-dq_)*ub_(i,j-1))/(2*dq_)));
  }}
  for(int j=2;j<=n;++j)for(int i=2;i<=n;++i)phiw_(i,j)=-od_(i,j)*omega_(i-1,j)+(1-oe_(i,j))*omega_(i,j)-of_(i,j)*omega_(i+1,j)+lorentz_(i,j)*dt/2;
  for(int i=2;i<=n;++i){for(int j=2;j<=n;++j){a[j-1]=oa_(i,j);b[j-1]=ob_(i,j)+1;c[j-1]=oc_(i,j);r[j-1]=phiw_(i,j);}r[1]=phiw_(i,2)-a[1]*omega_(i,1);r[n-1]=phiw_(i,n)-c[n-1]*omega_(i,n+1);tridag(a,b,c,r,x,n-1);for(int j=2;j<=n;++j){omega_(i,j)=x[j-1];phiw_(i,j)=-phiw_(i,j)+2*x[j-1];}}
  for(int i=2;i<=n;++i)for(int j=2;j<=n;++j){dror_(i,j)=(omega_(i+1,j)-omega_(i-1,j))/(2*dp_);drot_(i,j)=(omega_(i,j+1)-omega_(i,j-1))/(2*dq_);}
  for(int i=2;i<=n;++i)for(int j=2;j<=n;++j)phi_(i,j)=-md_(i,j)*u_(i-1,j)+(1-me_(i,j))*u_(i,j)-mf_(i,j)*u_(i+1,j)+parameters_.al0*alpha_(i,j)*ub_(i,j)*dt/2;
  for(int i=2;i<=n;++i){for(int j=2;j<=n;++j){a[j-1]=ma_(i,j);b[j-1]=mb_(i,j)+1;c[j-1]=mc_(i,j);r[j-1]=phi_(i,j);}r[1]=phi_(i,2)-a[1]*u_(i,1);r[n-1]=phi_(i,n)-c[n-1]*u_(i,n+1);tridag(a,b,c,r,x,n-1);for(int j=2;j<=n;++j){u_(i,j)=x[j-1];phi_(i,j)=-phi_(i,j)+2*x[j-1];}}
  for(int i=2;i<=n;++i){double p=radius_[i];for(int j=2;j<=n;++j){double q=theta_[j];double br=(u_(i,j+1)*std::sin(q+dq_)-u_(i,j-1)*std::sin(q-dq_))*dt/(4*dq_);double bt=(u_(i-1,j)*(p-dp_)-u_(i+1,j)*(p+dp_))*std::sin(q)*dt/(4*p*dp_);phib_(i,j)=-td_(i,j)*ub_(i-1,j)+(1-te_(i,j))*ub_(i,j)-tf_(i,j)*ub_(i+1,j)+br*dror_(i,j)+bt*drot_(i,j);}}
  for(int i=2;i<=n;++i){for(int j=2;j<=n;++j){a[j-1]=ta_(i,j);b[j-1]=tb_(i,j)+1;c[j-1]=tc_(i,j);r[j-1]=phib_(i,j);}r[1]=phib_(i,2)-a[1]*ub_(i,1);r[n-1]=phib_(i,n)-c[n-1]*ub_(i,n+1);tridag(a,b,c,r,x,n-1);for(int j=2;j<=n;++j){ub_(i,j)=x[j-1];phib_(i,j)=-phib_(i,j)+2*x[j-1];phi_(i,j)+=parameters_.al0*alpha_(i,j)*ub_(i,j)*dt/2;}}
  for(int i=2;i<=n;++i){double p=radius_[i];for(int j=2;j<=n;++j){double q=theta_[j];lorentz_(i,j)=1/(mu0*rho_(i,j)*std::sin(q)*p*p)*(
    ((std::sin(q+dq_)*u_(i,j+1)-std::sin(q-dq_)*u_(i,j-1))/(2*dq_))*(((p+dp_)*ub_(i+1,j)-(p-dp_)*ub_(i-1,j))/(2*dp_))
    -(((p+dp_)*u_(i+1,j)-(p-dp_)*u_(i-1,j))/(2*dp_))*((std::sin(q+dq_)*ub_(i,j+1)-std::sin(q-dq_)*ub_(i,j-1))/(2*dq_)));phiw_(i,j)+=lorentz_(i,j)*dt/2;}}
  for(int j=2;j<=n;++j){for(int i=2;i<=n;++i){a[i-1]=od_(i,j);b[i-1]=oe_(i,j)+1;c[i-1]=of_(i,j);r[i-1]=phiw_(i,j);}r[1]=phiw_(2,j)-a[1]*omega_(1,j);r[n-1]=phiw_(n,j)-c[n-1]*omega_(n+1,j);tridag(a,b,c,r,x,n-1);for(int i=2;i<=n;++i)omega_(i,j)=x[i-1];}
  for(int i=2;i<=n;++i)for(int j=2;j<=n;++j){dror_(i,j)=(omega_(i+1,j)-omega_(i-1,j))/(2*dp_);drot_(i,j)=(omega_(i,j+1)-omega_(i,j-1))/(2*dq_);}
  for(int i=2;i<=n;++i){double p=radius_[i];for(int j=2;j<=n;++j){double q=theta_[j];double br=(u_(i,j+1)*std::sin(q+dq_)-u_(i,j-1)*std::sin(q-dq_))*dt/(4*dq_);double bt=(u_(i-1,j)*(p-dp_)-u_(i+1,j)*(p+dp_))*std::sin(q)*dt/(4*p*dp_);phib_(i,j)+=br*dror_(i,j)+bt*drot_(i,j);}}
  for(int j=2;j<=n;++j){for(int i=2;i<=n;++i){a[i-1]=md_(i,j);b[i-1]=me_(i,j)+1;c[i-1]=mf_(i,j);r[i-1]=phi_(i,j);}r[1]=phi_(2,j)-a[1]*u_(1,j);r[n-1]=phi_(n,j)-c[n-1]*u_(n+1,j);tridag(a,b,c,r,x,n-1);for(int i=2;i<=n;++i)u_(i,j)=x[i-1];}
  for(int j=2;j<=n;++j){for(int i=2;i<=n;++i){a[i-1]=td_(i,j);b[i-1]=te_(i,j)+1;c[i-1]=tf_(i,j);r[i-1]=phib_(i,j);}r[1]=phib_(2,j)-a[1]*ub_(1,j);r[n-1]=phib_(n,j)-c[n-1]*ub_(n+1,j);tridag(a,b,c,r,x,n-1);for(int i=2;i<=n;++i)ub_(i,j)=x[i-1];}
  apply_boundaries();
}

void Model::apply_boundaries(){
  for(int j=1;j<=nmax;++j){omega_(1,j)=parameters_.omega0;omega_(nmax,j)=omega_(n,j);}
  for(int i=2;i<=n;++i){omega_(i,1)=omega_(i,2);omega_(i,nmax)=omega_(i,n);}
  for(int j=1;j<=n+1;++j){u_(1,j)=0;ub_(1,j)=0;}
  std::vector<double> ss1(lmax+1);for(int pass=1;pass<=11;++pass){for(int l=1;l<=lmax;++l)ss1[l]=ss1_p_[l]*simpson_ss(l);for(int j=1;j<=n+1;++j)u_(n+1,j)=u_(n,j)+upper_boundary_coefficient(ss1,j)*dp_;}
  for(int j=1;j<=n+1;++j)ub_(n+1,j)=0;
  for(int i=2;i<=n+1;++i){u_(i,n+1)=0;ub_(i,n+1)=0;u_(i,1)=0;ub_(i,1)=0;}
}

void Model::step(std::size_t count){if(!initialized_)throw std::logic_error("call initialize() first");for(std::size_t k=0;k<count;++k){advance_one();time_+=parameters_.dt;++step_number_;}}
void Model::ensure_directory(const std::filesystem::path& p){if(!p.empty())std::filesystem::create_directories(p);}
void Model::step_with_output(std::size_t count,const std::filesystem::path& dir){if(!initialized_)throw std::logic_error("call initialize() first");ensure_directory(dir);for(std::size_t k=0;k<count;++k){advance_one();const auto legacy_k=step_number_+1;if(legacy_k%200==0)write_omega_latitude(dir/"omg_lat.dat");if(legacy_k%50==0)write_diffrot_snapshot(dir/"diffrot.dat");time_+=parameters_.dt;++step_number_;}}
void Model::run(const std::filesystem::path& dir){std::size_t target=static_cast<std::size_t>(parameters_.tmax/parameters_.dt);if(target>step_number_)step_with_output(target-step_number_,dir);write_final_outputs(dir);}

void Model::write_omega_latitude(const std::filesystem::path& path)const{
  std::ofstream out(path,std::ios::app);if(!out)throw std::runtime_error("cannot write "+path.string());int ir=static_cast<int>(1+((.85*pm_)-pb_)*n/(pm_-pb_));
  out<<std::fixed<<std::setprecision(7);for(int j=1;j<=nmax;++j)out<<std::setw(13)<<theta_[j]<<' '<<std::setw(13)<<time_<<' '<<std::setw(13)<<omega_(ir,j)<<'\n';
}
void Model::write_diffrot_snapshot(const std::filesystem::path& path)const{
  std::ofstream out(path,std::ios::trunc);if(!out)throw std::runtime_error("cannot write "+path.string());out<<std::fixed<<std::setprecision(5);
  for(int i=1;i<=nmax;++i)for(int j=1;j<=nmax;++j)out<<std::setw(13)<<theta_[j]<<' '<<std::setw(13)<<radius_[i]<<' '<<std::setw(13)<<omega_(i,j)<<' '<<std::setw(13)<<time_<<'\n';
}
void Model::write_final_outputs(const std::filesystem::path& dir)const{
  if(!initialized_)throw std::logic_error("call initialize() first");
  ensure_directory(dir);
  {std::ofstream out(dir/"diffrot.dat",std::ios::trunc);if(!out)throw std::runtime_error("cannot write diffrot.dat");out<<std::fixed<<std::setprecision(7);for(int i=1;i<=nmax;++i)for(int j=1;j<=nmax;++j)out<<std::setw(13)<<theta_[j]<<' '<<std::setw(13)<<radius_[i]<<' '<<std::setw(13)<<omega_(i,j)<<'\n';}
  {std::ofstream out(dir/"final.dat",std::ios::trunc);if(!out)throw std::runtime_error("cannot write final.dat");out<<std::fixed<<std::setprecision(7);for(int i=1;i<=nmax;++i)for(int j=1;j<=nmax;++j)out<<std::setw(13)<<theta_[j]<<' '<<std::setw(13)<<radius_[i]<<' '<<std::setw(13)<<omega_(i,j)<<' '<<std::setw(13)<<radius_[i]*std::sin(theta_[j])*u_(i,j)<<' '<<std::setw(13)<<ub_(i,j)<<'\n';}
}

Snapshot Model::snapshot()const{
  Snapshot s;s.time=time_;s.step=step_number_;s.radius.reserve(nmax);s.theta.reserve(nmax);for(int i=1;i<=nmax;++i)s.radius.push_back(radius_[i]);for(int j=1;j<=nmax;++j)s.theta.push_back(theta_[j]);
  const auto size=nmax*nmax;s.poloidal_u.reserve(size);s.poloidal_output.reserve(size);s.toroidal.reserve(size);s.omega.reserve(size);
  for(int i=1;i<=nmax;++i)for(int j=1;j<=nmax;++j){s.poloidal_u.push_back(u_(i,j));s.poloidal_output.push_back(radius_[i]*std::sin(theta_[j])*u_(i,j));s.toroidal.push_back(ub_(i,j));s.omega.push_back(omega_(i,j));}
  return s;
}

} // namespace solar_dynamo
