#pragma once
#include <cstddef>
#include <filesystem>
#include <vector>

namespace solar_dynamo {
struct Parameters {
  double tmax = 10.0, v0 = -20.0, et0 = 3.0, et1 = 0.04, al0 = 30.0;
  double omega0 = 2.0 * 43.28 * 3.14159265358979323846, dt = 0.00005;
  bool relaxed_initial_state = true;
};
struct Snapshot {
  double time{}; std::size_t step{};
  std::vector<double> radius, theta, poloidal_u, poloidal_output, toroidal, omega;
};
class Model {
 public:
  static constexpr int nmax=257, n=256, lmax=96, halfmax=514;
  explicit Model(Parameters parameters = {});
  void initialize(const std::filesystem::path& init_file="init.dat");
  void step(std::size_t count=1);
  void run(const std::filesystem::path& output_directory=".");
  void step_with_output(std::size_t count, const std::filesystem::path& output_directory=".");
  void write_final_outputs(const std::filesystem::path& output_directory=".") const;
  Snapshot snapshot() const;
  double time() const noexcept{return time_;}
  std::size_t step_number() const noexcept{return step_number_;}
  bool initialized() const noexcept{return initialized_;}
  const Parameters& parameters() const noexcept{return parameters_;}
 private:
  class Grid {
   public:
    explicit Grid(int extent=nmax+1):extent_(extent),values_(extent*extent,0.0){}
    double& operator()(int i,int j){return values_[i*extent_+j];}
    double operator()(int i,int j)const{return values_[i*extent_+j];}
   private:int extent_;std::vector<double> values_;
  };
  void build_profiles(); void load_initial_state(const std::filesystem::path&);
  void build_coefficients(); void advance_one(); void apply_boundaries();
  void write_omega_latitude(const std::filesystem::path&)const;
  void write_diffrot_snapshot(const std::filesystem::path&)const;
  static void tridag(const std::vector<double>&,const std::vector<double>&,
                     const std::vector<double>&,const std::vector<double>&,
                     std::vector<double>&,int);
  double simpson_ss(int ell)const;
  double upper_boundary_coefficient(const std::vector<double>&,int)const;
  static double associated_legendre(int,int,double); static double legacy_erf(double);
  static void ensure_directory(const std::filesystem::path&);
  Parameters parameters_; double pi_,pm_,pb_,pk_,pw_,qm_,dp_,dq_,fac_;
  double time_=0.0;std::size_t step_number_=0;bool initialized_=false;
  Grid u_,ub_,omega_,alpha_,dom_,eta_,etab_,rho_,nu_d_,wr_,wt_;
  Grid dror_,drot_,lorentz_,phi_,phib_,phiw_;
  Grid ma_,mb_,mc_,md_,me_,mf_,ta_,tb_,tc_,td_,te_,tf_,oa_,ob_,oc_,od_,oe_,of_;
  Grid vp_,vq_,vp1_,vq1_,psi_,etab2_,dvp_,vpb_,deta_;
  Grid rho2_,nut2_,ft1_,gt1_,ft2_,gt2_,fr1_,gr1_,fr2_,gr2_,v00_,h1_;
  std::vector<double> radius_,theta_,ra_,smooth_,ss1_p_,sn_,pl_;
};
} // namespace solar_dynamo
