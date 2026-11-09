// Confronto UML SF/HSF sugli stessi eventi, con ROOT Minuit2.
// I sorgenti originali e i risultati dell'analisi non vengono modificati.
#include <TChain.h>
#include <Math/Functor.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "rgc_snapshot.h"

namespace fs = std::filesystem;
using Vec = std::array<double,7>;
using Mat = std::array<Vec,7>;
constexpr double is2 = 0.70710678118654752440;
const std::array<std::string,7> sf_names = {
    "AUL_sinPhi","AUL_sin2Phi","ALL","ALL_cosPhi","AUU_cosPhi","AUU_cos2Phi","ALU_sinPhi"};
const std::array<std::string,7> hsf_names = {
    "ReGpp_over_U0","ImGpp_over_U0","ReGmm_over_U0","ImGmm_over_U0",
    "ReG_TT_over_U0","ImG_TT_over_U0","DeltaG_T_over_U0"};

Vec to_hsf(const Vec& a) {
    return {-(a[4]+a[3])*is2, -(a[6]+a[0])*is2,
            -(a[4]-a[3])*is2, -(a[6]-a[0])*is2, -a[5], -a[1], 2*a[2]};
}
Vec to_sf(const Vec& q) {
    return {-(q[1]-q[3])*is2, -q[5], q[6]/2,
            -(q[0]-q[2])*is2, -(q[0]+q[2])*is2, -q[4], -(q[1]+q[3])*is2};
}
Mat inverse_matrix() {
    Mat b{};
    for (int j=0;j<7;++j) {
        Vec unit{}; unit[j]=1; auto a=to_sf(unit);
        for (int i=0;i<7;++i) b[i][j]=a[i];
    }
    return b;
}
Mat to_sf_covariance(const Mat& v) {
    const auto b=inverse_matrix(); Mat result{};
    for (int i=0;i<7;++i) for (int j=0;j<7;++j)
        for (int k=0;k<7;++k) for (int l=0;l<7;++l)
            result[i][j]+=b[i][k]*v[k][l]*b[j][l];
    return result;
}

struct Config {
    double pb=.8372, wpp=.934205, wmp=.935465, wpm=.998652, wmm=1;
};
Config configuration(const std::string& target,const std::string& campaign) {
    if (target=="NH3") {
        if(campaign=="fall22") return {};
        if(campaign=="sum22") return {.8384,.996267,1,.683416,.685976};
        if(campaign=="spring23") return {.8040,.545392,.545939,.99902,1};
    }
    if (target=="ND3") {
        if(campaign=="fall22") return {.8372,.998544,1,.894826,.896131};
        if(campaign=="sum22") return {.8384,.99734,1,.998677,.996021};
    }
    throw std::runtime_error("Configurazione campagna/target non supportata");
}
double weight(double lam,double target,const Config& cfg) {
    return target>0 ? (lam>0 ? cfg.wpp:cfg.wmp) : (lam>0 ? cfg.wpm:cfg.wmm);
}
struct Raw { double phi,eps,pol,lam,y; };
struct Row { Vec b; double w; };
struct Bin {
    std::vector<Row> sf;
    std::vector<Raw> reference;
    double sumw=0,x=0,q2=0,z=0,pt=0,eps=0,dilution=0,m1=0,m2=0;
};
Row prepare(const Raw& e,double dilution,const Config& cfg) {
    const double v=std::sqrt(2*e.eps*(1+e.eps));
    const double w=std::sqrt(2*e.eps*(1-e.eps));
    const double c=std::sqrt(1-e.eps*e.eps);
    const double p=dilution*e.pol, lp=e.lam*cfg.pb;
    return {{p*v*std::sin(e.phi),p*e.eps*std::sin(2*e.phi),lp*p*c,
             lp*p*w*std::cos(e.phi),v*std::cos(e.phi),e.eps*std::cos(2*e.phi),
             lp*w*std::sin(e.phi)},weight(e.lam,e.pol,cfg)};
}
Vec hsf_coefficients(const Vec& f) {
    // Coefficienti diretti HSF dell'intensita', non parametri SF da rifittare.
    return {-(f[4]+f[3])*is2,-(f[6]+f[0])*is2,
            -(f[4]-f[3])*is2,-(f[6]-f[0])*is2,-f[5],-f[1],f[2]/2};
}
std::vector<Row> prepare_hsf(const std::vector<Row>& sf) {
    std::vector<Row> h; h.reserve(sf.size());
    for(const auto& r:sf) h.push_back({hsf_coefficients(r.b),r.w});
    return h;
}
double dot(const Vec& a,const Vec& b) {
    double result=0; for(int j=0;j<7;++j) result+=a[j]*b[j]; return result;
}
struct Objective {
    const std::vector<Row>& rows;
    double sumw;
    Vec norm, bounds;
    bool hsf;
    double operator()(const double* parameters) const {
        Vec p; std::copy(parameters,parameters+7,p.begin());
        const Vec a=hsf ? to_sf(p):p;
        for(int j=0;j<7;++j)
            if(!std::isfinite(a[j]) || std::abs(a[j])>bounds[j]) return 1e12;
        const double n=1+dot(p,norm);
        if(!(n>0) || !std::isfinite(n)) return 1e12;
        double result=sumw*std::log(n);
        for(const auto& row:rows) {
            const double arg=1+dot(p,row.b);
            if(!(arg>0) || !std::isfinite(arg)) return 1e12;
            result-=row.w*std::log(arg);
        }
        return result;
    }
    Mat hessian(const Vec& p) const {
        // Derivate seconde esatte della stessa NLL, in coordinate esterne.
        std::array<std::array<long double,7>,7> accumulator{};
        const double n=1+dot(p,norm);
        for(int i=0;i<7;++i)for(int j=0;j<=i;++j)
            accumulator[i][j]=-sumw*norm[i]*norm[j]/(n*n);
        for(const auto& row:rows) {
            const double arg=1+dot(p,row.b);
            const double factor=row.w/(arg*arg);
            for(int i=0;i<7;++i)for(int j=0;j<=i;++j)
                accumulator[i][j]+=factor*row.b[i]*row.b[j];
        }
        Mat result{};
        for(int i=0;i<7;++i)for(int j=0;j<=i;++j)
            result[i][j]=result[j][i]=double(accumulator[i][j]);
        return result;
    }
};

Mat inverse_positive_definite(const Mat& matrix) {
    Mat l{},inverse{};
    for(int i=0;i<7;++i)for(int j=0;j<=i;++j) {
        double v=matrix[i][j];for(int k=0;k<j;++k)v-=l[i][k]*l[j][k];
        if(i==j) {
            if(!(v>0))throw std::runtime_error("Curvatura analitica non positiva definita");
            l[i][j]=std::sqrt(v);
        } else l[i][j]=v/l[j][j];
    }
    for(int col=0;col<7;++col) {
        Vec y{},x{};
        for(int i=0;i<7;++i) {
            double v=i==col ? 1:0;for(int k=0;k<i;++k)v-=l[i][k]*y[k];y[i]=v/l[i][i];
        }
        for(int i=6;i>=0;--i) {
            double v=y[i];for(int k=i+1;k<7;++k)v-=l[k][i]*x[k];x[i]=v/l[i][i];
        }
        for(int i=0;i<7;++i)inverse[i][col]=x[i];
    }
    return inverse;
}

struct Fit {
    Vec p{},err{}; Mat cov{};
    double f=0,edm=0; int status=-1,cov_status=-1,calls=0;
    bool ok=false,hesse=false;
};
Fit fit(const Objective& objective) {
    ROOT::Minuit2::Minuit2Minimizer m(ROOT::Minuit2::kMigrad);
    ROOT::Math::Functor f([&](const double* p){return objective(p);},7);
    m.SetFunction(f); m.SetPrintLevel(-1); m.SetMaxFunctionCalls(50000);
    m.SetMaxIterations(10000); m.SetTolerance(.001); m.SetStrategy(2); m.SetErrorDef(.5);
    Vec bounds=objective.bounds;
    if(objective.hsf) {
        const double r=(bounds[4]+bounds[3])*is2, i=(bounds[6]+bounds[0])*is2;
        bounds={r,i,r,i,bounds[5],bounds[1],2*bounds[2]};
    }
    for(int j=0;j<7;++j) {
        const auto& name=objective.hsf ? hsf_names[j]:sf_names[j];
        m.SetLimitedVariable(j,name,0,j==6&&objective.hsf ? .002:.001,-bounds[j],bounds[j]);
    }
    Fit result; result.ok=m.Minimize(); result.status=m.Status();
    result.hesse=result.ok && m.Hesse(); result.cov_status=m.CovMatrixStatus();
    result.f=m.MinValue(); result.edm=m.Edm(); result.calls=m.NCalls();
    for(int i=0;i<7;++i) {
        result.p[i]=m.X()[i]; result.err[i]=m.Errors()[i];
        for(int j=0;j<7;++j) result.cov[i][j]=m.CovMatrix(i,j);
    }
    return result;
}

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result; std::stringstream s(line); std::string field;
    while(std::getline(s,field,',')) result.push_back(field);
    return result;
}
void load_aux(const fs::path& project,const std::string& campaign,std::vector<Bin>& bins,
              const std::string& kind) {
    std::ifstream mc(project/"RECO_CSV"/("table_RGC_MC_"+campaign+"_AUU_"+kind+".csv"));
    if(!mc) throw std::runtime_error("Tabella momenti MC assente: "+kind);
    std::string line; std::getline(mc,line); std::vector<bool> found(bins.size(),false);
    while(std::getline(mc,line)) {
        if(line.find_first_not_of(" \r\t")==std::string::npos) continue;
        auto c=split(line); if(c.size()!=6) throw std::runtime_error("Formato momenti MC");
        int i=std::stoi(c[0])-1;
        if(i<0 || i>=int(bins.size()) || found[i]) throw std::runtime_error("Bin MC errato");
        bins[i].m1=std::stod(c[4]);bins[i].m2=std::stod(c[5]);found[i]=true;
    }
    if(std::find(found.begin(),found.end(),false)!=found.end()) throw std::runtime_error("Bin MC mancante");
    std::ifstream df(project/"dilution_factor"/("DF_Kplus_"+campaign+".csv"));
    if(!df) throw std::runtime_error("Tabella diluizione assente");
    std::getline(df,line);std::fill(found.begin(),found.end(),false);
    const std::string category=kind=="xQ2" ? "xB_Q2":"z_Pt";
    while(std::getline(df,line)) {
        auto c=split(line);if(c.size()<9 || c[1]!=category) continue;
        int i=std::stoi(c[2])-1;if(i<0 || i>=int(bins.size())) continue;
        bins[i].dilution=std::stod(c[7]);found[i]=true;
        if(!(bins[i].dilution>0)) throw std::runtime_error("Diluizione non valida");
    }
    if(std::find(found.begin(),found.end(),false)!=found.end()) throw std::runtime_error("Bin DF mancante");
}

void add(Bin& bin,const Raw& e,double xb,double q2,double z,double pt,const Config& cfg) {
    auto row=prepare(e,bin.dilution,cfg);bin.sf.push_back(row);bin.sumw+=row.w;
    bin.x+=xb;bin.q2+=q2;bin.z+=z;bin.pt+=pt;bin.eps+=e.eps;
    if(bin.reference.size()<64) bin.reference.push_back(e);
}
double check_snapshot(const Bin& bin,const Config& cfg,const std::string& target,
                      const std::string& campaign) {
    std::vector<double> phi,spin,y,eps,pol,lam,st;std::vector<Row> rows;double sumw=0;
    for(const auto& e:bin.reference) {
        phi.push_back(e.phi);spin.push_back(e.pol>0 ? 1:-1);y.push_back(e.y);
        eps.push_back(e.eps);pol.push_back(e.pol);lam.push_back(e.lam);st.push_back(0);
        rows.push_back(prepare(e,bin.dilution,cfg));sumw+=rows.back().w;
    }
    Vec norm{};norm[4]=bin.m1;norm[5]=bin.m2;Vec bounds;bounds.fill(.9);
    Objective fast{rows,sumw,norm,bounds,false};
    std::mt19937 random(39017);std::uniform_real_distribution<double> draw(-.08,.08);
    double maximum=0;
    for(int k=0;k<12;++k) {
        Vec a;for(auto& v:a) v=draw(random);
        const double original=rgc_snapshot::AUL_loglike_withLL(a.data(),phi,spin,y,eps,pol,
            lam,st,bin.dilution,bin.m1,bin.m2,true,target.c_str(),campaign.c_str());
        maximum=std::max(maximum,std::abs(original-fast(a.data())));
    }
    return maximum;
}

void self_test() {
    std::mt19937 random(38129); std::uniform_real_distribution<double> draw(-.08,.08);
    Config cfg;Bin b;b.dilution=.24;b.m1=.21;b.m2=-.07;
    for(int i=0;i<800;++i) {
        Raw e{-3.14+6.28*(i+.5)/800,.15+.7*(i%17)/17.,i%3 ? .7:-.65,i%2 ? 1.:-1.,.45};
        add(b,e,.2,2.,.4,.3,cfg);
    }
    auto h=prepare_hsf(b.sf);Vec norm{};norm[4]=b.m1;norm[5]=b.m2;Vec bounds;bounds.fill(.8);
    Objective sf{b.sf,b.sumw,norm,bounds,false};
    Objective hsf{h,b.sumw,hsf_coefficients(norm),bounds,true};
    for(int i=0;i<100;++i) {
        Vec a;for(auto& v:a)v=draw(random);auto q=to_hsf(a);auto back=to_sf(q);
        for(int j=0;j<7;++j)if(std::abs(a[j]-back[j])>1e-14)throw std::runtime_error("Inversione");
        if(std::abs(sf(a.data())-hsf(q.data()))>1e-10)throw std::runtime_error("Chiusura NLL");
    }
    Vec outside{};outside[0]=.9;auto qo=to_hsf(outside);
    if(sf(outside.data())!=hsf(qo.data()))throw std::runtime_error("Dominio diverso");
    if(check_snapshot(b,cfg,"NH3","fall22")>1e-10)throw std::runtime_error("Snapshot diverso");
    Mat covariance{};for(int j=0;j<7;++j)covariance[j][j]=1;
    if(std::abs(to_sf_covariance(covariance)[2][2]-.25)>1e-14)throw std::runtime_error("Covarianza ALL");
    Vec a{.01,.02,.03,.04,.05,.06,.07};auto q=to_hsf(a);
    auto cs=inverse_positive_definite(sf.hessian(a));
    auto ch=to_sf_covariance(inverse_positive_definite(hsf.hessian(q)));
    for(int i=0;i<7;++i)for(int j=0;j<7;++j)
        if(std::abs(cs[i][j]-ch[i][j])>1e-11)throw std::runtime_error("Curvatura analitica diversa");
    std::cout<<"Self-test: inversione, 100 chiusure NLL, dominio, snapshot, covarianza e Hessiana analitica OK\n";
}

void curvature_check(std::vector<Bin>& bins,const std::string& kind,const fs::path& out) {
    std::vector<Vec> sf_values(bins.size()),hsf_as_sf(bins.size()),sf_errors(bins.size()),hsf_errors(bins.size());
    std::vector<std::array<bool,7>> found(bins.size());
    std::ifstream input(out/("amplitudes_comparison_"+kind+".csv"));
    if(!input)throw std::runtime_error("Eseguire prima i due fit");
    std::string line;std::getline(input,line);
    while(std::getline(input,line)) {
        auto c=split(line);int k=std::stoi(c.at(0))-1;
        auto at=std::find(sf_names.begin(),sf_names.end(),c.at(7));
        if(k<0 || k>=int(bins.size()) || at==sf_names.end())throw std::runtime_error("Tabella fit errata");
        int j=int(at-sf_names.begin());sf_values[k][j]=std::stod(c.at(8));
        sf_errors[k][j]=std::stod(c.at(9));hsf_as_sf[k][j]=std::stod(c.at(10));
        hsf_errors[k][j]=std::stod(c.at(11));found[k][j]=true;
    }
    std::ofstream result(out/("analytic_curvature_"+kind+".csv"));
    result<<std::setprecision(17)<<"bin,parameter,sf_analytic_error,hsf_analytic_error,analytic_error_ratio,"
        "sf_minuit_reported_error,hsf_minuit_transformed_error,sf_minuit_over_analytic,hsf_minuit_over_analytic\n";
    for(int k=0;k<int(bins.size());++k) {
        if(std::find(found[k].begin(),found[k].end(),false)!=found[k].end())throw std::runtime_error("Fit incompleto");
        const auto& b=bins[k];auto h=prepare_hsf(b.sf);Vec norm{};norm[4]=b.m1;norm[5]=b.m2;
        Vec bounds;bounds.fill(.9);
        Objective sf{b.sf,b.sumw,norm,bounds,false};
        Objective hsf{h,b.sumw,hsf_coefficients(norm),bounds,true};
        auto vs=inverse_positive_definite(sf.hessian(sf_values[k]));
        auto vh=to_sf_covariance(inverse_positive_definite(hsf.hessian(to_hsf(hsf_as_sf[k]))));
        double maximum=0;
        for(int j=0;j<7;++j) {
            double es=std::sqrt(vs[j][j]),eh=std::sqrt(vh[j][j]);
            maximum=std::max(maximum,std::abs(eh/es-1));
            result<<k+1<<','<<sf_names[j]<<','<<es<<','<<eh<<','<<eh/es<<','
                  <<sf_errors[k][j]<<','<<hsf_errors[k][j]<<','<<sf_errors[k][j]/es<<','<<hsf_errors[k][j]/eh<<'\n';
        }
        std::cout<<kind<<" bin "<<k+1<<" max differenza relativa errori analitici="<<maximum<<std::endl;
    }
}

void compare(std::vector<Bin>& bins,const std::string& kind,const fs::path& out,
             const Config& cfg,const std::string& target,const std::string& campaign) {
    std::ofstream summary(out/("fit_comparison_"+kind+".csv"));
    std::ofstream amplitudes(out/("amplitudes_comparison_"+kind+".csv"));
    std::ofstream covariance(out/("covariances_"+kind+".csv"));
    for(auto* f:{&summary,&amplitudes,&covariance}) {
        if(!*f) throw std::runtime_error("Impossibile creare output");*f<<std::setprecision(17);
    }
    summary<<"bin,entries,mean_xB,mean_Q2,mean_z,mean_PhT,mean_epsilon,nll_sf,nll_hsf,delta_nll,"
             "sf_ok,hsf_ok,sf_status,hsf_status,sf_cov_status,hsf_cov_status,sf_hesse,hsf_hesse,"
             "sf_edm,hsf_edm,sf_calls,hsf_calls,max_abs_delta_over_sigma_sf,max_relative_error_difference,"
             "max_nll_identity_difference,snapshot_identity_difference,boundary_distance_sf,boundary_distance_hsf\n";
    amplitudes<<"bin,entries,mean_xB,mean_Q2,mean_z,mean_PhT,mean_epsilon,parameter,sf,sf_err,"
                "hsf_reconstructed,hsf_reconstructed_err,delta,delta_over_sigma_sf,error_ratio\n";
    covariance<<"bin,basis,row,column,value\n";
    for(int k=0;k<int(bins.size());++k) {
        auto& bin=bins[k];if(bin.sf.empty())throw std::runtime_error("Bin vuoto");
        auto h=prepare_hsf(bin.sf);
        Vec norm{};norm[4]=bin.m1;norm[5]=bin.m2;Vec bounds;bounds.fill(.8);
        if(kind=="xQ2")bounds[4]=bounds[5]=.9;
        Objective sf{bin.sf,bin.sumw,norm,bounds,false};
        Objective hsf{h,bin.sumw,hsf_coefficients(norm),bounds,true};
        double identity=0;
        std::mt19937 random(731+k);std::uniform_real_distribution<double> draw(-.06,.06);
        for(int t=0;t<3;++t) {
            Vec a;for(auto& v:a)v=draw(random);const auto q=to_hsf(a);
            identity=std::max(identity,std::abs(sf(a.data())-hsf(q.data())));
        }
        const double snapshot=check_snapshot(bin,cfg,target,campaign);
        if(identity>1e-6 || snapshot>1e-10)throw std::runtime_error("Mancata chiusura della funzione");
        const auto a=fit(sf);const auto q=fit(hsf);const auto ah=to_sf(q.p);const auto va=to_sf_covariance(q.cov);
        const double n=double(bin.sf.size());
        double maxdelta=0,maxerror=0,bd_sf=1e9,bd_hsf=1e9;
        for(int j=0;j<7;++j) {
            const double err=std::sqrt(std::max(0.,va[j][j]));
            const double delta=ah[j]-a.p[j];
            maxdelta=std::max(maxdelta,std::abs(delta)/std::max(a.err[j],1e-20));
            maxerror=std::max(maxerror,std::abs(err/std::max(a.err[j],1e-20)-1));
            bd_sf=std::min(bd_sf,bounds[j]-std::abs(a.p[j]));
            bd_hsf=std::min(bd_hsf,bounds[j]-std::abs(ah[j]));
            amplitudes<<k+1<<','<<bin.sf.size()<<','<<bin.x/n<<','<<bin.q2/n<<','<<bin.z/n<<','<<bin.pt/n<<','
                <<bin.eps/n<<','<<sf_names[j]<<','<<a.p[j]<<','<<a.err[j]<<','<<ah[j]<<','<<err<<','<<delta<<','
                <<delta/std::max(a.err[j],1e-20)<<','<<err/std::max(a.err[j],1e-20)<<'\n';
            for(int l=0;l<7;++l) {
                covariance<<k+1<<",SF,"<<sf_names[j]<<','<<sf_names[l]<<','<<a.cov[j][l]<<'\n';
                covariance<<k+1<<",HSF,"<<hsf_names[j]<<','<<hsf_names[l]<<','<<q.cov[j][l]<<'\n';
                covariance<<k+1<<",HSF_to_SF,"<<sf_names[j]<<','<<sf_names[l]<<','<<va[j][l]<<'\n';
            }
        }
        summary<<k+1<<','<<bin.sf.size()<<','<<bin.x/n<<','<<bin.q2/n<<','<<bin.z/n<<','<<bin.pt/n<<','<<bin.eps/n<<','
            <<a.f<<','<<q.f<<','<<q.f-a.f<<','<<a.ok<<','<<q.ok<<','<<a.status<<','<<q.status<<','
            <<a.cov_status<<','<<q.cov_status<<','<<a.hesse<<','<<q.hesse<<','<<a.edm<<','<<q.edm<<','
            <<a.calls<<','<<q.calls<<','<<maxdelta<<','<<maxerror<<','<<identity<<','<<snapshot<<','<<bd_sf<<','<<bd_hsf<<'\n';
        summary.flush();amplitudes.flush();covariance.flush();
        std::cout<<kind<<" bin "<<k+1<<" N="<<bin.sf.size()<<" status="<<a.status<<'/'<<q.status
            <<" cov="<<a.cov_status<<'/'<<q.cov_status<<" deltaNLL="<<q.f-a.f
            <<" max|delta/sigmaSF|="<<maxdelta<<std::endl;
    }
}

int main(int argc,char** argv) {
    try {
        if(argc==2 && std::string(argv[1])=="--self-test") {self_test();return 0;}
        fs::path project=".",out="sf_hsf_fit_results";std::string campaign="fall22",target="NH3";
        long long maximum=-1;bool curvature_only=false;
        for(int i=1;i<argc;++i) {
            std::string arg=argv[i];
            if(arg=="--curvature-only") {curvature_only=true;continue;}
            if(i+1>=argc)throw std::runtime_error("Argomento senza valore");
            std::string value=argv[++i];
            if(arg=="--project")project=value;else if(arg=="--outdir")out=value;
            else if(arg=="--campaign")campaign=value;else if(arg=="--target")target=value;
            else if(arg=="--max-entries")maximum=std::stoll(value);
            else throw std::runtime_error("Argomento sconosciuto: "+arg);
        }
        if(target!="NH3")throw std::runtime_error("Questa selezione usa DF_Kplus NH3; usare --target NH3");
        const auto cfg=configuration(target,campaign);
        std::vector<Bin> xq2(14),zpt(25);load_aux(project,campaign,xq2,"xQ2");load_aux(project,campaign,zpt,"zPt");
        TChain chain("Kaon+");std::vector<fs::path> files;
        for(const auto& entry:fs::directory_iterator(project/(campaign+"_"+target+"_data")))
            if(entry.path().extension()==".root")files.push_back(entry.path());
        std::sort(files.begin(),files.end());
        for(const auto& file:files)if(chain.Add(file.c_str())==0)throw std::runtime_error("File ROOT non letto");
        if(files.empty())throw std::runtime_error("Nessun file ROOT");
        chain.SetBranchStatus("*",0);
        double phi,eps,pol,lam,y,xb,q2,z,pt,chi2,pid,rl,rq,ntot,rpid,rchi,mom,labphi;
        auto bind=[&](const char* name,double* value) {
            chain.SetBranchStatus(name,1);if(chain.SetBranchAddress(name,value)<0)
                throw std::runtime_error(std::string("Branch non disponibile: ")+name);
        };
        bind("kaon_phi_h",&phi);bind("epsilon",&eps);bind("Polarization",&pol);bind("Helicity",&lam);
        bind("y",&y);bind("xB",&xb);bind("Q2",&q2);bind("z",&z);bind("kaon_Pt",&pt);
        bind("clas_chi2",&chi2);bind("clas_pid",&pid);bind("Rich_RL",&rl);bind("Rich_RQ",&rq);
        bind("Rich_nTot",&ntot);bind("Rich_PID",&rpid);bind("Rich_chi2",&rchi);
        bind("kaon_mom",&mom);bind("kaon_phi_lab",&labphi);
        const auto total=chain.GetEntries();const auto entries=maximum<0 ? total:std::min(total,maximum);
        std::cout<<"Lettura "<<entries<<"/"<<total<<" eventi da "<<files.size()<<" file ROOT"<<std::endl;
        long long accepted=0;
        for(long long i=0;i<entries;++i) {
            if(chain.GetEntry(i)<=0)throw std::runtime_error("Errore lettura evento");
            if(pid==-11 || y<=.2 || std::abs(pol)<.1 || lam==0 || !(std::abs(chi2)<3))continue;
            if(!rgc_snapshot::passRichSelection(rl,rq,int(ntot),int(rpid),int(pid),rchi,mom,labphi))continue;
            if(!std::isfinite(phi)||!std::isfinite(eps)||eps<0||eps>=1||!std::isfinite(pol))
                throw std::runtime_error("Evento selezionato non valido per la likelihood");
            Raw event{phi,eps,pol,lam,y};
            const int ix=rgc_snapshot::getBinIndex_xQ2(xb,q2),iz=rgc_snapshot::getBinIndex_zPt(z,pt);
            if(ix>0)add(xq2.at(ix-1),event,xb,q2,z,pt,cfg);
            if(iz>0)add(zpt.at(iz-1),event,xb,q2,z,pt,cfg);
            if(ix>0||iz>0)++accepted;
        }
        fs::create_directories(out);
        if(curvature_only) {
            curvature_check(xq2,"xQ2",out);curvature_check(zpt,"zPt",out);return 0;
        }
        std::ofstream manifest(out/"run_metadata.txt");
        manifest<<"project="<<fs::absolute(project)<<"\ncampaign="<<campaign<<"\ntarget="<<target
                <<"\nroot_files="<<files.size()<<"\ninput_entries="<<total<<"\nprocessed_entries="<<entries
                <<"\nselected_union_entries="<<accepted<<"\nfull_sample="<<(entries==total)
                <<"\nstart=zero_in_both_bases\nerror_def=0.5\nstrategy=2\ntolerance=0.001"
                <<"\nSF_bounds=0.8_except_xQ2_UU_0.9\nHSF_domain=exact_inverse_image_of_SF_box"
                <<"\ncovariance=Minuit_Hesse_curvature_not_weight_corrected\n";
        std::cout<<"Selezione completata: "<<accepted<<" eventi nell'unione delle due griglie"<<std::endl;
        compare(xq2,"xQ2",out,cfg,target,campaign);compare(zpt,"zPt",out,cfg,target,campaign);
        curvature_check(xq2,"xQ2",out);curvature_check(zpt,"zPt",out);
        return 0;
    } catch(const std::exception& error) {std::cerr<<"Errore: "<<error.what()<<'\n';return 1;}
}
