//Written by John G Baker NASA-GSFC (2014-16)
//This is jointly based on the gleam.cc ptmcmc-based gravitational-lens code
//and the bambi-based LISAinterface.c 

#include <valarray>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include "omp.h"
#include "options.hh"
#include "bayesian.hh"
#include "proposal_distribution.hh"
#include "ptmcmc.hh"
#include "LISAinference_common.h"

using namespace std;

typedef initializer_list<double> dlist;
typedef initializer_list<int> ilist;

shared_ptr<Random> globalRNG;//used for some debugging... 

//Global Control parameters (see more in main())
const int Npar=9;
int output_precision;

//First we define the likelihood function object
class flare_likelihood : public bayes_likelihood {
  sampleable_probability_function * prior;
  int count,nevery;
  double total_eval_time;
  double best_post;
  state best;
  void *context;
public:
  flare_likelihood(stateSpace *sp, void *context, sampleable_probability_function *prior=nullptr):context(context),prior(prior),bayes_likelihood(sp,nullptr,nullptr){
    best=state(sp,sp->size());
    reset();
    nevery=0;
  };
  void info_every(int n){nevery=n;};
  void reset(){
    best_post=-INFINITY;
    best=best.scalar_mult(0);
    count=0;
    total_eval_time=0;
  }
  state bestState(){return best;};
  double bestPost(){return best_post;};
  state transformDataState(const state &s)const{return s;};
  state transformSignalState(const state &s)const{return s;};
  double evaluate_log(state &s){
    valarray<double>params=s.get_params();
    double tstart=omp_get_wtime();
    
    //Here the param indices are hard coded, but we can use sp to find them by name
    LISAParams templateparams;
    templateparams.m1 = params[0];
    templateparams.m2 = params[1];
    templateparams.tRef = params[2];
    templateparams.distance = params[3];
    if(!priorParams->flat_distprior)//If not using flat prior on distance then we need to transform from s(D) to D.
      templateparams.distance= pow(params[3] * pow(priorParams->dist_max, 3) + (1.0 - params[3]) * pow(priorParams->dist_min, 3), 1.0 / 3.0);
    templateparams.phiRef = params[4];
    templateparams.inclination = params[5];
    templateparams.lambda = params[6];
    templateparams.beta = params[7];
    templateparams.polarization = params[8];
    templateparams.nbmode = globalparams->nbmodetemp; /* Using the global parameter for the number of modes in templates */
    
    double result=0;
    /* Note: context points to a LISAContext structure containing a LISASignal* */
    if(globalparams->tagint==0) {
      //LISAInjectionCAmpPhase* injection = ((LISAInjectionCAmpPhase*) context);
      LISAInjectionCAmpPhase injection = *((LISAInjectionCAmpPhase*) context);
      double like;
      like = CalculateLogLCAmpPhase(&templateparams, &injection);
      result = like - logZdata;
      //result = CalculateLogLCAmpPhase(&templateparams, injection) - logZdata;
    }
    else if(globalparams->tagint==1) {
      LISAInjectionReIm* injection = ((LISAInjectionReIm*) context);
      result = CalculateLogLReIm(&templateparams, injection) - logZdata;
    }
    
    double post=result;
    double lpriorval=0;
    if(prior)lpriorval=prior->evaluate_log(s);
    post+=lpriorval;
    double tend=omp_get_wtime();
    double eval_time = tend-tstart;
#pragma omp critical
    {     
      total_eval_time+=eval_time;
      count++;
      if(nevery>0&&0==count%nevery){
	cout<<"eval_time = "<<eval_time<<"  result="<<result<<" prior="<<lpriorval<<" logZdata="<<logZdata<<" post="<<post;
	print_info();
      }
      if(post>best_post){
	best_post=post;
	best=state(s);
      }
      //cout<<"loglike="<<result<<"<="<<maxLike<<endl;   
      if(!isfinite(result)){
	cout<<"Whoa dude, loglike is NAN! What's up with that?"<<endl;
	cout<<"params="<<s.get_string()<<endl;
	result=-INFINITY;
      }
    }
    return result;
  };

  ///from stateSpaceInterface
  virtual void defWorkingStateSpace(const stateSpace &sp){};

  virtual void write(ostream &out,state &st){cout<<"flare_likelihood::write: No write routine defined!"<<endl;};
  virtual void writeFine(ostream &out,state &st,int ns=-1, double ts=0, double te=0){cout<<"flare_likelihood::writeFine: No write routine defined!"<<endl;};
  virtual void getFineGrid(int & nfine, double &tstart, double &tend)const{cout<<"flare_likelihood::getFineGrid: No routine defined!"<<endl;};
  void print_info(){cout<<" mean = "<<total_eval_time/count<<" over "<<count<<" total evals"<<endl; };

};


//***************************************************************************************8
//main test program
int main(int argc, char*argv[]){
  //string datafile;
  const int NparRead=Npar; 

  //Create the sampler
  ptmcmc_sampler mcmc;
  bayes_sampler *s0=&mcmc;

  Options opt;

  s0->addOptions(opt,"");
  opt.add(Option("nchains","Number of consequtive chain runs. Default 1","1"));
  opt.add(Option("info_every","How often to dump likehood eval info to stdout. Default never","0"));
  opt.add(Option("rng_seed","Pseudo random number grenerator seed in [0,1). (Default=-1, use clock to seed.)","-1"));
  opt.add(Option("precision","Set output precision digits. (Default 13).","13"));
  opt.add(Option("help","Print help message."));
  //First we parse the ptmcmc-related parameters like un gleam. 
  opt.parse(argc,argv,false);
  if(opt.set("help")){
    cout<<"****** PTMCMC options *******\n"<<opt.print_usage()<<endl;
    static char help[]="--help";
    argv[argc]=help;
    argc++;
  } else cout<<"ptmcmc flags=\n"<<opt.report()<<endl;

  
  //Next we perform the initializations from LISAinference.c
  LISARunParams runParams={};
  int ndim=0,nPar=0;
  int *freeparamsmap = NULL;
  void *context = NULL;
  double logZtrue=0;
  addendum(argc,argv,&runParams,&ndim,&nPar,&freeparamsmap,&context,&logZtrue);

  //Post parse setup
  double seed;
  int Nchain,save_every,info_every;
  int Nsigma=1;
  int Nbest=10;
  istringstream(opt.value("nchains"))>>Nchain;
  istringstream(opt.value("rng_seed"))>>seed;
  istringstream(opt.value("info_every"))>>info_every;
  
  //if seed<0 set seed from clock
  if(seed<0)seed=fmod(time(NULL)/3.0e7,1);
  istringstream(opt.value("precision"))>>output_precision;
  ProbabilityDist::setSeed(seed);
  globalRNG.reset(ProbabilityDist::getPRNG());//just for safety to keep us from deleting main RNG in debugging.

  //output location
  string outname(runParams.outroot);
  ostringstream ss("");

  //report
  cout.precision(output_precision);
  cout<<"\noutname = '"<<outname<<"'"<<endl;
  cout<<"seed="<<seed<<endl; 
  cout<<"Running on "<<omp_get_max_threads()<<" thread"<<(omp_get_max_threads()>1?"s":"")<<"."<<endl;

  //Set up the parameter space
  stateSpace space(Npar);
  /* Note: here we output physical values in the cube (overwriting), and we keep the original order for physical parameters */ 
  string names[]={"m1","m2","tRef","dist","phase","inc","lambda","beta","pol"};
  if(!priorParams->flat_distprior)names[3]="s(D)";
  space.set_names(names);  
  space.set_bound(4,boundary(boundary::wrap,boundary::wrap,0,2*M_PI));//set 2-pi-wrapped space for phi.
  space.set_bound(6,boundary(boundary::wrap,boundary::wrap,0,2*M_PI));//set 2-pi-wrapped space for phi.
  space.set_bound(8,boundary(boundary::wrap,boundary::wrap,0,M_PI));//set pi-wrapped space for phi.
  cout<<"Parameter space:\n"<<space.show()<<endl;

  //Set the prior:
  //Eventually this should move to the relevant constitutent code elements where the params are given meaning.
  const int uni=mixed_dist_product::uniform, gauss=mixed_dist_product::gaussian, polar=mixed_dist_product::polar, copol=mixed_dist_product::copolar; 
  double lambda0=(priorParams->lambda_min+priorParams->lambda_max)/2,dlambda=(priorParams->lambda_max-priorParams->lambda_min)/2;
  double beta0=(priorParams->beta_min+priorParams->beta_max)/2,dbeta=(priorParams->beta_max-priorParams->beta_min)/2;//Note: these are ignored by co-polar prior
  double tref0=injectedparams->tRef,dtref=priorParams->deltaT;
  double phase0=(priorParams->phase_min+priorParams->phase_max)/2,dphase=(priorParams->phase_max-priorParams->phase_min)/2;
  double inc0=(priorParams->inc_min+priorParams->inc_max)/2,dinc=(priorParams->inc_max-priorParams->inc_min)/2;//Note: these are ignored by polar prior
  double dist0=(priorParams->dist_min+priorParams->dist_max)/2,ddist=(priorParams->dist_max-priorParams->dist_min)/2;
  double m10=(priorParams->comp_max+priorParams->comp_min)/2,m20=m10,dm1=(priorParams->comp_max-priorParams->comp_min)/2,dm2=dm1;
  double pol0=(priorParams->pol_min+priorParams->pol_max)/2,dpol=(priorParams->pol_max-priorParams->pol_min)/2;
  //If the parameters are "fixed" then we approximately realize that by restricting the range by 1e-5.
  if (!isnan(priorParams->fix_lambda))cout<<" ** Parameter fixing not supported.  Ignoring request! **"<<endl;
  if (!isnan(priorParams->fix_beta))cout<<" ** Parameter fixing not supported.  Ignoring request! **"<<endl;
  if (!isnan(priorParams->fix_time))cout<<" ** Parameter fixing not supported.  Ignoring request! **"<<endl;
  if (!isnan(priorParams->fix_phase))cout<<" ** Parameter fixing not supported.  Ignoring request! **"<<endl;
  if (!isnan(priorParams->fix_inc))cout<<" ** Parameter fixing not supported.  Ignoring request! **"<<endl;
  if (!isnan(priorParams->fix_pol))cout<<" ** Parameter fixing not supported.  Ignoring request! **"<<endl;
  if (!isnan(priorParams->fix_dist))cout<<" ** Parameter fixing not supported.  Ignoring request! **"<<endl;
  if (!isnan(priorParams->fix_m1))cout<<" ** Parameter fixing not supported.  Ignoring request! **"<<endl;
  if (!isnan(priorParams->fix_m2))cout<<" ** Parameter fixing not supported.  Ignoring request! **"<<endl;
  if(priorParams->logflat_massprior){
    dm1=sqrt((m10+dm1)/(m10-dm1));
    dm2=sqrt((m20+dm2)/(m20-dm2));
    m10=priorParams->comp_min*dm1;
    m20=priorParams->comp_min*dm2;
  }
  if(!priorParams->flat_distprior) {
    //In this case we transform the distance parameter to s(D)=(D^3-Dmin^3)/(Dmax^3-Dmin^3)
    double D0,dD;
    double Dmin3=pow(priorParams->dist_min,3);
    double denom=pow(priorParams->dist_max,3)-Dmin3;
    double smin= (pow(dist0-ddist,3) - Dmin3)/denom;
    double smax= (pow(dist0+ddist,3) - Dmin3)/denom;
    dist0=(smax+smin)/2;
    ddist=(smax-smin)/2;
  }
  //                                     m1      m2     tRef     dist   phase    inc    lambda    beta     pol  
  valarray<double>    centers((dlist){  m10,    m20,   tref0,  dist0,  phase0,  inc0, lambda0,  beta0,   pol0  });
  valarray<double> halfwidths((dlist){  dm1,    dm2,   dtref,  ddist,  dphase,  dinc, dlambda,  dbeta,   dpol  });
  valarray<int>         types((ilist){  uni,    uni,     uni,     uni,    uni, polar,     uni,  copol,    uni  });
  if (!isnan(priorParams->fix_beta))types[7]=uni;//"fixed" parameter set to narrow uniform range
  if (!isnan(priorParams->fix_inc)) types[5]=uni;//"fixed" parameter set to narrow uniform range
  if(priorParams->logflat_massprior)types[0]=types[1]=mixed_dist_product::log;
  sampleable_probability_function *prior;  
  prior=new mixed_dist_product(&space,types,centers,halfwidths);
  cout<<"Prior is:\n"<<prior->show()<<endl;
  
  //test the prior:
  //for(int i=0;i<5;i++){state s=prior->drawSample(*globalRNG);cout<<"test state "<<i<<"="<<s.show()<<endl;}
  
  //Set the likelihood
  bayes_likelihood *like=nullptr;
  flare_likelihood fl(&space, context, prior);
  fl.addOptions(opt,"");
  if(info_every>0)fl.info_every(info_every);
  like = &fl;
  like->defWorkingStateSpace(space);

  //Test the injection value:
  fl.info_every(1);
  valarray<double> inj_pars((dlist){injectedparams->m1,injectedparams->m2,injectedparams->tRef,injectedparams->distance,injectedparams->phiRef,injectedparams->inclination,injectedparams->lambda,injectedparams->beta,injectedparams->polarization});
  if(!priorParams->flat_distprior) {
    //In this case we transform the distance parameter to s(D)=(D^3-Dmin^3)/(Dmax^3-Dmin^3)
    double D0,dD;
    double Dmin3=pow(priorParams->dist_min,3);
    double denom=pow(priorParams->dist_max,3)-Dmin3;
    inj_pars[3] = (pow(injectedparams->distance,3) - Dmin3)/denom;
  }
  state injected_state(&space,inj_pars);
  cout<<"injected state="<<injected_state.show()<<endl;
  double injloglike=like->evaluate_log(injected_state);
  cout<<"injected log-like "<<injloglike<<" post="<<like->bestPost()<<endl;
  //restore likelihood states
  like->reset();
  if(info_every>0)fl.info_every(info_every);
  
  //assuming mcmc:
  //Set the proposal distribution 
  int Ninit;
  proposal_distribution *prop=ptmcmc_sampler::new_proposal_distribution(Npar,Ninit,opt,prior,&halfwidths);
  cout<<"Proposal distribution is:\n"<<prop->show()<<endl;
  //set up the mcmc sampler (assuming mcmc)
  mcmc.setup(Ninit,*like,*prior,*prop,output_precision);


  //Prepare for chain output
  //ss<<"gle_"<<outname;
  ss<<outname;
  string base=ss.str();
  
  //Loop over Nchains
  for(int ic=0;ic<Nchain;ic++){
    bayes_sampler *s=s0->clone();
    s->initialize();
    s->run(base,ic);
    s->analyze(base,ic,Nsigma,Nbest,*like);
    delete s;
  }
  
  //Dump summary info
  cout<<"best_post "<<like->bestPost()<<", state="<<like->bestState().get_string()<<endl;
  fl.print_info();
}


