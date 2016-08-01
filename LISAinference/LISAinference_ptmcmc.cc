//Written by John G Baker NASA-GSFC (2014-16)
//This is jointly based on the gleam.cc ptmcmc-based gravitational-lens code
//and the bambi-based LISAinterface.c 

#ifdef PARALLEL
#include "mpi.h"
#endif
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
#include "gsl/gsl_linalg.h"
using namespace std;

typedef initializer_list<double> dlist;
typedef initializer_list<int> ilist;

shared_ptr<Random> globalRNG;//used for some debugging... 

//Global Control parameters (see more in main())
const int Npar=9;
int output_precision;
double fisher_err_target=0.001;

//First we define the all-in-one likelihood function object
class flare_likelihood : public bayes_likelihood {
  sampleable_probability_function * prior;
  int count,nevery;
  double total_eval_time;
  double best_post;
  state best;
  void *context;
public:
  flare_likelihood(stateSpace *sp, void *context, sampleable_probability_function *prior=nullptr):context(context),prior(prior),bayes_likelihood(sp,nullptr,nullptr){
    setPrior(prior);
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
    double tstart=omp_get_wtime();
    LISAParams templateparams =state2LISAParams(s);
    double result=0;

    /* Note: context points to a LISAContext structure containing a LISASignal* */
    if(globalparams->tagint==0) {
      LISAInjectionCAmpPhase injection = *((LISAInjectionCAmpPhase*) context);
      double like;
      like = CalculateLogLCAmpPhase(&templateparams, &injection);
      result = like - logZdata;
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
  
  LISAParams state2LISAParams(const state &s){
    valarray<double>params=s.get_params();

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
    
    return templateparams;
  };

  ///from stateSpaceInterface
  virtual void defWorkingStateSpace(const stateSpace &sp){};

  //virtual void write(ostream &out,state &st){cout<<"flare_likelihood::write: No write routine defined!"<<endl;};
  //virtual void writeFine(ostream &out,state &st,int ns=-1, double ts=0, double te=0){cout<<"flare_likelihood::writeFine: No write routine defined!"<<endl;};
  //virtual void getFineGrid(int & nfine, double &tstart, double &tend)const{cout<<"flare_likelihood::getFineGrid: No routine defined!"<<endl;};
  void print_info(){cout<<" mean = "<<total_eval_time/count<<" through "<<count<<" total evals"<<endl; };
  double getFisher(const state &s0, vector<vector<double> >&fisher_matrix)override{
    //First we must set the injection context
    /* Initialize the data structure for the injection */
    //LISAInjectionCAmpPhase* injectedsignalCAmpPhase = NULL;
    LISAInjectionReIm* injectedsignalReIm = NULL;
    if(globalparams->tagint==0) {
      //LISAInjectionCAmpPhase_Init(&injectedsignalCAmpPhase);
      //LISAGenerateInjectionCAmpPhase(injectedparams, injectedsignalCAmpPhase);
      cout<<"Fisher matrix computation not yet implemented for AmpPhase polar wf representation (try --tagint==1)"<<endl;
      exit(-1);
    }
    else if(globalparams->tagint==1) {
      LISAInjectionReIm_Init(&injectedsignalReIm);
      LISAGenerateInjectionReIm(injectedparams, globalparams->minf, globalparams->nbptsoverlap, 1, injectedsignalReIm); /* Use here logarithmic sampling as a default */
    }
    int dim=s0.size();
    int maxFisherIter=3*dim;
    //double deltafactor=0.001;
    double deltafactor=fisher_err_target;
    double tol=5000*deltafactor*deltafactor*deltafactor;
    valarray<double> scales;
    nativePrior->getScales(scales);
    double err=1;
    int count=0;
    vector< vector<double> >last_fisher_matrix(dim,vector<double>(dim,0));
    while(err>tol&&count<maxFisherIter){
      //cout<<"\ncount="<<count<<"\nscales = ";for(int i=0;i<dim;i++)cout<<scales[i]<<"\t";cout<<endl;
      for(int i=0;i<dim;i++){
	double hi=scales[i]*deltafactor;
	state sPlusi=s0;
	sPlusi.set_param(i,s0.get_param(i)+hi);
	state sMinusi=s0;
	sMinusi.set_param(i,s0.get_param(i)-hi);
	for(int j=i;j<dim;j++){
	  //compute j derivative of model
	  double hj=scales[j]*deltafactor;
	  state sPlusj=s0;
	  sPlusj.set_param(j,s0.get_param(j)+hj);
	  state sMinusj=s0;
	  sMinusj.set_param(j,s0.get_param(j)-hj);
	  //for(int k=0;k<N;k++)dmudj[k]=(modelPlus[k]-modelMinus[k])/h/2.0;
	  //Compute fisher matrix element
	  //cout<<"Fisher i,j="<<i<<","<<j<<endl;
	  //cout<<"            hi,hj="<<hi<<","<<hj<<endl;
	  fisher_matrix[i][j]
	    = CalculateOverlapReIm(state2LISAParams( sPlusi), state2LISAParams( sPlusj), injectedsignalReIm)
	    - CalculateOverlapReIm(state2LISAParams( sPlusi), state2LISAParams(sMinusj), injectedsignalReIm)
	    - CalculateOverlapReIm(state2LISAParams(sMinusi), state2LISAParams( sPlusj), injectedsignalReIm)
	    + CalculateOverlapReIm(state2LISAParams(sMinusi), state2LISAParams(sMinusj), injectedsignalReIm);	 
	  fisher_matrix[i][j]/=4*hi*hj;
	  
	  fisher_matrix[j][i] = fisher_matrix[i][j];
	}
      }
      
      //estimate error
      err=0;
      double square=0;
      for(int i=0;i<dim;i++)for(int j=0;j<dim;j++){
	  double delta=(fisher_matrix[i][j]-last_fisher_matrix[i][j]);///scales[i]/scales[j];
	  square+=fisher_matrix[i][j]*fisher_matrix[i][j];
	  err+=delta*delta;
	  //if(isnan(delta)||delta*delta>tol/10)cout<<"delta["<<i<<","<<j<<"]="<<delta*delta<<endl;
	}
      if(isnan(err)){
	for(int i=0;i<dim;i++)scales[i]/=3;
	err=1;
      }else{
	//set scale estimate based on result
	for(int i=0;i<dim;i++)scales[i]=1.0/sqrt(1/scales[i]+fisher_matrix[i][i]);
	//prep for next version of fisher calc;
	for(int i=0;i<dim;i++)for(int j=0;j<dim;j++)last_fisher_matrix[i][j]=fisher_matrix[i][j];
      }
      count++;
      //cout<<"err->"<<err<<"( goal="<<tol<<")"<<endl;
    }
    err=sqrt(err);
    //cout<<"err="<<err<<endl;
    //cout<<"tol="<<tol<<endl;
    if(err<tol)return tol;
    return err; 
  };

  double CalculateOverlapReIm(LISAParams params1, LISAParams params2, LISAInjectionReIm * injection)
  {
    double overlap = -DBL_MAX;
    int ret;
    
    /* Frequency vector - assumes common to A,E,T, i.e. identical fLow, fHigh in all channels */
    gsl_vector* freq = injection->freq;
    
    /* Generating the signal in the three detectors for the input parameters */
    LISASignalReIm* signal1 = NULL;
    LISASignalReIm* signal2 = NULL;
    LISASignalReIm_Init(&signal1);
    LISASignalReIm_Init(&signal2);
    ret = LISAGenerateSignalReIm(&params1, freq, signal1);
    if(ret==SUCCESS){
      ret = LISAGenerateSignalReIm(&params2, freq, signal2);
    }
    /* If LISAGenerateSignal failed (e.g. parameters out of bound), silently return -Infinity logL */
    if(ret==FAILURE) {
      overlap = -DBL_MAX;
    }
    else if(ret==SUCCESS) {
      /* Computing the likelihood for each TDI channel - fstartobs has already been taken into account */
      double loglikelihoodTDI1 = FDLogLikelihoodReIm(signal1->TDI1Signal, signal2->TDI1Signal, injection->noisevalues1);
      double loglikelihoodTDI2 = FDLogLikelihoodReIm(signal1->TDI2Signal, signal2->TDI2Signal, injection->noisevalues2);
      double loglikelihoodTDI3 = FDLogLikelihoodReIm(signal1->TDI3Signal, signal2->TDI3Signal, injection->noisevalues3);
      overlap = loglikelihoodTDI1 + loglikelihoodTDI2 + loglikelihoodTDI3;
    }
    
    /* Clean up */
    LISASignalReIm_Cleanup(signal1);
    LISASignalReIm_Cleanup(signal2);

    //cout<<" overlap="<<overlap<<endl;
    return overlap;
  }

};

//***************************************************************************************8
//main test program
int main(int argc, char*argv[]){
  int myid = 0;

  //This tell LISAinterface_common not to call MPI stuff;
  noMPI=1;

  //string datafile;
  const int NparRead=Npar; 

  //Create the sampler
  ptmcmc_sampler mcmc;
  bayes_sampler *s0=&mcmc;

  Options opt;

  s0->addOptions(opt,"");
  opt.add(Option("nchains","Number of consequtive chain runs. Default 1","1"));
  opt.add(Option("info_every","How often to dump likehood eval info to stdout. Default never","10000"));
  opt.add(Option("rng_seed","Pseudo random number grenerator seed in [0,1). (Default=-1, use clock to seed.)","-1"));
  opt.add(Option("precision","Set output precision digits. (Default 13).","13"));
  opt.add(Option("Fisher_err_target","Set target for Fisher error measure. (Default 0.001).","0.001"));
  opt.add(Option("help","Print help message."));
  //First we parse the ptmcmc-related parameters like un gleam. 
  opt.parse(argc,argv,false);
  
  //Next we perform the initializations from LISAinference.c
  LISARunParams runParams={};
  int ndim=0,nPar=0;
  int *freeparamsmap = NULL;
  void *context = NULL;
  double logZtrue=0;
  if(opt.set("help")){
    cout<<"****** PTMCMC options *******\n"<<opt.print_usage()<<endl;
    //The rest of this is just to make the code generate the usage message.
    LISAParams arg3;
    LISAGlobalParams arg4;
    LISAPrior arg5;
    LISARunParams arg6;
    char* tmpargv[2]={argv[0],(char*)"--help"};
    parse_args_LISA(2, tmpargv, &arg3, &arg4, &arg5, &arg6);
  } else cout<<"ptmcmc flags=\n"<<opt.report()<<endl;

  addendum(argc,argv,&runParams,&ndim,&nPar,&freeparamsmap,&context,&logZtrue);

  //Post parse setup
  double seed;
  int Nchain,save_every,info_every;
  int Nsigma=1;
  int Nbest=10;
  istringstream(opt.value("nchains"))>>Nchain;
  istringstream(opt.value("rng_seed"))>>seed;
  istringstream(opt.value("info_every"))>>info_every;
  istringstream(opt.value("Fisher_err_target"))>>fisher_err_target;
  
  //if seed<0 set seed from clock
  if(seed<0)seed=fmod(time(NULL)/3.0e7,1);
  istringstream(opt.value("precision"))>>output_precision;
  ProbabilityDist::setSeed(seed);
  globalRNG.reset(ProbabilityDist::getPRNG());//just for safety to keep us from deleting main RNG in debugging.

  //output location
  string outname(runParams.outroot);
  ostringstream ss("");
  ss<<outname;
  string base=ss.str();
  

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
  if (!std::isnan(priorParams->fix_lambda))cout<<" ** Parameter fixing not supported.  Ignoring request! **"<<endl;
  if (!std::isnan(priorParams->fix_beta))cout<<" ** Parameter fixing not supported.  Ignoring request! **"<<endl;
  if (!std::isnan(priorParams->fix_time))cout<<" ** Parameter fixing not supported.  Ignoring request! **"<<endl;
  if (!std::isnan(priorParams->fix_phase))cout<<" ** Parameter fixing not supported.  Ignoring request! **"<<endl;
  if (!std::isnan(priorParams->fix_inc))cout<<" ** Parameter fixing not supported.  Ignoring request! **"<<endl;
  if (!std::isnan(priorParams->fix_pol))cout<<" ** Parameter fixing not supported.  Ignoring request! **"<<endl;
  if (!std::isnan(priorParams->fix_dist))cout<<" ** Parameter fixing not supported.  Ignoring request! **"<<endl;
  if (!std::isnan(priorParams->fix_m1))cout<<" ** Parameter fixing not supported.  Ignoring request! **"<<endl;
  if (!std::isnan(priorParams->fix_m2))cout<<" ** Parameter fixing not supported.  Ignoring request! **"<<endl;
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
  //valarray<int>         types((ilist){  uni,    uni,     uni,     uni,    uni, polar,     uni,  copol,    uni  });
  valarray<int>         types{  uni,    uni,     uni,     uni,    uni, polar,     uni,  copol,    uni  };
  if (!std::isnan(priorParams->fix_beta))types[7]=uni;//"fixed" parameter set to narrow uniform range
  if (!std::isnan(priorParams->fix_inc)) types[5]=uni;//"fixed" parameter set to narrow uniform range
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
  
  //Next compute the Fisher matrix at the injected state.
  ss.str("");ss<<base<<"_fishcov.dat";
  ofstream outfish(ss.str().c_str());
  vector< vector<double> > fim(Npar,vector<double>(Npar));
  double err=like->getFisher(injected_state,fim);
  cout<<"Fisher err ~"<<err<<endl;
  cout<<"At pars:\n"<<endl;
  for(int i=0;i<Npar;i++)cout<<names[i]<<"\t";
  cout<<endl;
  cout<<"\nFisher at injection:"<<endl;
  outfish<<"#At pars:\n#";
  for(int i=0;i<Npar;i++)outfish<<names[i]<<"\t";
  outfish<<endl;
  for(int i=0;i<Npar;i++)outfish<<injected_state.get_param(i)<<"\t";
  outfish<<endl;
  outfish<<"\n#Fisher:"<<endl;
  for(int i=0;i<Npar;i++){
    cout<<names[i]<<"\t";
    for(int j=0;j<=i;j++)cout<<"\t"<<fim[i][j];
    cout<<endl;
    for(int j=0;j<=i;j++)outfish<<"\t"<<fim[i][j];
    for(int j=i+1;j<Npar;j++)outfish<<"\t"<<fim[i][j];
    outfish<<endl;
  }
  //Next dump Fisher scaled by diagonals;
  cout<<"Fisher scales:"<<endl;
  for(int i=0;i<Npar;i++)cout<<sqrt(fim[i][i])<<"\t";
  cout<<endl;
  outfish<<"\n#Fisher scales:"<<endl;
  for(int i=0;i<Npar;i++)outfish<<sqrt(fim[i][i])<<"\t";
  outfish<<endl;
  cout<<"\nScaled Fisher:"<<endl;
  outfish<<"\n#Scaled Fisher:"<<endl;
  for(int i=0;i<Npar;i++){
    cout<<names[i]<<"\t";
    for(int j=0;j<=i;j++)cout<<"\t"<<fim[i][j]/sqrt(fim[i][i]*fim[j][j]);
    cout<<endl;
    for(int j=0;j<=i;j++)outfish<<"\t"<<fim[i][j]/sqrt(fim[i][i]*fim[j][j]);
    for(int j=i+1;j<Npar;j++)outfish<<"\t"<<fim[i][j]/sqrt(fim[i][i]*fim[j][j]);
    outfish<<endl;
  }
  //Now invert it to get the covariance matrix
  //We use GSLs Cholesky decomposition for this.
  gsl_set_error_handler_off();
  gsl_matrix * fishcov=gsl_matrix_alloc(Npar,Npar);
  for(int i=0;i<Npar;i++)for(int j=0;j<=i;j++)gsl_matrix_set(fishcov,i,j,fim[i][j]);//Don't need upper triangle for GSLs routine
  if(gsl_linalg_cholesky_decomp(fishcov))cout<<"Fisher matrix Cholesky decomposition failed."<<endl;
  else if(gsl_linalg_cholesky_invert(fishcov))cout<<"Fisher matrix Cholesky inverse failed."<<endl;
  else {
    cout<<"\nCovariance:"<<endl;
    outfish<<"\n#Covariance:"<<endl;
    for(int i=0;i<Npar;i++){
      cout<<names[i]<<"\t";
      for(int j=0;j<=i;j++)cout<<"\t"<<gsl_matrix_get(fishcov,i,j);
      cout<<endl;
      for(int j=0;j<=i;j++)outfish<<"\t"<<gsl_matrix_get(fishcov,i,j);
      for(int j=i+1;j<Npar;j++)outfish<<"\t"<<gsl_matrix_get(fishcov,i,j);
      outfish<<endl;
    }
    //Now the correlation matrix:
    cout<<"\nCovar scales:"<<endl;
    for(int i=0;i<Npar;i++)cout<<sqrt(gsl_matrix_get(fishcov,i,i))<<"\t";
    outfish<<"\n#Covar scales:"<<endl;
    for(int i=0;i<Npar;i++)outfish<<sqrt(gsl_matrix_get(fishcov,i,i))<<"\t";
    cout<<"\n\nCorrelation:"<<endl;
    outfish<<"\n\n#Coorrelation:"<<endl;
    for(int i=0;i<Npar;i++){
      cout<<names[i]<<"\t";
      for(int j=0;j<=i;j++)cout<<"\t"<<gsl_matrix_get(fishcov,i,j)/sqrt(gsl_matrix_get(fishcov,i,i)*gsl_matrix_get(fishcov,j,j));
      cout<<endl;
      for(int j=0;j<=i;j++)outfish<<"\t"<<gsl_matrix_get(fishcov,i,j)/sqrt(gsl_matrix_get(fishcov,i,i)*gsl_matrix_get(fishcov,j,j));
  for(int j=i+1;j<Npar;j++)outfish<<"\t"<<gsl_matrix_get(fishcov,i,j)/sqrt(gsl_matrix_get(fishcov,i,i)*gsl_matrix_get(fishcov,j,j));
      outfish<<endl;
    }
  }
  //Close-out Fisher/Covar
  gsl_matrix_free(fishcov);
  outfish.close();
  if(stoi(opt.value("nsteps"))<=0){
    cout<<"No MCMC steps requested"<<endl;
    exit(0);
  }

  //assuming mcmc:
  //Set the proposal distribution 
  int Ninit;
  proposal_distribution *prop=ptmcmc_sampler::new_proposal_distribution(Npar,Ninit,opt,prior,&halfwidths);
  cout<<"Proposal distribution is:\n"<<prop->show()<<endl;
  //set up the mcmc sampler (assuming mcmc)
  mcmc.setup(Ninit,*like,*prior,*prop,output_precision);

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


