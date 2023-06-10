#include "BioMiCo.h"

//unsigned int myRmultinom( NumericVector probs ){                                                                                                                                  
//  unsigned int rmultinomDraw = 0;                                                                                                                                                 
//                                                                                                                                                                                  
//  return rmultinomDraw ;                                                                                                                                                          
//}                                                                                                                                                                                 


#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/normal.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <vector>
#include <algorithm>
#include <numeric>

#include <iostream>
#include <fstream>

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_errno.h>

#include <math.h>


using namespace boost::numeric::ublas;
using namespace boost::numeric;

boost::mt19937 gen(static_cast<unsigned int>(std::time(0)));;


/*`
  Now define a function that simulates rolling a loaded die.
  Note that the C++0x library contains a `discrete_distribution`
  class which would be a better way to do this.
*/
int roll_weighted_die( const std::vector<double> probabilities) {
  std::vector<double> cumulative;
  std::partial_sum(probabilities.begin(), probabilities.end(),
		   std::back_inserter(cumulative));
  boost::uniform_real<> dist(0, cumulative.back());
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > die(gen, dist);

  return (std::lower_bound(cumulative.begin(), cumulative.end(), die()) - cumulative.begin()) ;
}

matrix<double> Alpha_theta_Moment_Matching_update(const ublas::vector< matrix<unsigned> > &Sample_Count_G_V, const matrix<unsigned int> &Count_S_V, unsigned effective_V, unsigned effective_G ){
  unsigned N = Count_S_V.size1();
  unsigned V = Count_S_V.size2();
  unsigned G = Sample_Count_G_V(0).size1();
  matrix<double> RET_alpha_theta_mat(V, G ); RET_alpha_theta_mat.clear();
  matrix<double> var_mat(V, G); var_mat.clear();
  matrix<double> mean_mat(V, G); mean_mat.clear();
  matrix<double> m_mat(V, G); m_mat.clear();

  std::vector<double> sum_mean_mat_vec(V, 0);
  /*
  cout << "N=" << N << ", V=" << V << ", G=" << G << endl;  
  for ( unsigned n=0; n< N; n++){
    cout << "Sample_Count_G_V for n=" << n << endl;
    for (unsigned v=0; v < V; v++){
      for (unsigned g=0; g < G; g++)
	cout << Sample_Count_G_V(n)(g, v) << " ";
      cout << endl;
    }
  }
  */
  std::vector<unsigned> Samples_from_V(V,0);
  for ( unsigned n=0; n< N; n++)
    for (unsigned v=0; v < effective_V; v++)
      	if (Count_S_V(n,v) > 0)
	  Samples_from_V[v] ++;
  
  for ( unsigned n=0; n< N; n++)
    for (unsigned v=0; v < effective_V; v++)
      	if (Count_S_V(n,v) > 0)
	  for (unsigned g=0; g < effective_G; g++){
	    mean_mat(v,g) += ( Sample_Count_G_V(n)(g, v) + 0.01) / (Samples_from_V[v] * (Count_S_V(n,v) + effective_G * 0.01) );
	    sum_mean_mat_vec[v] += mean_mat(v,g);
	  }
  /*
  cout << "mean_mat" << endl;
  for (unsigned v=0; v < V; v++){
    for (unsigned g=0; g < G; g++)
      cout << mean_mat(v,g) << " " ;
    cout << endl;
  }
  */
  for ( unsigned n=0; n< N; n++)
    for (unsigned v=0; v < effective_V; v++)
      if (Count_S_V(n,v) > 0)
	for (unsigned g=0; g < effective_G; g++)
	  var_mat(v,g) += pow( (Sample_Count_G_V(n)(g, v) + 0.01) / ( Count_S_V(n,v) + effective_G * 0.01) - mean_mat(v,g)  , 2) / Samples_from_V[v];
  /*
  cout << "var_mat" << endl;
  for (unsigned v=0; v < V; v++){
    for (unsigned g=0; g < G; g++)
      cout << var_mat(v,g) << " " ;
    cout << endl;
  }
  */
  for (unsigned v=0; v < effective_V; v++)
    for (unsigned g=0; g < effective_G; g++)
      m_mat(v,g) = mean_mat(v,g) * (1 - mean_mat(v,g)) / var_mat(v,g) -1;
  /*
  cout << "m_mat" << endl;
  for (unsigned v=0; v < V; v++){
    for (unsigned g=0; g < G; g++)
      cout << m_mat(v,g) << " " ;
    cout << endl;
  }
  */

  std::vector<double> sum_alpha_theta_vec(V, 0);
  for (unsigned v=0; v < effective_V; v++){
    double sum_log_m_mat_row = 0;
    for (unsigned g=0; g < effective_G; g++)
      sum_log_m_mat_row += log(m_mat(v,g));
    sum_alpha_theta_vec[v] = exp(sum_log_m_mat_row / (effective_G-1));
    for (unsigned g=0; g < effective_G; g++)
      RET_alpha_theta_mat(v,g) = mean_mat(v,g) * sum_alpha_theta_vec[v] / sum_mean_mat_vec[v];
  }
  /*
  cout << "RET_alpha_theta_mat" << endl;
  for (unsigned v=0; v < V; v++){
    for (unsigned g=0; g < G; g++)
      cout << RET_alpha_theta_mat(v,g) << " " ;
    cout << endl;
  }
  */
  return RET_alpha_theta_mat;
}

matrix<double> Alpha_theta_Metropolis_Hastings_update(const matrix<double>& alpha_theta_mat, std::vector<double> sum_alpha_theta_vec, const matrix<unsigned int> Count_G_V, const std::vector<unsigned int> sum_Count_G_V, unsigned effective_V, unsigned effective_G){
  //RNGScope scope;
  boost::normal_distribution<> nd(0.0, 1.0);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(gen, nd);
  //var_nor.engine().seed(static_cast<unsigned int>(std::time(NULL) + getpid()));
  //var_nor.distribution().reset();
  boost::uniform_real<> dist(0, 1);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > unif_rand(gen, dist);
  using boost::math::normal;
  normal s;

  double SearchWidth = 0.1;
  matrix<double> RET_alpha_theta_mat( alpha_theta_mat );
  matrix<double> NEW_alpha_theta_mat( alpha_theta_mat );
  unsigned G = Count_G_V.size1();
  unsigned V = Count_G_V.size2();  
  double pAccept = 0;
  bool Accept = 0;
  double AcceptanceRate = 0;
  double RET_alpha_theta_val = RET_alpha_theta_mat(0, 0);
  double NEW_alpha_theta_val = RET_alpha_theta_val + var_nor() * RET_alpha_theta_val*SearchWidth;
  for (unsigned v=0; v < effective_V; v++)
    for (unsigned g=0; g < effective_G; g++)
      NEW_alpha_theta_mat( v,g ) = NEW_alpha_theta_val;

  double RET_sum_alpha_theta = RET_alpha_theta_val * effective_G;
  double NEW_sum_alpha_theta = NEW_alpha_theta_val * effective_G;

  pAccept += effective_V*( lgamma(NEW_sum_alpha_theta) - effective_G*lgamma(NEW_alpha_theta_val) );
  pAccept -= effective_V*( lgamma(RET_sum_alpha_theta) - effective_G*lgamma(RET_alpha_theta_val) );
  for (unsigned v=0; v < effective_V; v++){
    pAccept -= lgamma(sum_Count_G_V[v] + NEW_sum_alpha_theta);
    pAccept += lgamma(sum_Count_G_V[v] + RET_sum_alpha_theta);
    for (unsigned g=0; g < effective_G; g++){
      pAccept += ( lgamma(NEW_alpha_theta_val + Count_G_V(g, v)) );
      pAccept -= ( lgamma(RET_alpha_theta_val + Count_G_V(g, v)) );
    }
  }
  pAccept = exp(pAccept);
  pAccept *= ( pdf(s, (RET_alpha_theta_val - NEW_alpha_theta_val) /(SearchWidth*NEW_alpha_theta_val) ) ) / ( pdf(s, (NEW_alpha_theta_val - RET_alpha_theta_val) /(SearchWidth*RET_alpha_theta_val) ) );
  if (pAccept >= 1)
    Accept = 1;
  else
    Accept = (unif_rand() < pAccept);

  if (Accept)
    return NEW_alpha_theta_mat;
  else
    return RET_alpha_theta_mat;
  
  
}



matrix<double> Alpha_theta_Mat_Metropolis_Hastings_update(const matrix<double>& alpha_theta_mat, std::vector<double> sum_alpha_theta_vec, const matrix<unsigned int> Count_G_V, const std::vector<unsigned int> sum_Count_G_V, unsigned effective_V, unsigned effective_G){
  //RNGScope scope;
  boost::normal_distribution<> nd(0.0, 1.0);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(gen, nd);
  //var_nor.engine().seed(static_cast<unsigned int>(std::time(NULL) + getpid()));
  //var_nor.distribution().reset();
  boost::uniform_real<> dist(0, 1);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > unif_rand(gen, dist);
  using boost::math::normal;
  normal s;

  double SearchWidth = 0.1;
  matrix<double> RET_alpha_theta_mat( alpha_theta_mat );
  std::vector<double> RET_sum_alpha_theta_vec( sum_alpha_theta_vec );
  matrix<double> NEW_alpha_theta_mat( alpha_theta_mat );
  std::vector<double> NEW_sum_alpha_theta_vec( sum_alpha_theta_vec );
  unsigned G = Count_G_V.size1();
  unsigned V = Count_G_V.size2();
  
  double pAccept = 0;
  bool Accept = 0;
  double AcceptanceRate = 0;

  for (unsigned v=0; v < effective_V; v++)
    for (unsigned g=0; g < effective_G; g++){
      NEW_alpha_theta_mat( v,g ) = RET_alpha_theta_mat(v, g) + var_nor() * RET_alpha_theta_mat(v, g)*SearchWidth;
      NEW_sum_alpha_theta_vec[v] += (NEW_alpha_theta_mat( v,g ) - RET_alpha_theta_mat( v,g ) );
      pAccept += ( lgamma( NEW_sum_alpha_theta_vec[v]) + lgamma( Count_G_V(g,v)+NEW_alpha_theta_mat(v,g) ) - lgamma( NEW_alpha_theta_mat(v,g) ) - lgamma( sum_Count_G_V[v]+NEW_sum_alpha_theta_vec[v]) );
      pAccept -= ( lgamma( RET_sum_alpha_theta_vec[v]) + lgamma( Count_G_V(g,v)+RET_alpha_theta_mat(v,g) ) - lgamma( RET_alpha_theta_mat(v,g) ) - lgamma( sum_Count_G_V[v]+RET_sum_alpha_theta_vec[v]) );
      /*
	double p_NEW_alpha = 0;
	for (unsigned i=0; i < Count_G_V(g,v); i++)
	p_NEW_alpha += log( NEW_alpha_theta_mat( v,g ) + i);
	for (unsigned i=0; i < sum_Count_G_V[v]; i++)
	p_NEW_alpha -= log( NEW_sum_alpha_theta_vec[v] + i);

	double p_RET_alpha = 0;
	for (unsigned i=0; i < Count_G_V(g,v); i++)
	p_RET_alpha += log( RET_alpha_theta_mat( v,g ) + i);
	for (unsigned i=0; i < sum_Count_G_V[v]; i++)
	p_RET_alpha -= log( RET_sum_alpha_theta_vec[v] + i);

	pAccept = p_NEW_alpha - p_RET_alpha;
      */
      pAccept = exp(pAccept);
      pAccept *= ( pdf(s, (RET_alpha_theta_mat(v, g) - NEW_alpha_theta_mat( v,g )) /(SearchWidth*NEW_alpha_theta_mat( v,g )) ) ) / ( pdf(s, (NEW_alpha_theta_mat(v, g) - RET_alpha_theta_mat( v,g )) /(SearchWidth*RET_alpha_theta_mat( v,g )) ) );

      if (pAccept >= 1)
        Accept = 1;
      else
        Accept = (unif_rand() < pAccept);

      if (Accept){
	//cout << "old value = " <<  RET_alpha_theta_mat(v, g) << " new value = " << NEW_alpha_theta_mat(v, g) << endl;
	AcceptanceRate += 1/(V*G);
	RET_alpha_theta_mat( v,g ) = NEW_alpha_theta_mat( v,g );
	RET_sum_alpha_theta_vec[v] = NEW_sum_alpha_theta_vec[v];
      }

    }
  
  return RET_alpha_theta_mat;
}

double Alpha_phi_unknown_Metropolis_Hastings_update(const double alpha_phi_unknown, const std::vector<unsigned int> Count_T_G_unknown, unsigned int sum_Count_T_G_unknown){
  boost::normal_distribution<> nd(0.0, 1.0);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(gen, nd);
  //var_nor.engine().seed(static_cast<unsigned int>(std::time(NULL) + getpid()));
  //var_nor.distribution().reset();
  boost::uniform_real<> dist(0, 1);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > unif_rand(gen, dist);
  using boost::math::normal;
  normal s;

  double SearchWidth = 0.1;
  double RET_alpha_phi = alpha_phi_unknown;
  unsigned T = Count_T_G_unknown.size();
  double RET_sum_alpha_phi = alpha_phi_unknown * T;

  double pAccept = 0;
  bool Accept = 0;
  double NEW_alpha_phi = RET_alpha_phi + var_nor() * RET_alpha_phi*SearchWidth;
  double NEW_sum_alpha_phi = NEW_alpha_phi * T;

  pAccept +=  lgamma(NEW_sum_alpha_phi) - T*lgamma(NEW_alpha_phi);
  pAccept -=  lgamma(RET_sum_alpha_phi) - T*lgamma(RET_alpha_phi);
  pAccept -= lgamma(sum_Count_T_G_unknown + NEW_sum_alpha_phi);
  pAccept += lgamma(sum_Count_T_G_unknown + RET_sum_alpha_phi);
  for (unsigned t=0; t < T; t++){
    pAccept += ( lgamma(NEW_alpha_phi + Count_T_G_unknown[t]) );
    pAccept -= ( lgamma(RET_alpha_phi + Count_T_G_unknown[t]) );
  }


  pAccept = exp(pAccept);
  pAccept *= ( pdf(s, (RET_alpha_phi - NEW_alpha_phi) /(SearchWidth*NEW_alpha_phi) ) ) / ( pdf(s, (NEW_alpha_phi - RET_alpha_phi) /(SearchWidth*RET_alpha_phi) ) );

  if (pAccept >= 1)
    Accept = 1;
  else
    Accept = (unif_rand() < pAccept);
  
  if (Accept)
    RET_alpha_phi = NEW_alpha_phi;
  
  return RET_alpha_phi;

}


double Alpha_phi_Metropolis_Hastings_update(const double alpha_phi, const matrix<unsigned int> Count_T_G, const std::vector<unsigned int> sum_Count_T_G, unsigned effective_G){
  boost::normal_distribution<> nd(0.0, 1.0);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(gen, nd);
  //var_nor.engine().seed(static_cast<unsigned int>(std::time(NULL) + getpid()));
  //var_nor.distribution().reset();
  boost::uniform_real<> dist(0, 1);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > unif_rand(gen, dist);
  using boost::math::normal;
  normal s;

  double SearchWidth = 0.1;
  double RET_alpha_phi = alpha_phi;
  unsigned T = Count_T_G.size1();
  unsigned G = Count_T_G.size2();
  double RET_sum_alpha_phi = alpha_phi * T;

  double pAccept = 0;
  bool Accept = 0;
  double NEW_alpha_phi = RET_alpha_phi + var_nor() * RET_alpha_phi*SearchWidth;
  double NEW_sum_alpha_phi = NEW_alpha_phi * T;

  pAccept += effective_G*( lgamma(NEW_sum_alpha_phi) - T*lgamma(NEW_alpha_phi) );
  pAccept -= effective_G*( lgamma(RET_sum_alpha_phi) - T*lgamma(RET_alpha_phi) );
  for (unsigned g=0; g < effective_G; g++){
    pAccept -= lgamma(sum_Count_T_G[g] + NEW_sum_alpha_phi);
    pAccept += lgamma(sum_Count_T_G[g] + RET_sum_alpha_phi);
    for (unsigned t=0; t < T; t++){
      pAccept += ( lgamma(NEW_alpha_phi + Count_T_G(t, g)) );
      pAccept -= ( lgamma(RET_alpha_phi + Count_T_G(t, g)) );
    }
  }

  pAccept = exp(pAccept);
  pAccept *= ( pdf(s, (RET_alpha_phi - NEW_alpha_phi) /(SearchWidth*NEW_alpha_phi) ) ) / ( pdf(s, (NEW_alpha_phi - RET_alpha_phi) /(SearchWidth*RET_alpha_phi) ) );
  if (pAccept >= 1)
    Accept = 1;
  else
    Accept = (unif_rand() < pAccept);
  
  if (Accept)
    RET_alpha_phi = NEW_alpha_phi;
  
  return RET_alpha_phi;
}

double Alpha_pi_Metropolis_Hastings_update(const double alpha_pi, const matrix<unsigned int> Count_S_V, const std::vector<unsigned int> sum_Count_S_V, Rcpp::List data_envs, std::vector<unsigned int> num_envs){
  boost::normal_distribution<> nd(0.0, 1.0);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(gen, nd);
  //var_nor.engine().seed(static_cast<unsigned int>(std::time(NULL) + getpid()));
  //var_nor.distribution().reset();
  boost::uniform_real<> dist(0, 1);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > unif_rand(gen, dist);
  using boost::math::normal;
  normal s;

  double SearchWidth = 0.1;
  double RET_alpha_pi = alpha_pi;
  unsigned N = Count_S_V.size1();
  unsigned V = Count_S_V.size2();

  double pAccept = 0;
  bool Accept = 0;
  double NEW_alpha_pi = RET_alpha_pi + var_nor() * RET_alpha_pi*SearchWidth;
  
  for (unsigned n=0; n<N; n++){
    pAccept += ( lgamma( num_envs[n] * NEW_alpha_pi) - num_envs[n] * lgamma( NEW_alpha_pi) );
    pAccept -= ( lgamma( num_envs[n] * RET_alpha_pi) - num_envs[n] * lgamma( RET_alpha_pi) );

    pAccept -= lgamma(sum_Count_S_V[n] + num_envs[n] * NEW_alpha_pi);
    pAccept += lgamma(sum_Count_S_V[n] + num_envs[n] * RET_alpha_pi);

    Rcpp::IntegerVector cur_envs_lst = as<Rcpp::IntegerVector>(data_envs(n));
    for (unsigned x=0; x<num_envs[n]; x++){
      pAccept += lgamma( Count_S_V(n, cur_envs_lst(x)) + NEW_alpha_pi);
      pAccept -= lgamma( Count_S_V(n, cur_envs_lst(x)) + RET_alpha_pi);
    }
  }

  pAccept = exp(pAccept);
  pAccept *= ( pdf(s, (RET_alpha_pi - NEW_alpha_pi) /(SearchWidth*NEW_alpha_pi) ) ) / ( pdf(s, (NEW_alpha_pi - RET_alpha_pi) /(SearchWidth*RET_alpha_pi) ) );
  if (pAccept >= 1)
    Accept = 1;
  else
    Accept = (unif_rand() < pAccept);
  
  if (Accept)
    RET_alpha_pi = NEW_alpha_pi;
  
  return RET_alpha_pi;
  
}


RcppExport SEXP runGibbs(SEXP Rdata_samples, SEXP Rdata_envs, SEXP RG, SEXP RV, SEXP RT, SEXP RN, SEXP RunknownEnv,//"run.gibbs" <- function(data.samples, data.envs, G, V, T, N,
			 SEXP RZ, SEXP RX,   //      Z=NULL, X=NULL, 
			 SEXP Rburnin, SEXP Rnrestarts, SEXP Rndraws_per_restart, SEXP Rdelay,  //burnin=100, nrestarts=10, ndraws.per.restart=10, delay=10,
			 SEXP Ralpha_pi, SEXP Ralpha_phi, SEXP Ralpha_theta, SEXP Rmaxdepth,  //      alpha.pi=1e-2, alpha.phi=1e-2, alpha.theta=1e-2, maxdepth=NULL,
			 SEXP Rverbosity, SEXP Rprinting_index, SEXP Rprinting_total){  //      verbosity=1, printing.index=NULL, printing.total=NULL)

  cout << "Entering runGibbs" << endl;
  Rcpp::IntegerMatrix data_samples = Rcpp::as<Rcpp::IntegerMatrix>(Rdata_samples);
  cout << "data_samples has " << data_samples.nrow() << " rows and " << data_samples.ncol() << " columns." << endl;
  Rcpp::List data_envs(Rdata_envs);
  unsigned int G = Rcpp::as<int>(RG);
  cout << "#communities: " << G << endl;

  unsigned int V = Rcpp::as<int>(RV);
  cout << "#environments: " << V << endl;
  unsigned int T = Rcpp::as<int>(RT);
  cout << "#taxa: " << T << endl;
  unsigned int N = Rcpp::as<int>(RN);
  cout << "#samples: " << data_envs.size() << endl;

  bool unknownEnv = Rcpp::as<bool>(RunknownEnv);
  //if (unknownEnv){
  //  V++;
  //  G++
  //}

  Rcpp::List Z(RZ);
  Rcpp::List X(RX);
  unsigned int burnin = Rcpp::as<int>(Rburnin);
  unsigned int nrestarts = Rcpp::as<int>(Rnrestarts);
  unsigned int ndraws_per_restart = Rcpp::as<int>(Rndraws_per_restart);
  unsigned int delay = Rcpp::as<int>(Rdelay);

  double alpha_pi = Rcpp::as<double>(Ralpha_pi);
  std::vector<double> alpha_pi_vec(V, alpha_pi);
  double alpha_pi_unknown = alpha_pi;
  if (unknownEnv){
    alpha_pi_unknown = alpha_pi / 10; 
    alpha_pi_vec[V-1] = alpha_pi_unknown;
  }
  double alpha_pi_sum = (V-1)*alpha_pi + alpha_pi_vec[V-1];

  double alpha_phi = Rcpp::as<double>(Ralpha_phi);
  double alpha_phi_unknown = alpha_phi;// / 50;

  double alpha_theta = Rcpp::as<double>(Ralpha_theta);
  //matrix<double> alpha_theta_mat(V, G);
 Rcpp:NumericMatrix alpha_theta_mat(V, G);
  std::vector<double> sum_alpha_theta_vec(V, G*alpha_theta);
  for (unsigned i=0; i < V; i++)
    for (unsigned j=0; j < G; j++)
      alpha_theta_mat(i,j) = alpha_theta;
  if (unknownEnv){
    for (unsigned j=0; j < G-1; j++)
      alpha_theta_mat(V-1,j) = 0;
    for (unsigned i=0; i < V-1; i++)
      alpha_theta_mat(i, G-1) = 0;
    alpha_theta_mat(V-1, G-1) = alpha_theta; //1;
  }
    

  unsigned int maxdepth = Rcpp::as<int>(Rmaxdepth);
  unsigned int verbosity = Rcpp::as<int>(Rverbosity);
  unsigned int printing_index = Rcpp::as<int>(Rprinting_index);
  unsigned int printing_total = Rcpp::as<int>(Rprinting_total);

  cout << "Parsing arguments finished" << endl;


  gen.seed(static_cast<unsigned int>(std::time(NULL) + getpid()));


  matrix<unsigned int> Count_T_G(T, G); //  Count.T.G = matrix(0, nrow=T, ncol=G) # counts of taxas assigned to groups, Count.T.G[t,g]: number of times taxa t is assigned to group g
  for (unsigned t=0; t < T; t++)
    for (unsigned g=0; g < G; g++)
      Count_T_G(t, g) =0;
  std::vector<unsigned int> sum_Count_T_G(G,0); //  sum.Count.T.G = rep(0,G)

  matrix<unsigned int> Count_G_V(G, V); //  Count.G.V = matrix(0, nrow=G, ncol=V) # counts of taxa groups being observed in environments, Count.G.V[g,v]: number of times taxa group g was observed in environment v
  for (unsigned x=0; x < V; x++)
    for (unsigned g=0; g < G; g++)
      Count_G_V(g, x) =0;
  std::vector<unsigned int> sum_Count_G_V(V,0); //  sum.Count.G.V = rep(0,V)

  ublas::vector< matrix<unsigned> > Sample_Count_G_V;
  Sample_Count_G_V.resize(N);
  for (unsigned n=0; n < N; n++){
    Sample_Count_G_V(n).resize(G, V);
    Sample_Count_G_V(n).clear();
  }

  matrix<unsigned int> Count_S_V(N, V); //  Count.S.V = matrix(0, nrow=N, ncol=V) # counts of taxa observed in a sample from an environment, Count.S.V[s,v]: number of times a taxa came from environment v in sample s
  for (unsigned n=0; n < N; n++)
    for (unsigned x=0; x < V; x++)
      Count_S_V(n, x) =0;
  std::vector<unsigned int> sum_Count_S_V(N,0); //  sum.Count.S.V = rep(0,N)

  
  std::vector<unsigned int> num_envs(N, 0);
  for  (unsigned n=0; n<N; n++){
    num_envs[n] = as<Rcpp::IntegerVector>(data_envs(n)).size();
    if (unknownEnv && ( as<Rcpp::IntegerVector>(data_envs(n))(num_envs[n]-1) != V-1 ))
      cout << "Something wrong, assuming unknown env exists but the last env is not unknown!" << endl;
    if (n==1 || n==4){
      cout << "envs list for " << n+1 << " : ";
      for (unsigned x=0; x < num_envs[n]; x++)
	cout << as<Rcpp::IntegerVector>(data_envs(n))(x) << "\t";
      cout << endl;
    }
    if (num_envs[n] < 1)
      cout << "Something wrong, no environment for sample " << n+1 << endl;
  }

  if (Z.size() == 0){ //  if (is.null(Z)){
    cout << "Initializing counter variables" << endl; //    cat("Initializing counter variables\n")
    Z = Rcpp::List(N*T); //    Z = list()  //    length(Z) = N*T
    X = Rcpp::List(N*T); //    X = list()  //    length(X) = N*T
    for (unsigned n=0; n<N; n++){  //    for (n in 1:N){ # looping over all samples

      std::vector<double> probZ;
      std::vector<double> probX;
      if (unknownEnv){
	probZ = std::vector<double>(G-1, 1.0/(G-1));
	probX = std::vector<double>(num_envs[n], 1.0/num_envs[n]);
	//probX = std::vector<double>(num_envs[n]-1, 0.9/(num_envs[n]-1)); // 0.1 for the unknown env and 0.9 for all other envs
	//probX.resize(num_envs[n], 0.1);
      } else {
	probZ = std::vector<double>(G, 1.0/G);
	probX = std::vector<double>(num_envs[n], 1.0/num_envs[n]);
      }

      for (unsigned t=0; t<T; t++){ //      for (t in 1:T){# looping over all taxa
	if ( data_samples(n,t) > 0){ //        if (data.samples[n,t] > 0){
	  Z(n*T + t) = rep(0, data_samples(n,t)); // Z.insert( n*T + t, rep(0, data_samples(n,t) )); //          Z[[(n-1)*T + t]] = rep(0, data.samples[n,t])
	  X(n*T + t) = rep(0, data_samples(n,t)); // X.insert( n*T + t, rep(0, data_samples(n,t) )); //          X[[(n-1)*T + t]] = rep(0, data.samples[n,t])
	  for (unsigned i=0; i < data_samples(n,t); i++){ //          for (i in 1:data.samples[n,t]){
	    unsigned cur_envs = roll_weighted_die(probX); //            cur.envs = sample( data.envs[[n]], 1, replace=TRUE)
	    as<Rcpp::IntegerVector>(X( n*T + t ))(i) = as<Rcpp::IntegerVector>(data_envs(n))(cur_envs); //            X[[(n-1)*T + t]][i] = cur.envs
	    unsigned cur_G = roll_weighted_die(probZ);
	    if (unknownEnv && ( cur_envs == (num_envs[n]-1) ) ) // if cur_envs is not unknown, then taxa comes from ones of the first communities
	      cur_G = G-1;
	    as<Rcpp::IntegerVector>(Z( n*T + t ))(i) =  cur_G; //  cur.G = sample( 1:G, 1, replace=TRUE)  //            Z[[(n-1)*T + t]][i] = cur.G
	  } //          }
	} //        }
      } //      }
    } //    }      
  } //  }

  cout << "Filling counter variables" << endl; //  cat("Filling counter variables\n")
  for (unsigned n=0; n<N; n++){ //  for (n in 1:N){ # looping over all samples
    for (unsigned t=0; t<T; t++){ //    for (t in 1:T){# looping over all taxa
      if ( data_samples(n,t) > 0){ //      if (data.samples[n,t] > 0){
	for (unsigned i=0; i < data_samples(n,t); i++){ //        for (i in 1:data.samples[n,t]){
	  int cur_envs = as<Rcpp::IntegerVector>(X( n*T + t ))(i); //          cur.envs = X[[(n-1)*T + t]][i]
	  Count_S_V(n, cur_envs) = Count_S_V(n, cur_envs) +1; //          Count.S.V[n, cur.envs] = Count.S.V[n, cur.envs] +1
	  sum_Count_S_V[n] = sum_Count_S_V[n] +1; //          sum.Count.S.V[n] = sum.Count.S.V[n] +1
	  int cur_G = as<Rcpp::IntegerVector>(Z( n*T + t ))(i); //          cur.G = Z[[(n-1)*T + t]][i]
	  Count_G_V(cur_G, cur_envs) = Count_G_V(cur_G, cur_envs) + 1; //          Count.G.V[cur.G, cur.envs] = Count.G.V[cur.G, cur.envs] + 1
	  sum_Count_G_V[cur_envs] = sum_Count_G_V[cur_envs] +1; //          sum.Count.G.V[cur.envs] = sum.Count.G.V[cur.envs] +1
	  Sample_Count_G_V(n)(cur_G, cur_envs) = Sample_Count_G_V(n)(cur_G, cur_envs) +1;
	  Count_T_G(t, cur_G) = Count_T_G(t, cur_G) +1; //          Count.T.G[t, cur.G] = Count.T.G[t, cur.G] +1
	  sum_Count_T_G[cur_G] = sum_Count_T_G[cur_G] + 1; //          sum.Count.T.G[cur.G] = sum.Count.T.G[cur.G] + 1
	} //        }
      } //      }
    } //    }
  } //  }  

  cout << endl;
  for (unsigned x=0; x < V; x++){
    for (unsigned g=0; g < G; g++)
      cout << Count_G_V(g, x) << "\t";
    cout << endl;
  }

  //  # start Gibbs sampling
  //

  unsigned nonzero_counts = 0;
  std::vector<double> T_counts(T, 0.0);
  for (unsigned n=0; n<N; n++) 
    for (unsigned t=0; t<T; t++) 
      if ( data_samples(n,t) > 0){
  	nonzero_counts += data_samples(n,t);
	T_counts[t] += data_samples(n,t);
      }

  std::vector< std::vector<unsigned> >  data_iterator(nonzero_counts, std::vector<unsigned>(3,0));//  data.iterator = list();
  //data_iterator.resize(nonzero_counts); //  length(data.iterator) = sum(data.samples > 0)
  unsigned data_iterator_idx = 0; //  data.iterator.idx = 1
  for (unsigned n=0; n<N; n++) //  for (n in 1:N)
    for (unsigned t=0; t<T; t++) //    for (t in 1:T)
      if ( data_samples(n,t) > 0) //      if (data.samples[n,t] > 0)
	for (unsigned i=0; i < data_samples(n,t); i++){ //        for (i in 1:data.samples[n,t]){
	  data_iterator[data_iterator_idx][0] = n; 
	  data_iterator[data_iterator_idx][1] = t; 
	  data_iterator[data_iterator_idx][2] = i; //          data.iterator[[data.iterator.idx]] = c(n,t,i) #c(data.iterator, list(c(n,t,i)))
	  data_iterator_idx ++; //          data.iterator.idx = data.iterator.idx +1
	} //        }

  unsigned ndraws = nrestarts * ndraws_per_restart; //  ndraws <- nrestarts * ndraws.per.restart # total number of draws
  unsigned npasses = burnin + (ndraws_per_restart-1) * delay + 1; //  npasses <- burnin + (ndraws.per.restart-1) * delay + 1 # passes per restart
  unsigned drawcount = 1; //  drawcount <- 1 # keeps running count of draws for this sample
  Rcpp::List savedGibbsIterations(ndraws); //  savedGibbsIterations = list()
  //  length(savedGibbsIterations) = ndraws
  cout << "starting Gibbs sampler" << endl;//  cat("starting Gibbs sampler\n")  
  for (unsigned restart_counter=1;  restart_counter <= nrestarts; restart_counter++){ //  for (restart.counter in 1:nrestarts){
    if (verbosity>=1) //    if(verbosity>=1) cat('.')
      cout << endl << "chain " <<  restart_counter << " of " << nrestarts << endl;
    //    #options(warn=-1)
    //
    for (unsigned itr=1; itr <= npasses; itr++){ //    for (itr in 1:npasses){
      //cout << "iteration: " << itr << " of chain: " << restart_counter << endl; //      cat(paste("iteration: ", itr, " of chain: ", restart.counter, "\n", sep=""))
      if (verbosity>=1) //    if(verbosity>=1) cat('.')
	cout << "." << flush;
      unsigned effective_G = G;
      unsigned effective_X;
      std::vector< std::vector<unsigned> > data_iterator_permuted = data_iterator;
      std::random_shuffle(data_iterator_permuted.begin(), data_iterator_permuted.end() ); //      data.iterator.permuted = data.iterator[sample(1:length(data.iterator))]
      for (std::vector< std::vector<unsigned> >::iterator taxa_bit=data_iterator_permuted.begin(); taxa_bit != data_iterator_permuted.end(); taxa_bit++ ){ //      for (taxa.bit in data.iterator.permuted){ # looping over all samples
	unsigned n = (*taxa_bit)[0]; //        n = taxa.bit[1]
	unsigned t = (*taxa_bit)[1]; //        t = taxa.bit[2]
	unsigned i = (*taxa_bit)[2]; //        i = taxa.bit[3]
	      
	unsigned int cur_G = as<Rcpp::IntegerVector>(Z( n*T + t ))(i); //        cur.G = Z[[(n-1)*T + t]][i]
	unsigned int cur_envs = as<Rcpp::IntegerVector>(X( n*T + t ))(i); //        cur.envs = X[[(n-1)*T + t]][i]
	Count_S_V(n, cur_envs) = Count_S_V(n, cur_envs) -1; //        Count.S.V[n, cur.envs] = Count.S.V[n, cur.envs] -1
	sum_Count_S_V[n] = sum_Count_S_V[n] -1; //        sum.Count.S.V[n] = sum.Count.S.V[n] -1
	Count_G_V(cur_G, cur_envs) = Count_G_V(cur_G, cur_envs) -1; //        Count.G.V[cur.G, cur.envs] = Count.G.V[cur.G, cur.envs] -1
	sum_Count_G_V[cur_envs] = sum_Count_G_V[cur_envs] -1; //        sum.Count.G.V[cur.envs] = sum.Count.G.V[cur.envs] -1
	Sample_Count_G_V(n)(cur_G, cur_envs) = Sample_Count_G_V(n)(cur_G, cur_envs) -1;
	Count_T_G(t, cur_G) = Count_T_G(t, cur_G) -1; //        Count.T.G[t, cur.G] = Count.T.G[t, cur.G] -1
	sum_Count_T_G[cur_G] = sum_Count_T_G[cur_G] -1; //        sum.Count.T.G[cur.G] = sum.Count.T.G[cur.G] -1
		
	effective_G = G;
	effective_X = num_envs[n];
	if (unknownEnv){  
	  alpha_pi_sum = (num_envs[n]-1) * alpha_pi + alpha_pi_unknown;
	  effective_G = G-1;
	  effective_X = num_envs[n]-1;
	}
	else
	  alpha_pi_sum = num_envs[n] * alpha_pi;	
	
	std::vector<double> x_z_prob(effective_X*effective_G, 0); //        x.z.prob = rep(0, length(data.envs[[n]])*G)
	if (unknownEnv)
	  x_z_prob.resize(effective_X*effective_G+1, 0);
	
	//double sum_x_z_prob = 0;
	Rcpp::IntegerVector cur_envs_lst = as<Rcpp::IntegerVector>(data_envs(n));

	for (unsigned xidx=0; xidx < effective_X; xidx++ ) //        for (x.idx in 1:length(data.envs[[n]]) )
	  for (unsigned zidx=0; zidx < effective_G; zidx++){ //          for (z.idx in 1:G)   
	    x_z_prob[ xidx*effective_G + zidx ] = ( (Count_T_G(t, zidx) + alpha_phi) / ( sum_Count_T_G[zidx] + T*alpha_phi) ) 
	      * ( (Count_S_V(n, cur_envs_lst(xidx)) + alpha_pi_vec[cur_envs_lst(xidx)] ) ); // / (sum_Count_S_V[n] + alpha_pi_sum)  );
	    
	    //if (!unknownEnv || ( xidx != num_envs[n]-1) )
	    x_z_prob[ xidx*effective_G + zidx ] *= ( (Count_G_V(zidx, cur_envs_lst(xidx) ) + alpha_theta_mat( cur_envs_lst(xidx), zidx ) ) / ( sum_Count_G_V[ cur_envs_lst(xidx)] + sum_alpha_theta_vec[cur_envs_lst(xidx)] ) ) ;

	    //sum_x_z_prob += x_z_prob[ xidx*G + zidx ];
	  }
	
	if (unknownEnv){
	  //x_z_prob[ effective_X*effective_G ] = ((T_counts[t] + 0.1)/(nonzero_counts + 0.1*T)) * ( (Count_S_V(n, cur_envs_lst(effective_X)) + alpha_pi_vec[cur_envs_lst(effective_X)] ) ); // / (sum_Count_S_V[n] + alpha_pi_sum)  );
	  //x_z_prob[ effective_X*effective_G ] = 1.0/T * ( (Count_S_V(n, cur_envs_lst(effective_X)) + alpha_pi_vec[cur_envs_lst(effective_X)] ) ); // / (sum_Count_S_V[n] + alpha_pi_sum)  );
	  
	  x_z_prob[ effective_X*effective_G ] = ( (Count_T_G(t, effective_G) + alpha_phi_unknown) / ( sum_Count_T_G[effective_G] + T*alpha_phi_unknown) ) 
	      * ( (Count_S_V(n, cur_envs_lst(effective_X)) + alpha_pi_vec[cur_envs_lst(effective_X)] ) ); // / (sum_Count_S_V[n] + alpha_pi_sum)  );
	  //x_z_prob[ effective_X*effective_G ] *= 
	  //  ( (Count_G_V(effective_G, cur_envs_lst(effective_X) ) + alpha_theta_mat( cur_envs_lst(effective_X), effective_G ) ) / ( sum_Count_G_V[ cur_envs_lst(effective_X)] + sum_alpha_theta_vec[cur_envs_lst(effective_X)] ) ) ;
	  
	}
	
	
	//for (unsigned xidx=0; xidx < num_envs[n]; xidx++ )
	//  for (unsigned zidx=0; zidx < G; zidx++)
	//    x_z_prob[ xidx*G + zidx ] /= sum_x_z_prob;
	
	unsigned XZ = roll_weighted_die( x_z_prob ); //        XZ = which(rmultinom(1,1, x.z.prob)==1);
	if (unknownEnv && (XZ == effective_X*effective_G)){
	  cur_G = effective_G;
	  cur_envs = effective_X;
	}else{
	  cur_G = XZ % effective_G; //        cur.G = XZ %% G;
	  cur_envs = (XZ - cur_G)/effective_G; //        cur.envs = (XZ - cur.G)/G;
	}

	cur_envs = cur_envs_lst(cur_envs); //        cur.envs = data.envs[[n]][cur.envs]
	        
	as<Rcpp::IntegerVector>(X( n*T + t ))(i) = cur_envs;//        X[[(n-1)*T + t]][i] = cur.envs
	as<Rcpp::IntegerVector>(Z( n*T + t ))(i) = cur_G;//        Z[[(n-1)*T + t]][i] = cur.G
	Count_S_V(n, cur_envs) = Count_S_V(n, cur_envs) +1; //        Count.S.V[n, cur.envs] = Count.S.V[n, cur.envs] +1
	sum_Count_S_V[n] = sum_Count_S_V[n] +1; //        sum.Count.S.V[n] = sum.Count.S.V[n] +1
	Count_G_V(cur_G, cur_envs) = Count_G_V(cur_G, cur_envs) + 1; //        Count.G.V[cur.G, cur.envs] = Count.G.V[cur.G, cur.envs] + 1
	sum_Count_G_V[cur_envs] = sum_Count_G_V[cur_envs] +1; //        sum.Count.G.V[cur.envs] = sum.Count.G.V[cur.envs] +1
	Sample_Count_G_V(n)(cur_G, cur_envs) = Sample_Count_G_V(n)(cur_G, cur_envs) +1;
	Count_T_G(t, cur_G) = Count_T_G(t, cur_G) +1; //        Count.T.G[t, cur.G] = Count.T.G[t, cur.G] +1
	sum_Count_T_G[cur_G] = sum_Count_T_G[cur_G] + 1; //        sum.Count.T.G[cur.G] = sum.Count.T.G[cur.G] + 1
	
      } //      }
      //      

      // Update alpha_theta prior hyperparameters
      for (unsigned MH_itr=0; MH_itr < 1; MH_itr++){
	matrix<double> new_alpha_theta_mat(V, G);
	for (unsigned v=0; v < V; v++)
	  for (unsigned g=0; g < G; g++)
	    new_alpha_theta_mat(v,g) = alpha_theta_mat(v,g);
	new_alpha_theta_mat = Alpha_theta_Metropolis_Hastings_update( new_alpha_theta_mat, sum_alpha_theta_vec, Count_G_V, sum_Count_G_V, V, G);
	//new_alpha_theta_mat = Alpha_theta_Moment_Matching_update(Sample_Count_G_V, Count_S_V, V, G );
	for (unsigned v=0; v < V; v++){
	  sum_alpha_theta_vec[v] = 0;
	  for (unsigned g=0; g < G; g++){
	    //cout << alpha_theta_mat(v,g) << " ";
	    alpha_theta_mat(v,g) = new_alpha_theta_mat(v,g);
	    sum_alpha_theta_vec[v] += alpha_theta_mat(v,g);
	  }
	  //cout << endl;
	}
	
	alpha_phi = Alpha_phi_Metropolis_Hastings_update(alpha_phi, Count_T_G, sum_Count_T_G, effective_G);
	if (unknownEnv){
	  std::vector<unsigned int> Count_T_G_unknown(T,0);
	  for (unsigned t=0; t<T; t++)
	    Count_T_G_unknown[t] = Count_T_G(t, G-1);
	  alpha_phi_unknown = Alpha_phi_unknown_Metropolis_Hastings_update(alpha_phi_unknown, Count_T_G_unknown, sum_Count_T_G[G-1]);
	}
	
	alpha_pi = Alpha_pi_Metropolis_Hastings_update( alpha_pi, Count_S_V, sum_Count_S_V, data_envs, num_envs);
	for (unsigned xidx=0; xidx<V; xidx++) alpha_pi_vec[xidx] = alpha_pi;
      }

      if(itr > burnin && (((itr - burnin) % delay)==1 || delay<=1)){ //      if(itr > burnin && (((itr - burnin) %% delay)==1 || delay<=1)){                        
	//        # save current Gibbs iteration
	cout << endl << "saving sample " << drawcount << endl; //        cat(paste("saving sample ", drawcount, "\n", sep=""))
	savedGibbsIterations(drawcount-1) =  Rcpp::List::create(Rcpp::Named("Z")=clone(Z), Rcpp::Named("X")=clone(X), Rcpp::Named("ALPHA_THETA_MAT")=clone(alpha_theta_mat), Rcpp::Named("ALPHA_PHI")=alpha_phi) ;
	//savedGibbsIterations.insert(drawcount-1, Rcpp::List::create(Rcpp::Named("Z")=Z, Rcpp::Named("X")=X) ); //        savedGibbsIterations[[drawcount]] = list(Z=Z, X=X)
	drawcount ++; //        drawcount <- drawcount + 1
      } //      }
      //      
    } //    }
  } //  }
  if (verbosity>=1) //    if(verbosity>=1) cat('.')
      cout << endl ;  
  
  //
  
  //for (unsigned n=0; n<N; n++) 
  //  for (unsigned t=0; t<T; t++) 
  //    if ( data_samples(n,t) > 0)
  //	for (unsigned i=0; i < data_samples(n,t); i++){ 
  //	  unsigned int cur_G = as<Rcpp::IntegerVector>( as<Rcpp::List>( as<Rcpp::List>(savedGibbsIterations(ndraws-1))["Z"]) ( n*T + t ) )(i);
  //	  cout << n << " , " << t << " , " << i << " , " << cur_G << endl;
  //	}

  return savedGibbsIterations; //  return(invisible(savedGibbsIterations))
  //

}










// function for sampling community composition and environment assignments for aaaa set of samples

RcppExport SEXP runGibbsOnTest(SEXP Rdata_samples, SEXP RG, SEXP RV, SEXP RT, SEXP RN, SEXP RunknownEnv,
			 SEXP RCount_T_G, SEXP RCount_G_V,   
			 SEXP Rburnin, SEXP Rnrestarts, SEXP Rndraws_per_restart, SEXP Rdelay,  
			 SEXP Ralpha_pi, SEXP Ralpha_phi, SEXP Ralpha_theta, SEXP Rmaxdepth,  
			 SEXP Rverbosity, SEXP Rprinting_index, SEXP Rprinting_total){  


  cout << "Entering runGibbsOnTest" << endl;
  Rcpp::IntegerMatrix data_samples = Rcpp::as<Rcpp::IntegerMatrix>(Rdata_samples);
  cout << "data_samples has " << data_samples.nrow() << " rows and " << data_samples.ncol() << " columns." << endl;
  unsigned int G = Rcpp::as<int>(RG);
  cout << "#communities: " << G << endl;
  unsigned int V = Rcpp::as<int>(RV);
  cout << "#environments: " << V << endl;
  unsigned int T = Rcpp::as<int>(RT);
  cout << "#taxa: " << T << endl;
  unsigned int N = Rcpp::as<int>(RN);
  cout << "#samples: " << data_samples.nrow() << endl;
  bool unknownEnv = Rcpp::as<bool>(RunknownEnv);
  if (unknownEnv){
    V++;
    G++;
  }

  Rcpp::IntegerMatrix Count_T_G_train = Rcpp::as<Rcpp::IntegerMatrix>(RCount_T_G); // counts of taxas assigned to groups, Count.T.G[t,g]: number of times taxa t is assigned to group g
  Rcpp::IntegerMatrix Count_G_V_train = Rcpp::as<Rcpp::IntegerMatrix>(RCount_G_V); // counts of taxa groups being observed in environments, Count.G.V[g,v]: number of times taxa group g was observed in environment v

  matrix<unsigned int> Count_G_V(G, V);
  matrix<unsigned int> Count_T_G(T, G);
  for (unsigned t=0; t < T; t++)
    for (unsigned g=0; g < G; g++)
      Count_T_G(t, g) =0;
  for (unsigned x=0; x < V; x++)
    for (unsigned g=0; g < G; g++)
      Count_G_V(g, x) =0;

  unsigned train_G = Count_G_V_train.nrow();
  unsigned train_V = Count_G_V_train.ncol();
  for (unsigned g_idx=0; g_idx< train_G; g_idx++){
    for (unsigned v_idx=0; v_idx< train_V; v_idx++)
      Count_G_V(g_idx, v_idx) = Count_G_V_train(g_idx, v_idx);
    for (unsigned t=0; t<T; t++)
      Count_T_G(t, g_idx) = Count_T_G_train(t, g_idx);
  }
    
  unsigned int burnin = Rcpp::as<int>(Rburnin);
  unsigned int nrestarts = Rcpp::as<int>(Rnrestarts);
  unsigned int ndraws_per_restart = Rcpp::as<int>(Rndraws_per_restart);
  unsigned int delay = Rcpp::as<int>(Rdelay);
  double alpha_pi = Rcpp::as<double>(Ralpha_pi);
  std::vector<double> alpha_pi_vec(V, alpha_pi);
  double alpha_pi_unknown = alpha_pi;
  if (unknownEnv)
    alpha_pi_unknown = alpha_pi ; // / 10; 
  alpha_pi_vec[V-1] = alpha_pi_unknown;
  double alpha_pi_sum = (V-1)*alpha_pi + alpha_pi_vec[V-1];

  double alpha_phi = Rcpp::as<double>(Ralpha_phi);
  double alpha_phi_unknown = alpha_phi;// / 50;
  double alpha_theta = Rcpp::as<double>(Ralpha_theta);
  matrix<double> alpha_theta_mat(V, G);
  std::vector<double> sum_alpha_theta_vec(V, G*alpha_theta);
  for (unsigned i=0; i < V; i++)
    for (unsigned j=0; j < G; j++)
      alpha_theta_mat(i,j) = alpha_theta;
  if (unknownEnv){
    for (unsigned j=0; j < G-1; j++)
      alpha_theta_mat(V-1,j) = 0;
    for (unsigned i=0; i < V-1; i++){
      alpha_theta_mat(i, G-1) = 0;
      sum_alpha_theta_vec[i] = (G-1)*alpha_theta;
    }
    alpha_theta_mat(V-1, G-1) = alpha_theta; //1;
  }

  unsigned int maxdepth = Rcpp::as<int>(Rmaxdepth);
  unsigned int verbosity = Rcpp::as<int>(Rverbosity);
  unsigned int printing_index = Rcpp::as<int>(Rprinting_index);
  unsigned int printing_total = Rcpp::as<int>(Rprinting_total);

  cout << "Parsing arguments finished" << endl;

  gen.seed(static_cast<unsigned int>(std::time(NULL) + getpid()));

  std::vector<unsigned int> sum_Count_T_G(G,0); 
  for  (unsigned g=0; g<G; g++)  for ( unsigned t=0; t<T; t++) sum_Count_T_G[g] += Count_T_G(t,g);

  std::vector<unsigned int> sum_Count_G_V(V,0); 
  for  (unsigned g=0; g<G; g++)  for ( unsigned v=0; v<V; v++) sum_Count_G_V[v] += Count_G_V(g,v);

  ublas::vector< matrix<unsigned> > Sample_Count_G_V;
  Sample_Count_G_V.resize(N);
  for (unsigned n=0; n < N; n++){
    Sample_Count_G_V(n).resize(G, V);
    Sample_Count_G_V(n).clear();
  }


  matrix<unsigned int> Count_S_V(N, V); //  counts of taxa observed in a sample from an environment, Count.S.V[s,v]: number of times a taxa came from environment v in sample s
  for (unsigned n=0; n < N; n++)
    for (unsigned x=0; x < V; x++)
      Count_S_V(n, x) =0;
  std::vector<unsigned int> sum_Count_S_V(N,0); 

  Rcpp::List data_envs(N);
  std::vector<unsigned int> num_envs(N, 0);
  Rcpp::IntegerVector all_envs_idx(V);
  for (unsigned i=0; i<V; i++) all_envs_idx(i) = i;
  for  (unsigned n=0; n<N; n++){
    data_envs(n) = all_envs_idx;
    num_envs[n] = V; 
  }

  Rcpp::List Z;
  Rcpp::List X;

  /*  if (Z.size() == 0){ //  if (is.null(Z)){
    cout << "Initializing counter variables" << endl; //    cat("Initializing counter variables\n")
    Z = Rcpp::List(N*T); //    Z = list()  //    length(Z) = N*T
    X = Rcpp::List(N*T); //    X = list()  //    length(X) = N*T
    for (unsigned n=0; n<N; n++){  //    for (n in 1:N){ # looping over all samples
      std::vector<double> probZ(G, 1.0/G);
      std::vector<double> probX(num_envs[n], 1.0/num_envs[n]);
      for (unsigned t=0; t<T; t++){ //      for (t in 1:T){# looping over all taxa
	if ( data_samples(n,t) > 0){ //        if (data.samples[n,t] > 0){
	  Z(n*T + t) = rep(0, data_samples(n,t)); // Z.insert( n*T + t, rep(0, data_samples(n,t) )); //          Z[[(n-1)*T + t]] = rep(0, data.samples[n,t])
	  X(n*T + t) = rep(0, data_samples(n,t)); // X.insert( n*T + t, rep(0, data_samples(n,t) )); //          X[[(n-1)*T + t]] = rep(0, data.samples[n,t])
	  for (unsigned i=0; i < data_samples(n,t); i++){ //          for (i in 1:data.samples[n,t]){
	    unsigned cur_envs = roll_weighted_die(probX); //            cur.envs = sample( data.envs[[n]], 1, replace=TRUE)
	    as<Rcpp::IntegerVector>(X( n*T + t ))(i) = as<Rcpp::IntegerVector>(data_envs(n))(cur_envs); //            X[[(n-1)*T + t]][i] = cur.envs
	    as<Rcpp::IntegerVector>(Z( n*T + t ))(i) = roll_weighted_die(probZ); //  cur.G = sample( 1:G, 1, replace=TRUE)  //            Z[[(n-1)*T + t]][i] = cur.G
	  } //          }
	} //        }
      } //      }
    } //    }      
    }*/ //  }

  /*  if (Z.size() == 0){ 
    cout << "Initializing counter variables" << endl; 
    Z = Rcpp::List(N*T); 
    X = Rcpp::List(N*T); 
    for (unsigned n=0; n<N; n++){  
      std::vector<double> probZ(G, 1.0/G);
      std::vector<double> probX(num_envs[n], 1.0/num_envs[n]);
      for (unsigned t=0; t<T; t++){
	if ( data_samples(n,t) > 0){
	  Z(n*T + t) = rep(0, data_samples(n,t)); 
	  X(n*T + t) = rep(0, data_samples(n,t)); 
	  for (unsigned i=0; i < data_samples(n,t); i++){ 
	    std::vector<double> x_z_prob(num_envs[n]*G, 0); 
	    for (unsigned xidx=0; xidx < num_envs[n]; xidx++ ) 
	      for (unsigned zidx=0; zidx < G; zidx++) 
		x_z_prob[ xidx*G + zidx ] = ( (Count_T_G(t, zidx) + alpha_phi) / ( sum_Count_T_G[zidx] + T*alpha_phi) ) 
		  * ( (Count_G_V(zidx, xidx ) + alpha_theta) / ( sum_Count_G_V[ xidx] + G*alpha_theta ) ) 
		  * ( 1.0 / num_envs[n]   );

	    unsigned XZ = roll_weighted_die( x_z_prob ); 
	    unsigned cur_G = XZ % G; 
	    unsigned cur_envs = (XZ - cur_G)/G; 

	    as<Rcpp::IntegerVector>(X( n*T + t ))(i) = cur_envs; 
	    as<Rcpp::IntegerVector>(Z( n*T + t ))(i) = cur_G;
	  } 
	} 
      } 
    } 
  } 
  */

  if (Z.size() == 0){ //  if (is.null(Z)){
    cout << "Initializing counter variables" << endl; //    cat("Initializing counter variables\n")
    Z = Rcpp::List(N*T); 
    X = Rcpp::List(N*T); 
    for (unsigned n=0; n<N; n++){  
      std::vector<double> probZ;
      std::vector<double> probX;
      if (unknownEnv){
	probZ = std::vector<double>(G-1, 1.0/(G-1));
	probX = std::vector<double>(num_envs[n], 1.0/num_envs[n]);
      } else {
	probZ = std::vector<double>(G, 1.0/G);
	probX = std::vector<double>(num_envs[n], 1.0/num_envs[n]);
      }
      for (unsigned t=0; t<T; t++){ //      for (t in 1:T){# looping over all taxa
	if ( data_samples(n,t) > 0){ //        if (data.samples[n,t] > 0){
	  Z(n*T + t) = rep(0, data_samples(n,t)); // Z.insert( n*T + t, rep(0, data_samples(n,t) )); //          Z[[(n-1)*T + t]] = rep(0, data.samples[n,t])
	  X(n*T + t) = rep(0, data_samples(n,t)); // X.insert( n*T + t, rep(0, data_samples(n,t) )); //          X[[(n-1)*T + t]] = rep(0, data.samples[n,t])
	  for (unsigned i=0; i < data_samples(n,t); i++){ //          for (i in 1:data.samples[n,t]){
	    unsigned cur_envs = roll_weighted_die(probX); //            cur.envs = sample( data.envs[[n]], 1, replace=TRUE)
	    as<Rcpp::IntegerVector>(X( n*T + t ))(i) = as<Rcpp::IntegerVector>(data_envs(n))(cur_envs); //            X[[(n-1)*T + t]][i] = cur.envs
	    unsigned cur_G = roll_weighted_die(probZ);
	    if (unknownEnv && ( cur_envs == (num_envs[n]-1) ) ) // if cur_envs is not unknown, then taxa comes from ones of the first communities
	      cur_G = G-1;
	    as<Rcpp::IntegerVector>(Z( n*T + t ))(i) =  cur_G; //  cur.G = sample( 1:G, 1, replace=TRUE)  //            Z[[(n-1)*T + t]][i] = cur.G
	  } //          }
	} //        }
      } //      }
    } //    }      
  } //  }


  cout << "Filling counter variables" << endl; 
  for (unsigned n=0; n<N; n++){ 
    for (unsigned t=0; t<T; t++){ 
      if ( data_samples(n,t) > 0){
	for (unsigned i=0; i < data_samples(n,t); i++){ 
	  int cur_envs = as<Rcpp::IntegerVector>(X( n*T + t ))(i); 
	  Count_S_V(n, cur_envs) = Count_S_V(n, cur_envs) +1; 
	  sum_Count_S_V[n] = sum_Count_S_V[n] +1; 
	  int cur_G = as<Rcpp::IntegerVector>(Z( n*T + t ))(i);  
	  Count_G_V(cur_G, cur_envs) = Count_G_V(cur_G, cur_envs) + 1; 
	  sum_Count_G_V[cur_envs] = sum_Count_G_V[cur_envs] +1; 
	  Sample_Count_G_V(n)(cur_G, cur_envs) = Sample_Count_G_V(n)(cur_G, cur_envs) +1;
	  Count_T_G(t, cur_G) = Count_T_G(t, cur_G) +1; 
	  sum_Count_T_G[cur_G] = sum_Count_T_G[cur_G] + 1; 
	} 
      } 
    } 
  } 

  //  # start Gibbs sampling

  unsigned nonzero_counts = 0;
  for (unsigned n=0; n<N; n++) 
    for (unsigned t=0; t<T; t++) 
      if ( data_samples(n,t) > 0)
  	nonzero_counts += data_samples(n,t);

  std::vector< std::vector<unsigned> >  data_iterator(nonzero_counts, std::vector<unsigned>(3,0));
  unsigned data_iterator_idx = 0; 
  for (unsigned n=0; n<N; n++) 
    for (unsigned t=0; t<T; t++) 
      if ( data_samples(n,t) > 0) 
	for (unsigned i=0; i < data_samples(n,t); i++){ 
	  data_iterator[data_iterator_idx][0] = n; 
	  data_iterator[data_iterator_idx][1] = t; 
	  data_iterator[data_iterator_idx][2] = i; 
	  data_iterator_idx ++; 
	} 

  unsigned ndraws = nrestarts * ndraws_per_restart; 
  unsigned npasses = burnin + (ndraws_per_restart-1) * delay + 1; 
  unsigned drawcount = 1; 
  Rcpp::List savedGibbsIterations(ndraws); 

  cout << "starting Gibbs sampler" << endl;
  for (unsigned restart_counter=1;  restart_counter <= nrestarts; restart_counter++){ 
    if (verbosity>=1) 
      cout << endl << "chain " <<  restart_counter << " of " << nrestarts << endl;

    for (unsigned itr=1; itr <= npasses; itr++){ 
      if (verbosity>=1) 
	cout << "." << flush;

      unsigned effective_G = G;
      unsigned effective_X;

      std::vector< std::vector<unsigned> > data_iterator_permuted = data_iterator;
      std::random_shuffle(data_iterator_permuted.begin(), data_iterator_permuted.end() ); 
      for (std::vector< std::vector<unsigned> >::iterator taxa_bit=data_iterator_permuted.begin(); taxa_bit != data_iterator_permuted.end(); taxa_bit++ ){
	unsigned n = (*taxa_bit)[0]; 
	unsigned t = (*taxa_bit)[1]; 
	unsigned i = (*taxa_bit)[2]; 
	      
	unsigned int cur_G = as<Rcpp::IntegerVector>(Z( n*T + t ))(i); 
	unsigned int cur_envs = as<Rcpp::IntegerVector>(X( n*T + t ))(i); 
	Count_S_V(n, cur_envs) = Count_S_V(n, cur_envs) -1; 
	sum_Count_S_V[n] = sum_Count_S_V[n] -1; 
	Count_G_V(cur_G, cur_envs) = Count_G_V(cur_G, cur_envs) -1; 
	sum_Count_G_V[cur_envs] = sum_Count_G_V[cur_envs] -1; 
	Sample_Count_G_V(n)(cur_G, cur_envs) = Sample_Count_G_V(n)(cur_G, cur_envs) -1;
	Count_T_G(t, cur_G) = Count_T_G(t, cur_G) -1; 
	sum_Count_T_G[cur_G] = sum_Count_T_G[cur_G] -1; 

	effective_G = G;
	effective_X = num_envs[n];
	if (unknownEnv){  
	  alpha_pi_sum = (num_envs[n]-1) * alpha_pi + alpha_pi_unknown;
	  effective_G = G-1;
	  effective_X = num_envs[n]-1;
	}
	else
	  alpha_pi_sum = num_envs[n] * alpha_pi;	
	
	std::vector<double> x_z_prob(effective_X*effective_G, 0); //        x.z.prob = rep(0, length(data.envs[[n]])*G)
	if (unknownEnv)
	  x_z_prob.resize(effective_X*effective_G+1, 0);
	
	Rcpp::IntegerVector cur_envs_lst = as<Rcpp::IntegerVector>(data_envs(n));

	for (unsigned xidx=0; xidx < effective_X; xidx++ ) //        for (x.idx in 1:length(data.envs[[n]]) )
	  for (unsigned zidx=0; zidx < effective_G; zidx++){ //          for (z.idx in 1:G)   
	    x_z_prob[ xidx*effective_G + zidx ] = ( (Count_T_G(t, zidx) + alpha_phi) / ( sum_Count_T_G[zidx] + T*alpha_phi) ) 
	      * ( (Count_S_V(n, cur_envs_lst(xidx)) + alpha_pi_vec[cur_envs_lst(xidx)] ) ); // / (sum_Count_S_V[n] + alpha_pi_sum)  );
	    
	    x_z_prob[ xidx*effective_G + zidx ] *= ( (Count_G_V(zidx, cur_envs_lst(xidx) ) + alpha_theta_mat( cur_envs_lst(xidx), zidx ) ) / ( sum_Count_G_V[ cur_envs_lst(xidx)] + sum_alpha_theta_vec[cur_envs_lst(xidx)] ) ) ;

	  }
	
	if (unknownEnv){ 
	  x_z_prob[ effective_X*effective_G ] = ( (Count_T_G(t, effective_G) + alpha_phi_unknown) / ( sum_Count_T_G[effective_G] + T*alpha_phi_unknown) ) 
	      * ( (Count_S_V(n, cur_envs_lst(effective_X)) + alpha_pi_vec[cur_envs_lst(effective_X)] ) ); // / (sum_Count_S_V[n] + alpha_pi_sum)  );
	}



	/*
	std::vector<double> x_z_prob(num_envs[n]*G, 0); 
	Rcpp::IntegerVector cur_envs_lst = as<Rcpp::IntegerVector>(data_envs(n));
	if (unknownEnv)   
	  alpha_pi_sum = (num_envs[n]-1) * alpha_pi + alpha_pi_unknown;
	else
	  alpha_pi_sum = num_envs[n] * alpha_pi;
	for (unsigned xidx=0; xidx < num_envs[n]; xidx++ ) 
	  for (unsigned zidx=0; zidx < G; zidx++) 
	    x_z_prob[ xidx*G + zidx ] = ( (Count_T_G(t, zidx) + alpha_phi) / ( sum_Count_T_G[zidx] + T*alpha_phi) ) 
	      * ( (Count_G_V(zidx, cur_envs_lst(xidx) ) + alpha_theta) / ( sum_Count_G_V[ cur_envs_lst(xidx)] + G*alpha_theta ) ) 
	      * ( (Count_S_V(n, cur_envs_lst(xidx)) + alpha_pi_vec[cur_envs_lst(xidx)] ) / (sum_Count_S_V[n] + alpha_pi_sum)  );
	
	if (unknownEnv){
	  for (unsigned xidx=0; xidx < (num_envs[n]-1); xidx++ )
	    x_z_prob[ xidx*G + G-1 ] = 0;
	  for (unsigned zidx=0; zidx < G-1; zidx++) 
	    x_z_prob[ (num_envs[n]-1)*G + zidx ] = 0;
	}
	*/

	unsigned XZ = roll_weighted_die( x_z_prob ); 
	if (unknownEnv && (XZ == effective_X*effective_G)){
	  cur_G = effective_G;
	  cur_envs = effective_X;
	}else{
	  cur_G = XZ % effective_G; //        cur.G = XZ %% G;
	  cur_envs = (XZ - cur_G)/effective_G; //        cur.envs = (XZ - cur.G)/G;
	}

	/*
	cur_G = XZ % G; 
	cur_envs = (XZ - cur_G)/G; 
	*/
	cur_envs = cur_envs_lst[cur_envs]; 
  
	as<Rcpp::IntegerVector>(X( n*T + t ))(i) = cur_envs;
	as<Rcpp::IntegerVector>(Z( n*T + t ))(i) = cur_G;
	Count_S_V(n, cur_envs) = Count_S_V(n, cur_envs) +1; 
	sum_Count_S_V[n] = sum_Count_S_V[n] +1; 
	Count_G_V(cur_G, cur_envs) = Count_G_V(cur_G, cur_envs) + 1; 
	sum_Count_G_V[cur_envs] = sum_Count_G_V[cur_envs] +1; 
	Sample_Count_G_V(n)(cur_G, cur_envs) = Sample_Count_G_V(n)(cur_G, cur_envs) +1;
	Count_T_G(t, cur_G) = Count_T_G(t, cur_G) +1; 
	sum_Count_T_G[cur_G] = sum_Count_T_G[cur_G] + 1; 
	
      } 
      
      // Update alpha_theta prior hyperparameters
      //
      for (unsigned MH_itr=0; MH_itr < 1; MH_itr++){
	if (unknownEnv){
	  //alpha_theta_mat = Alpha_theta_Moment_Matching_update(Sample_Count_G_V, Count_S_V, V-1, G-1 );
	  alpha_theta_mat = Alpha_theta_Metropolis_Hastings_update( alpha_theta_mat, sum_alpha_theta_vec, Count_G_V, sum_Count_G_V, V-1, G-1);
	  alpha_phi = Alpha_phi_Metropolis_Hastings_update(alpha_phi, Count_T_G, sum_Count_T_G, G-1);
	  std::vector<unsigned int> Count_T_G_unknown(T,0);
	  for (unsigned t=0; t<T; t++)
	    Count_T_G_unknown[t] = Count_T_G(t, G-1);
	  alpha_phi_unknown = Alpha_phi_unknown_Metropolis_Hastings_update(alpha_phi_unknown, Count_T_G_unknown, sum_Count_T_G[G-1]);
	}else{
	  //alpha_theta_mat = Alpha_theta_Moment_Matching_update(Sample_Count_G_V, Count_S_V, V, G );
	  alpha_theta_mat = Alpha_theta_Metropolis_Hastings_update( alpha_theta_mat, sum_alpha_theta_vec, Count_G_V, sum_Count_G_V, V, G);
	  alpha_phi = Alpha_phi_Metropolis_Hastings_update(alpha_phi, Count_T_G, sum_Count_T_G, G);
	}
	for (unsigned v=0; v < V; v++){
	  sum_alpha_theta_vec[v] = 0;
	  for (unsigned g=0; g < G; g++)
	    sum_alpha_theta_vec[v] += alpha_theta_mat(v,g);
	}
	alpha_pi = Alpha_pi_Metropolis_Hastings_update( alpha_pi, Count_S_V, sum_Count_S_V, data_envs, num_envs);
	for (unsigned xidx=0; xidx<V; xidx++) alpha_pi_vec[xidx] = alpha_pi;

      }
      
      if(itr > burnin && (((itr - burnin) % delay)==1 || delay<=1)){ 
	// save current Gibbs iteration
	cout << endl << "saving sample " << drawcount << endl; 
	//savedGibbsIterations(drawcount-1) =  Rcpp::List::create(Rcpp::Named("Z")=Z, Rcpp::Named("X")=X) ;
	savedGibbsIterations(drawcount-1) =  Rcpp::List::create(Rcpp::Named("Z")=clone(Z), Rcpp::Named("X")=clone(X) ) ;
	drawcount ++; 
      }       
    } 
  } 
  if (verbosity>=1) 
      cout << endl ;  
  
  return savedGibbsIterations; 

}
