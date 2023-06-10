#ifndef _BIOMICO_H
#define _BIOMICO_H

#include <Rcpp.h>
#include <vector>
#include <list>

using namespace Rcpp;
using namespace std;

/*                                                                                                                                                                                  
 * note : RcppExport is an alias to `extern "C"` defined by Rcpp.                                                                                                                   
 *                                                                                                                                                                                  
 * It gives C calling convention to the rcpp_hello_world function so that                                                                                                           
 * it can be called from .Call in R. Otherwise, the C++ compiler mangles the                                                                                                        
 * name of the function and .Call can't find it.                                                                                                                                    
 *                                                                                                                                                                                  
 * It is only useful to use RcppExport when the function is intended to be called                                                                                                   
 * by .Call. See the thread http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672                                                                                            
 * on Rcpp-devel for a misuse of RcppExport                                                                                                                                         
 */
RcppExport SEXP runGibbs(SEXP Rdata_samples, SEXP Rdata_envs, SEXP RG, SEXP RV, SEXP RT, SEXP RN, SEXP RunknownEnv,
			 SEXP RZ, SEXP RX, 
			 SEXP Rburnin, SEXP Rnrestarts, SEXP Rndraws_per_restart, SEXP Rdelay, 
			 SEXP Ralpha_pi, SEXP Ralpha_phi, SEXP Ralpha_theta, 
			 SEXP Rmaxdepth, SEXP Rverbosity, SEXP Rprinting_index, SEXP Rprinting_total);

RcppExport SEXP runGibbsOnTest(SEXP Rdata_samples, SEXP RG, SEXP RV, SEXP RT, SEXP RN, SEXP RunknownEnv,
			 SEXP RCount_T_G, SEXP RCount_G_V,   
			 SEXP Rburnin, SEXP Rnrestarts, SEXP Rndraws_per_restart, SEXP Rdelay,  
			 SEXP Ralpha_pi, SEXP Ralpha_phi, SEXP Ralpha_theta, SEXP Rmaxdepth,  
			       SEXP Rverbosity, SEXP Rprinting_index, SEXP Rprinting_total);

#endif
