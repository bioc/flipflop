/* exporting methods */
#if (__GNUC__ >= 4) || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
#  ifndef GCC_HASCLASSVISIBILITY
#    define GCC_HASCLASSVISIBILITY
#  endif
#endif

#ifndef SWIGEXPORT
# if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
#   if defined(STATIC_LINKED)
#     define SWIGEXPORT
#   else
#     define SWIGEXPORT __declspec(dllexport)
#   endif
# else
#   if defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
#     define SWIGEXPORT __attribute__ ((visibility("default")))
#   else
#     define SWIGEXPORT
#   endif
# endif
#endif


/* attribute recognised by some compilers to avoid 'unused' warnings */
/*
#ifndef SWIGUNUSED
# if defined(__GNUC__)
#   if !(defined(__cplusplus)) || (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4))
#     define SWIGUNUSED __attribute__ ((__unused__)) 
#   else
#     define SWIGUNUSED
#   endif
# elif defined(__ICC)
#   define SWIGUNUSED __attribute__ ((__unused__)) 
# else
#   define SWIGUNUSED 
# endif
#endif

#ifndef SWIGINTERN
//# define SWIGINTERN static SWIGUNUSED
# define SWIGINTERN SWIGUNUSED
#endif

#ifndef SWIGRUNTIME
# define SWIGRUNTIME SWIGINTERN
#endif
*/

int ffProcesssam(char * in, char * pre, char * an);

//SWIGRUNTIME SEXP R_SWIG_debug_getCallbackFunctionData();

//SWIGRUNTIME SEXP R_SWIG_R_pushCallbackFunctionData(SEXP fun, SEXP userData);

SWIGEXPORT SEXP R_swig_multLeftDiag ( SEXP X, SEXP Y, SEXP XY);

SWIGEXPORT SEXP R_swig_fistaFlat ( SEXP X, SEXP D, SEXP alpha0, SEXP alpha, SEXP num_threads, SEXP max_it, SEXP L0, SEXP fixed_step, SEXP gamma, SEXP s_lambda, SEXP delta, SEXP lambda2, SEXP lambda3, SEXP a, SEXP b, SEXP c, SEXP tol, SEXP it0, SEXP max_iter_backtracking, SEXP compute_gram, SEXP lin_admm, SEXP admm, SEXP intercept, SEXP resetflow, SEXP name_regul, SEXP name_loss, SEXP verbose, SEXP pos, SEXP clever, SEXP log, SEXP ista, SEXP subgrad, SEXP logName, SEXP is_inner_weights, SEXP inner_weights, SEXP size_group, SEXP sqrt_step, SEXP transpose);

SWIGEXPORT SEXP R_swig_evalPathCoding ( SEXP alpha0, SEXP dual_val, SEXP precision, SEXP weights, SEXP ir, SEXP jc, SEXP start_weights, SEXP stop_weights, SEXP num_threads, SEXP lambda1, SEXP lambda2, SEXP intercept, SEXP resetflow, SEXP name_regul, SEXP verbose, SEXP pos, SEXP clever, SEXP eval, SEXP eval_dual, SEXP size_group, SEXP transpose);

SWIGEXPORT SEXP R_swig_sepCostsPathCoding ( SEXP alpha0, SEXP alpha, SEXP weights, SEXP ir, SEXP jc, SEXP start_weights, SEXP stop_weights, SEXP max_capacity, SEXP epsilon_flow, SEXP prices, SEXP num_threads, SEXP lambda, SEXP tol, SEXP delta, SEXP loss_weights, SEXP name_regul, SEXP name_loss, SEXP pos, SEXP mode_decomposition);

SWIGEXPORT SEXP R_swig_solverPoisson ( SEXP y, SEXP X, SEXP beta0, SEXP beta, SEXP weights, SEXP delta, SEXP max_iter, SEXP tol);

SWIGEXPORT SEXP R_swig_solverPoissonFull ( SEXP y, SEXP X, SEXP beta0, SEXP beta, SEXP weights, SEXP delta, SEXP max_iter, SEXP tol);
