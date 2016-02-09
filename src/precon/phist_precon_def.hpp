// this file maps stuff implemented in a C++ traits class (e.g. PreconTraits<ST,IFPACK>::Apply) to plain C 
// (e.g. void (*apply)(...) in a linearOp_t struct. It therefore instantiates all the supported PreconTraits
// templates and selects them with the following macro
#ifdef SELECT_PT_MEMBER
#undef SELECT_PT_MEMBER
#endif

#ifdef CALL_PT_MEMBER
#undef CALL_PT_MEMBER
#endif


#define SELECT_PT_MEMBER(PRECON_TYPE,MEMBER) \
        PRECON_TYPE==NONE? phist::PreconTraits<_ST_,NONE>:: ## MEMBER: \
        PRECON_TYPE==IFPACK? phist::PreconTraits<_ST_,IFPACK>:: ## MEMBER: \
        PRECON_TYPE==ML? phist::PreconTraits<_ST_,ML>:: ## MEMBER: \
        PRECON_TYPE==IFPACK2? phist::PreconTraits<_ST_,IFPACK2>:: ## MEMBER: \
        PRECON_TYPE==MUELU? phist::PreconTraits<_ST_,MUELU>:: ## MEMBER: \
        PRECON_TYPE==AMESOS2? phist::PreconTraits<_ST_,AMESOS2>:: ## MEMBER: \
        NULL;

#define CALL_PT_MEMBER(PRECON_TYPE,MEMBER,...) \
        if (PRECON_TYPE==NONE) PHIST_CHK_IERR(phist::PreconTraits<_ST_,NONE>:: ## MEMBER(__VA_ARGS__),*iflag); \
        else if (PRECON_TYPE==IFPACK) PHIST_CHK_IERR(phist::PreconTraits<_ST_,IFPACK>:: ## MEMBER(__VA_ARGS__),*iflag): \
        else if (PRECON_TYPE==ML) PHIST_CHK_IERR(phist::PreconTraits<_ST_,ML>:: ## MEMBER(__VA_ARGS__),*iflag): \
        else if (PRECON_TYPE==IFPACK2) PHIST_CHK_IERR(phist::PreconTraits<_ST_,IFPACK2>:: ## MEMBER(__VA_ARGS__),*iflag): \
        else if (PRECON_TYPE==MUELU) PHIST_CHK_IERR(phist::PreconTraits<_ST_,MUELU>:: ## MEMBER(__VA_ARGS__),*iflag): \
        else if(PRECON_TYPE==AMESOS2) PHIST_CHK_IERR(phist::PreconTraits<_ST_,AMESOS2>:: ## MEMBER(__VA_ARGS__),*iflag): \
        else PHIST_CHK_IERR(*iflag=PHIST_INVALID_INPUT,*iflag);

// create a preconditioner for an iterative linear solver

// this function can be used to create an operator that can be used to precondition linear systems     
// with A-sigma*B. If sigma is 0, B is not touched. If sigma!=0 but B==NULL, B=I (identity matrix)     
// is assumed.                                                                                         
//
// Applying the preconditioner is done via the usual apply member functions of the linearOp struct.    
//
// For some preconditioners it may be useful to provide (an approximation of) the kernel of            
// A-sigma*B, this can be done via Vkern and BVkern. If they are NULL, they are not used anyway.       
//
// The preconditioning method is selected based on the string <method>. Options are                    
// passed via the options string. The methods available depend on the kernel lib used and the          
// optional libraries found while building phist. Passing in method="usage". The option string         
// in turn depends on the input <method>. Passing in a supported method and options="usage" prints     
// usage information for that particular preconditioner.                                               
//                                                                                                     
// The special string method="update" can be used if a preconditioner has already                      
// been created and only the values (but not the pattern) of A, B and sigma changed.                   
// method="update shift" indicates that only the value of sigma changed. In both cases, the options    
// string is ignored.                                                                                  
//
// Example:                                                                                            
//                                                                                                     
// If you use the Epetra kernel library and the Trilinos library Ifpack (incomplete factorization      
// preconditioners) is available, you could use                                                        
//                                                                                                     
// DlinearOp_t P;                                                                                      
// phist_Dprecon_create(&P, &A, 0, NULL, NULL, NULL, "ifpack", "ifpack_params.xml", &iflag);           
//                                                                                                     
// which will read the preconditioner settings from an XML file compatible with the Teuchos            
// ParameterList.                                                                                      
//                                                                                                     
extern "C" void SUBR(precon_create)(TYPE(linearOp_ptr) op, TYPE(const_sparseMat_ptr) A, 
                         _ST_ sigma, TYPE(const_sparseMat_ptr) B,
                         TYPE(const_mvec_ptr) Vkern, TYPE(const_mvec_ptr) BVkern,
                         const char* method, const char* options, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  
  if (!strcasecmp(method,"usage"))
  {
    PHIST_SOUT(PHIST_ERROR,"Your PHIST installation supports the following preconditioners:\n");
#ifdef PHIST_KERNEL_LIB_EPETRA
# ifdef PHIST_HAVE_IFPACK
        PHIST_SOUT(PHIST_ERROR,"\t'%s'\n",precon2str(IFPACK));
# endif
# ifdef PHIST_HAVE_ML
        PHIST_SOUT(PHIST_ERROR,"\t'%s'\n",precon2str(ML));
# endif
# ifdef PHIST_HAVE_MUELU
        PHIST_SOUT(PHIST_ERROR,"\t'%s'\n",precon2str(MUELU));
# endif
#elif defined(PHIST_KERNEL_LIB_TPETRA
# ifdef PHIST_HAVE_IFPACK2
        PHIST_SOUT(PHIST_ERROR,"\t'%s'\n",precon2str(IFPACK2));
# endif
# ifdef PHIST_HAVE_MUELU
        PHIST_SOUT(PHIST_ERROR,"\t'%s'\n",precon2str(MUELU));
# endif
# ifdef PHIST_HAVE_AMESOS2
        PHIST_SOUT(PHIST_ERROR,"\t'%s'\n",precon2str(AMESOS2));
# endif
#else
        PHIST_SOUT(PHIST_ERROR,"\t'%s'\n",p);
#endif
    *iflag=0;
    return;
  }
  
  phist_internal_precon_t* pt = new phist_internal_precon_t;
  
  precon_t precType=str2precon(method);
  pt->type_ = precType;
  
  if (precType==INVALID_PRECON_T)
  {
    PHIST_SOUT(PHIST_ERROR,"your given precon type '%s' was not recognized,\n"
                           "please check the spelling and if the required TPLs are\n"
                           "available in PHIST (cf. phist_config.h)\n",method);
    SUBR(precon_create)(NULL,NULL,sigma,NULL,NULL,NULL,"usage",NULL,iflag);
    *iflag=PHIST_INVALID_INPUT;
    return;
  }
  
  PHIST_CHK_IERR( CALL_PT_MEMBER(precType, Create,&pt->P_,A,sigma,B,Vkern,BVkern,options,iflag), *iflag);

  op->A_=pt;
  op->apply = SELECT_PT_MEMBER(precType,Apply);
  op->applyT = SELECT_PT_MEMBER(precType,ApplyT);
  op->apply_shifted = SELECT_PT_MEMBER(precType,ApplyShifted);

}

// destroy preconditioner
extern "C" void SUBR(precon_delete)(TYPE(linearOp_ptr) op, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  CAST_PTR_FROM_VOID(phist_internal_precon_t, pt, op->A_,*iflag);
  precon_t precType=pt->type_;
  PHIST_CHK_IERR( CALL_PT_MEMBER(precType,Delete,pt,iflag), *iflag);
}

//@}



