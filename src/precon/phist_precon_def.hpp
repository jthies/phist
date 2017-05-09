/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
namespace phist 
{

template<> void PreconTraits<_ST_,phist_NO_PRECON>::Apply
                                         (_ST_ alpha, void const* P, TYPE(const_mvec_ptr) X, 
                                          _ST_ beta,                 TYPE(mvec_ptr)       Y, 
                                            int* iflag)
  {
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(alpha,X,beta,Y,iflag),*iflag);
    return;
  }
  template<> void PreconTraits<_ST_,phist_NO_PRECON>::ApplyShifted
         (_ST_ alpha, const void* P, _ST_ const * sigma, TYPE(const_mvec_ptr) X, 
          _ST_ beta,                                     TYPE(mvec_ptr) Y, int* iflag)
  {
    int nvec;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&nvec,iflag),*iflag);
    _ST_ alphas[nvec];
    for (int i=0; i<nvec;i++) alphas[i]=alpha*(st::one()-sigma[i]);
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(alphas,X,beta,Y,iflag),*iflag);
    return;
  }
}//namespace phist
  
// this file maps stuff implemented in a C++ traits class (e.g. PreconTraits<ST,phist_IFPACK>::Apply) to plain C 
// (e.g. void (*apply)(...) in a linearOp_t struct. It therefore instantiates all the supported PreconTraits
// templates and selects them with the following macro
#ifdef SELECT_PT_MEMBER
#undef SELECT_PT_MEMBER
#endif

#ifdef CALL_PT_MEMBER
#undef CALL_PT_MEMBER
#endif


#define SELECT_PT_MEMBER(PRECON_TYPE,MEMBER) \
        PRECON_TYPE==phist_NO_PRECON? phist::PreconTraits<_ST_,phist_NO_PRECON>::MEMBER: \
        PRECON_TYPE==phist_IFPACK? phist::PreconTraits<_ST_,phist_IFPACK>::MEMBER: \
        PRECON_TYPE==phist_ML? phist::PreconTraits<_ST_,phist_ML>::MEMBER: \
        PRECON_TYPE==phist_MUELU? phist::PreconTraits<_ST_,phist_MUELU>::MEMBER: \
        PRECON_TYPE==phist_AMESOS2? phist::PreconTraits<_ST_,phist_AMESOS2>::MEMBER: \
        PRECON_TYPE==phist_USER_PRECON? phist::PreconTraits<_ST_,phist_USER_PRECON>::MEMBER: \
        NULL;

#define CALL_PT_MEMBER(PRECON_TYPE,MEMBER,...) \
        if (PRECON_TYPE==phist_NO_PRECON) PHIST_CHK_IERR((phist::PreconTraits<_ST_,phist_NO_PRECON>::MEMBER)(__VA_ARGS__),*iflag) \
        else if (PRECON_TYPE==phist_IFPACK) PHIST_CHK_IERR((phist::PreconTraits<_ST_,phist_IFPACK>::MEMBER)(__VA_ARGS__),*iflag) \
        else if (PRECON_TYPE==phist_ML) PHIST_CHK_IERR((phist::PreconTraits<_ST_,phist_ML>::MEMBER)(__VA_ARGS__),*iflag) \
        else if (PRECON_TYPE==phist_MUELU) PHIST_CHK_IERR((phist::PreconTraits<_ST_,phist_MUELU>::MEMBER)(__VA_ARGS__),*iflag) \
        else if(PRECON_TYPE==phist_AMESOS2) PHIST_CHK_IERR((phist::PreconTraits<_ST_,phist_AMESOS2>::MEMBER)(__VA_ARGS__),*iflag) \
        else if(PRECON_TYPE==phist_USER_PRECON) PHIST_CHK_IERR((phist::PreconTraits<_ST_,phist_USER_PRECON>::MEMBER)(__VA_ARGS__),*iflag) \
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
// The string method="usage" will lead to a list of available preconditioners being printed.
// The string method="user_defined" can be used if an application provides its own precondi-           
// tioning. In that case, the class phist::PreconTraits<_ST_,phist_USER_PRECON> should be implemented  
// by the application. The "last_arg" pointer can be used to provide the actual preconditioner object  
// to the PreconTraits::Create member function.                                                        
//
// Example:                                                                                            
//                                                                                                     
// If you use the Epetra kernel library and the Trilinos library Ifpack (incomplete factorization      
// preconditioners) is available, you could use                                                        
//                                                                                                     
// DlinearOp_t P;                                                                                      
// phist_Dprecon_create(&P, &A, 0, NULL, NULL, NULL, "ifpack", "ifpack_params.xml", NULL, &iflag);     
//                                                                                                     
// which will read the preconditioner settings from an XML file compatible with the Teuchos            
// ParameterList.                                                                                      
//                                                                                                     
extern "C" void SUBR(precon_create)(TYPE(linearOp_ptr) op, TYPE(const_sparseMat_ptr) A, 
                         _ST_ sigma, TYPE(const_sparseMat_ptr) B,
                         TYPE(const_mvec_ptr) Vkern, TYPE(const_mvec_ptr) BVkern,
                         const char* method, const char* options, 
                         void* last_arg, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  
  if (!strcasecmp(method,"usage"))
  {
    PHIST_SOUT(PHIST_ERROR,"Your PHIST installation supports the following preconditioners:\n");
    PHIST_SOUT(PHIST_ERROR,"\t'%s'\n",precon2str(phist_NO_PRECON));
#ifdef PHIST_KERNEL_LIB_EPETRA
# ifdef PHIST_HAVE_IFPACK
        PHIST_SOUT(PHIST_ERROR,"\t'%s'\n",precon2str(phist_IFPACK));
# endif
# ifdef PHIST_HAVE_ML
        PHIST_SOUT(PHIST_ERROR,"\t'%s'\n",precon2str(phist_ML));
# endif
# ifdef PHIST_HAVE_MUELU
        PHIST_SOUT(PHIST_ERROR,"\t'%s'\n",precon2str(phist_MUELU));
# endif
#elif defined(PHIST_KERNEL_LIB_TPETRA)
# ifdef PHIST_HAVE_IFPACK2
        PHIST_SOUT(PHIST_ERROR,"\t'%s'\n",precon2str(phist_IFPACK));
# endif
# ifdef PHIST_HAVE_MUELU
        PHIST_SOUT(PHIST_ERROR,"\t'%s'\n",precon2str(phist_MUELU));
# endif
# ifdef PHIST_HAVE_AMESOS2
        PHIST_SOUT(PHIST_ERROR,"\t'%s'\n",precon2str(phist_AMESOS2));
# endif
// TODO
//#else
//        PHIST_SOUT(PHIST_ERROR,"\t'%s'\n",p);
#endif
    PHIST_SOUT(PHIST_ERROR,"In order to wrap an existing preconditioner of the corresponding type, \n"
                           "you can pass in options=NULL and a pointer to the existing object via last_arg\n");
    *iflag=0;
    return;
  }
  
  phist_Eprecon precType=str2precon(method);

  if (options!=NULL && !strcasecmp(options,"usage"))
  {
    CALL_PT_MEMBER(precType,Usage);
    return;
  }

  phist_internal_precon* pt = new phist_internal_precon;
  
  pt->type_ = precType;
  
  if (precType==phist_INVALID_PRECON)
  {
    PHIST_SOUT(PHIST_ERROR,"your given precon type '%s' was not recognized,\n"
                           "please check the spelling and if the required TPLs are\n"
                           "available in PHIST (cf. phist_config.h)\n",method);
    SUBR(precon_create)(NULL,NULL,sigma,NULL,NULL,NULL,"usage",NULL,NULL,iflag);
    *iflag=PHIST_INVALID_INPUT;
    return;
  }
  
  if (options==NULL && last_arg!=NULL)
  {
    CALL_PT_MEMBER(precType,Wrap,&pt->P_,A,sigma,B,Vkern,BVkern,last_arg,iflag);
  }
  else
  {
    CALL_PT_MEMBER(precType,Create,&pt->P_,A,sigma,B,Vkern,BVkern,options,last_arg,iflag);
  }
  
  pt->A_=A;
  pt->B_=B;
  pt->Vkern_=Vkern;
  pt->BVkern_=BVkern;
  
  op->A=pt;
  op->apply = SUBR(precon_apply);
  op->applyT = SUBR(precon_applyT);
  op->apply_shifted = SUBR(precon_apply_shifted);
  PHIST_CHK_IERR(SUBR(sparseMat_get_range_map)(A,&op->range_map,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sparseMat_get_domain_map)(A,&op->domain_map,iflag),*iflag);

}

// destroy preconditioner
extern "C" void SUBR(precon_delete)(TYPE(linearOp_ptr) op, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(phist_internal_precon, pt, op->A,*iflag);
  phist_Eprecon precType=pt->type_;
  CALL_PT_MEMBER(precType,Delete,pt->P_,iflag);
  delete pt;
}

// given an existing preconditioner, recompute it for a new shift sigma and (near) kernel Vkern.     

// The matrices A and B and preconditioning options remain unchanged from the cann SUBR(precon_create).
// This means that the preconditioner needs to store pointers to A and B at create() time and can assume
// those matrices are still there and unchanged. If they aren't (there or unchanged), a new preconditioner
// must be created instead.
extern "C" void SUBR(precon_update)(TYPE(linearOp_ptr) op, _ST_ sigma,
                         TYPE(const_mvec_ptr) Vkern,
                         TYPE(const_mvec_ptr) BVkern,
                         int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(phist_internal_precon, pt, op->A,*iflag);
  phist_Eprecon precType=pt->type_;
  pt->Vkern_=Vkern;
  pt->BVkern_=BVkern;
  CALL_PT_MEMBER(precType,Update,pt->P_,pt->A_,sigma,pt->B_,Vkern,BVkern,iflag);
}

// apply preconditioner
extern "C" void SUBR(precon_apply)(_ST_ alpha, void const* vP, TYPE(const_mvec_ptr) X, 
                                   _ST_ beta,  TYPE(mvec_ptr) Y,int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(phist_internal_precon, pt, vP,*iflag);
  phist_Eprecon precType=pt->type_;
  CALL_PT_MEMBER(precType,Apply,alpha,pt->P_,X,beta,Y,iflag);
}

// apply preconditioner
extern "C" void SUBR(precon_applyT)(_ST_ alpha, void const* vP, TYPE(const_mvec_ptr) X, 
                                   _ST_ beta,  TYPE(mvec_ptr) Y,int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(phist_internal_precon, pt, vP,*iflag);
  phist_Eprecon precType=pt->type_;
  CALL_PT_MEMBER(precType,ApplyT,alpha,pt->P_,X,beta,Y,iflag);
}

// apply preconditioner
extern "C" void SUBR(precon_apply_shifted)(_ST_ alpha, void const* vP, _ST_ const* sigma, TYPE(const_mvec_ptr) X, 
                                   _ST_ beta,  TYPE(mvec_ptr) Y,int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(phist_internal_precon, pt, vP,*iflag);
  phist_Eprecon precType=pt->type_;
  CALL_PT_MEMBER(precType,ApplyShifted,alpha,pt->P_,sigma,X,beta,Y,iflag);
}

//@}



