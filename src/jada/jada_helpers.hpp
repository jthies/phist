/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! \file jada_helpers.hpp
//! \brief compare functors with tolerance

#include "phist_config.h"

#ifndef DOXYGEN

#include <cstdlib>
#include <algorithm>
#include <vector>
#include <complex>
#include "phist_enums.h"
#include "phist_ScalarTraits.hpp"

#endif

// compare functors with tolerance

template<typename T>
class SelectLM
{
  public:
  
    typedef typename phist::ScalarTraits<T>::magn_t MT;
  
    SelectLM(MT tol) : tol_(tol) {}

    bool operator()(std::pair<T,int> a_, std::pair<T,int> b_)
    {
      T& a = a_.first;
      T& b = b_.first; 
      if ( std::abs(std::abs(a)-std::abs(b)) < tol_ )
      {
        if ( std::abs(std::real(a)-std::real(b)) < tol_ )
        {
          // possibly conj pair, positive imag first
          return std::imag(a)>std::imag(b)+tol_;
        }
        else
        {
          return std::real(a)>std::real(b);
        }
      }
      return std::abs(a)>std::abs(b);
    }

  private:
    MT tol_;
};

template<typename T>
class SelectSM
{
  public:

    typedef typename phist::ScalarTraits<T>::magn_t MT;

    SelectSM(MT tol) : tol_(tol) {}

    bool operator()(std::pair<T,int> a_, std::pair<T,int> b_)
    {
      T& a = a_.first;
      T& b = b_.first; 
      if ( std::abs(std::abs(a)-std::abs(b)) < tol_ )
      {
        if ( std::abs(std::real(a)-std::real(b)) < tol_ )
        {
          // conj pair, positive imag first
          return std::imag(a)>std::imag(b)+tol_;
        }
        else
        {
          return std::real(a)<std::real(b);
        }
      }
      return std::abs(a)<std::abs(b);
    }

  private:
    MT tol_;
};

template<typename T>
class SelectLR
{
  public:

    typedef typename phist::ScalarTraits<T>::magn_t MT;

    SelectLR(MT tol) : tol_(tol) {}

    bool operator()(std::pair<T,int> a_, std::pair<T,int> b_)
    {
      T& a = a_.first;
      T& b = b_.first; 
      return std::real(a)>std::real(b)+tol_;
    }

  private:
    MT tol_;
};

template<typename T>
class SelectSR
{
  public:

    typedef typename phist::ScalarTraits<T>::magn_t MT;

    SelectSR(MT tol) : tol_(tol) {}

    bool operator()(std::pair<T,int> a_, std::pair<T,int> b_)
    {
      T& a = a_.first;
      T& b = b_.first; 
      return std::real(a)<std::real(b)-tol_;
    }

  private:
    MT tol_;
};



template<typename MT>
void SortEig(std::complex<MT>* ev,int n,int* idx,phist_EeigSort which, MT tol, int* iflag)
{
  typedef std::complex<MT> ST;
  typedef std::pair<ST,int> PT;
  std::vector<PT> v(n);
  for (int i=0;i<n;i++)
    v[i]=PT(ev[i],i);

  // we need to use tol/n because all fp values would be conisidered equal in (0 tol, 2*tol, 3*tol) and others
  // and thats problably not what we want!
  tol = tol/n;

  if (which==phist_LM)
    std::stable_sort(v.begin(),v.end(),SelectLM<ST>(tol));
  else if (which==phist_SM)
    std::stable_sort(v.begin(),v.end(),SelectSM<ST>(tol));
  else if (which==phist_LR)
    std::stable_sort(v.begin(),v.end(),SelectLR<ST>(tol));
  else if (which==phist_SR)
    std::stable_sort(v.begin(),v.end(),SelectSR<ST>(tol));
  else if (which!=phist_NO_EIGSORT)
  {
    // sort type not implemented
    *iflag=PHIST_NOT_IMPLEMENTED;
  }

  for (int i=0;i<n;i++)
  {
    ev[i]=v[i].first;
    idx[i]=v[i].second;
    if (std::imag(ev[i])!=(MT)0.0)
      idx[i]*=-1;
  }
}

