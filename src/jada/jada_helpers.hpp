#include "phist_config.h"

#ifndef DOXYGEN

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <complex>
#include "phist_enums.h"

#endif

// compare functors with tolerance

template<typename T>
class SelectLM
{
  public:
    SelectLM(_MT_ tol) : tol_(tol) {}

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
    _MT_ tol_;
};

template<typename T>
class SelectSM
{
  public:
    SelectSM(_MT_ tol) : tol_(tol) {}

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
    _MT_ tol_;
};

template<typename T>
class SelectLR
{
  public:
    SelectLR(_MT_ tol) : tol_(tol) {}

    bool operator()(std::pair<T,int> a_, std::pair<T,int> b_)
    {
      T& a = a_.first;
      T& b = b_.first; 
      return std::real(a)>std::real(b)+tol_;
    }

  private:
    _MT_ tol_;
};

template<typename T>
class SelectSR
{
  public:
    SelectSR(_MT_ tol) : tol_(tol) {}

    bool operator()(std::pair<T,int> a_, std::pair<T,int> b_)
    {
      T& a = a_.first;
      T& b = b_.first; 
      return std::real(a)<std::real(b)-tol_;
    }

  private:
    _MT_ tol_;
};



template<typename MT>
void SortEig(std::complex<MT>* ev,int n,int* idx,phist_EeigSort which, _MT_ tol, int* iflag)
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

