#include <cstdlib>
#include <algorithm>
#include <vector>
#include <complex>
#include "phist_enums.h"

template<typename T>
bool SelectLM(std::pair<T,int> a_, std::pair<T,int> b_)
  {
  T& a = a_.first;
  T& b = b_.first; 
  if (std::abs(a)==std::abs(b))
    {
    if (std::real(a)==std::real(b))
      {
      // conj pair, positive imag first
      return std::imag(a)>std::imag(b);
      }
    else
      {
      return std::real(a)>std::real(b);
      }
    }
  return std::abs(a)>std::abs(b);
  }

template<typename T>
bool SelectSM(std::pair<T,int> a_, std::pair<T,int> b_)
  {
  T& a = a_.first;
  T& b = b_.first; 
  if (std::abs(a)==std::abs(b))
    {
    if (std::real(a)==std::real(b))
      {
      // conj pair, positive imag first
      return std::imag(a)>std::imag(b);
      }
    else
      {
      return std::real(a)<std::real(b);
      }
    }
  return std::abs(a)<std::abs(b);
  }

template<typename T>
bool SelectLR(std::pair<T,int> a_, std::pair<T,int> b_)
  {
  T& a = a_.first;
  T& b = b_.first; 
  return std::real(a)>std::real(b);
  }

template<typename T>
bool SelectSR(std::pair<T,int> a_, std::pair<T,int> b_)
  {
  T& a = a_.first;
  T& b = b_.first; 
  return std::real(a)<std::real(b);
  }



template<typename MT>
void SortEig(std::complex<MT>* ev,int n,int* idx,eigSort_t which,int* ierr)
  {
  typedef std::complex<MT> ST;
  typedef std::pair<ST,int> PT;
  std::vector<PT> v(n);
  for (int i=0;i<n;i++)
    {
    v[i]=PT(ev[i],i);
    }

  if (which==LM)
    {
    std::sort(v.begin(),v.end(),SelectLM<ST>);
    }
  else if (which==SM)
    {
    std::sort(v.begin(),v.end(),SelectSM<ST>);
    }
  else if (which==LR)
    {
    std::sort(v.begin(),v.end(),SelectLR<ST>);
    }
  else if (which==SR)
    {
    std::sort(v.begin(),v.end(),SelectSR<ST>);
    }
  else if (which!=NONE)
    {
    // sort type not implemented
    *ierr=-99;
    }

  for (int i=0;i<n;i++)
    {
    ev[i]=v[i].first;
    idx[i]=v[i].second;
    if (std::imag(ev[i])!=(MT)0.0) idx[i]*=-1;
    }
  }

