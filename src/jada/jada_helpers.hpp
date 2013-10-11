// templated helper functions to sort arrays of eigenvalues

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
void SortEig(std::complex<MT>* ev,int n,int* idx,sortEig_t which,int* ierr)
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
    if (std::imag(ev[i])!=mt::zero()) idx[i]*=-1;
    }
  }

// swap two rows and columns in the Schur form T, taking the Schur vectors in S along.
// The indices i1 and i2 should be determined using SortEig, that subroutine marks 
// complex eigenvalues by negative indices which allows this subroutine to function
// even for complex pairs of real matrices.
template<typename ST>
void SwapSchur(int m, ST* T, lidx_t ldT, ST* S, lidx_t ldS,
        int i1, int i2)
  {
  *ierr=-99;
  }

// specialization for the complex case, which is simpler because we
// store all eigenvalues in the same type as the matrices T and S.
template<typename ST>
void SwapSchur(int m, ST* T, lidx_t ldT, ST* S, lidx_t ldS,
        int i1, int i2)
  {
  // we don't mind if it is a real or complex eigenvalue:
  i1=std::abs(i1); i2=std::abs(i2);
 
 //TODO - change this MATLAB code to C/LAPACK 
  
  kk=min(length(I),size(S,1)-1);
  j=1; 
while (j<=kk & j==I(j))
  {
  j=j+1; 
  while (j<=kk)
    {
    i=I(j);
    for (k=i-1:-1:j)
      {
      q = [S(k,k)-S(k+1,k+1),S(k,k+1)];
      if (q(1) != 0)
        {
        q = q/norm(q);
        G = [[q(2);-q(1)],q'];
        J = [k,k+1];
        Q(:,J) = Q(:,J)*G;
        S(:,J) = S(:,J)*G;
        S(J,:) = G'*S(J,:);
        }
      S(k+1,k) = 0;
      }
    I=I+(I<i);
    j=j+1; 
    while (j<=kk & j==I(j))
      {
      j=j+1; 
      }
    }
  }
}
