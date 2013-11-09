function [R,Q,resnorm,lambda_history,resnorm_history,m,restarts,lambda_iter,totalMatVecs,correctionRatio]=melven_subspace_jada(A,B,v0,nEig,arnoldiIter,maxIter,minBas,maxBas,res_eps,maxDeflationVecs)

totalMatVecs = 0;
n = size(A,1);                                  % matrix dimension
k = size(v0,2);                                 % blocksize, e.g. number of eigenvalues wanted
totalCorrections = 0;
totalNewCorrections = 0;

V = zeros(n,maxBas);                         % search space empty initially
BV = zeros(n,maxBas);
nV = 0;
W = zeros(n,maxBas);
t = v0;                                         % initial new direction
H = zeros(maxBas);                           % preallocate H
restarts = 0;

%data for already converged eigenvalues
resnorm = inf(nEig,1);
resnorm_history = inf(nEig,arnoldiIter+maxIter);
lambda_history = nan(nEig,arnoldiIter+maxIter);
lambda_iter = zeros(nEig,1);
nConv = 0;
R = zeros(nEig,nEig);
Q = zeros(n,nEig);
BQ = zeros(n,nEig);
nNewConv = 0;
firstDeflationVec = 0;

for m = 1:maxIter+arnoldiIter                               % main iteration loop
    nV_old = nV;
    if( m > 1 )
        if( firstDeflationVec > 1 )
            t = t - Q(:,1:firstDeflationVec-1) * (BQ(:,1:firstDeflationVec-1)' * t);
        end
        t = t - V(:,1:nV_old) * (BV(:,1:nV_old)' * t);     % orthogonalize t wrt V
    end
totalCorrections = totalCorrections + size(t,2);
    [t,rank_t] = orth(t,B,res_eps/100); % need B-orthogonality
totalNewCorrections = totalNewCorrections + rank_t;
correctionRatio = totalNewCorrections/totalCorrections;
    fprintf(' rank of correction: %d (previously converged EV: %d)\n', rank_t, nNewConv);
    if( rank_t == 0 && nNewConv == 0 )
        warning('using new random vector!');
        t = rand(n,1);
        t = t - Q(:,1:nConv)*BQ(:,1:nConv)'*t;
        t = t - V(:,1:nV)*BV(:,1:nV)'*t;
        [t,rank_t] = orth(t,B,res_eps/100); % need B-orthogonality
        fprintf(' rank of random correction: %d\n', rank_t);
    end
    if( nV_old + rank_t > maxBas )
      restarts = restarts + 1;
      fprintf('restarting with %d basis vectors after subspace dimension of %d was reached.\n', minBas, nV_old + rank_t);
      % we need to get the "best" vectors out of V
      V(:,1:minBas) = V(:,1:nV_old)*Q_H(:,1:minBas);
      % also update BV, W and H
      BV(:,1:minBas) = BV(:,1:nV_old)*Q_H(:,1:minBas);
      W(:,1:minBas) = W(:,1:nV_old)*Q_H(:,1:minBas);
      H(1:minBas,1:minBas) = Q_H(:,1:minBas)' * H(1:nV_old,1:nV_old) * Q_H(:,1:minBas);
      nV_old = minBas;
    end
    nV = nV_old + rank_t;
    V(:,nV_old+1:nV) = t;                     % update subspace V
    BV(:,nV_old+1:nV) = B*t;
V_orth = max(max(abs(eye(nV) - BV(:,1:nV)' * V(:,1:nV))));
if( V_orth > res_eps/100 )
    %% reorthogonalize
    warning('V not orthogonal (%g), reorthogonalizing', V_orth);
    V(:,1:nV) = V(:,1:nV) - Q(:,1:nConv)*BQ(:,1:nConv)' * V(:,1:nV);
    [t,nV] = orth(V(:,1:nV),B,res_eps);
    nV_old = 0;
    V(:,1:nV) = t;
    BV(:,1:nV) = B*t;
    V_orth = max(max(abs(eye(nV) - BV(:,1:nV)' * t)));
    if( V_orth > res_eps )
        error('could not reorthogonalize, error: %g', V_orth);
    end
end
    W(:,nV_old+1:nV) = A*t;                   % update W = A*V
totalMatVecs = totalMatVecs + rank_t;


    % update H = V' * W
    H(1:nV_old,    nV_old+1:nV)  = V(:,1:nV_old)'    * W(:,nV_old+1:nV);
    H(nV_old+1:nV, 1:nV_old)     = V(:,nV_old+1:nV)' * W(:,1:nV_old);
    H(nV_old+1:nV, nV_old+1:nV)  = V(:,nV_old+1:nV)' * W(:,nV_old+1:nV);

    % compute the schur form of H
    [Q_H,R_H] = SortSchur(H(1:nV,1:nV),'LM',1,nV);

%test_H = max(max(abs(H(1:nV,1:nV)*Q_H - Q_H*R_H)));
%if( test_H > res_eps )
%    error('test Q_H R_H wrong: %e', test_H);
%end

    
    % check for converged eigenvalues
    r = R_H(1:k,1:k);
    q = V(:,1:nV) * Q_H(:,1:k);
    Bq = BV(:,1:nV) * Q_H(:,1:k);
    % calculate residuum
    Aq = W(:,1:nV) * Q_H(:,1:k);
    res = Aq - Bq*r;
    a_ = BQ(:,1:nConv)'*res;
    res = res - Q(:,1:nConv)*a_;
    for i = 1:k
        resnorm(nConv+i) = norm(res(:,i),2);
    end
    lambda_history(1:nConv,m) = diag(R(1:nConv,1:nConv));
    lambda_history(nConv+1:nConv+k,m) = diag(r);
    resnorm_history(:,m) = resnorm;
    nNewConv = 0;
    for i = 1:k
        if( resnorm(nConv+i) < res_eps )
            fprintf('\nconverged eigenvalue %d: %s\n', nConv+i, num2str(r(i,i),'%8.4g'));
            nNewConv = nNewConv + 1;
        else
            %% TODO: also handle converged eigenvalues in the middle of the block q
            break
        end
    end
    if( nNewConv > 0 )
        % add it to converged data for deflation
        Q(:,                     nConv+1:nConv+nNewConv) = q(:,         1:nNewConv);
        R(1:nConv,               nConv+1:nConv+nNewConv) = a_(1:nConv,  1:nNewConv);
        R(nConv+1:nConv+nNewConv,nConv+1:nConv+nNewConv) = r(1:nNewConv,1:nNewConv);
        BQ(:,                    nConv+1:nConv+nNewConv) = Bq(:,        1:nNewConv);
        % remove it from V!
        V (:,1:nV-nNewConv) = V (:,1:nV) * Q_H(:,nNewConv+1:nV);
        BV(:,1:nV-nNewConv) = BV(:,1:nV) * Q_H(:,nNewConv+1:nV);
        W (:,1:nV-nNewConv) = W (:,1:nV) * Q_H(:,nNewConv+1:nV);
        H (1:nV-nNewConv, 1:nV-nNewConv) = Q_H(:,nNewConv+1:nV)' * H(1:nV,1:nV) * Q_H(:,nNewConv+1:nV);
        V_orth = max(max(abs(eye(nV-nNewConv) - BV(:,1:nV-nNewConv)' * V(:,1:nV-nNewConv))));
        % remove it from Q_H, R_H
        Q_H = eye(nV-nNewConv);
        R_H = R_H(nNewConv+1:nV,nNewConv+1:nV);
        nV = nV - nNewConv;
        nConv = nConv + nNewConv;
        if( nConv == nEig )
            return;
        end
        k = min(k,nEig-nConv);

%test_H = V(:,1:nV)' * W(:,1:nV);
%if( H(1:nV,1:nV) - test_H > res_eps )
%    error('H wrong!');
%end
%test_H = max(max(abs(H(1:nV,1:nV)*Q_H - Q_H*R_H)));
%if( test_H > res_eps )
%    error('Q_H R_H wrong: %e', test_H);
%end
        
        % recalculate r,q,...
        r = R_H(1:k,1:k);
        q = V(:,1:nV) * Q_H(:,1:k);
        Bq = BV(:,1:nV) * Q_H(:,1:k);
        % calculate residuum
        Aq = W(:,1:nV) * Q_H(:,1:k);
        res = Aq - Bq*r;
        a_ = BQ(:,1:nConv)'*res;
        res = res - Q(:,1:nConv)*a_;
        for i = 1:k
            resnorm(nConv+i) = norm(res(:,i),2);
        end
        lambda_history(1:nConv,m) = diag(R(1:nConv,1:nConv));
        lambda_history(nConv+1:nConv+k,m) = diag(r);
        resnorm_history(:,m) = resnorm;
    end
    
    
    % display stuff
    if( m > arnoldiIter )
        fprintf('\nBlockJaDa-iteration %d:',m - arnoldiIter);
    else
        fprintf('\nBlockArnoldi-iteration %d:', m);
    end
    fprintf('\nconverged eigenvalues: %d (of %d)', nConv, nEig);
    fprintf('\nsubspace dimension %d:',nV);
    fprintf('\n V-orthog. %e:',V_orth);
    fprintf('\n app. eigenvalues: %s', num2str([diag(R(1:nConv,1:nConv)); diag(r)]','%8.4g'));
    fprintf('\n   schur residuum:');  fprintf(' %6.2g', resnorm);
    fprintf('\n');

    
    if( m <= arnoldiIter )
        t = res;
        continue;
    end

    % construct q_ and Bq_ which deflate already converged eigenvalues
    firstDeflationVec = max(1,nConv-maxDeflationVecs);
    q_ = [Q(:,firstDeflationVec:nConv) q];
    Bq_ = [BQ(:,firstDeflationVec:nConv) Bq];
%if( max(max(abs(q_' * Bq_ - eye(nConv+k)))) > res_eps )
%    error('q_Bq_')
%end
%r_ = [R(1:nConv,1:nConv) a_(1:nConv,1:k); zeros(k,nConv) r];
%q_Orth = max(max(abs(q_' * Bq_ - eye(nConv+k)))); 
%if( q_Orth > res_eps )
%    error('q_/Bq_ not orthogonal: %e', q_Orth);
%end
    t = zeros(n,k);
    if( true )
        % solve approximately
        % (I-Bqq')A(I-qBq')*t - (I-Bqq')*B*t*r = -res
        % as r is an upper triangular matrix, we can do this column for column
        % (I-Bqq')(A-r(i,i)*B)(I-qBq')*t(:,i) = -res(:,i) + (I-Bqq')*B*t(:,1:i-1)*r(1:i-1,i)
        for i = 1:k
            lambda_iter(nConv+k) = lambda_iter(nConv+k) + 1;
            m_ = lambda_iter(nConv+k);
            % i-th residuum
            res_i = -res(:,i);
            % analyse behaviour if we omit the coupling here...
            if( false )
                for j = 1:1:i-1
                    res_i = res_i - r(j,i)*B*t(:,j);
                end
            end
            % solve i-th equation
            %A_proj_i = (eye(n)-q*q')*(A-r(i,i)*eye(n))*(eye(n)-q*q');
            %t(:,i) = res_i\A_proj_i;
            shift = r(i,i);
            [t(:,i)] = gmres(@A_proj, res_i, ...
                           min(40,10*m_), ...
                           min(0.5,max(res_eps,1/m_^2)), ...
                           max(floor(4/m_),1));
            % use the projection of t onto (I-qBq')
            t(:,i) = t(:,i) - q_*(Bq_'*t(:,i));
        end
    else % standard block_jada?
        % solve approximately
        %(I-Bqq')(A-lambda(i) B)(I-qBq')t(:,i) = -res_ev(:,i)
        for i = 1:k
            if( resnorm(i) < res_eps )
                t(:,i) = 0;
                continue;
            end
            %A_proj_i = (eye(n)-q*Bq')*(A-r(i,i)*B)*(eye(n)-q*Bq');
            %t(:,i) = -res_ev(:,i)\A_proj_i;
            res_ev(:,i) = res_ev(:,i) - q*(Bq'*res_ev(:,i));
            shift = r(i,i);
            [t(:,i)] = gmres(@A_proj, -res_ev(:,i), ...
                           min(40,10*m), ...
                           min(0.5,max(res_eps,1/m^2)), ...
                           max(floor(4/m),1));
            % use the projection of t onto (I-qBq')
            t(:,i) = t(:,i) - q_*(Bq_'*t(:,i));
        end
    end
end
warning('reached maxIter');


function y = A_proj(x)
    y = x-q_*(Bq_'*x);
    y = A*y - shift*(B*y);
    y = y - Bq_*(q_'*y);
    totalMatVecs = totalMatVecs + 1;
end

end


function [y,rank_y] = orth(x,D,eps)
    n = size(x,1);
    k = size(x,2);
    y = zeros(n,k);
    l = 1;
    for i = 1:k;
        y(:,l) = x(:,i);
        for j = 1:i-1
            alpha_jl = y(:,j)'*D*y(:,l);
            y(:,l) = y(:,l) - y(:,j)*alpha_jl;
        end
        alpha_ii = sqrt(y(:,l)'*D*y(:,l));
        if( alpha_ii > eps )
            y(:,l) = y(:,l)/alpha_ii;
            l = l + 1;
        end
    end
    rank_y = l-1;
    y = y(:,1:rank_y);
end
