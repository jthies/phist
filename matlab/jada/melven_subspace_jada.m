function [r,q,resnorm,resnorm_history,m]=melven_subspace_jada(A,B,v0,maxIter,res_eps)

n = size(A,1);                                  % matrix dimension
k = size(v0,2);                                 % blocksize, e.g. number of eigenvalues wanted

V = zeros(n,maxIter*k);                         % search space empty initially
BV = zeros(n,maxIter*k);
nV = 0;
W = zeros(n,maxIter*k);
t = v0;                                         % initial new direction
H = zeros(maxIter*k);                           % preallocate H
resnorm = zeros(k,1);
resnorm_history = zeros(k,maxIter);
resnorm_ev = zeros(k,1);
lambda = zeros(k,1);


for m = 1:maxIter                               % main iteration loop
    nV_old = nV;
    if( m > 1 )
        t = t - V(:,1:nV_old) * (BV(:,1:nV_old)' * t);     % orthogonalize t wrt V
    end
    [t,rank_t] = orth(t,B,res_eps); % need B-orthogonality
    %inv_norm_t = diag(1./sqrt(sum(t'*t)));
    %t = t*inv_norm_t;
    %rank_t = rank(t,res_eps);
    fprintf(' rank of correction: %d\n', rank_t);
    if( rank(t) == 0 )
        return
    end
    nV = nV_old + rank_t;
    V(:,nV_old+1:nV) = t;                     % update subspace V
    BV(:,nV_old+1:nV) = B*t;
    V_orth = max(max(abs(eye(nV) - BV(:,1:nV)' * V(:,1:nV))));
    if( V_orth > res_eps )
        error('V not orthogonal: %g', V_orth);
    end
    W(:,nV_old+1:nV) = A*t;                   % update W = A*V

    % update H
    H(1:nV_old,nV_old+1:nV)     = V(:,1:nV_old)'     * W(:,nV_old+1:nV);
    H(nV_old+1:nV,1:nV_old)     = V(:,nV_old+1:nV)' * W(:,1:nV_old);
    H(nV_old+1:nV,nV_old+1:nV) = V(:,nV_old+1:nV)' * W(:,nV_old+1:nV);

    % compute the schur form of H
    [Q_H,R_H] = schur(H(1:nV,1:nV),'complex');
    % extract largest eigenvalues (in modulus)
    DH = ordeig(R_H);
    [~,sort_index] = sort(abs(DH),'descend');
    select_ev = false(nV,1);
    select_ev(sort_index(1:k)) = true;
    % sort schur form
    [Q_H,R_H] = ordschur(Q_H,R_H,select_ev);

    % get our current approximation of  A*Q = Q*R
    r = R_H(1:k,1:k);
    q = V(:,1:nV) * Q_H(:,1:k);
    Bq = BV(:,1:nV) * Q_H(:,1:k);
    % calculate residuum
    Aq = W(:,1:nV) * Q_H(:,1:k);
    res = Aq - Bq*r;
    for i = 1:k
        resnorm(i) = norm(res(:,i),2);
    end
    resnorm_history(:,m) = resnorm;
    
    
    % explicitly calculate eigenvalues and -vectors 
    v = q;
    for i = 1:k
        for j = 1:1:i-1
            v(:,i) = v(:,i) + r(j,i)/(r(i,i)-r(j,j))*q(:,j);
        end
    end
    delta_lambda = lambda - diag(r);
    lambda = diag(r);
    res_ev = A*v - B*v*diag(lambda);
    for i = 1:k
        resnorm_ev(i) = norm(res_ev(:,i),2);
    end
    
    % display stuff
    fprintf('\niteration %d:',m);
    fprintf('\n V-orthog. %e:',V_orth);
    fprintf('\n app. eigenvalues: %s', num2str(diag(r)','%8.4g'));
    fprintf('\n   schur residuum:');  fprintf(' %6.2g', resnorm);
    %fprintf('\n eigenv. residuum:');  fprintf(' %10.4g', resnorm_ev);
    if( m > 1 )
        fprintf('\n    eigenv. delta:');  fprintf(' %6.2g', abs(delta_lambda));
    end
    fprintf('\n');


    if( true )
        % solve approximately
        % (I-Bqq')A(I-qBq')*t - (I-Bqq')*B*t*r = -res
        % as r is an upper triangular matrix, we can do this column for column
        % (I-Bqq')(A-r(i,i)*B)(I-qBq')*t(:,i) = -res(:,i) + (I-Bqq')*B*t(:,1:i-1)*r(1:i-1,i)
        for i = 1:k
            % don't correct already converged eigenvalues
            if( resnorm(i) < res_eps )
                t(:,i) = 0;
                continue
            end
            % i-th residuum
            res_i = -res(:,i);
            % analyse behaviour if we omit the coupling here...
            if( true )
                for j = 1:1:i-1
                    res_i = res_i - r(j,i)*B*t(:,j);
                end
            end
            res_i = res_i - q*(Bq'*res_i);
            % solve i-th equation
            %A_proj_i = (eye(n)-q*q')*(A-r(i,i)*eye(n))*(eye(n)-q*q');
            %t(:,i) = res_i\A_proj_i;
            shift = r(i,i);
            t(:,i) = gmres(@A_proj, res_i, ...
                           min(40,10*m), ...
                           min(0.5,max(res_eps,1/m^2)), ...
                           max(floor(4/m),1));
            % use the projection of t onto (I-qBq')
            t(:,i) = t(:,i) - q*(Bq'*t(:,i));
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
            t(:,i) = gmres(@A_proj, -res_ev(:,i), ...
                           min(40,10*m), ...
                           min(0.5,max(res_eps,1/m^2)), ...
                           max(floor(4/m),1));
            t(:,i) = t(:,i) - q*(Bq'*t(:,i));
        end
    end
end

function y = A_proj(x)
    y = x-q*(Bq'*x);
    y = A*y - shift*(B*y);
    y = y - Bq*(q'*y);
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