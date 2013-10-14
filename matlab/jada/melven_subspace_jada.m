function [r,q,resnorm,resnorm_history,m]=melven_subspace_jada(A,v0,maxIter,res_eps)

n = size(A,1);                                  % matrix dimension
k = size(v0,2);                                 % blocksize, e.g. number of eigenvalues wanted

V = zeros(n,maxIter*k);                         % search space empty initially
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
        t = t - V(:,1:nV_old) * (V(:,1:nV_old)' * t);     % orthogonalize t wrt V
    end
    t = orth(t);
    inv_norm_t = diag(1./sqrt(sum(t'*t)));
    t = t*inv_norm_t;
    rank_t = rank(t,res_eps);
    fprintf(' rank of correction: %d\n', rank_t);
    if( rank(t) == 0 )
        return
    end
    nV = nV_old + rank_t;
    V(:,nV_old+1:nV) = t;                     % update subspace V
    V_orth = max(max(abs(eye(nV) - V(:,1:nV)' * V(:,1:nV))));
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
    % calculate residuum
    Aq = W(:,1:nV) * Q_H(:,1:k);
    res = Aq - q*r;
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
    res_ev = A*v - v*diag(lambda);
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
        % (I-qq')A(I-qq')*t - (I-qq')*t*r = -res
        % as r is an upper triangular matrix, we can do this column for column
        % (I-qq')(A-r(i,i)*I)(I-qq')*t(:,i) = -res(:,i) + (I-qq')*t(:,1:i-1)*r(1:i-1,i)
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
                    res_i = res_i + r(j,i)*t(:,j);
                end
            else
                res_i = res_i - q*(q'*res_i);
            end
            % solve i-th equation
            %A_proj_i = (eye(n)-q*q')*(A-r(i,i)*eye(n))*(eye(n)-q*q');
            %t(:,i) = res_i\A_proj_i;
            shift = r(i,i);
            t(:,i) = gmres(@A_proj, res_i);
            % use the projection of t onto (I-qq')
            t(:,i) = t(:,i) - q*(q'*t(:,i));
        end
    else % standard block_jada?
        % solve approximately
        %(I-qq')(A-lambda(i) I)(I-qq')t(:,i) = -res_ev(:,i)
        for i = 1:k
            if( resnorm(i) < res_eps )
                t(:,i) = 0;
                continue;
            end
            %A_proj_i = (eye(n)-q*q')*(A-r(i,i)*eye(n))*(eye(n)-q*q');
            %t(:,i) = -res_ev(:,i)\A_proj_i;
            res_ev(:,i) = res_ev(:,i) - q*(q'*res_ev(:,i));
            shift = r(i,i);
            t(:,i) = gmres(@A_proj, -res_ev(:,i));
            t(:,i) = t(:,i) - q*(q'*t(:,i));
        end
    end
end

function y = A_proj(x)
    y = x-q*(q'*x);
    y = A*y - shift*y;
    y = y - q*(q'*y);
end

end
