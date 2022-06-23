function X = FastPFP(X0, A1, A2, B1, B2, alpha, lambda, threshold1, threshold2, max_iter1, max_iter2)
N1 = size(A1, 1);
N2 = size(A2, 1);

if N2 > N1
    error('Note N2 must be less than or equal to N1');
end

ones1 = ones(N1, 1);
ones2 = ones(N2, 1);


if ~isempty(X0)
    X = X0;
else
    X = ones1 *  ones2.' / sqrt(N1 * N2);
end

X_diff = Inf;
while X_diff > 1e-6
    X_new = X ./ sum(X, 1);
    X_new = X_new ./ sum(X_new, 2);
    X_diff = max(abs(X - X_new), [], 'all');
    X = X_new;
end

Y = zeros(N1, N1);

K = B1 * B2.';

epsilon1 = Inf;
iter1 = 0;

loss_val = loss_function(X, A1, A2, B1, B2, lambda);
fprintf('iter1 = %d, loss = %e, eps1 = %e\n', iter1, loss_val, epsilon1);
while epsilon1 > threshold1 && iter1 < max_iter1
    Y(1:N1,1:N2) = A1 * X * A2 + lambda * K;
    epsilon2 = Inf;
    iter2 = 0;
    
    while epsilon2 > threshold2 && iter2 < max_iter2
        tmp = eye(N1, N1) / N1;
        tmp = tmp + ones1.' * Y * ones1 / (N1 * N1) * eye(N1, N1);
        tmp = tmp - Y / N1;
        tmp = tmp * (ones1 * ones1.');
        Y_new = Y + tmp - (ones1 * ones1.') * Y / N1;
        Y_new = (Y_new + abs(Y_new)) / 2;
        epsilon2 = max(abs(Y_new - Y), [], 'all');
        Y = Y_new;
        iter2 = iter2 + 1;
    end
    
    X_new = (1 - alpha) * X + alpha * Y(1:N1, 1:N2);
    X_new = X_new / max(X_new, [], 'all');
    epsilon1 = max(abs(X_new - X), [], 'all');
    X = X_new;
    
    iter1 = iter1 + 1;
    loss_val = loss_function(X, A1, A2, B1, B2, lambda);
    fprintf('iter1 = %d, loss = %e, eps1 = %e\n', iter1, loss_val, epsilon1);
end

end

