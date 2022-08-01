function P = ConstrainedMinimization_Equal(X0, A1, A2, B1, B2, lambda)
% division_constraints [p x 1] are the indices of nodes in A1 that are NOT allowed to divide

% lineage_constraints [q x 2] is a list of matches which must be enforced. Row is [ind1, ind2] where ind1 
% is the index of a node in A1 and ind2 is the index of a node in A2.


N1 = size(A1, 1);
N2 = size(A2, 1);

ones1 = ones(N1, 1);
ones2 = ones(N2, 1);

if isempty(X0)
    X0 = ones1 *  ones2.' / sqrt(N1 * N2); % uniform prior
end

X0_diff = Inf;
while X0_diff > 1e-6
    X0_new = X0 ./ sum(X0, 1);
    X0_new = X0_new ./ sum(X0_new, 2);
    X0_diff = max(abs(X0 - X0_new), [], 'all');
    X0 = X0_new;
end

K1 = B1 * B2.';
K2 = B2 * B2.';

% Future graph is A2 and previous graph is A1
Aeq_1 = kron(ones(1, N2), eye(N1));
Aeq_2 = kron(eye(N2), ones(1, N1));

Aeq = [Aeq_1; Aeq_2];
beq = ones(size(Aeq, 1), 1);

A = [];
b = [];

lb = zeros(N1*N2, 1);
ub = [];
nonlcon = [];
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter-detailed', ...
    'SpecifyObjectiveGradient', true);

p_vec = fmincon(@(X) objective(X, A1, A2, B1, B2, K1, K2, lambda), X0(:), A, b, Aeq, beq, lb, ub, nonlcon, options);
P = reshape(p_vec, N1, N2);

end

function x0 = fminpenalty(obj, x0, A,b, Aeq, beq, lb,up, nonlcon, options)

error('not working');
A = [A; Aeq];
b = [b; beq];

nx = size(x0, 1);
nc = size(A, 1);
lambda = zeros(nc, 1);

x_k = x0;
[f_k, df_dx_k] = obj(x_k);

s_k = -df_dx_k;
x_k_plus1 = x_k + s_k;

[f_k_plus1, df_dx_k_plus1] = obj(x_k_plus1);
y_k = df_dx_k_plus1 - df_dx_k;

H = eye(nx);
H = H + y_k * y_k.' / (y_k.' * s_k) - H * (s_k * s_k.') * H.' / (s_k.' * H  * s_k);

zero_mat = zeros(nc);
step = 1;
dx_norm = Inf;
while dx_norm > 1e-6  
    f_k = f_k_plus1;
    x_k = x_k_plus1;


    [f, df_dx] = obj(x0);

    c  = A * x0 - b;
    dc_dx = A;
    
    R = -[df_dx - dc_dx.' * lambda; c];

    Hfull = [H, dc_dx.'; dc_dx, zero_mat];
    
    d = Hfull \ R;
    dx = d(1:nx);
    dlambda = d(nx+1:end);
    
    dx_norm = norm(dx);
    
    if write_output
        fprintf('Step %d: Fun = %e, d_norm = %e\n', step, f, dx_norm);
    end
    step = step + 1;
end

end

function [f, grad] = objective(x_vec, A1, A2, B1, B2, K1, K2, lambda)
N1 = size(A1, 1);
N2 = size(A2, 1);
X = reshape(x_vec, N1, N2);

D1 = A1 - X * A2 * X.';
D2 = B1 - X * B2;

f = 0.5 * trace(D1.' * D1) + 0.5 * lambda * trace(D2.' * D2);
grad = (X * A2 * (X.' * X) * A2.' + X * A2.' * (X.' * X) * A2) - (A1.' * X * A2 + A1 * X * A2.') + ...
       lambda * (0.5 * X * (K2.' + K2) - K1);

grad = grad(:);

end