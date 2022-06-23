function P = ConstrainedMinimization_General(X0, A1, A2, B1, B2, lambda, division_constraints, lineage_constraints)
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

if ~isempty(lineage_constraints)
    % make sure initial condition satisfies the constraints

    for ii = 1:size(lineage_constraints, 1)
        X0(:,lineage_constraints(ii,2)) = 0;
        X0(lineage_constraints(ii,1), lineage_constraints(ii,2)) = 1;
    end
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
temp_A_N2 = ones(1, N2);
if ~isempty(lineage_constraints)
    temp_A_N2(lineage_constraints(:,2)) = 0;
end
A_base1 = kron(temp_A_N2, eye(N1));
A_base2 = -A_base1;

if ~isempty(lineage_constraints)
    union_div_lin_constraints = union(division_constraints, lineage_constraints(:,1));
else
    union_div_lin_constraints = division_constraints;
end

for ii = 1:size(lineage_constraints, 1)
    A_base1(lineage_constraints(ii,1),lineage_constraints(ii,1) + (lineage_constraints(ii,2) - 1) * N1) = 1;
end

Aeq_extra = A_base1(division_constraints, :);

A_base1(division_constraints,:) = [];
A_base2(union_div_lin_constraints,:) = [];

A = [A_base1; A_base2];
b = [2 * ones(size(A_base1, 1), 1); -ones(size(A_base2, 1), 1)];

temp_Aeq_N1 = ones(1, N1);
temp_Aeq_N2 = ones(N2, 1);
if ~isempty(lineage_constraints)
    temp_Aeq_N2(lineage_constraints(:,2)) = 0;
end
Aeq = kron(diag(temp_Aeq_N2), temp_Aeq_N1);
for ii = 1:size(lineage_constraints, 1)
    Aeq(lineage_constraints(ii,2),lineage_constraints(ii,1) + (lineage_constraints(ii,2) - 1) * N1) = 1;
end

Aeq = [Aeq; Aeq_extra; ones(1, N1 * N2)];
beq = [ones(N2, 1); ones(size(Aeq_extra, 1), 1)];
beq = [beq; N2];

lb = zeros(N1*N2, 1);
ub = [];
nonlcon = [];
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter-detailed', 'SpecifyObjectiveGradient', true);

p_vec = fmincon(@(X) objective(X, A1, A2, B1, B2, K1, K2, lambda), X0(:), A, b, Aeq, beq, lb, ub, nonlcon, options);
P = reshape(p_vec, N1, N2);

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
