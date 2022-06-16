function val = loss_function(P, A1, A2, B1, B2, lambda)
val = 0.5 * norm(A1 - P * A2 * P.', 'fro')^2 + lambda * norm(B1 - P * B2, 'fro')^2;
end