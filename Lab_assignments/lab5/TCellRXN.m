function [S,P,K] = TCellRXN(params)
% takes in parameters from t-cell reaction model
% returns the S, P, K matrices
% 6 reactions, 6 elements, 6 x 6 matrices
% S - substractes
% P - products
% K - reaction rates
% S --> P (K)

% A1 --> E + A1         (k1)
% A2 --> S + A2         (k2)
% E + A1 --> E + A1_p   (k3)
% A1_p --> E + A1_p     (k4)
% S + A1 --> S + A1_i   (k5)
% S --> 0; E --> 0      (kD)

% order: A1, A1_p, A1_i, A2, E, S

S = [
    1 0 0 0 0 0
    0 0 0 1 0 0
    1 0 0 0 1 0
    0 1 0 0 0 0
    1 0 0 0 0 1
    0 0 0 0 1 1
    ];

P = [
    1 0 0 0 1 0
    0 0 0 1 0 1
    0 1 0 0 1 0
    0 1 0 0 1 0
    0 0 1 0 0 1
    0 0 0 0 0 0
    ];

K = params';

end