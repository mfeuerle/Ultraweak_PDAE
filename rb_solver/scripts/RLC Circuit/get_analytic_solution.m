function [x_analy, u_analy, dx_analy] = get_analytic_solution(R,L,C,U_hat,phi_U,omega)
%Calculate analytic solution for RLC circuit with cosine as control via 
%phase vector representation

% ----------------------------------------------------------------------
% REFERENCE:
%  E.Beurer, M.Feuerle, N.Reich, K.Urban
%  "An ultraweak variational method for parmeterized linear 
%  differential-algebraic equations"
%  Ulm University, 2022
%  https://doi.org/10.48550/arXiv.2202.12834
%  https://github.com/mfeuerle/Ultraweak_PDAE
% ----------------------------------------------------------------------

if nargin == 5 || omega == 0
    omega = 1/sqrt(L*C);
end

u_analy =@(t) U_hat*cos(omega*t+phi_U);

%% Analytic solution
%Impedances
Z_R = R;
Z_L = 1i*omega*L;
Z_C = -1i/(omega*C);
Z_Ges = Z_R + Z_L + Z_C;
%Volatge from the voltage source
U_ = U_hat*exp(1i*phi_U);
%Results as phase vectors
I_ = U_ / Z_Ges;
U_R = Z_R*I_;
U_L = Z_L*I_;
U_C = Z_C*I_;
%Transformation to cosine functions
phi_I = angle(I_);
phi_U_R = angle(U_R);
phi_U_L = angle(U_L);
phi_U_C = angle(U_C);
I_hat = abs(I_);
U_R_hat = abs(U_R);
U_L_hat = abs(U_L);
U_C_hat = abs(U_C);

x_analy =@(t) [  I_hat   * cos(omega*t + phi_I);
    U_C_hat * cos(omega*t + phi_U_C);
    U_L_hat * cos(omega*t + phi_U_L);
    U_R_hat * cos(omega*t + phi_U_R)];

dx_analy =@(t) [ -I_hat   * omega * sin(omega*t + phi_I);
    -U_C_hat * omega * sin(omega*t + phi_U_C);
    -U_L_hat * omega * sin(omega*t + phi_U_L);
    -U_R_hat * omega * sin(omega*t + phi_U_R)];
end

