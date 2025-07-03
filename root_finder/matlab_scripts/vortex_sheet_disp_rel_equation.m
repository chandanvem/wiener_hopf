function f = vortex_sheet_disp_rel_equation(s,omega,k_T,M_j,M_o,m,pole_zero_string)

scale = 1;
k_rho = 1/k_T^2;
lambda_j_plus =  sqrt(1-(s*(M_j + 1)));
lambda_j_minus = sqrt(1-(s*(M_j - 1)));
lambda_1 = lambda_j_plus.*lambda_j_minus;

lambda_o_plus  =  sqrt(k_T-(s*((k_T*M_o) + 1)));
lambda_o_minus =  sqrt(k_T-(s*((k_T*M_o) - 1)));
lambda_2 = lambda_o_plus.*lambda_o_minus;


if strcmp(pole_zero_string,'zero_mode')

    term_1_num = k_rho.*((1 - (s*M_j)).^2).*besselj(m,lambda_1.*omega,scale);
    term_1_den = 0.5*(besselj(m-1,lambda_1*omega,scale)-besselj(m+1,lambda_1*omega,scale));
    term_1_den = lambda_1.*term_1_den;

    term_1 = term_1_num./term_1_den;

    term_2_num = ((1 - (s*M_o)).^2).*besselh(m,1,lambda_2*omega,scale);
    term_2_den = 0.5*(besselh(m-1,1,lambda_2*omega,scale)-besselh(m+1,1,lambda_2*omega,scale));
    term_2_den = lambda_2.*term_2_den;
    term_2 = term_2_num./term_2_den;

    f = term_1 - term_2;

elseif strcmp(pole_zero_string,'pole_mode')

    term_1 = lambda_1;
    j_prime = 0.5*(besselj(m-1,lambda_1*omega,scale)-besselj(m+1,lambda_1*omega,scale));
    term_1 = term_1*j_prime;

    term_2 = lambda_2;
    h_prime = 0.5*(besselh(m-1,1,lambda_2*omega,scale)-besselh(m+1,1,lambda_2*omega,scale));
    term_2 = term_2*h_prime;

    f = term_1 * term_2;

end

end