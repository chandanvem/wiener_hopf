function f = vortex_sheet_disp_rel_equation(s,omega,k_T,M_j,M_o,m)

      scale = 1;
      k_rho = 1/k_T^2;
      
      lambda_j_plus =  sqrt(1-(s*(M_j + 1)));
      lambda_j_minus = sqrt(1-(s*(M_j - 1)));
      lambda_j = lambda_j_plus.*lambda_j_minus;      

      lambda_o_plus  =  sqrt(k_T-(s*((k_T*M_o) + 1)));
      lambda_o_minus =  sqrt(k_T-(s*((k_T*M_o) - 1)));
      lambda_o = lambda_o_plus.*lambda_o_minus;

      term_1_num = k_rho.*((1 - (s*M_j)).^2).*besselj(m,lambda_j.*omega,scale);
      term_1_den = 0.5*(besselj(m-1,lambda_j*omega,scale)-besselj(m+1,lambda_j*omega,scale));
      term_1_den = lambda_j.*term_1_den;

      term_1 = term_1_num./term_1_den;

      term_2_num = ((1 - (s*M_o)).^2).*besselh(m,1,lambda_o*omega,scale);
%       term_2_den = (m*besselh(m,1,lambda_o*omega,scale)/lambda_o*omega)-besselh(m+1,1,lambda_o*omega,scale);
%       if lambda_o == 0 
%           term_2_den = 0;
%       else
       term_2_den = 0.5*(besselh(m-1,1,lambda_o*omega,scale)-besselh(m+1,1,lambda_o*omega,scale));
%       end

      term_2_den = lambda_o.*term_2_den;

      term_2 = term_2_num./term_2_den;

      f = term_1 - term_2;


end