function f = plot_pressure_eigenmode(s_r,s_i,omega,k_T,M_j,M_o,azim_wave_number,fig_id,scale)

Nr = 501;
r  = linspace(1e-6,4,Nr);

s = s_r + s_i*1i;

lambda_j_plus =  sqrt(1-(s*(M_j + 1)));
lambda_j_minus = sqrt(1-(s*(M_j - 1)));
lambda_j = lambda_j_plus*lambda_j_minus;

lambda_o_plus  =  sqrt(k_T-(s*((k_T*M_o) + 1)));
lambda_o_minus =  sqrt(k_T-(s*((k_T*M_o) - 1)));
lambda_o = lambda_o_plus*lambda_o_minus;

c_o = ((1-(s*M_j))/(1-(s*M_o)))*...
      besselj(azim_wave_number,lambda_j*omega,scale)/besselh(azim_wave_number,1,lambda_o*omega,scale);

p_hat =              (1-(s*M_j))*besselj(azim_wave_number,lambda_j*omega*r,scale).*heaviside(1-r);
p_hat = p_hat +  c_o* (1-(s*M_o))*besselh(azim_wave_number,1,lambda_o*omega*r,scale).*heaviside(r-1);
% p_hat = p_hat +  c_o*besselk(azim_wave_number,lambda_o*omega*r,scale).*heaviside(r-1);
% p_hat = c_o*besselh(azim_wave_number,1,lambda_o*omega*r,scale).*heaviside(r-1);

figure(fig_id);
clf;
w1 = plot(r,abs(p_hat)/max(abs(p_hat)));
set(w1,'marker','none','LineStyle','-','LineWidth',3);
grid on;
% xlim([-4 2]);
% ylim([-40 40]);
xlim([0 2]);
ylim([0 1]);
box on;
xlabel('$r/D$','interpreter','latex','fontsize',36);
ylabel('$p$','interpreter','latex','fontsize',36,'rotation',0);
set(gca,'fontsize',32,'TickLabelInterpreter','latex');
axis square;

end