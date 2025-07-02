clear all;
%close all;

%%
op_dir = '/work/home/chandan/Chandan_case_files/trapped_acoustic_modes/AS_disp_relation/validation/omega_corrected/for_tracking_co_flow';

if exist(op_dir,'dir')== 0
    mkdir(op_dir);
end

TOL = 1e-15;
ERR_TOL = 1e-6;

k_T = sqrt(1);
M_j = 0.9;
M_o = 1e-10;

St_list = [0.4];
St_min = St_list(1);
St_max = St_list(end);
d_St = 0; 

num_frequencies = length(omega_list);
s_solutions_vs_freq = cell(num_frequencies,1);

azim_wave_number = 0;

s_real_max =   2;
s_real_min =  -8;
N_real = 250;
% kD_real_stencil = [-8:2:8]; %
s_real_stencil = linspace(s_real_min,s_real_max,N_real);
% N_real = length(kD_r5eal_stencil);

s_imag_max =   0.5;
s_imag_min =  -1.5;
N_imag = 2e2;
s_imag_stencil = linspace(s_imag_min,s_imag_max,N_imag);

[s_real,s_imag] = meshgrid(s_real_stencil,s_imag_stencil);
k_solution = zeros(size(s_real));

options = optimset('TolFun',TOL);
% options = optimset('TolFun',TOL,'Algorithm','trust-region','FinDiffType','central');
% options = optimset('TolFun',ERR_TOL,'FinDiffType','central');
% options = optimset('TolFun',TOL,'Algorithm','levenberg-marquardt');
% k_t = kD_real ;
%%

for freq_idx = 1 : num_frequencies
    count = 1;
    result_count = 1;
    omega = omega_list(freq_idx);
    St_target = St_list(freq_idx);
    op_mat_file =  sprintf('%s/s_data_M_%1.3g_k_T_%d_m_%d_St_%1.3g.mat',...
                                     op_dir,M_j,k_T,azim_wave_number,St_target);
    
    if exist(op_mat_file,'file') == 2
        delete(op_mat_file);
    end

    for real_idx = 1 : N_real
        for imag_idx = 1 : N_imag
            s_guess = s_real(imag_idx,real_idx) + 1i*s_imag(imag_idx,real_idx);

            disp_rel = @(s) vortex_sheet_disp_rel_equation(s,omega,k_T,M_j,M_o,azim_wave_number);
            [s_solution_iter,~,exit_flag] = fsolve(disp_rel,s_guess,options);
             s_solution_iter_residue= ...
                        vortex_sheet_disp_rel_equation(s_solution_iter,omega,k_T,M_j,M_o,azim_wave_number);
            if exit_flag <=4 && exit_flag >= 1 && abs(s_solution_iter_residue) < ERR_TOL          
                if result_count == 1
                    s_solution_list(result_count) = s_solution_iter;
                    s_solution_residue(result_count) = ...
                        vortex_sheet_disp_rel_equation(s_solution_iter,omega,k_T,M_j,M_o,azim_wave_number);
                    result_count = result_count + 1;
                else
                    abs_err = abs(s_solution_list - s_solution_iter);
                    if min(abs_err) > ERR_TOL
                        s_solution_list(result_count) = s_solution_iter;
                        s_solution_residue(result_count) = ...
                            vortex_sheet_disp_rel_equation(s_solution_iter,omega,k_T,M_j,M_o,azim_wave_number);
                        result_count = result_count + 1;
                    end
                end
            end
            count = count + 1;
            fprintf('Count = %d/%d for freq idx %d/%d\n',count-1,N_real*N_imag,freq_idx,num_frequencies);
        end
    end
    save(op_mat_file,"s_solution_list");
%     s_solutions_vs_freq{freq_idx} = s_solution_list;
end
