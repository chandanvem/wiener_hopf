clear all; close all; clc;

%%
num_of_threads = 24;

%%
op_dir = '/work/home/chandan/Chandan_case_files/wiener_hopf_results/two_stream_subsonic_runs/results_St_1.7/roots_delta_vary/matlab_results_del_0.001';
data_op_dir = sprintf('%s/DataDump',op_dir);
log_op_dir = sprintf('%s/log',op_dir);

if exist(op_dir,'dir')== 0
    mkdir(op_dir);
end

if exist(data_op_dir,'dir')== 0
    mkdir(data_op_dir);
end

if exist(log_op_dir,'dir')== 0
    mkdir(log_op_dir);
end

%%
TOL = 1e-15;
ERR_TOL = 1e-6;

k_T = sqrt(1);
M_j = 0.9;
M_o = 1e-10;

delta_degrees = 0.001;

St_list = 0.1:0.01:2.0;
St_min = St_list(1);
St_max = St_list(end);

if length(St_list) == 1
    d_St = 0;
else
    d_St = St_list(2)-St_list(1);
end

num_frequencies = length(St_list);
omega_list = pi*M_j*St_list*exp(1i*delta_degrees*pi/180);
s_solutions_vs_freq = cell(num_frequencies,1);

azim_wave_number = 0;

pole_zero_string = 'zero_mode';

s_real_max =   2;
s_real_min =  -2;
N_real = 100;

s_real_stencil = linspace(s_real_min,s_real_max,N_real);
% N_real = length(kD_r5eal_stencil);

s_imag_max =   2;
s_imag_min =  -2;
N_imag = 2e2;
s_imag_stencil = linspace(s_imag_min,s_imag_max,N_imag);

num_of_solutions = 100;

options = optimset('TolFun',TOL,'Display','off');
%%
parpool('local', num_of_threads);
parfor freq_idx = 1 : num_frequencies
    [s_real,s_imag] = meshgrid(s_real_stencil,s_imag_stencil);
    s_solution_list = zeros(1,num_of_solutions);
    s_solution_residue = s_solution_list;

    log_file_name = sprintf('%s/log_%d.dat',log_op_dir,freq_idx);
    fileID_log = fopen(log_file_name, 'w');

    count = 1;
    result_count = 1;
    omega = omega_list(freq_idx);
    St_target = St_list(freq_idx);

    for real_idx = 1 : N_real
        for imag_idx = 1 : N_imag
            s_guess = s_real(imag_idx,real_idx) + 1i*s_imag(imag_idx,real_idx);

            disp_rel = @(s) vortex_sheet_disp_rel_equation(s,omega,k_T,M_j,M_o,...
                azim_wave_number,pole_zero_string);
            [s_solution_iter,~,exit_flag] = fsolve(disp_rel,s_guess,options);
            s_solution_iter_residue= ...
                vortex_sheet_disp_rel_equation(s_solution_iter,omega,k_T,M_j,M_o,...
                azim_wave_number,pole_zero_string);
            if exit_flag <=4 && exit_flag >= 1 && abs(s_solution_iter_residue) < ERR_TOL
                if result_count == 1
                    s_solution_list(result_count) = s_solution_iter;
                    s_solution_residue(result_count) = ...
                        vortex_sheet_disp_rel_equation(s_solution_iter,omega,k_T,M_j,M_o,...
                        azim_wave_number,pole_zero_string);
                    result_count = result_count + 1;
                else
                    abs_err = abs(s_solution_list - s_solution_iter);
                    if min(abs_err) > ERR_TOL
                        s_solution_list(result_count) = s_solution_iter;
                        s_solution_residue(result_count) = ...
                            vortex_sheet_disp_rel_equation(s_solution_iter,omega,k_T,M_j,M_o,...
                            azim_wave_number,pole_zero_string);
                        result_count = result_count + 1;
                    end
                end
            end
            count = count + 1;
            fprintf(fileID_log,...
                'Grid point = %d/%d for freq idx %d/%d\n',count-1,N_real*N_imag,freq_idx,num_frequencies);
        end
    end
    fclose(fileID_log);
    if result_count > 1
        if strcmp(pole_zero_string,'zero_mode')
            op_data_file =  sprintf('%s/zerolist_m_%d_St_%1.3g.dat',...
                data_op_dir,azim_wave_number,St_target);
        else
            op_data_file =  sprintf('%s/polelist_m_%d_St_%1.3g.dat',...
                data_op_dir,azim_wave_number,St_target);    
        end
        op_dat_file_ID = fopen(op_data_file, 'w');
        for point_idx = 1 : result_count-1
            fprintf(op_dat_file_ID, '%12.6f %12.6f %18.12f %18.12f\n', ...
                real(s_solution_list(point_idx)),...
                imag(s_solution_list(point_idx)),...
                real(s_solution_residue(point_idx)),...
                imag(s_solution_residue(point_idx)));
        end
        fclose(op_dat_file_ID);
    end
end
