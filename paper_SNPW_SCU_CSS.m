%% Modified discerete SNPW.  Match with SZKTW.
%% Init parallel
% clear all;
numberOfWorkers = 100;                       % Maximum number of cores.
parpool('local', numberOfWorkers);
%% Geometry
t0 = 0;
tMax = 50;
tstepsNPW = 1e5+1;
tarrayNPW = linspace(t0, tMax, tstepsNPW);
dtNPW = tarrayNPW(2) - tarrayNPW(1);
t_fb_on_NPW = 1;
t_fb_off_NPW = tstepsNPW;

% Measurement params
dt_meas_int = 0.8;                                                                                                                       % interval between measurements in non-dim time.
% dt_meas_int = 0.4;
dt_meas_int_timesteps = dt_meas_int/dtNPW;
tarrayNPW_meas = tarrayNPW(dt_meas_int_timesteps:dt_meas_int_timesteps:end);                       % Array of measurement times.
num_meas = length(tarrayNPW_meas);                                                                                           % total number of measurements.

tstepsNPW_samp = 5000+1;                                       % number of points to sample.  Pick divisor.
tarrayNPW_samp = linspace(t0, tMax, tstepsNPW_samp);
dt_samp_timesteps = (tstepsNPW-1)/(tstepsNPW_samp-1);
dtNPW_samp = tarrayNPW_samp(2)-tarrayNPW_samp(1);

%check sampling is set up correctly; directly after measurement. 
disp([tarrayNPW(dt_meas_int_timesteps), tarrayNPW_meas(1), tarrayNPW_samp(dt_meas_int_timesteps/dt_samp_timesteps+1), dtNPW, num_meas])
%% Atom-light parameters and state init.
% SKZTW params
N_atoms = 100;
N_light = 1e7;
b1_0 = sqrt(N_light);
chi = 0.01;
% chi = 0;
kappa = 0.09;  
% kappa = 0;
k_fb = 0.1;                                                                       % feedback gain parameter;
lambda = 0.6e-4;                                                              % SKZTW meas param.
% lambda = 0.8e-4;                                                              % SKZTW meas param.

% NPW params
t_p = dtNPW;                                                                    % instantaneous pulse.
gamma = lambda^2*b1_0^2/t_p;
tsteps_meas = 500;
dtNPW_meas = t_p/tsteps_meas;           

% Setting up cond. and fic. trajectories.
N_real_trajec = 500;               % Number of real noises (measurement)
N_fic_trajec = 500;                     % Number of fictitious noises (unravelling trajecs); 

% Weight initialisation
weights_init = 1/N_fic_trajec*ones(1, N_fic_trajec);            % initialise weights equal NORMAL SPACE
weights_old_log = log(weights_init);                            % initial weights equal in log space.
resamp_thresh = 1/2;                                            % resampling threshold between 0 and 1.
% resamp_thresh = 0.7;

% State initialisation.
n10 = 20;
n20 = 80;
phi10 = pi/2;
phi20 = 0;
specify = "CSS_custom";
normalise = 1;
[n1_old, n2_old, phi1_old, phi2_old, a1_old, a2_old] = npw_init_TM(n10, n20, phi10, phi20, N_fic_trajec, specify, normalise, 1);
[jx_npw_cond_avg, jy_npw_cond_avg, jz_npw_cond_avg, jx2_npw_cond_avg, jy2_npw_cond_avg, jz2_npw_cond_avg] = NPW_cond_avg(a1_old, a2_old, log(weights_init));
disp([(jx2_npw_cond_avg - jx_npw_cond_avg^2)/N_atoms, (jy2_npw_cond_avg - jy_npw_cond_avg^2)/N_atoms, (jz2_npw_cond_avg - jz_npw_cond_avg^2)/N_atoms]);
%% Initialise storage
% store (un)conditioned averages
store_jx_npw_cond_avg = zeros(N_real_trajec, tstepsNPW_samp); 
store_jy_npw_cond_avg = zeros(N_real_trajec, tstepsNPW_samp);
store_jz_npw_cond_avg = zeros(N_real_trajec, tstepsNPW_samp); 
store_jx2_npw_cond_avg = zeros(N_real_trajec, tstepsNPW_samp);
store_jy2_npw_cond_avg = zeros(N_real_trajec, tstepsNPW_samp);
store_jz2_npw_cond_avg = zeros(N_real_trajec, tstepsNPW_samp);

store_jx3_npw_cond_avg = zeros(N_real_trajec, tstepsNPW_samp);
store_jy3_npw_cond_avg = zeros(N_real_trajec, tstepsNPW_samp);
store_jz3_npw_cond_avg = zeros(N_real_trajec, tstepsNPW_samp);

store_OBDM_npw_cond_avg = zeros(N_real_trajec, tstepsNPW_samp, 2, 2);
store_energy_npw_cond_avg = zeros(N_real_trajec, tstepsNPW_samp);

store_fb = zeros(N_real_trajec, num_meas);                                   % store feedback signal derived from fic. average.
store_jz_est = zeros(N_real_trajec, num_meas);                             % store estimate of Jz from measurement signal.
store_jz_est_perf_diff = zeros(N_real_trajec, num_meas);
SNPW_flag = 1;
%% NPW Block
sim_NPW = 1;
if sim_NPW == 1
    parfor nn = 1:N_real_trajec

        resamp_counter = 0;                         % count number of times we resamped.
        resamp_flag = 0;
        fb = 0;                                              % initialise feedback.
        y_meas_sig = 0;
        
        % initialisation for parfor loop.
        [n1_old, n2_old, phi1_old, phi2_old, a1_old, a2_old] = npw_init_TM(n10, n20, phi10, phi20, N_fic_trajec, specify, normalise, 1);
        weights_init = 1/N_fic_trajec*ones(1, N_fic_trajec);            % initialise weights equal NORMAL SPACE
        weights_old_log = log(weights_init);                                    % initial weights equal in log space.
        weights_new_log = weights_old_log;

        tarrayNPW = linspace(t0, tMax, tstepsNPW);

        %% Initialise temp storage for parfor loop.
        store_fb_temp = zeros(1, num_meas);
        store_jz_est_temp = zeros(1, num_meas);
        store_jz_perf_temp = zeros(1, num_meas);
    
        store_jx_npw_cond_avg_temp = zeros(1, tstepsNPW_samp); 
        store_jy_npw_cond_avg_temp = zeros(1, tstepsNPW_samp);
        store_jz_npw_cond_avg_temp = zeros(1, tstepsNPW_samp); 

        store_jx2_npw_cond_avg_temp = zeros(1, tstepsNPW_samp);
        store_jy2_npw_cond_avg_temp = zeros(1, tstepsNPW_samp);
        store_jz2_npw_cond_avg_temp = zeros(1, tstepsNPW_samp);

        store_jx3_npw_cond_avg_temp = zeros(1, tstepsNPW_samp);
        store_jy3_npw_cond_avg_temp = zeros(1, tstepsNPW_samp);
        store_jz3_npw_cond_avg_temp = zeros(1, tstepsNPW_samp);

        store_energy_npw_cond_avg_temp = zeros(1, tstepsNPW_samp);
        store_OBDM_npw_cond_temp = zeros(tstepsNPW_samp, 2, 2);

        samp_counter = 1;               % counter for sampling array.
        meas_counter = 1;               % counter for measurement.
        for ii = 1:tstepsNPW         
            %% Block 1:  Meas. and Feedback turned off
            if ii < t_fb_on_NPW || ii >= t_fb_off_NPW
                fb = 0;
                [a1_new, a2_new] = TWA_CAHO(a1_old, a2_old, chi, kappa, dtNPW, fb, "jz_damp");          % standard evolution
                weights_new_log = weights_old_log;                                                                                 % weights don't change.

            % Block 2:  Make measurement, swap to new time array for measurement.
            elseif (ii >= t_fb_on_NPW) && (mod(ii,  dt_meas_int_timesteps) == 0)
               tarrayNPW_meas = linspace(tarrayNPW(ii), tarrayNPW(ii) + t_p, tsteps_meas);                    % new time array for measurement only.
               dtNPW_meas = tarrayNPW_meas(2) - tarrayNPW_meas(1);          
               y_meas_sig = 0;
               store_meas_sig_temp = zeros(1, tsteps_meas);                                                                   % store meas signal.
               
               % implement continuous time NPW in instantaneous time.
               for jj =1:tsteps_meas                                                                                                         
                    % generate noise trace records.
                    dW = sqrt(dtNPW_meas)*randn(1, 1);
                    dV_array = sqrt(dtNPW_meas)*randn(1, N_fic_trajec);
                    
                    % get perfect estimate from stochastic avg.
                    [a1_new, a2_new] = SIEM_TM_atomic(a1_old, a2_old, dtNPW_meas, fb, gamma, kappa, chi, dV_array, SNPW_flag);
                    % [a1_new, a2_new] = SIEM_TM_atomic_SNPW(a1_old, a2_old, dtNPW_meas, fb, gamma, kappa, chi, dV_array);
                    [weights_new_log] = SIEM_TM_weights(a1_old, a2_old, a1_new, a2_new, weights_old_log, gamma, dtNPW_meas, dW);

                    % check resampling criteria
                    [resamp, ess] = resamp_ess(weights_new_log, resamp_thresh);
                    resamp_flag = resamp;
                    % resample fields using Kitagawa deterministic sorted algorithm.
                    if resamp_flag == 1
                        % fprintf('\nI resamped at timestep %d out of %d during measurement number %d as the ESS ratio was %.2f.', jj, tsteps_meas, meas_counter, ess/N_fic_trajec);
                        [a1_resamp, a2_resamp, weights_resamp_log] = NPW_TM_resamp(a1_new, a2_new, weights_new_log);
                        resamp_counter = resamp_counter + 1;

                        %check the resamped have the same moments and then update.
                        a1_new = a1_resamp;                                                                   % resample fields.
                        a2_new = a2_resamp;
                        weights_new_log = weights_resamp_log;                                       % this resets weights.
                    end
               
                    % Integrate measurement signal       
                    [jx_npw_cond_avg, jy_npw_cond_avg, jz_npw_cond_avg, jx2_npw_cond_avg, jy2_npw_cond_avg, jz2_npw_cond_avg] = NPW_cond_avg(a1_old, a2_old, weights_old_log);
                    y_meas_sig_new = y_meas_sig + jz_npw_cond_avg*dtNPW_meas + dW/sqrt(2*2*gamma);                      % put in factor of 2 in signal.
                    store_meas_sig_temp(1, jj) = y_meas_sig;
                    y_meas_sig = y_meas_sig_new;                                  % update measurement signal;       

                    a1_old = a1_new;                                                      % update state and weights.
                    a2_old = a2_new;
                    weights_old_log = weights_new_log;
               end
               
               % Integrate over entire measurement pulse and make est.
                fit_slope = fit(tarrayNPW_meas', store_meas_sig_temp', 'poly1');   % fit to get slope, over dt prop to y(t).
                slope_fit = coeffvalues(fit_slope);
                slope_est = slope_fit(1);
                % jz_est = slope_est;
                
                % Using simon's two-point estimate.
                jz_est = (store_meas_sig_temp(end)-store_meas_sig_temp(1))/t_p;
                store_jz_est_temp(1, meas_counter) = jz_est;                % we make one estimate over entire meas. pulse
                
                % Using perfect measurement instead. (cond. average at end of seq.)
                % jz_perf = (sum(exp(weights_old_log).*(abs(a1_old).^2 - abs(a2_old).^2)/2))/(sum(exp(weights_old_log)));
                % store_jz_est_temp(1, meas_counter) = jz_perf;
                % store_jz_est_perf_diff(nn, meas_counter) = jz_est-jz_perf;
               if meas_counter ~= 1                                                                   % Compute jz derivative finite differences.
                     djz_dt = (store_jz_est_temp(1, meas_counter) - store_jz_est_temp(1, meas_counter-1))/(dt_meas_int);
                     fb = k_fb*djz_dt;
                     store_fb_temp(1, meas_counter) = fb;                                      % Only start applying feedback after the second measurement.
               end
               meas_counter = meas_counter+1;
               % fprintf('\nI finished measurement number %d out of %d.', meas_counter-1, num_meas);
               
        
               % evolve at the end of measurement step.
               [a1_new, a2_new] = TWA_CAHO(a1_old, a2_old, chi, kappa, dtNPW, fb, "jz_damp");          % standard evolution
               weights_new_log = weights_old_log;   
        
            % Block 3:  Hamiltonian evolution between measurement pulses.
            else
                % Hamiltonian evolution between measurement pulses. Weights  don't change but NOT equal.
                [a1_new, a2_new] = TWA_CAHO(a1_old, a2_old, chi, kappa, dtNPW, fb, "jz_damp");          % standard evolution
                weights_new_log = weights_old_log;                                                                                % weights don't change.
            end
  
            %% Compute moments computed using previous timestep values (beginning of timestep).
            if mod(ii-1, dt_samp_timesteps) == 0 
                [jx_npw_cond_avg, jy_npw_cond_avg, jz_npw_cond_avg, jx2_npw_cond_avg, jy2_npw_cond_avg, jz2_npw_cond_avg, jx3_npw_cond_avg, jy3_npw_cond_avg, jz3_npw_cond_avg] = NPW_cond_avg(a1_old, a2_old, weights_old_log);
                store_jx_npw_cond_avg_temp(1, samp_counter) = jx_npw_cond_avg;     
                store_jy_npw_cond_avg_temp(1, samp_counter) = jy_npw_cond_avg;
                store_jz_npw_cond_avg_temp(1, samp_counter) = jz_npw_cond_avg;
                store_jx2_npw_cond_avg_temp(1, samp_counter) = jx2_npw_cond_avg;
                store_jy2_npw_cond_avg_temp(1, samp_counter) = jy2_npw_cond_avg;
                store_jz2_npw_cond_avg_temp(1, samp_counter) = jz2_npw_cond_avg;
                store_jx3_npw_cond_avg_temp(1, samp_counter) = jx3_npw_cond_avg;
                store_jy3_npw_cond_avg_temp(1, samp_counter) = jy3_npw_cond_avg;
                store_jz3_npw_cond_avg_temp(1, samp_counter) = jz3_npw_cond_avg;

                % compute energy and OBDM.
                [energy_npw_cond_avg, OBDM_npw_cond] = npw_energy_OBDM_cond_avg(a1_old, a2_old, weights_old_log, chi, kappa);
                store_energy_npw_cond_avg_temp(1, samp_counter) = energy_npw_cond_avg;
                store_OBDM_npw_cond_temp(samp_counter, :, :) = OBDM_npw_cond;
                samp_counter = samp_counter + 1;
            end

            % reassign variables for next timestep.
            weights_old_log = weights_new_log;
            resamp_flag = 0;
            a1_old = a1_new;
            a2_old = a2_new;
        end
        % assign temp storage to global storage
        store_fb(nn, :) = store_fb_temp;   
        store_jz_est(nn, :) = store_jz_est_temp;
        
        store_jx_npw_cond_avg(nn, :) = store_jx_npw_cond_avg_temp;
        store_jy_npw_cond_avg(nn, :) = store_jy_npw_cond_avg_temp;
        store_jz_npw_cond_avg(nn, :) = store_jz_npw_cond_avg_temp;

        store_jx2_npw_cond_avg(nn, :) = store_jx2_npw_cond_avg_temp; 
        store_jy2_npw_cond_avg(nn, :) = store_jy2_npw_cond_avg_temp;
        store_jz2_npw_cond_avg(nn, :) = store_jz2_npw_cond_avg_temp;    

        store_jx3_npw_cond_avg(nn, :) = store_jx3_npw_cond_avg_temp; 
        store_jy3_npw_cond_avg(nn, :) = store_jy3_npw_cond_avg_temp;
        store_jz3_npw_cond_avg(nn, :) = store_jz3_npw_cond_avg_temp;    
        
        store_energy_npw_cond_avg(nn, :) = store_energy_npw_cond_avg_temp;
        store_OBDM_npw_cond_avg(nn, :, :, :) = store_OBDM_npw_cond_temp;

        fprintf('\nI resampled a total of %d times for real trajectory %d.', resamp_counter, nn)
    end
end

%% Do some averaging to get uncond. state.

store_jx_npw_uncond_avg = mean(store_jx_npw_cond_avg, 1);
store_jy_npw_uncond_avg = mean(store_jy_npw_cond_avg, 1);
store_jz_npw_uncond_avg = mean(store_jz_npw_cond_avg, 1);
store_jx2_npw_uncond_avg = mean(store_jx2_npw_cond_avg, 1);
store_jy2_npw_uncond_avg = mean(store_jy2_npw_cond_avg, 1);
store_jz2_npw_uncond_avg = mean(store_jz2_npw_cond_avg, 1);
store_jx3_npw_uncond_avg = mean(store_jx3_npw_cond_avg, 1);
store_jy3_npw_uncond_avg = mean(store_jy3_npw_cond_avg, 1);
store_jz3_npw_uncond_avg = mean(store_jz3_npw_cond_avg, 1);

store_energy_npw_uncond_avg = mean(store_energy_npw_cond_avg, 1);
store_OBDM_npw_uncond_avg = squeeze(mean(store_OBDM_npw_cond_avg, 1));
%% saving files.
% save("paper_npw_discrete_NCI.mat", '-v7.3', '-regexp', '^(?!(store_fb|store_jz_est|store_weights_log|fic_noise_array|meas_noise_array|store_meas_sig)$).')
%save("SNPW_CSS_300R_300F_new_params_all_trajecs.mat", '-v7.3', '-regexp', ['^(?!(store_fb|store_jz_est|store_weights_log|fic_noise_array|meas_noise_array|store_meas_sig' ...
 %  '|store_jx_npw_cond_avg|store_jy_npw_cond_avg|store_jz_npw_cond_avg|store_jx2_npw_cond_avg|store_jy2_npw_cond_avg|store_jz2_npw_cond_avg|store_jx3_npw_cond_avg|store_jy3_npw_cond_avg|store_jz3_npw_cond_avg|store_energy_npw_cond_avg|store_OBDM_npw_cond_avg|store_OBDM_npw_cond_avg)$).'])
% save("KNPW_mixed_200R_200F_all.mat", '-v7.3', '-regexp', ['^(?!(store_fb|store_jz_est|store_weights_log|fic_noise_array|meas_noise_array|store_meas_sig' ...
%     '|store_jx_npw_cond_avg|store_jy_npw_cond_avg|store_jz_npw_cond_avg|store_jx2_npw_cond_avg|store_jy2_npw_cond_avg|store_jz2_npw_cond_avg|store_jx3_npw_cond_avg|store_jy3_npw_cond_avg|store_jz3_npw_cond_avg|store_energy_npw_cond_avg|store_frac_npw_cond_avg)$).'])
% save("KNPW_01_sing_trajec.mat", '-v7.3', '-regexp', '^(?!(store_weights_log|fic_noise_array|meas_noise_array)$).')

save("SNPW_CSS_500R_500F_new_final.mat", '-v7.3')