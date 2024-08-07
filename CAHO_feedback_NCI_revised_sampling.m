% numberOfWorkers = 12;                       % Maximum number of cores.
% parpool('local', numberOfWorkers);
%% Geometry
t0 = 0;
tMax = 50;
tstepsTW = 1e5+1;
tarrayTW = linspace(t0, tMax, tstepsTW);
dtTW = tarrayTW(2)-tarrayTW(1);
t_fb_on_TW = tstepsTW;
t_fb_off_TW = tstepsTW;

% Measurement params
% dt_meas_int = 0.4;
dt_meas_int = 0.8;                                                                                                                       % interval between measurements in non-dim time.
dt_meas_int_timesteps = dt_meas_int/dtTW;
tarrayTW_meas = tarrayTW(dt_meas_int_timesteps:dt_meas_int_timesteps:end);                       % Array of measurement times.
num_meas = length(tarrayTW_meas);                                                                                           % total number of measurements.

tstepsTW_samp = 5000+1;                                       % number of points to sample.  Pick divisor.
tarrayTW_samp = linspace(t0, tMax, tstepsTW_samp);
dt_samp_timestepsTW = (tstepsTW-1)/(tstepsTW_samp-1);
dtTW_samp = tarrayTW_samp(2)-tarrayTW_samp(1);

% Set up SU2 timesteps = number of samples.                                                                              
tstepsSU2 = tstepsTW_samp;
tarraySU2 = linspace(t0, tMax, tstepsSU2);                                                                              % No need for explicit sampling array.
dtSU2 = tarraySU2(2)-tarraySU2(1);
dt_meas_timestepsSU2 = dt_meas_int/dtSU2;                                                 
tarraySU2_meas = tarraySU2(dt_meas_timestepsSU2:dt_meas_timestepsSU2:end);
t_fb_on_SU2 = 1;
t_fb_off_SU2 = tstepsTW_samp;

% disp([tarrayNPW(dt_meas_int_timesteps), tarrayNPW_meas(1), tarrayNPW_samp(dt_meas_int_timesteps/dt_samp_timesteps+1), dtNPW, num_meas])
%% Misc. number and int. strengths, Initial state
% SKZTW params
N_atoms = 100;
N_light = 1e7;
b1_0 = sqrt(N_light);
chi = 0.01;
kappa = 0.09;  
% chi = 0;
% kappa = 0;
k_fb = 0.10;                                                                       % feedback gain parameter;
% lambda = 0.8e-4;                                                              % SKZTW meas param.
lambda = 0.6e-4;                                                              % SKZTW meas param.
threshold = (2*pi/N_atoms/2)/(k_fb*N_atoms*dtTW);                      % set threshold for strength of feedback term.  implements rotation.  

% Setting up cond. and fic. trajectories.
N_trajecTW = 2000;
N_trajecSU2 = 5000;

% State initilisation
state = "CSS_custom";
fb_type = "jz_damp";
n10 = 20;
n20 = 80;
phi10 = pi/2;
phi20 = 0;
a10 = sqrt(20);
a20 = sqrt(80);
normalise = 1;
[~, ~, b1_old, jx0, jy0, jz0, varjx, varjy, varjz, rho_initial] = TWSU2init(N_atoms, N_light, N_trajecTW, state, a10, a20, normalise);
[n1_old, n2_old, phi1_old, phi2_old, a1_old, a2_old] = npw_init_TM(n10, n20, phi10, phi20, N_trajecTW, state, normalise, 1);
[mean_jx, mean_jy, mean_jz, var_jx, var_jy, var_jz] = computeJv2(a1_old, a2_old);

% SU2 Evolution
[Jx, Jy, Jz, Jp, Jm, J2, Na, Nb, Jx2, Jy2, Jz2] = su2matrices(N_atoms);
H2M_free = chi*(2*Jz2 + Na*Nb) + 2*kappa*Jx;
Ut_free = expm(-1i*H2M_free*dtSU2);
Ut_free_t = Ut_free';

% debugging flags
backaction_on = 1;                  % There is no backaction without noise!  
noise_on = 1;
perf_meas = 0;

store_index = 1;
fb = zeros(1, N_trajecTW);      % Initialise feedback.
%% Bulk Storage
% =============================================================== TW storage
store_fb_TW = zeros(num_meas, N_trajecTW);                                      % feedback is the same between measurements.

store_jx_TW_mean = zeros(1, tstepsTW_samp);
store_jy_TW_mean = zeros(1, tstepsTW_samp);
store_jz_TW_mean = zeros(1, tstepsTW_samp);
store_jx_TW_var = zeros(1, tstepsTW_samp);
store_jy_TW_var = zeros(1, tstepsTW_samp);
store_jz_TW_var = zeros(1, tstepsTW_samp);

store_jx_TW_skew = zeros(1, tstepsTW_samp);
store_jy_TW_skew = zeros(1, tstepsTW_samp);
store_jz_TW_skew = zeros(1, tstepsTW_samp);

% storing estimates
store_jz_est_TW = zeros(num_meas, N_trajecTW);
store_OBDM_TW_mean = zeros(tstepsTW_samp, 2, 2);
store_OBDM_TW_trajec = zeros(tstepsTW_samp, N_trajecTW, 2, 2);                  % store for bootstrapping.

% single trajectory storage for errors etc.
store_jx_TW_trajec = zeros(tstepsTW_samp, N_trajecTW);
store_jy_TW_trajec = zeros(tstepsTW_samp, N_trajecTW);
store_jz_TW_trajec = zeros(tstepsTW_samp, N_trajecTW);

store_jx2_TW_trajec = zeros(tstepsTW_samp, N_trajecTW);
store_jy2_TW_trajec = zeros(tstepsTW_samp, N_trajecTW);
store_jz2_TW_trajec = zeros(tstepsTW_samp, N_trajecTW);

store_jx3_TW_trajec = zeros(tstepsTW_samp, N_trajecTW);
store_jy3_TW_trajec = zeros(tstepsTW_samp, N_trajecTW);
store_jz3_TW_trajec = zeros(tstepsTW_samp, N_trajecTW);

% ================================ SU2 Storage ====================================
store_jx_SU2_trajec = zeros(tstepsSU2, N_trajecSU2);
store_jy_SU2_trajec = zeros(tstepsSU2, N_trajecSU2);
store_jz_SU2_trajec = zeros(tstepsSU2, N_trajecSU2);
store_jx2_SU2_trajec = zeros(tstepsSU2, N_trajecSU2);
store_jy2_SU2_trajec = zeros(tstepsSU2, N_trajecSU2);
store_jz2_SU2_trajec = zeros(tstepsSU2, N_trajecSU2);          

store_jx3_SU2_trajec = zeros(tstepsSU2, N_trajecSU2);
store_jy3_SU2_trajec = zeros(tstepsSU2, N_trajecSU2);
store_jz3_SU2_trajec = zeros(tstepsSU2, N_trajecSU2);

%OBDM stuff
store_rho_SU2_uncond = zeros(tstepsSU2, N_atoms+1, N_atoms+1);        %% Just store density matrix at all times really.
store_rho_SU2_mean_final = zeros(N_atoms+1, N_atoms+1);     % store final density matrix.
store_OBDM_SU2_mean = zeros(tstepsSU2, 2, 2);                    % two modes, literally just two.

% Storing energies.
store_TW_energies = zeros(1, tstepsTW_samp);
store_TW_energies_trajec = zeros(tstepsTW_samp, N_trajecTW);
store_SU2_energies_mean = zeros(1, tstepsSU2);
%% SZKTW Block
sim_SKZTW = 1;
store_index_TW = 1;
meas_counter_TW = 1;
if sim_SKZTW == 1
    for ii = 1:tstepsTW
                % Block 1:  When feedback/meas is off.
                if (ii < t_fb_on_TW) || (ii > t_fb_off_TW)
                    fb = zeros(1, N_trajecTW);
                    [a1_new, a2_new] = TWA_CAHO(a1_old, a2_old, chi, kappa, dtTW, fb, fb_type);
                 
                % Block 2: Make measurement.
                elseif (ii >= t_fb_on_TW) && (mod(ii, dt_meas_int_timesteps) == 0) 
                    % Entangle and Estimate
                    b1_old = b1_0 + noise_on*1/2*(randn(1, N_trajecTW) + 1i*randn(1, N_trajecTW));
                    [a1_entang, a2_entang, b1_entang] = entanglement(a1_old, a2_old, b1_old, lambda, N_light, backaction_on);  
                    
                    % Case for perfect measurement.
                    if perf_meas == 0
                        jz_est = measure_jz(b1_entang, N_light, lambda);
                    elseif perf_meas == 1      % magical perfect estimate.
                        jz_est = 1/2*(abs(a1_old).^2 - abs(a2_old.^2));
                    end
                    store_jz_est_TW(meas_counter_TW, :) = jz_est;

                    if meas_counter_TW ~= 1
                          djz_dt = (store_jz_est_TW(meas_counter_TW, :) - store_jz_est_TW(meas_counter_TW-1, :))/(dt_meas_int);
                          fb = k_fb*djz_dt;
                          fb = 0*min(abs(fb), threshold).*sign(fb);                      % capping threshold on potential.  this *dtTW << single rotation.
                          store_fb_TW(meas_counter_TW, :) = fb;                 % convention is period after each measurement.
                          % disp([ii, max(fb-fb1)])
                    end             
                    % Evolve at end of measurement step.
                    meas_counter_TW = meas_counter_TW+1;                    % update counter.
                    [a1_new, a2_new] = TWA_CAHO(a1_entang, a2_entang, chi, kappa, dtTW, fb, fb_type); 
                
                %Block 3:  Hamiltonian evolution between measurements.
                else     
                    [a1_new, a2_new] = TWA_CAHO(a1_old, a2_old, chi, kappa, dtTW, fb, fb_type);     
                end
                
                %% Compute moments computed using previous timestep values (beginning of timestep).
                if (mod(ii-1, dt_samp_timestepsTW) == 0) 
                    [sing_jx, sing_jy, sing_jz, jx2, jy2, jz2, jx3, jy3, jz3, jx4, jy4, jz4] = computeJ(a1_old, a2_old);
                    store_jx_TW_trajec(store_index_TW, :) = sing_jx;
                    store_jy_TW_trajec(store_index_TW, :) = sing_jy;
                    store_jz_TW_trajec(store_index_TW, :) = sing_jz;

                    store_jx2_TW_trajec(store_index_TW, :) = jx2;                % need these for binning.
                    store_jy2_TW_trajec(store_index_TW, :) = jy2;
                    store_jz2_TW_trajec(store_index_TW, :) = jz2;
                
                    store_jx3_TW_trajec(store_index_TW, :) = jx3;
                    store_jy3_TW_trajec(store_index_TW, :) = jy3;
                    store_jz3_TW_trajec(store_index_TW, :) = jz3;
                 
                    % Compute using trajec at initial step.
                    [mean_jx, mean_jy, mean_jz, var_jx, var_jy, var_jz] = computeJv2(a1_old, a2_old);
                    store_jx_TW_mean(1, store_index_TW) = mean_jx;
                    store_jy_TW_mean(1, store_index_TW) = mean_jy;
                    store_jz_TW_mean(1, store_index_TW) = mean_jz;
                    store_jx_TW_var(1, store_index_TW) = var_jx;
                    store_jy_TW_var(1, store_index_TW) = var_jy;
                    store_jz_TW_var(1, store_index_TW) = var_jz;
                    store_jx_TW_skew(1, store_index_TW) = (mean((1/2*(a1_old.*conj(a2_old) + a2_old.*conj(a1_old))).^3)-3*mean_jx*var_jx - mean_jx^3)/(var_jx^(3/2));
                    store_jy_TW_skew(1, store_index_TW) = (mean((1i/2*(conj(a1_old).*a2_old - conj(a2_old).*a1_old)).^3)-3*mean_jy*var_jy - mean_jy^3)/(var_jy^(3/2));
                    store_jz_TW_skew(1, store_index_TW) = (mean((1/2*(a1_old.*conj(a1_old) - a2_old.*conj(a2_old))).^3)-3*mean_jz*var_jz - mean_jz^3)/(var_jz^(3/2));
           
                    % COMPUTE OBDM and energies.
                    OBDM = twomode_OBDM_TW(a1_old, a2_old);             % this is actual OBDM averaged over trajec.        
                    store_OBDM_TW_mean(store_index_TW, :, :) = OBDM;
                    OBDM_single = OBDM_TW_trajec(a1_old, a2_old);       % this is OBDM from single trajec.  for bootstrapping.
                    store_OBDM_TW_trajec(store_index_TW, :, :, :) = OBDM_single;
                    [energy_mean, energy_trajecs] = tw_energy(a1_old, a2_old, chi, kappa);
                    store_TW_energies(1, store_index_TW) = energy_mean;
                    store_TW_energies_trajec(store_index_TW, :) = energy_trajecs;
                    store_index_TW = store_index_TW + 1;        
                end
   
                % Reassign variables for next step.
                a1_old = a1_new;
                a2_old = a2_new;
    end
end
%% SSE Exact Block
sim_SSE = 0;
if sim_SSE == 1
    parfor nn = 1:N_trajecSU2
        % Initialise state and temp storages for parfor.
        rho_init = rho_initial;
        fb = 0;
        meas_counter_SU2 = 1;
        % temp storage.
        store_energy_SU2_temp = zeros(1, tstepsSU2);
        store_OBDM_SU2_temp = zeros(tstepsSU2, 2, 2);
        store_jz_est_SU2_temp = zeros(1, num_meas);
        store_rho_temp = zeros(1, N_atoms+1, N_atoms+1);
        Ut = Ut_free;
        Ut_t = Ut';

        for ii = 1:tstepsSU2
            % disp(ii)
            % Block 1: no meas or feedback.
            if (ii < t_fb_on_SU2) || (ii > t_fb_off_SU2)
               rho_init_new = Ut_free*rho_init*Ut_free_t;       

            % Block 2:  Measurement + evolve.
            elseif (ii >= t_fb_on_SU2) && (mod(ii, dt_meas_timestepsSU2) == 0)
                % Measure sample and feedback
                rho_diag = sqrt(diag(rho_init));     % coefficients are square root of diagonal.
                [y_meas, y_array, prob_y, jz_y_est, prob_jz_est] = samplerY(rho_diag, 1, N_atoms, N_light, lambda);
                rho_measure = measure_rho(rho_init, y_meas, N_light, lambda);  
                jz_est = 1/lambda*asin(y_meas/(2*sqrt(N_light)));
                store_jz_est_SU2_temp(1, meas_counter_SU2) = jz_est;
                % disp('meased')
                % construct feedback.
                 if meas_counter_SU2 ~=1
                       djz_dt = (store_jz_est_SU2_temp(1, meas_counter_SU2) - store_jz_est_SU2_temp(1, meas_counter_SU2-1))/(dt_meas_int);
                       fb = k_fb*djz_dt;
                       H2M = H2M_free + (-fb*Jz);                                      % reconstruct Hamiltonian.
                       Ut = expm(-1i*H2M*dtSU2);
                       Ut_t = Ut';          
                 end
                % Evolve at end of meas step.
                meas_counter_SU2 = meas_counter_SU2+1;
                rho_init_new = Ut*rho_measure*Ut_t;
            
            % Block 3: Hamiltonian evolution.
            else
                rho_init_new = Ut*rho_init*Ut_t;
            end
            %% Compute moments at start of every timestep.
            [jx, jy, jz, j2, jx2, jy2, jz2, jx3, jy3, jz3, jx4, jy4, jz4] = computeJ_SU2rho(rho_init, Jx, Jy, Jz, J2, Jx2, Jy2, Jz2);
            store_jx_SU2_trajec(ii, nn) = jx; 
            store_jy_SU2_trajec(ii, nn) = jy;
            store_jz_SU2_trajec(ii, nn) = jz;
            store_jx2_SU2_trajec(ii, nn) = jx2;
            store_jy2_SU2_trajec(ii, nn) = jy2;
            store_jz2_SU2_trajec(ii, nn) = jz2;
            store_jx3_SU2_trajec(ii, nn) = jx3;
            store_jy3_SU2_trajec(ii, nn) = jy3;
            store_jz3_SU2_trajec(ii, nn) = jz3;

            % OBDM and Energies
            OBDM_SU2 = twomode_OBDM(rho_init, Na, Nb, Jp, Jm);
            store_OBDM_SU2_temp(ii, :, :) = OBDM_SU2;
            store_energy_SU2_temp(1, ii) = expval(H2M_free, rho_init);                           % Evaluate with respect to free Hamiltonian.

            % Reassign state at end of timestep.
            rho_init = rho_init_new;         
            % disp(ii)     
        end
        % Assign temp storages to global.
        store_OBDM_SU2_mean = store_OBDM_SU2_mean + 1/N_trajecSU2*store_OBDM_SU2_temp;
        store_SU2_energies_mean = store_SU2_energies_mean + 1/N_trajecSU2*store_energy_SU2_temp;
        store_rho_SU2_uncond = store_rho_SU2_uncond+1/N_trajecSU2*store_rho_temp;
        disp([nn])
    end
end

%% Compute some means for SU2 and also standard sampling errors.
store_jx_SU2_mean = mean(store_jx_SU2_trajec, 2);
store_jy_SU2_mean = mean(store_jy_SU2_trajec, 2);
store_jz_SU2_mean = mean(store_jz_SU2_trajec, 2);
store_jx_SU2_var = mean(store_jx2_SU2_trajec, 2) - store_jx_SU2_mean.^2;
store_jy_SU2_var = mean(store_jy2_SU2_trajec, 2) - store_jy_SU2_mean.^2;
store_jz_SU2_var = mean(store_jz2_SU2_trajec, 2) - store_jz_SU2_mean.^2;

store_jx3_SU2_skew = (mean(store_jx3_SU2_trajec, 2) - 3*store_jx_SU2_mean.*store_jx_SU2_var - store_jx_SU2_mean.^3)./(store_jx_SU2_var.^(3/2));
store_jy3_SU2_skew = (mean(store_jy3_SU2_trajec, 2) - 3*store_jy_SU2_mean.*store_jy_SU2_var - store_jy_SU2_mean.^3)./(store_jy_SU2_var.^(3/2));
store_jz3_SU2_skew = (mean(store_jz3_SU2_trajec, 2) - 3*store_jz_SU2_mean.*store_jz_SU2_var - store_jz_SU2_mean.^3)./(store_jz_SU2_var.^(3/2));

% standard errors in mean and variance
[tw_jx_std, tw_varjx_std, tw_skewjx_std] = bootstrapperTW(store_jx_TW_trajec, store_jx2_TW_trajec, store_jx3_TW_trajec, N_atoms);
[tw_jy_std, tw_varjy_std, tw_skewjy_std] = bootstrapperTW(store_jy_TW_trajec, store_jy2_TW_trajec, store_jy3_TW_trajec, N_atoms);
[tw_jz_std, tw_varjz_std, tw_skewjz_std] = bootstrapperTW(store_jz_TW_trajec, store_jz2_TW_trajec, store_jz3_TW_trajec, N_atoms);

%% Save stuff
% save("SKZTW_TW_only_meas_only.mat", '-v7.3', '-regexp', '^(?!(store_jx2_SU2_trajec|store_jy2_SU2_trajec|store_jz2_SU2_trajec|store_jx2_TW_trajec|store_jy2_TW_trajec|store_jz2_TW_trajec|store_jx_SU2_trajec|store_jy_SU2_trajec|store_jz_SU2_trajec|store_jx3_TW_trajec|store_jy3_TW_trajec|store_jz3_TW_trajec|store_jx4_TW_trajec|store_jy4_TW_trajec|store_jz4_TW_trajec)$).')
% save("SKZTW_CSS_v2_TW_only", '-v7.3');
% save("SKZTW_CSS_v2_TW_only", '-v7.3', '-regexp', '^(?!(store_fb_TW)$).')
