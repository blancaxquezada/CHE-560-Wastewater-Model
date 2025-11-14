%% ===== Simplified Random Steps Generator =====
% Generates random step heights and durations over a 14-day period
% Outputs: KLa5_step, QiZi_step (From Workspace, ZOH, time in DAYS)

%% --- User Settings ---
Tend_days = 14;              % Total simulation time (days)
Nsteps = 24;                 % Number of steps to generate

% Random duration settings (days)
min_dur = 0.2;               % Minimum duration per step
max_dur = 1.0;               % Maximum duration per step

% Random magnitude settings (relative to nominal, as fraction)
kla5_max_dev = 0.10;         % KLa5 varies ±10% from nominal
qizi_max_dev = 0.08;         % QiZi varies ±8% from nominal

% Nominal values
kla5_ss = 84;
qizi_ss = 55338;

% Reproducibility
rng_seed = [];               % [] = different random steps each run (for training data variety)
if ~isempty(rng_seed), rng(rng_seed); end

%% --- Generate Random Durations ---
% Random durations between min and max
dur = min_dur + (max_dur - min_dur) * rand(1, Nsteps);

% Calculate step start times (cumulative sum of durations)
% Start at time 0
time_days = [0; cumsum(dur(:))];

% If we exceed Tend_days, truncate and adjust last duration
if time_days(end) > Tend_days
    % Find how many step start times fit (excluding initial 0)
    idx = find(time_days(2:end) <= Tend_days, 1, 'last');  % idx is number of steps that fit
    if isempty(idx) || idx == 0
        % No steps fit, just use initial and final
        time_days = [0; Tend_days];
        Nsteps = 0;
        dur = [];
    else
        % Keep idx steps, add Tend_days at end
        time_days = [0; time_days(2:idx+1); Tend_days];  % [0, step1, step2, ..., step_idx, Tend_days]
        Nsteps = idx;
        dur = dur(1:Nsteps);
    end
else
    % Add final time point
    time_days = [time_days; Tend_days];
end

%% --- Generate Random Step Heights ---
% Random deviations from nominal (uniform distribution)
% Use the actual Nsteps (may have been truncated)
if Nsteps > 0
    dK = kla5_max_dev * (2*rand(1, Nsteps) - 1);  % [-kla5_max_dev, +kla5_max_dev]
    dQ = qizi_max_dev * (2*rand(1, Nsteps) - 1);  % [-qizi_max_dev, +qizi_max_dev]
    
    % Calculate actual values for each step
    valK_steps = kla5_ss * (1 + dK);   % KLa5 values for each step
    valQ_steps = qizi_ss * (1 + dQ);   % QiZi values for each step
    
    % Build complete value arrays to match time_days structure
    % time_days has: [0, step1_time, step2_time, ..., stepN_time, Tend_days] = Nsteps+2 elements
    % valK/valQ need: [initial_val, step1_val, step2_val, ..., stepN_val, stepN_val] = Nsteps+2 elements
    valK = [kla5_ss; valK_steps(:); valK_steps(end)];
    valQ = [qizi_ss; valQ_steps(:); valQ_steps(end)];
else
    % No steps, just initial and final values (both at nominal)
    valK = [kla5_ss; kla5_ss];
    valQ = [qizi_ss; qizi_ss];
end

% Ensure dimensions match
assert(length(time_days) == length(valK), 'Dimension mismatch: time_days and valK must have same length');
assert(length(time_days) == length(valQ), 'Dimension mismatch: time_days and valQ must have same length');

%% --- Create Simulink Structures ---
KLa5_step.time               = time_days;
KLa5_step.signals.values     = valK;
KLa5_step.signals.dimensions = 1;

QiZi_step.time               = time_days;
QiZi_step.signals.values     = valQ;
QiZi_step.signals.dimensions = 1;

%% --- Plots ---
% KLa5 step plot
figure(21); clf
stairs(KLa5_step.time, KLa5_step.signals.values, 'LineWidth', 1.2); grid on
xlabel('Days'); ylabel('KLa5'); 
title('KLa5 random steps (random durations)');
xlim([0 Tend_days])

% QiZi step plot
figure(22); clf
stairs(QiZi_step.time, QiZi_step.signals.values, 'LineWidth', 1.2); grid on
xlabel('Days'); ylabel('QiZi'); 
title('QiZi random steps (random durations)');
xlim([0 Tend_days])

% Duration bar plot (show actual durations used)
if Nsteps > 0
    % Calculate actual durations from time_days (excluding final time point)
    actual_dur = diff(time_days(1:end-1));  % Differences between consecutive time points
    figure(23); clf
    bar(actual_dur); grid on
    xlabel('Step #'); ylabel('Duration [days]');
    title(sprintf('Random plateau durations (sum = %.3f d)', sum(actual_dur)));
end

