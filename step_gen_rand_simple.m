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
    % Find how many steps fit
    idx = find(time_days <= Tend_days, 1, 'last');
    time_days = time_days(1:idx);
    time_days(end) = Tend_days;  % Set last time to Tend_days
    Nsteps = length(time_days) - 1;  % Update actual number of steps
else
    % Add final time point
    time_days = [time_days; Tend_days];
end

%% --- Generate Random Step Heights ---
% Random deviations from nominal (uniform distribution)
dK = kla5_max_dev * (2*rand(1, Nsteps) - 1);  % [-kla5_max_dev, +kla5_max_dev]
dQ = qizi_max_dev * (2*rand(1, Nsteps) - 1);  % [-qizi_max_dev, +qizi_max_dev]

% Calculate actual values
valK = kla5_ss * (1 + dK);   % KLa5 values for each step
valQ = qizi_ss * (1 + dQ);   % QiZi values for each step

% Add initial value at time 0 (nominal)
valK = [kla5_ss; valK(:)];
valQ = [qizi_ss; valQ(:)];

% Duplicate last value at final time (ZOH hold)
valK = [valK; valK(end)];
valQ = [valQ; valQ(end)];

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

