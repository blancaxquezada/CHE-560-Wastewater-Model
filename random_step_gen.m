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
    % valK/valQ need: [initial_val, step1_val, step2_val, ..., stepN_val, steady_state] = Nsteps+2 elements
    % Final value returns to steady state at Tend_days
    valK = [kla5_ss; valK_steps(:); kla5_ss];
    valQ = [qizi_ss; valQ_steps(:); qizi_ss];
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

%% --- Run Simulink ---
olasm1init;
olsettler1dinit;

load constinfluent;
load dryinfluent;
load raininfluent;
load storminfluent;

sim('openloop_steps14.slx')

%% --- Plot Data from Simulink

% Input Steps
    % Qint
    subplot(2,2,1);
    stairs(Qint_log.time, Qint_log.signals.values); grid on
    xlabel('Time [days]'); ylabel('Q_{intr} [m^3/d]'); 
    title('Random Signals for Q_{intr}');
    % KLa5
    subplot(2,2,2);
    stairs(Kla5_log.time, Kla5_log.signals.values); grid on
    xlabel('Time [days]'); ylabel('K_{La,5} [d^{-1}]'); 
    title('Random Signals for K_{La,5}');

% Output Responses
    % SNO2
    subplot(2,2,3);
    plot(SNO2.time, SNO2.signals.values); grid on; title('SNO2')
    xlabel('Time [days]'); ylabel('S_NO_2 [mg N/l]'); title('S_NO in Reactor 2 Response');
    % SO5
    subplot(2,2,4);
    plot(SO5.time,  SO5.signals.values);  grid on; title('SO5')
    xlabel('Time [days]'); ylabel('S_O_5 [mg -COD/l]'); title('S_O in Reactor 5 Response');

% Data Random
KLa5_rand = [Kla5_log.time, Kla5_log.signals.values];
SO5_rand = [SO5.time, SO5.signals.values];

Qint_rand = [Qint_log.time, Qint_log.signals.values];
SNO2_rand = [SNO2.time, SNO2.signals.values];

%% --- Build Model ---
% Extract sampling time (assuming uniform sampling)
Ts = KLa5_rand(3,1) - KLa5_rand(2,1);

% Compute deviations from initial steady-state (as instructed by professor)
u1 = KLa5_rand(2:end-1,2) - KLa5_rand(2,2);  % KLa5 deviation
u2 = Qint_rand(2:end-1,2) - Qint_rand(2,2);   % Qint deviation
y1 = SNO2_rand(2:end-1,2) - SNO2_rand(2,2);   % SNO2 deviation
y2 = SO5_rand(2:end-1,2) - SO5_rand(2,2);     % SO5 deviation

% Combine inputs and outputs
U = [u1, u2];  % Input matrix: [u1, u2] - 2 inputs
Y = [y1, y2];  % Output matrix: [y1, y2] - 2 outputs

% Create iddata object for System Identification Toolbox
% Format: iddata(outputs, inputs, sampling_time)
data = iddata(Y, U, Ts, 'TimeUnit', 'days');

% Display data info
fprintf('Data prepared:\n');
fprintf('  Sampling time: %.6f days\n', Ts);
fprintf('  Number of samples: %d\n', length(u1));
fprintf('  Inputs: KLa5, Qint\n');
fprintf('  Outputs: SNO2, SO5\n\n');

%% --- State-Space Model Identification ---
% Method 1: Using ssest() - State-space estimation (recommended)
% This estimates a state-space model of specified order
nx = 4;  % State order (number of states) - adjust as needed
% Try different orders: 2, 3, 4, 5, 6, etc. to find best fit

fprintf('Identifying state-space model with order = %d...\n', nx);
sys_ss = ssest(data, nx, 'Ts', Ts);

% Display model information
fprintf('\nIdentified State-Space Model:\n');
fprintf('  Order: %d\n', size(sys_ss.A, 1));
fprintf('  Inputs: %d (KLa5, Qint)\n', size(sys_ss.B, 2));
fprintf('  Outputs: %d (SNO2, SO5)\n', size(sys_ss.C, 1));

% Display fit percentages for each output
% FitPercent is a vector: [fit_output1, fit_output2, ...]
fprintf('  Fit to data:\n');
fprintf('    Output 1 (SNO2): %.2f%%\n', sys_ss.Report.Fit.FitPercent(1));
fprintf('    Output 2 (SO5):  %.2f%%\n', sys_ss.Report.Fit.FitPercent(2));
fprintf('    Average fit:     %.2f%%\n', mean(sys_ss.Report.Fit.FitPercent));

% Note: The fit percentage is calculated as:
% Fit = 100 * (1 - norm(y_measured - y_predicted) / norm(y_measured - mean(y_measured)))
% This measures how much better the model is compared to just using the mean value.
% The compare() plot may show slightly different values due to different normalization
% or calculation methods used in the visualization.

%% --- Model Validation ---
% Compare model output with measured data
figure;
compare(data, sys_ss);
title('Model Validation: Measured vs. Model Output');

% Plot step response
figure;
step(sys_ss);
title('Step Response of Identified State-Space Model');
grid on;

% Display state-space matrices
fprintf('\nState-Space Matrices:\n');
fprintf('A matrix (state dynamics):\n');
disp(sys_ss.A);
fprintf('B matrix (input to state):\n');
disp(sys_ss.B);
fprintf('C matrix (state to output):\n');
disp(sys_ss.C);
fprintf('D matrix (direct feedthrough):\n');
disp(sys_ss.D);

% Save the model
save('identified_ss_model.mat', 'sys_ss', 'data', 'Ts');
fprintf('\nModel saved to: identified_ss_model.mat\n');



