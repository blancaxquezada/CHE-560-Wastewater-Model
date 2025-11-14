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
Ts = KLa5_rand(3,1) - KLa5_rand(2,1);
u1 = KLa5_rand(2:end-1,2) - KLa5_rand(2,2);
u2 = Qint_rand(2:end-1,2) - Qint_rand(2,2);
y1 = SNO2_rand(2:end-1,2) - SNO2_rand(2,2);
y2 = SO5_rand(2:end-1,2) - SO5_rand(2,2);



