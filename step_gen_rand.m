%% ===== Random step magnitudes + random durations (ALL IN DAYS) =====
% Steps occur only in [Warmup_days .. StepsEnd_days], then hold to Tend_days.
% Outputs: KLa5_step, QiZi_step (From Workspace, ZOH, time in DAYS)

%% --- User knobs ---
Tend_days      = 22;      % total horizon
StepsEnd_days  = 14;      % latest time a step may start
Warmup_days    = 5;       % flat at nominal before steps

% How many steps to place inside the window:
Nsteps_target  = 24;      % e.g., 24 steps between Warmup and StepsEnd

% Duration randomness (each plateau/step's length in days)
min_plateau_days = 0.20;  % hard minimum for any plateau
jitter_pct       = 0.30;  % ±30% jitter around the mean duration before rescale
% (durations are rescaled to exactly fill the usable window)

% Magnitude randomness (relative to nominal)
modeK = "uniform";   A_kla_max  = 0.10;    % uniform in [-0.10, +0.10]
modeQ = "uniform";   A_qizi_max = 0.08;    % uniform in [-0.08, +0.08]
% If using normal mode, set std and clamp to ±A_*_max:
A_kla_std  = 0.06;
A_qizi_std = 0.05;

% K & Q correlation for step direction
%  "independent" (default) or "synchronous" (share same sign per step)
corr_mode = "independent";

% Reproducibility
rng_seed = 42;        % [] to use current RNG state

%% --- Nominal values
kla5_ss = 84;
qizi_ss = 55338;

%% --- Timing window sanity
assert(StepsEnd_days <= Tend_days, 'StepsEnd_days must be <= Tend_days');
assert(Warmup_days   <  StepsEnd_days, 'Warmup_days must be < StepsEnd_days');

usable = StepsEnd_days - Warmup_days;     % window dedicated to steps
assert(usable > 0, 'No room for steps with current settings');

%% --- Draw random durations that exactly fill the usable window
if ~isempty(rng_seed), rng(rng_seed); end
Nsteps = max(0, floor(Nsteps_target));
if Nsteps == 0
    p = [];                                 % no steps after warmup
else
    mean_dur     = usable / Nsteps;
    raw_factors  = 1 + jitter_pct*(2*rand(1,Nsteps) - 1);  % in [1-ζ, 1+ζ]
    raw_durs     = max(min_plateau_days, mean_dur * raw_factors);

    % Rescale so total = usable (keeps min constraint, adjusts proportionally)
    raw_sum      = sum(raw_durs);
    dur          = raw_durs * (usable / raw_sum);

    % Edge case: enforce min and re-normalize again if needed
    below_min = dur < min_plateau_days;
    if any(below_min)
        % lift below-min to min; shrink others proportionally to keep sum usable
        deficit = sum(min_plateau_days - dur(below_min));
        dur(below_min) = min_plateau_days;
        if deficit > 0 && any(~below_min)
            dur(~below_min) = dur(~below_min) * ( (sum(dur(~below_min)) - deficit) / sum(dur(~below_min)) );
        end
        % final normalization
        dur = dur * (usable / sum(dur));
    end

    % Plateau start times (in days) after warmup
    p = Warmup_days + [0, cumsum(dur(1:end-1))];
end

% Time stamps for held values (in DAYS): 0 (warmup), each plateau start p, and Tend_days
time_days = [0; p(:); Tend_days];          % column, length = Nsteps+2

%% --- Random magnitudes (non-cumulative, around nominal) --------------------
draw_levels = @(mode, Amax, N, sigma) ...
    (mode=="uniform") .* (2*Amax*rand(1,N) - Amax) + ...
    (mode=="normal" ) .* max(min(sigma*randn(1,N), +Amax), -Amax);

dK = draw_levels(modeK, A_kla_max,  Nsteps, A_kla_std);
dQ = draw_levels(modeQ, A_qizi_max, Nsteps, A_qizi_std);

% Optionally synchronize signs step-by-step
if Nsteps > 0 && corr_mode == "synchronous"
    sgn = sign(2*rand(1,Nsteps)-1); sgn(sgn==0) = 1;   % ±1
    dK  = abs(dK) .* sgn;
    dQ  = abs(dQ) .* sgn;
end

% Build level sequences (length Nsteps+1 including warmup)
levK_seq = [kla5_ss,  kla5_ss*(1 + dK)];
levQ_seq = [qizi_ss,  qizi_ss*(1 + dQ)];

% Values at defined times (duplicate last value at Tend_days)
valK = [levK_seq(:); levK_seq(end)];       % column, length = Nsteps+2
valQ = [levQ_seq(:); levQ_seq(end)];

%% --- Checks
fprintf('Nsteps=%d, mean_dur=%.4g d, min=%.4g d, window=%.3f d\n', ...
        Nsteps, (usable/max(Nsteps,1)), min_plateau_days, usable);
assert(numel(time_days)==numel(valK) && numel(valK)==numel(valQ), 'Size mismatch');

%% --- From-Workspace structs (ZOH, time in DAYS)
KLa5_step.time               = time_days;
KLa5_step.signals.values     = valK;
KLa5_step.signals.dimensions = 1;

QiZi_step.time               = time_days;
QiZi_step.signals.values     = valQ;
QiZi_step.signals.dimensions = 1;

assert(abs(KLa5_step.time(end)-Tend_days)<1e-12, 'KLa5.time end != Tend_days');
assert(abs(QiZi_step.time(end)-Tend_days)<1e-12, 'QiZi.time end != Tend_days');

% %% --- Quick plots
% figure(21); clf
% stairs(KLa5_step.time, KLa5_step.signals.values,'LineWidth',1.2); grid on
% xlabel('Days'); ylabel('KLa5'); title('KLa5 random steps (random durations)');
% xline(StepsEnd_days,'--'); xlim([0 Tend_days])
% 
% figure(22); clf
% stairs(QiZi_step.time, QiZi_step.signals.values,'LineWidth',1.2); grid on
% xlabel('Days'); ylabel('QiZi'); title('QiZi random steps (random durations)');
% xline(StepsEnd_days,'--'); xlim([0 Tend_days])
% 
% % (Optional) visualize the random durations:
% if Nsteps>0
%     figure(23); clf
%     bar(dur); grid on; xlabel('Step #'); ylabel('Duration [days]');
%     title(sprintf('Random plateau durations (sum = %.3f d)', sum(dur)));
% end

% What Simulink will ingest (in DAYS):
disp([max(KLa5_step.time), numel(KLa5_step.time), max(QiZi_step.time), numel(QiZi_step.time)])
