
% Input Steps
    % Qint
    subplot(2,2,1);
    stairs(Qint_log.time, Qint_log.signals.values); grid on
    xlabel('Time [days]'); ylabel('Q_{intr} [m^3/d]'); 
    %title('PRBS for Q_{intr}'); %PRBS
    title('Random Signals for Q_{intr}'); %Random
    % KLa5
    subplot(2,2,2);
    stairs(Kla5_log.time, Kla5_log.signals.values); grid on
    xlabel('Time [days]'); ylabel('K_{La,5} [d^{-1}]'); 
    %title('PRBS for K_{La,5}'); %PRBS
    title('Random Signals for K_{La,5}'); %PRBS


% Output Responses
    % SNO2
    subplot(2,2,3);
    plot(SNO2.time, SNO2.signals.values); grid on; title('SNO2')
    xlabel('Time [days]'); ylabel('S_NO_2 [mg N/l]'); title('S_NO in Reactor 2 Response');
    % SO5
    subplot(2,2,4);
    plot(SO5.time,  SO5.signals.values);  grid on; title('SO5')
    xlabel('Time [days]'); ylabel('S_O_5 [mg -COD/l]'); title('S_O in Reactor 5 Response');


% % Data PRBS
% KLa5_prbs = [Kla5_log.time, Kla5_log.signals.values];
% SO5_prbs = [SO5.time, SO5.signals.values];
% 
% Qint_prbs = [Qint_log.time, Qint_log.signals.values];
% SNO2_prbs = [SNO2.time, SNO2.signals.values];

% Data Random
KLa5_rand = [Kla5_log.time, Kla5_log.signals.values];
SO5_rand = [SO5.time, SO5.signals.values];

Qint_rand = [Qint_log.time, Qint_log.signals.values];
SNO2_rand = [SNO2.time, SNO2.signals.values];

