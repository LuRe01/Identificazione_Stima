%% TEST DI BIANCHEZZA DI ANDERSON SULL'INNOVAZIONE
% Test basato sul correlogramma normalizzato:
% rho_hat(tau) = gamma_hat(tau)/gamma_hat(0)
%
% Se epsilon(t) è rumore bianco:
% sqrt(N)*rho_hat(tau) ~ N(0,1)
%
% Intervallo di confidenza:
% rho_hat(tau) in (-beta, beta), beta = z_(1-alpha/2)/sqrt(N)

clc

alpha_test = 0.05;      % livello di significatività
maxLag = 50;            % numero di ritardi del correlogramma

%% Lettura segnali da Simulink
e_acc_raw = out.e_acc.Data;

flag_acc_raw = out.flag_acc.Data;

%% Conversione in formati comodi

e_acc = convert_scalar_signal(e_acc_raw);

flag_acc = convert_scalar_signal(flag_acc_raw);

%% Selezione dei campioni effettivamente usati nella correzione

idx_acc = flag_acc > 0.5;

e_acc_valid = e_acc(idx_acc);

%% Test sulle tre componenti di innovazione

res_acc = anderson_whiteness_test(e_acc_valid, maxLag, alpha_test, ...
    'Innovazione accelerometro');

%% Tabella riassuntiva

Nome = {'acc'};

N = res_acc.N;
Beta = res_acc.beta;
NumFuori = res_acc.num_out;
PercFuori = res_acc.perc_out;
TestSuperato = res_acc.test_passed;

T = table(Nome, N, Beta, NumFuori, PercFuori, TestSuperato);

disp(' ');
disp('RISULTATI TEST DI BIANCHEZZA DI ANDERSON');
disp(T);

function result = anderson_whiteness_test(epsilon, maxLag, alpha, signal_name)

    epsilon = epsilon(:);
    epsilon = epsilon(~isnan(epsilon));

    % Rimozione della media campionaria
    epsilon = epsilon - mean(epsilon);

    N = length(epsilon);

    if N <= maxLag + 1
        error('Numero di campioni insufficiente per il test su %s.', signal_name);
    end

    % gamma_hat(0): varianza stimata
    gamma0 = (1/N) * sum(epsilon.^2);

    if gamma0 <= eps
        error('Varianza nulla o troppo piccola nel segnale %s.', signal_name);
    end

    % Calcolo del correlogramma normalizzato
    rho = zeros(maxLag,1);

    for tau = 1:maxLag
        gamma_tau = (1/N) * sum(epsilon(1:N-tau) .* epsilon(1+tau:N));
        rho(tau) = gamma_tau / gamma0;
    end

    % Intervallo di confidenza
    % Per alpha = 0.05, z = 1.96
    z = sqrt(2) * erfinv(1 - alpha);
    beta = z / sqrt(N);

    % Conteggio dei valori fuori banda
    out_of_bounds = abs(rho) > beta;

    num_out = sum(out_of_bounds);
    perc_out = num_out / maxLag;

    % Criterio del test:
    % superato se la frazione di valori fuori banda è minore di alpha
    test_passed = perc_out < alpha;

    %% Plot correlogramma

    figure;
    stem(1:maxLag, rho, 'filled');
    hold on;
    yline(beta, '--', 'Limite superiore');
    yline(-beta, '--', 'Limite inferiore');
    grid on;
    xlabel('\tau');
    ylabel('\hat{\rho}(\tau)');
    title(['Test di bianchezza Anderson - ', signal_name]);

    %% Stampa risultati

    fprintf('\n%s\n', signal_name);
    fprintf('N = %d\n', N);
    fprintf('maxLag = %d\n', maxLag);
    fprintf('beta = %.6f\n', beta);
    fprintf('Numero valori fuori banda = %d su %d\n', num_out, maxLag);
    fprintf('Percentuale fuori banda = %.2f %%\n', 100*perc_out);

    if test_passed
        fprintf('Esito: test superato. Innovazione compatibile con rumore bianco.\n');
    else
        fprintf('Esito: test NON superato. Innovazione autocorrelata.\n');
    end

    result.N = N;
    result.rho = rho;
    result.beta = beta;
    result.num_out = num_out;
    result.perc_out = perc_out;
    result.test_passed = test_passed;

end


function x = convert_scalar_signal(xraw)

    x = squeeze(xraw);

    if isrow(x)
        x = x';
    end

end


function X = convert_vector_signal(Xraw, nv)

    sz = size(Xraw);

    % Caso N x nv
    if ismatrix(Xraw) && sz(2) == nv
        X = Xraw;
        return;
    end

    % Caso nv x N
    if ismatrix(Xraw) && sz(1) == nv
        X = Xraw';
        return;
    end

    % Caso 1 x nv x N
    if numel(sz) == 3 && sz(1) == 1 && sz(2) == nv
        X = squeeze(Xraw);
        if size(X,1) == nv
            X = X';
        end
        return;
    end

    % Caso nv x 1 x N
    if numel(sz) == 3 && sz(1) == nv && sz(2) == 1
        X = squeeze(Xraw)';
        return;
    end

    % Caso N x 1 x nv
    if numel(sz) == 3 && sz(2) == 1 && sz(3) == nv
        X = squeeze(Xraw);
        return;
    end

    error('Formato vettoriale non riconosciuto. Size: [%s]', num2str(sz));

end

