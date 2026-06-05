%% TEST DI BIANCHEZZA - INNOVAZIONE MAGNETOMETRICA ANGOLARE

clc

alpha_test = 0.05;
maxLag = 50;

%% Lettura segnali

% Misura magnetometrica grezza [m_x; m_y]
y_magn_raw = out.y_magn_raw.Data;

% Stima predetta EKF
X_hat_pred_raw = out.X_hat_pred.Data;

% Flag misura magnetometrica usata nella correzione
flag_magn_raw = out.flag_magn.Data;

%% Conversione formati

y_magn = convert_vector_signal(y_magn_raw, 2);       % N x 2
X_hat_pred = convert_state_signal(X_hat_pred_raw, 4); % N x 4
flag_magn = convert_scalar_signal(flag_magn_raw);     % N x 1

%% Riallineamento temporale, se necessario

t_y_magn = out.y_magn_raw.Time;
t_x_pred = out.X_hat_pred.Time;
t_flag_magn = out.flag_magn.Time;

% Uso come riferimento il tempo di y_magn
t_ref = t_y_magn;

% Allineo X_hat_pred ai tempi di y_magn
X_hat_pred_aligned = align_matrix_signal_to_time(X_hat_pred, t_x_pred, t_ref);

% Allineo flag_magn ai tempi di y_magn
flag_magn_aligned = align_scalar_signal_to_time(flag_magn, t_flag_magn, t_ref);

%% Costruzione misura angolare equivalente

mx = y_magn(:,1);
my = y_magn(:,2);

beta_magn = atan2(-my, mx);

beta_pred = X_hat_pred_aligned(:,3);

e_beta_magn = wrapToPi_local(beta_magn - beta_pred);

%% Uso solo campioni magnetometrici validi/accettati

idx_valid = flag_magn_aligned > 0.5;

e_beta_magn_valid = e_beta_magn(idx_valid);

%% Rimozione transitorio iniziale

perc_transitorio = 0.10;

N_valid = length(e_beta_magn_valid);
idx0 = floor(perc_transitorio*N_valid) + 1;

e_beta_magn_valid = e_beta_magn_valid(idx0:end);

%% Test di bianchezza Anderson

res_beta_magn = anderson_whiteness_test(e_beta_magn_valid, maxLag, alpha_test, ...
    'Innovazione magnetometrica equivalente angolare');

function result = anderson_whiteness_test(epsilon, maxLag, alpha, signal_name)

    epsilon = epsilon(:);
    epsilon = epsilon(~isnan(epsilon));

    % Rimozione media campionaria
    epsilon = epsilon - mean(epsilon);

    N = length(epsilon);

    if N <= maxLag + 1
        error('Numero di campioni insufficiente per il test su %s.', signal_name);
    end

    gamma0 = (1/N) * sum(epsilon.^2);

    if gamma0 <= eps
        error('Varianza nulla o troppo piccola nel segnale %s.', signal_name);
    end

    rho = zeros(maxLag,1);

    for tau = 1:maxLag
        gamma_tau = (1/N) * sum(epsilon(1:N-tau) .* epsilon(1+tau:N));
        rho(tau) = gamma_tau / gamma0;
    end

    z = sqrt(2) * erfinv(1 - alpha);
    beta = z / sqrt(N);

    out_of_bounds = abs(rho) > beta;

    num_out = sum(out_of_bounds);
    perc_out = num_out / maxLag;

    test_passed = perc_out < alpha;

    figure;
    stem(1:maxLag, rho, 'filled');
    hold on;
    yline(beta, '--', 'Limite superiore');
    yline(-beta, '--', 'Limite inferiore');
    grid on;
    xlabel('\tau');
    ylabel('\hat{\rho}(\tau)');
    title(['Test di bianchezza Anderson - ', signal_name]);

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

    if ismatrix(Xraw) && sz(2) == nv
        X = Xraw;
        return;
    end

    if ismatrix(Xraw) && sz(1) == nv
        X = Xraw';
        return;
    end

    if numel(sz) == 3 && sz(1) == 1 && sz(2) == nv
        X = squeeze(Xraw);
        if size(X,1) == nv
            X = X';
        end
        return;
    end

    if numel(sz) == 3 && sz(1) == nv && sz(2) == 1
        X = squeeze(Xraw)';
        return;
    end

    if numel(sz) == 3 && sz(2) == 1 && sz(3) == nv
        X = squeeze(Xraw);
        return;
    end

    error('Formato vettoriale non riconosciuto. Size: [%s]', num2str(sz));

end


function X = convert_state_signal(Xraw, nx)

    sz = size(Xraw);

    if ismatrix(Xraw) && sz(2) == nx
        X = Xraw;
        return;
    end

    if ismatrix(Xraw) && sz(1) == nx
        X = Xraw';
        return;
    end

    if numel(sz) == 3 && sz(1) == 1 && sz(2) == nx
        X = squeeze(Xraw);
        if size(X,1) == nx
            X = X';
        end
        return;
    end

    if numel(sz) == 3 && sz(1) == nx && sz(2) == 1
        X = squeeze(Xraw)';
        return;
    end

    if numel(sz) == 3 && sz(2) == 1 && sz(3) == nx
        X = squeeze(Xraw);
        return;
    end

    error('Formato stato non riconosciuto. Size: [%s]', num2str(sz));

end


function y_aligned = align_scalar_signal_to_time(y, t_y, t_ref)

    y = y(:);
    t_y = t_y(:);
    t_ref = t_ref(:);

    y_aligned = zeros(length(t_ref),1);

    for k = 1:length(t_ref)
        [~, idx] = min(abs(t_y - t_ref(k)));
        y_aligned(k) = y(idx);
    end

end


function X_aligned = align_matrix_signal_to_time(X, t_X, t_ref)

    t_X = t_X(:);
    t_ref = t_ref(:);

    Nref = length(t_ref);
    nx = size(X,2);

    X_aligned = zeros(Nref,nx);

    for k = 1:Nref
        [~, idx] = min(abs(t_X - t_ref(k)));
        X_aligned(k,:) = X(idx,:);
    end

end


function a = wrapToPi_local(a)

    a = mod(a + pi, 2*pi) - pi;

end