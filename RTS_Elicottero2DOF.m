%% RTS SMOOTHER - Elicottero 2DOF
% Segnali richiesti dentro out:
% out.X_hat_pred
% out.X_hat
% out.P_pred
% out.P
% out.F

clc

%% Estrazione dati dalle timeseries

x_pred_raw = out.X_hat_pred.Data;      % x_{k|k-1}
x_corr_raw = out.X_hat.Data;           % x_{k|k}

P_pred_raw = out.P_pred.Data;          % P_{k|k-1}
P_corr_raw = out.P.Data;               % P_{k|k}
F_raw      = out.F.Data;               % F_k

t = out.X_hat.Time;

nx = 4;

%% Conversione stati in formato N x 4

x_pred = convert_state_timeseries(x_pred_raw, nx);
x_corr = convert_state_timeseries(x_corr_raw, nx);

N = size(x_corr, 1);

%% Conversione matrici in formato 4 x 4 x N

P_pred_3D = convert_to_3D_matrices(P_pred_raw, nx);
P_corr_3D = convert_to_3D_matrices(P_corr_raw, nx);
F_3D      = convert_to_3D_matrices(F_raw, nx);

%% Controlli dimensionali

fprintf('Numero campioni N = %d\n', N);
fprintf('size x_pred = [%s]\n', num2str(size(x_pred)));
fprintf('size x_corr = [%s]\n', num2str(size(x_corr)));
fprintf('size P_pred_3D = [%s]\n', num2str(size(P_pred_3D)));
fprintf('size P_corr_3D = [%s]\n', num2str(size(P_corr_3D)));
fprintf('size F_3D = [%s]\n', num2str(size(F_3D)));

if size(P_pred_3D,3) ~= N || size(P_corr_3D,3) ~= N || size(F_3D,3) ~= N
    error('Numero di campioni non coerente tra stati, covarianze e F.');
end

%% Inizializzazione RTS

x_hat_smoothed = zeros(N, nx);
P_smoothed = zeros(nx, nx, N);

% All'ultimo istante, stima smoothed = stima filtrata EKF
x_hat_smoothed(N, :) = x_corr(N, :);
P_smoothed(:, :, N) = P_corr_3D(:, :, N);

%% Esecuzione backward dello smoother RTS

for k = N-1:-1:1

    % Covarianza corretta al tempo k
    Pkk = P_corr_3D(:, :, k);

    % Grandezze predette al tempo k+1
    Pkp1_pred = P_pred_3D(:, :, k+1);
    Fkp1 = F_3D(:, :, k+1);

    % Guadagno RTS
    Ck = Pkk * Fkp1' / Pkp1_pred;

    % Stato smoothed
    x_hat_smoothed(k, :) = x_corr(k, :)' + ...
        Ck * (x_hat_smoothed(k+1, :)' - x_pred(k+1, :)');

    % Covarianza smoothed
    P_smoothed(:, :, k) = Pkk + ...
        Ck * (P_smoothed(:, :, k+1) - Pkp1_pred) * Ck';

    % Simmetrizzazione numerica
    P_smoothed(:, :, k) = 0.5 * (P_smoothed(:, :, k) + P_smoothed(:, :, k)');
end

%% Grafici confronto EKF - RTS

figure;
subplot(4,1,1)
plot(t, x_corr(:,1), 'LineWidth', 1.2); hold on;
plot(t, x_hat_smoothed(:,1), 'LineWidth', 1.2);
grid on;
ylabel('\alpha [rad]');
legend('EKF', 'RTS');
title('Confronto EKF - RTS');

subplot(4,1,2)
plot(t, x_corr(:,2), 'LineWidth', 1.2); hold on;
plot(t, x_hat_smoothed(:,2), 'LineWidth', 1.2);
grid on;
ylabel('\dot{\alpha} [rad/s]');
legend('EKF', 'RTS');

subplot(4,1,3)
plot(t, x_corr(:,3), 'LineWidth', 1.2); hold on;
plot(t, x_hat_smoothed(:,3), 'LineWidth', 1.2);
grid on;
ylabel('\beta [rad]');
legend('EKF', 'RTS');

subplot(4,1,4)
plot(t, x_corr(:,4), 'LineWidth', 1.2); hold on;
plot(t, x_hat_smoothed(:,4), 'LineWidth', 1.2);
grid on;
ylabel('\dot{\beta} [rad/s]');
xlabel('Tempo [s]');
legend('EKF', 'RTS');

%% Grafici diagonali delle covarianze

P_corr_diag = zeros(N, nx);
P_smooth_diag = zeros(N, nx);

for k = 1:N
    P_corr_diag(k,:) = diag(P_corr_3D(:,:,k))';
    P_smooth_diag(k,:) = diag(P_smoothed(:,:,k))';
end

figure;
subplot(4,1,1)
plot(t, P_corr_diag(:,1), 'LineWidth', 1.2); hold on;
plot(t, P_smooth_diag(:,1), 'LineWidth', 1.2);
grid on;
ylabel('var(\alpha)');
legend('EKF', 'RTS');
title('Confronto covarianze EKF - RTS');

subplot(4,1,2)
plot(t, P_corr_diag(:,2), 'LineWidth', 1.2); hold on;
plot(t, P_smooth_diag(:,2), 'LineWidth', 1.2);
grid on;
ylabel('var(\dot{\alpha})');
legend('EKF', 'RTS');

subplot(4,1,3)
plot(t, P_corr_diag(:,3), 'LineWidth', 1.2); hold on;
plot(t, P_smooth_diag(:,3), 'LineWidth', 1.2);
grid on;
ylabel('var(\beta)');
legend('EKF', 'RTS');

subplot(4,1,4)
plot(t, P_corr_diag(:,4), 'LineWidth', 1.2); hold on;
plot(t, P_smooth_diag(:,4), 'LineWidth', 1.2);
grid on;
ylabel('var(\dot{\beta})');
xlabel('Tempo [s]');
legend('EKF', 'RTS');

%% Salvataggio risultati

RTS.x_hat_smoothed = x_hat_smoothed;
RTS.P_smoothed = P_smoothed;
RTS.time = t;

disp('RTS completato correttamente.');

%% Funzioni locali

function X = convert_state_timeseries(Xraw, nx)

    sz = size(Xraw);

    % Caso desiderato: N x nx
    if ismatrix(Xraw) && sz(2) == nx
        X = Xraw;
        return;
    end

    % Caso nx x N
    if ismatrix(Xraw) && sz(1) == nx
        X = Xraw';
        return;
    end

    % Caso 1 x nx x N
    if numel(sz) == 3 && sz(1) == 1 && sz(2) == nx
        X = squeeze(Xraw);      % nx x N oppure N x nx a seconda di squeeze
        if size(X,1) == nx
            X = X';
        end
        return;
    end

    % Caso nx x 1 x N
    if numel(sz) == 3 && sz(1) == nx && sz(2) == 1
        X = squeeze(Xraw);      % nx x N
        X = X';
        return;
    end

    % Caso N x 1 x nx
    if numel(sz) == 3 && sz(2) == 1 && sz(3) == nx
        X = squeeze(Xraw);      % N x nx
        return;
    end

    error('Formato stato non riconosciuto. Size trovata: [%s]', num2str(sz));

end


function P3 = convert_to_3D_matrices(Praw, nx)

    sz = size(Praw);

    if numel(sz) ~= 3
        error('Il segnale matrice non ha formato 3D. Size trovata: [%s]', num2str(sz));
    end

    if sz(1) == nx && sz(2) == nx
        % Formato nx x nx x N
        P3 = Praw;

    elseif sz(2) == nx && sz(3) == nx
        % Formato N x nx x nx
        N = sz(1);
        P3 = zeros(nx, nx, N);

        for k = 1:N
            P3(:,:,k) = squeeze(Praw(k,:,:));
        end

    elseif sz(1) == 1 && sz(2) == nx && sz(3) == nx
        % Caso singola matrice 1 x nx x nx, poco probabile
        P3 = squeeze(Praw);
        P3 = reshape(P3, nx, nx, 1);

    else
        error('Formato matrice 3D non riconosciuto. Size trovata: [%s]', num2str(sz));
    end

end