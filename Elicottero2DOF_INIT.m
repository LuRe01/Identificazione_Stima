%% FILE DI INIZIALIZZAZIONE PER STIMA ELICOTTERO 2DOF
% Autori: Edoardo Di Scalzo, Luigi Realacci, Francesco Bonciani
clear 
close all
clc

%-------------------------------------------------------
% Parametri del modello
%-------------------------------------------------------

% Parametri Inerziali e di massa
params.m = 0.2;             % [kg] Massa del sistema 
params.J_alpha = 0.012;     % [kg*m^2] Mom. di Inerzia rispetto all'asse di pitch 
params.J_y = 2.3e-4;        % [kg*m^2] Componente di Inerzia Trasversale 
params.J_z = 36.4e-4;       % [kg*m^2] Componente di Inerzia Verticale 
params.I_b = 2.3e-4;        % [kg*m^2] Componente di Inerzia della base 


% Parametri geometrici e meccanici
params.g = 9.81;            % [m/s^2] Accelerazione di gravità 
params.l = 0.2;             % [m] Distanza tra perno di rotazione e applicazione forze rotori 
params.c_alpha = 0.01;      % [N*m*s/rad] Coefficiente di attrito viscoso asse di pitch 
params.c_beta = 0.01;       % [N*m*s/rad] Coefficiente di attrito viscoso asse di yaw

% Parametri Aerodinamici
params.epsilon_p = 0.1;     % Coefficiente di cross-thrust del rotore di coda sull'asse del pitch
params.epsilon_y = 0.1;     % Coefficiente di cross-thrust del rotore principale sull'asse dello yaw

% Sample time
params.st = 0.00435;        % Pari al sample time del sensore piu veloce (accelerometro)

Tf = 20;                    % Tempo di simulazione
t = (0:params.st:Tf)';

%-------------------------------------------------------
% Condizioni Iniziali
%-------------------------------------------------------

alpha_0 = deg2rad(10);            % [rad] Angolo di pitch iniziale
alpha_dot_0 = 0;        % [rad/s] Derivata angolo di pitch iniziale
beta_0 = 0;             % [rad] Angolo di yaw iniziale
beta_dot_0 = 0;         % [rad/s] Derivata angolo di yaw iniziale


% Vettori condizioni iniziali (Vettori colonna!)
q_0 = [alpha_0 beta_0]';                 % Vettore condizioi iniziali per alpha e beta
q_dot_0 = [alpha_dot_0 beta_dot_0]';     % Vettore condizioni iniziali per alpha_dot e beta_dot

%-------------------------------------------------------
% Input del sistema (forze F1 ed F2)
%-------------------------------------------------------
F1 = 0.341 + 0.04 * sin(0.6 * t) + 0.02 * sin(1.4 * t);                % [N] Forza del rotore principale
F2 = -0.006 + 0.010 * sin(0.8 * t + 0.3) + 0.006 * sin(1.7 * t);       % [N] Forza del rotore di coda
u = [t, F1, F2];

%% SENSORI: Accelerometro e Magnetometro (Vectornav vn-100)
%-------------------------------------------------------
% Parametri Sensori
%-------------------------------------------------------
seme_rumore = 12345;        % Seme per rumore additivo dei sensori

% Accelerometro
sigma_acc = 0.0137;         % [m/s^2] Dev. Std. Accelerometro
st_acc = 0.00435;           % Sample time Accelerometro (230 Hz)

% Magnetometro
M_0 = 46.5;                 % [uT] Intensità campo magnetico locale (Es. Livorno)
sigma_magn = 0.14;          % [muT] Dev. Std. Magnetometro
st_magn = 0.005;            % Sample time Magnetometro (200 Hz)


%% PARAMETRI GESTIONE OUTLIER E MAHALANOBIS

% Attiva/Disattiva Outlyer (1 attivi, 0 disattivati)
Outlier_flag = 0;

% Flag per attivare/Disattivare Distanza di Mahalanobis
MAHALANOBIS = false;

% Treshold per distanza di mahalanobis
m_treshold = 3;

%% Parametri EKF

%-------------------------------------------------------
% Condizioni iniziali
%-------------------------------------------------------
P_0_alpha = deg2rad(5)^2;
P_0_alpha_dot = 0.20^2;
P_0_beta = deg2rad(8)^2;
P_0_beta_dot = 0.20^2;

x_0 = [(deg2rad(10) + deg2rad(3)) 0.05 deg2rad(5) -0.05]';           % Stato iniziale

% Matrice delle covarianze iniziale
P0 = diag([P_0_alpha P_0_alpha_dot P_0_beta P_0_beta_dot]);


% Covarianza rumore di processo
% params.sigma_alpha = 1e-5;                  % Dev. Std. alpha
% params.sigma_alpha_dot = 1e-4;              % Dev. Std. alpha_dot
% params.sigma_beta = 1e-5;                   % Dev. Std. beta
% params.sigma_beta_dot = 1e-4;               % Dev. Std. beta_dot
params.sigma_F = 1.8 * 0.05;                  % [N] Dev.Std. Ingresso

params.R_acc = sigma_acc^2;
params.R_magn = sigma_magn^2;

%% Parametri UKF

% SUKF (Filtro di Kalman Unscented Scalato)
params.Alpha_UKF = 1e-2;       % Parametro Alpha UKF (Dispersione sigma-points)
params.Beta_UKF = 2;        % Parametro Beta UKF (Ottimizza termini di ordine superiore)
params.Kappa_UKF = 0;       % Parametro Kappa UKF (Parametro di scaling secondario)

params.n_state = length(x_0);                   % Dimesnione stato
params.n_input = 2;                             % Dimensione input

% Caso con rumore di attuazione
params.n = params.n_state + params.n_input;     
params.N_sigma_points = 2 * params.n + 1;       % Numero di sigma points (2*6 + 1 = 13)

% Parametro Lambda UKF 
% Per fase di predizione
params.lambda = params.Alpha_UKF^2 * (params.n + params.Kappa_UKF) - params.n;
% Per fase di correzione
params.lambda_c = params.Alpha_UKF^2 * (params.n_state + params.Kappa_UKF) - params.n_state;

% Assegnamento pesi per media e covarianza
n = params.n;
params.Wm = zeros(2*n + 1, 1);          % Pesi media
params.Wc = zeros(2*n + 1, 1);          % Pesi covarianza
% Parte di predizione
for i =-n:n
    idx = i + n + 1;
    if i == 0  
        % Peso centrale
        params.Wm(idx) = params.lambda / (n  + params.lambda);
        params.Wc(idx) = (params.lambda / (n  + params.lambda)) + 1 - params.Alpha_UKF^2 + params.Beta_UKF;
    else
        % Pesi restanti
        params.Wm(idx) = 1 / (2 * (n + params.lambda));
        params.Wc(idx) = 1 / (2 * (n + params.lambda));
    end
end

% Parte di correzione (i sigma points sono 9 poichè si sfrutta rumore
% additivo)
n_state = params.n_state;
params.Wm_c = zeros(2*n_state + 1, 1);          % Pesi media
params.Wc_c = zeros(2*n_state + 1, 1);          % Pesi covarianza
for i =-n_state:n_state
    idx = i + n_state + 1;
    if i == 0  
        % Peso centrale
        params.Wm_c(idx) = params.lambda_c / (n_state  + params.lambda_c);
        params.Wc_c(idx) = (params.lambda_c / (n_state  + params.lambda_c)) + 1 - params.Alpha_UKF^2 + params.Beta_UKF;
    else
        % Pesi restanti
        params.Wm_c(idx) = 1 / (2 * (n_state + params.lambda_c));
        params.Wc_c(idx) = 1 / (2 * (n_state + params.lambda_c));
    end
end

% Chiamata al simulatore su Simulink
out = sim('Elicottero2DOF_Sim.slx');