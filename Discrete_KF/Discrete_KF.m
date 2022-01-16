clear; clc; close all;

%% Параметры системы
b22 = -0.3;
b23 = 5;
i1 = 0.8;
i2 = 1;

%% Параметры моделирования и вектор времени t
h = 0.01;
t_end = 10;
n = t_end / h + 1;
t = 0 : h : t_end;

%% СКО шумов
sig_syst = 0.001;
sig_meas = 0.01;

%% Матрицы
A = [ 0  1;
      0 b22 ];
      
B = [ 0;
     b23 ];

G = B;

H = [0 1];

Q = sig_syst^2;
R = sig_meas^2;

F = eye(2) + A*h;
psi = B*h;
ge = G*h;

%% Н.У.
% начальное значение корреляционной матрицы ошибок
P = [ sig_syst^2 0 ;
       0  sig_syst^2 ];

% начальный вектор состояния
z_real(:, 1) = [ 0.2; 0 ];

% начальное измерение угловой скорости крена
hi(1) = z_real(2, 1) + normrnd(0, sig_meas);

% начальная оценка вектора состояния
z_est(:, 1) = z_real(:, 1);

%% Решение
for i = 2 : 1 : n
         
    % модель реального вектора состояния z
    u = -i1*z_real(1, i-1) - i2*z_real(2, i-1);
    w = normrnd(0, sig_syst);
    z_real(:, i) = F * z_real(:, i-1)  +  psi * u  +  ge * w;
    
    % Нахождение экстраполированной оценки вектора состояния
    u_est = -i1*z_est(1, i-1) - i2*z_est(2, i-1);
    z_extrap = F * z_est(:, i-1)  +  psi * u_est;

    % Моделирование текущих измерений
    hi(i) = z_real(2, i) + normrnd(0, sig_meas);

    %% Нахождение опт. оценки вектора состояния и ковариационной матрицы P
    P_extrap = F * P * F' + ge * Q * ge';

    v = hi(i) - H * z_extrap;

    S = H * P_extrap * H' + R;

    K = P_extrap * H' * S^(-1);

    % Опт. оценка
    z_est(:, i) = z_extrap + K * v;

    % Ковариационная матрица
    P = P_extrap - K * H * P_extrap;
end

%% Графики
figure; % крен
title('\gamma(t)')
hold, grid on
plot( t, z_real(1, :), 'm', 'DisplayName', 'реальность'  );
plot( t, z_est(1, :), 'b', 'markersize', 3, 'DisplayName', 'оценка' );
legend show

figure; % угл. скорость крена
title('\omega(t)')
hold, grid on
plot( t, hi, 'c', 'DisplayName', 'измерение'  );
plot( t, z_real(2, :), 'm', 'DisplayName', 'реальность'  );
plot( t, z_est(2, :), 'b', 'markersize', 3, 'DisplayName', 'оценка' );
legend show

