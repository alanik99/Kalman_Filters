clear; clc; close all;

global i1 i2 sig_syst sig_meas A B G H R Q h

%% ��������� �������
b22 = -0.55;
b23 = 3;
i1 = 0.8;
i2 = 1;

%% ��� �����
sig_syst = 0.02;
sig_meas = 0.5;

%% �������
A = [ 0  1;
      0 b22 ];
      
B = [ 0;
     b23 ];

G = B;

H = [0 1];

Q = sig_syst^2;
R = sig_meas^2;

%% ��������� ������������� � ������ ������� t
h = 0.1;
t_end = 20;
n = t_end / h + 1;
t = 0 : h : t_end;

%% �.�.
% ��������� �������� �������������� ������� ������
P = [ sig_syst^2 0 ;
       0  sig_syst^2 ];

% ��������� ������ ���������
z_real(:, 1) = [ 5; 0 ];

% ��������� ��������� ���� �����
hi(1) = z_real(2, 1) + normrnd(0, sig_meas);

% ��������� ������ ������� ���������: ��������� ��� ����� - ������ ���. ���������, ��� ���. �������� - ����
% z_est(:, 1) = [ 0; hi(1) ];
z_est(:, 1) = z_real(:, 1);

%% �������
for i = 2 : 1 : n
    
    % ������ ������ ������� ��������� z
     z_est(:, i) = rungeKuttaStep(  @dzEst_dt,  z_est(:, i-1),   P,  hi(i-1)  );
         
    % ������ ��������� ������� ��������� z
    z_real(:, i) = rungeKuttaStep( @dzReal_dt,  z_real(:, i-1),  0,   0 );
    
    % ������� �� �������
               P = rungeKuttaStep(   @dP_dt,      P,             0,   0 );
      
    % ������� ��������� ���� �����
    hi(i) = z_real(2, i) + normrnd(0, sig_meas);
end

%% �������
figure; % ����
title('\gamma(t)')
hold, grid on
plot( t, z_real(1, :), 'm', 'DisplayName', '����������'  );
plot( t, z_est(1, :), 'b', 'markersize', 3, 'DisplayName', '������' );
legend show

figure; % ���. �������� �����
title('\omega(t)')
hold, grid on
plot( t, hi, 'c', 'DisplayName', '���������'  );
plot( t, z_real(2, :), 'm', 'DisplayName', '����������'  );
plot( t, z_est(2, :), 'b', 'markersize', 3, 'DisplayName', '������' );
legend show


%% �������
function RETURN = dzEst_dt( z_est, P, hi )
    
    global R A H B i1 i2

    K = P * H' * R^(-1); 
    
    u_est = -i1*z_est(1) - i2*z_est(2);
    
    RETURN = A*z_est + B*u_est + K*( hi - H*z_est );
end

function RETURN = dP_dt( P, unused1, unused2 )

    global R Q A G H

    RETURN = A*P + P*A'  -  P * H' * R^(-1) * H * P  +  G * Q * G';
end

function RETURN = dzReal_dt( z, unused1, unused2 )
    
    global i1 i2 sig_syst A B G
    
    u = -i1*z(1) - i2*z(2);
    
    RETURN = A*z + B*u + G*normrnd(0, sig_syst); 
end

function RETURN = rungeKuttaStep(f, z_curr, param1, param2)
    % params may be unused or used in some of the functions
    global h

    k1 = f(z_curr, param1, param2);
    k2 = f(z_curr + h*k1/2, param1, param2); 
    k3 = f(z_curr + h*k2/2, param1, param2); 
    k4 = f(z_curr + h*k3, param1, param2);

    RETURN = z_curr + h*(k1 + 2*k2 + 2*k3 + k4) / 6;
end
