clc;
syms theta1 theta2 theta3 d a t 
%% Configurazione iniziale manipolatore
theta1_deg = 90;
theta2_deg = 0;
theta3_deg = 0;
theta1_rad = deg2rad(theta1_deg);
theta2_rad = deg2rad(theta2_deg);
theta3_rad = deg2rad(theta3_deg);
d=0.5;
a=1;
%%------Parametri dell'algo di inversione------%%
ti = 0;
tf = 5;
deltaT = 0.002;
K = [500,0;0,500];
k0 = 5;
I = eye(4,4);
%%---------------------------------------------%%
q=[theta1_rad,theta2_rad,theta3_rad,d] %vettore spazio dei giunti
%task = 1 indica che il manipolatore deve eseguire solo il compito primario
%task = 2 indica che il manipolatore deve eseguire anche il compito secondario
task = 2; 
%% Cinematica diretta
[A10,A20,A30,A40] = cinematicaDiretta(a,q); 
%% Spazio operativo x = k(q) con solo la posizione
x_pos = [A40(1,4),A40(2,4)]
%% JACOBIANO geometrico = analitico
J = jacobianoGeometrico(q,A10,A20,A30,A40);
J = J([1:2],:) %Elimino le righe nulle J e non considero l'orientamento = Jp
k = rank(J) %calcolo rango 
Jinv = pinv(J)
%%Variabili per i grafici
pos_q1=[];
pos_q2=[];
pos_q3=[];
pos_q4=[];
vel_q1=[];
vel_q2=[];
vel_q3=[];
vel_q4=[];
err_pos=[];
manip=[];
PosEEx=[];
PosEEy=[];
Pos_des_EEx=[];
Pos_des_EEy=[];
%%--------------------%%
%% Algoritmo di inversione
fprintf('*******Start inverse kinematics********\n')
i = 0;
%% Inizio ciclo
for t = ti : deltaT : tf
fprintf('|---- Iter : %d\n', i);
fprintf('|---- Time : %f\n', t);
%%-------------------Specificazione della traiettoria--------------------%%
%% PEx e PEy approssimati
% PEx = 0.5 + (( -0.022 * (t^3) + 0.169 * (t^2)) / sqrt(2)); %Posizione lungo x dell'ee
% PEy = 2 + (( 0.022 * (t^3) - 0.169 * (t^2)) / sqrt(2)); %Posizione lungo y dell'ee
% PEx_der = -0.046 * (t^2) + 0.239 * t; %derivata della posizione lungo x dell'ee
% PEy_der = 0.046 * (t^2) - 0.239 * t; %derivata della posizione lungo y dell'ee
%% PEx e PEy reali
PEx = 0.5 + (( -2 * (t^3) + 15 * (t^2)) / 125); 
PEy = 2 - (( -2 * (t^3) + 15 * (t^2)) / 125); 
PEx_der = -(6/125) * (t^2) + (6/25) * t; 
PEy_der = (6/125) * (t^2) - (6/25) * t; 

Pd = [PEx,PEy]; %Vettore posizione desiderata
Pd_der = [PEx_der,PEy_der]; %Vettore derivata posizione desiderata
%%-----------------------------------------------------------------------%%
%%----------Calcolo errore--------%%
e = Pd-x_pos
%%--------------------------------%%
%%-----------------Definizione del compito secondario sigma--------------%%
sigma = (1/2)*(((sin(q(2)))^2)+((sin(q(3)))^2)) %misura di manipolabilità
grad_sigma =[0,sin(q(2))*cos(q(2)),sin(q(3))*cos(q(3)),0] %gradiente della misura di manipolabilità
q_der_0 = k0*(grad_sigma)'
%%-----------------------------------------------------------------------%%
%% Calcolo delle velocita' e posizioni di q
%%----------------Compito primario con fattore di correzione-------------%%
if(task == 1)
    q_der = Jinv*((Pd_der')+(K*e'))
end
%%-----------------------------------------------------------------------%%
%%-----------------------Con compito secondario--------------------------%%
if(task == 2)
    q_der = (Jinv*((Pd_der')+(K*e')))+((I-Jinv*J)*q_der_0)
end
%%-----------------------------------------------------------------------%%
%%--------------------------Calcolo posizioni di q-----------------------%%
q = q + (q_der * deltaT)'
%%-----------------------------------------------------------------------%%
[A10,A20,A30,A40] = cinematicaDiretta(a,q);
x_pos = [A40(1,4),A40(2,4)]
J = jacobianoGeometrico(q,A10,A20,A30,A40);
J = J([1:2],:)
Jinv = pinv(J)

pos_q1 = [pos_q1,q(1)];
pos_q2 = [pos_q2,q(2)];
pos_q3 = [pos_q3,q(3)];
pos_q4 = [pos_q4,q(4)];
vel_q1 = [vel_q1,q_der(1)];
vel_q2 = [vel_q2,q_der(2)];
vel_q3 = [vel_q3,q_der(3)];
vel_q4 = [vel_q4,q_der(4)];
err_pos = [err_pos,norm(e)];
manip = [manip,sigma];
PosEEx = [PosEEx , x_pos(1)];
PosEEy = [PosEEy , x_pos(2)];
Pos_des_EEx = [Pos_des_EEx , PEx];
Pos_des_EEy = [Pos_des_EEy , PEy];
i=i+1;
end

q_deg = [rad2deg(q(1)),rad2deg(q(2)),rad2deg(q(3)),q(4)] %posizione dei giunti in gradi

%% Stampa dei grafici
t=ti:deltaT:tf;

figure
subplot(2,3,1)
hold on; 
plot( t, pos_q1, 'r', 'LineWidth', 2);
plot( t, pos_q2, 'g', 'LineWidth', 2);
plot( t, pos_q3, 'b', 'LineWidth', 2);
title('Position of joints: q1,q2,q3')
legend('q1', 'q2','q3');
xlabel('[s]')
ylabel('[rad]')
hold off;

subplot(2,3,2)
plot(t,pos_q4,'c','LineWidth', 2);
title('Position of joint q4')
xlabel('[s]')
ylabel('[m]')

subplot(2,3,3)
plot(t,err_pos,'b', 'LineWidth', 2);
title('Norm of position error')
xlabel('[s]')
ylabel('[m]')

subplot(2,3,4)
hold on; 
plot( t, vel_q1,'r' ,'LineWidth', 2);
plot( t, vel_q2, 'g', 'LineWidth', 2);
plot( t, vel_q3, 'b', 'LineWidth', 2);
title('Velocity of joints: q1,q2,q3')
legend('q1', 'q2','q3');
xlabel('[s]')
ylabel('[rad/s]')
hold off;

subplot(2,3,5)
plot(t,vel_q4,'c','LineWidth', 2);
title('Velocity of joint q4')
xlabel('[s]')
ylabel('[m/s]')

subplot(2,3,6)
plot(t,manip,'r', 'LineWidth', 2);
title('Manip')
xlabel('[s]')
ylabel('[rad]')

% figure
% plot(PosEEx,PosEEy,'b','LineWidth', 3)
% title('Posizione dell organo terminale')
% xlabel('x')
% ylabel('y')
% hold on;
% plot(Pos_des_EEx,Pos_des_EEy,'r','LineWidth', 2,'LineStyle', ':')


