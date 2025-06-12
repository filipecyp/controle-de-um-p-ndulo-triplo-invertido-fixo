%Modelagem; linearização; função de transferência; controlabilidade; observabilidade
% clc;
clear;
close all;
format shortG

% 1. Declarar variáveis simbólicas
syms t theta1(t) theta2(t) theta3(t) 
syms L1 L2 L3 m g T1 T2 m1 m2 m3 W


G1x=L1*sin(theta1(t))/2;
G1y=L1*cos(theta1(t))/2;
G2x=(L1*sin(theta1(t))+L2*sin(theta2(t))/2); 
G2y=(L1*cos(theta1(t))+L2*cos(theta2(t))/2);
G3x=L1*sin(theta1(t))+L2*sin(theta2(t))+L3*sin(theta3(t))/2; 
G3y=(L1*cos(theta1(t))+L2*cos(theta2(t))+L3*cos(theta3(t))/2);


vg1x=diff(G1x,t);
vg1y=diff(G1y,t);

vg2x=diff(G2x,t);
vg2y=diff(G2y,t);

vg3x=diff(G3x,t);
vg3y=diff(G3y,t);

vg1=sqrt(vg1x^2+vg1y^2);
vg2=sqrt(vg2x^2+vg2y^2);
vg3=sqrt(vg3x^2+vg3y^2);

dtheta1 = diff(theta1(t),t);
dtheta2 = diff(theta2(t),t);
dtheta3 = diff(theta3(t),t);
Ecin1=m1*vg1^2/2+(m1*L1^2*(dtheta1)^2)/6;
Ecin2=m2*vg2^2/2+(m2*L2^2*(dtheta2)^2)/6;
Ecin3=m3*vg3^2/2+m3*L3^2*(dtheta3)^2/6;

V1=m1*g*G1y;
V2=m2*g*G2y;
V3=m3*g*G3y;

Ecin=Ecin1+Ecin2+Ecin3;
V=V1+V2+V3;
L=Ecin-V;

syms c
R1 = (1/2)*c*(dtheta1)^2;
R2 = (1/2)*c*(dtheta1 - dtheta2)^2;
R3 = (1/2)*c*(dtheta2 - dtheta3)^2;
R=R1+R2+R3;


dL_dt1=diff(L,theta1(t));

dL_Dt1=diff(L,dtheta1);
dL_Dt1_dt=diff(dL_Dt1,t);

eq1= -T1+T2+dL_Dt1_dt-dL_dt1 + diff(R, dtheta1)==0;


dL_dt2=diff(L,theta2(t));
dL_Dt2=diff(L,dtheta2);
dL_Dt2_dt=diff(dL_Dt2,t);
eq2=T1-T2+dL_Dt2_dt-dL_dt2 + diff(R, dtheta2)==0;


dL_dt3=diff(L,theta3(t));
dL_Dt3=diff(L,dtheta3);
dL_Dt3_dt=diff(dL_Dt3,t);
eq3=dL_Dt3_dt-dL_dt3-W == 0;

syms d2theta1 d2theta2 d2theta3

%dfthetas=[diff(theta1(t), t, t),diff(theta2(t), t, t),diff(theta3(t), t, t),diff(theta1(t), t),diff(theta2(t), t),diff(theta3(t), t)];
%dthetas=[]

eq1_subs = subs(eq1, diff(dtheta1, t), d2theta1);
eq1_subs = subs(eq1_subs, diff(dtheta2, t), d2theta2);
eq1_subs = subs(eq1_subs, diff(dtheta3, t), d2theta3);

eq2_subs = subs(eq2, diff(dtheta1, t), d2theta1);
eq2_subs = subs(eq2_subs, diff(dtheta2, t), d2theta2);
eq2_subs = subs(eq2_subs, diff(dtheta3, t), d2theta3);

eq3_subs = subs(eq3, diff(dtheta1, t), d2theta1);
eq3_subs = subs(eq3_subs, diff(dtheta2, t), d2theta2);
eq3_subs = subs(eq3_subs, diff(dtheta3, t), d2theta3);


[L, U] = equationsToMatrix([eq1_subs, eq2_subs, eq3_subs], [d2theta1, d2theta2, d2theta3]);
X=linsolve(L,U);


f1=dtheta1;
f2=X(1);
f3=dtheta2;
f4=X(2);
f5=dtheta3;
f6=X(3);



a21=diff(f2,theta1);
a22=diff(f2,theta2);
a23=diff(f2,theta3);
a24=diff(f2,dtheta1);
a25=diff(f2,dtheta2);
a26=diff(f2,dtheta3);

a41=diff(f4,theta1);
a42=diff(f4,theta2);
a43=diff(f4,theta3);
a44=diff(f4,dtheta1);
a45=diff(f4,dtheta2);
a46=diff(f4,dtheta3);

a61=diff(f6,theta1);
a62=diff(f6,theta2);
a63=diff(f6,theta3);
a64=diff(f6,dtheta1);
a65=diff(f6,dtheta2);
a66=diff(f6,dtheta3);




valores = [L1, L2, L3, g ,m1, m2, m3, dtheta1, dtheta2, dtheta3, theta1(t), theta2(t), theta3(t), T1,T2, c, W];
subs_valores = [0.5, 0.3, 0.2, 9.81, 1, 0.8, 0.6, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0];

a21_val = double(subs(a21, valores, subs_valores));
a22_val = double(subs(a22, valores, subs_valores));
a23_val = double(subs(a23, valores, subs_valores));
a24_val = double(subs(a24, valores, subs_valores));
a25_val = double(subs(a25, valores, subs_valores));
a26_val = double(subs(a26, valores, subs_valores));

a41_val = double(subs(a41, valores, subs_valores));
a42_val = double(subs(a42, valores, subs_valores));
a43_val = double(subs(a43, valores, subs_valores));
a44_val = double(subs(a44, valores, subs_valores));
a45_val = double(subs(a45, valores, subs_valores));
a46_val = double(subs(a46, valores, subs_valores));

a61_val = double(subs(a61, valores, subs_valores));
a62_val = double(subs(a62, valores, subs_valores));
a63_val = double(subs(a63, valores, subs_valores));
a64_val = double(subs(a64, valores, subs_valores));
a65_val = double(subs(a65, valores, subs_valores));
a66_val = double(subs(a66, valores, subs_valores));



A=[0 0 0 1 0 0
    0 0 0 0 1 0
    0 0 0 0 0 1 
    a21_val a22_val a23_val a24_val a25_val a26_val
    a41_val a42_val a43_val a44_val a45_val a46_val 
    a61_val a62_val a63_val a64_val a65_val a66_val];


disp(A);

b21=diff(f2,T1);
b41=diff(f4,T1);
b61=diff(f6,T1);

b22=diff(f2,T2);
b42=diff(f4,T2);
b62=diff(f6,T2);


b21_val=double(subs(b21,valores,subs_valores));
b22_val=double(subs(b22,valores,subs_valores));
b41_val=double(subs(b41,valores,subs_valores));
b42_val=double(subs(b42,valores,subs_valores));
b61_val=double(subs(b61,valores,subs_valores));
b62_val=double(subs(b62,valores,subs_valores));

B=[0 0
   0 0
   0 0
   b21_val b22_val
   b41_val b42_val
   b61_val b62_val
   ];

disp(B);


C = [1 0 0 0 0 0;
     0 1 0 0 0 0;
     0 0 1 0 0 0];


D = zeros(3, 2);  

e1= diff(f2,W);
e2= diff(f4,W);
e3= diff(f6,W);

e1v= double(subs(e1,valores,subs_valores));
e2v= double(subs(e2,valores,subs_valores));
e3v= double(subs(e3,valores,subs_valores));

E= [0;
    0;
    0;
    e1v;
    e2v;
    e3v];

disp(E);

syms s

I=eye(size(A));
G=C*(inv(I*s-A))*B;
%disp(simplify(G))

% Polinômio característico

polyA = poly(A);
disp('Polinômio característico:');
disp(polyA);

% Polos
poles = roots(polyA);
disp('Polos do sistema:');
disp(poles);

%controlabilidade

Co = ctrb(A, B);


% Verifica o posto (rank)
if rank(Co) == size(A,1)
    disp('O sistema é controlável.');
else
    disp('O sistema NÃO é controlável.');
end

%observabilidade

O = obsv(A, C); % matriz de observabilidade

if rank(O) == size(A, 1)
    disp('O sistema é observável.');
else
    disp('O sistema NÃO é observável.');
end
%%
%3
%Parte2 

sma=ss(A,B,C,D);

% 1 ALOCAÇÃO DE POLOS =====================================

% Polos para alocação
paloc=[-3, -3.5, -4, -4.5, -5, -5.5];

% Cálculo da matriz de ganho por alocação de polos
Kaloc=place(A,B,paloc);
disp(Kaloc);

% 2 - LQR ========================================

% Matrizes de ponderação
Q = 10^2*eye(6); R=10^0*eye(2);

disp(Q);
disp(R)

% Calculo da matriz de ganho por LQR
[Klqr,sre, polqr]=lqr(sma,Q,R);
disp(Klqr)

% 4 - SISTEMA CONTROLADO POR ALOCAÇÃO
Aaloc=A-B*Kaloc;
smfaloc=ss(Aaloc,B,C,D);

% 5 - SISTEMA CONTROLADO POR LQR
Alqr = A-B*Klqr;
smfklqr = ss(Alqr,B,C,D);
disp(eig(Alqr));

% 6 - SIMULAÇÕES DOS SISTEMAS COM ALOC. E LQR

% vetor de tempo:
ts=0:0.01:10;
% Condição inicial: 
x0=[0.1; -0.1; 0.1; 0; 0; 0];
% Função de entrada (degrau em t=0)
u1 = zeros(length(ts), 1);
u2= zeros(length(ts), 1);
ut = [u1(:), u2(:)];

% SIMULAÇÃO COM POLOS ALOCADOS
[yaloc,t,xaloc] = lsim(smfaloc,ut,ts,x0);
% SIMULACAO COM LQR
[ylqr,t,xlqr] = lsim(smfklqr,ut,ts,x0);

% GRÁFICOS
% DESLOCAMENTO ANGULAR DO PÊNDULO
figure(1);
plot(t,xaloc(:,1),'Linewidth',3); hold;
%comet(t,xaloc(:,3)); hold;
plot(t,xlqr(:,1),'Linewidth',3);
%comet(t,xlqr(:,3));
%plot(t,xaloc(:,3),t,xlqr(:,3));
legend('aloc','lqr'); title('Deslocamento angular do pêndulo (\theta1)');
xlabel('tempo (s)',FontSize=14);
ylabel('\theta1 (rad)',FontSize=14); grid;
% 
hold off;

%DESLOCAMENTO ANGULAR DA HASTE (ENTRADA DE CONTROLE)
figure(2);
plot(t,xaloc(:,2),'Linewidth',3); hold;
plot(t,xlqr(:,2),'Linewidth',3);
%plot(t,xaloc(:,3),t,xlqr(:,3));
legend('aloc','lqr'); title('Deslocamento angular da haste (\theta2)');
xlabel('tempo (s)',FontSize=14);
ylabel('\theta2 (rad)',FontSize=14); grid;
% 
hold off;

figure(3);
plot(t,xaloc(:,3),'Linewidth',3); hold;
plot(t,xlqr(:,3),'Linewidth',3);
%plot(t,xaloc(:,3),t,xlqr(:,3));
legend('aloc','lqr'); title('Deslocamento angular da haste (\theta3)');
xlabel('tempo (s)',FontSize=14);
ylabel('\theta3 (rad)',FontSize=14); grid;
% 
hold off;

%%
%observador

Polos_control_AmenosBKaloc = eig(A-B*Kaloc);

Polos_obs = [-10, -11, -12, -13, -14, -15];
L_obs = place(A', C', Polos_obs)';
disp('Matriz de ganho do Observador (L_obs):'); disp(L_obs);

polos_obs_AmenosLC = eig(A-L_obs*C);
disp('Polos efetivos do observador (autovalores de A-L_obs*C):');
disp(sort(polos_obs_AmenosLC));

% Tabela comparando polos
fprintf('\nComparação de Polos:\n');
fprintf('-------------------------------------------------\n');
fprintf('Polos Controlador (A-BK_aloc) | Polos Observador (A-LC)\n');
fprintf('-------------------------------------------------\n');
for k_pole = 1:size(A,1)
    fprintf('%10.4f + %10.4fi   |  %10.4f + %10.4fi\n', ...
            real(Polos_control_AmenosBKaloc(k_pole)), imag(Polos_control_AmenosBKaloc(k_pole)), ...
            real(polos_obs_AmenosLC(k_pole)), imag(polos_obs_AmenosLC(k_pole)));
end
fprintf('-------------------------------------------------\n');

% Gráfico comparando polos (opcional, pode ser complexo de visualizar bem)
figure;
plot(real(Polos_control_AmenosBKaloc), imag(Polos_control_AmenosBKaloc), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
plot(real(polos_obs_AmenosLC), imag(polos_obs_AmenosLC), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
grid on;
xlabel('Real Part');
ylabel('Imaginary Part');
title('Comparação: Polos do Controlador (x) vs. Polos do Observador (o)');
legend('Controlador (A-BK_{aloc})', 'Observador (A-LC)');
axis equal;


%4.2 Simulações do sistema com regulador e observador
% Regulador escolhido: Kaloc
K_ctrl = Kaloc; % Ou Klqr se preferir

% Sistema aumentado para simulação: [x; x_hat]
% dx/dt     = A*x - B*K_ctrl*x_hat + B*u_ext
% dx_hat/dt = L*C*x + (A - L*C - B*K_ctrl)*x_hat + B*u_ext

A_aug = [A,              -B*K_ctrl;
         L_obs*C,        A - L_obs*C - B*K_ctrl];

B_aug = [B;
         B];

% Para observar x e x_hat
C_aug = [eye(size(A,1)), zeros(size(A,1));
         zeros(size(A,1)), eye(size(A,1))];

D_aug = [zeros(size(A,1),size(B,2));
         zeros(size(A,1),size(B,2))];

sys_aug = ss(A_aug, B_aug, C_aug, D_aug);

% Condições iniciais
x0_plant = [0.1; -0.1; 0.1; 0; 0; 0]; % Condição inicial da planta
x0_hat   = zeros(size(A,1), 1);      % Condição inicial do observador (diferente para ver convergência)
z0_aug   = [x0_plant; x0_hat];

% Vetor de tempo e entrada externa u_ext (usando ut do código original)
ts_sim = 0:0.01:3; % Tempo de simulação reduzido para melhor visualização da convergência

% Entrada externa ut (igual à usada nas simulações de controle)
% Se for um regulador puro (sem u_ext), ut_sim should be zeros.
% Using the 'ut' from user's previous simulations.
u1_sim = zeros(length(ts_sim), 1);
u2_sim= zeros(length(ts_sim), 1);
ut_sim = [u1_sim, u2_sim];

disp('Simulando sistema com controlador e observador...');
[Y_aug, t_sim, X_aug_full] = lsim(sys_aug, ut_sim, ts_sim, z0_aug);

% Extrair estados reais e estimados
x_plant_sim = Y_aug(:, 1:size(A,1));
x_hat_sim   = Y_aug(:, size(A,1)+1:end);

% Calcular erro de estimação
error_est = x_plant_sim - x_hat_sim;

% Plotar resultados
state_names = {'\alpha (rad)', '\beta (rad)', 'x_0 (m)', 'd\alpha/dt (rad/s)', 'd\beta/dt (rad/s)', 'dx_0/dt (m/s)'};

% Plot dos três primeiros estados (posição) e seus estimados
for i = 1:3
    figure;
    plot(t_sim, x_plant_sim(:,i), 'b-', 'LineWidth', 2);
    hold on;
    plot(t_sim, x_hat_sim(:,i), 'r--', 'LineWidth', 2);
    xlabel('Tempo (s)');
    ylabel(state_names{i});
    title(['Estado Real vs. Estimado: ', state_names{i}]);
    legend('Real x_t', 'Estimado x_{hat}');
    grid on;
    hold off;
end

% Plot do erro de estimação para os três primeiros estados
for i = 1:3
    figure;
    plot(t_sim, error_est(:,i), 'k-', 'LineWidth', 2);
    xlabel('Tempo (s)');
    ylabel(['Erro e_',num2str(i),' = x_',num2str(i),' - x_{hat',num2str(i),'}']);
    title(['Erro de Estimação para ', state_names{i}]);
    grid on;
end

% Análise à luz do Princípio da Separação
poles_A_aug = eig(A_aug);
disp('Polos do sistema aumentado (controlador + observador):');
disp(sort(poles_A_aug));
disp('Polos do controlador (A-BK_ctrl):');
disp(sort(Polos_control_AmenosBKaloc));
disp('Polos do observador (A-L_obs*C):');
disp(sort(polos_obs_AmenosLC));
fprintf('\nO Princípio da Separação afirma que os polos do sistema combinado\n');
fprintf('são a união dos polos do controlador e dos polos do observador.\n');
fprintf('Verifique se os polos combinados são aproximadamente a união dos polos individuais.\n');

% Simulação do sistema original com Kaloc (SEM OBSERVADOR) para comparação de desempenho
% usando x0_plant como condição inicial
Aaloc_original = A - B*Kaloc;
smfaloc_original = ss(Aaloc_original, B, C, D); % Note: B is for u_ext here.
[yaloc_orig, t_orig_sim, xaloc_orig] = lsim(smfaloc_original, ut_sim, ts_sim, x0_plant);

figure;
subplot(3,1,1);
plot(t_sim, x_plant_sim(:,1), 'b-', t_orig_sim, yaloc_orig(:,1), 'r-');
title('Comparação \theta1: com Observador (azul) vs. Estado Real (vermelho)');
ylabel('\theta1 (rad)'); legend('Com Observador', 'Feedback Estado Real'); grid on;
subplot(3,1,2);
plot(t_sim, x_plant_sim(:,2), 'b-', t_orig_sim, yaloc_orig(:,2), 'r-');
title('\theta2');
ylabel('\theta2 (rad)'); grid on;
subplot(3,1,3);
plot(t_sim, x_plant_sim(:,3), 'b-', t_orig_sim, yaloc_orig(:,3), 'r-');
title('theta3');
ylabel('theta3'); xlabel('Tempo (s)'); grid on;
sgtitle('Desempenho do Sistema: Regulador Kaloc com Observador vs. Feedback de Estado Real');
