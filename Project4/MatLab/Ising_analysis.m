%% Ising plot maker
%% Ising 2x2
clear all
load 'Ising2x2.mat'

% compute model values
temp2 = 1.0;
model_chi_ass=(1./temp2).*(6*exp(8./temp2)+2*exp(-8./temp2)+6)./(19+cosh(16./temp2)+12*cosh(8./temp2));
model_cv=64*(3*cosh(8./temp2)+1)./((6+2*cosh(8./temp2)).^2.*temp2.^2);
model_energy=(16*exp(-8./temp2)-16*exp(8./temp2))./(12+2*exp(-8./temp2)+2*exp(8./temp2))/4;
model_magn=(8*exp(8./temp2)+16)./(12+2*exp(-8./temp2)+2*exp(8./temp2))/4;
model_chi=(4./temp2).*(exp(8./temp2)+1)./(6+2*cosh(8./temp2));

figure();
plot(cycles, e_m);
hold on
plot([0,1e6],[model_energy, model_energy])
xlabel('MC cycles','Interpreter','Latex','Fontsize',14);
ylabel('$ <E>$ $[J]$','Interpreter','Latex','Fontsize',14);
title('Mean energy for MC runs','Interpreter','Latex','Fontsize',18);
set(gca,'TickLabelInterpreter','latex');
grid on;

figure();
plot(cycles, e_var);
hold on
plot([0,1e6],[model_cv, model_cv])
xlabel('MC cycles','Interpreter','Latex','Fontsize',14);
ylabel('$c_V$ $[k_B]$','Interpreter','Latex','Fontsize',14);
title('Specific heat for MC runs','Interpreter','Latex','Fontsize',18);
set(gca,'TickLabelInterpreter','latex');
grid on;

figure();
plot(cycles, m_abs_m);
hold on
plot([0,1e6],[model_magn, model_magn])
xlabel('MC cycles','Interpreter','Latex','Fontsize',14);
ylabel('$|M|$','Interpreter','Latex','Fontsize',14);
title('Magnetization for MC runs','Interpreter','Latex','Fontsize',18);
set(gca,'TickLabelInterpreter','latex');
grid on;

figure();
plot(cycles, m_var);
hold on
plot([0,1e6],[model_chi, model_chi])
xlabel('MC cycles','Interpreter','Latex','Fontsize',14);
ylabel('$\chi$ $[J^{-1}]$','Interpreter','Latex','Fontsize',14);
title('Susceptibility for MC runs','Interpreter','Latex','Fontsize',18);
set(gca,'TickLabelInterpreter','latex');
ylim([3.6 4.1])
grid on;

%% Ising 20x20
clear all
close all 
clc
% Data for all spins up
load 'Ising20_up.mat'
% T=1.0
figure();
plot(cycles1,energy24,'color','r','LineWidth',1)
xlabel('MC cycles','Interpreter','Latex','Fontsize',14);
ylabel('$<E>$ $[J]$','Interpreter','Latex','Fontsize',14);
title('Energy for MC cycles - ordered lattice','Interpreter','Latex','Fontsize',18);
set(gca,'TickLabelInterpreter','latex');
xlim([0,1e5]);
%ylim([-2.5,-0.6])
grid on;

figure()
plot(cycles1,magn24,'color','m','LineWidth',1)
xlabel('MC cycles','Interpreter','Latex','Fontsize',14);
ylabel('$|M|$','Interpreter','Latex','Fontsize',14);
title('Magnetization for MC cycles - ordered lattice','Interpreter','Latex','Fontsize',18);
set(gca,'TickLabelInterpreter','latex');
xlim([0,1e5]);
%ylim([0,1.5]);
grid on;

%%
load 'Ising20_rand.mat'
load 'Ising20_up.mat'
figure();
loglog(cycles1, choices24)
hold on
loglog(cycles1, choices1)
xlabel('MC cycles','Interpreter','Latex','Fontsize',14);
ylabel('$N_{Moves}$','Interpreter','Latex','Fontsize',14);
title('Accepted moves for MC cycles - ordered lattice','Interpreter','Latex','Fontsize',18);
ll=legend(gca,'show','$T=2.4\,Jk_B^{-1}$','$T=1.0\,Jk_B^{-1}$');
set(ll,'Interpreter','Latex','Location','northwest')
set(gca,'TickLabelInterpreter','latex');
%xlim([0,1e5]);
%ylim([0,1.5]);
grid on;

figure();
loglog(cycles1, choices24_r)
hold on
loglog(cycles1, choices1_r)
xlabel('MC cycles','Interpreter','Latex','Fontsize',14);
ylabel('$N_{Moves}$','Interpreter','Latex','Fontsize',14);
title('Accepted moves for MC cycles - random lattice','Interpreter','Latex','Fontsize',18);
ll=legend(gca,'show','$T=2.4\,Jk_B^{-1}$','$T=1.0\,Jk_B^{-1}$');
set(ll,'Interpreter','Latex','Location','northwest')
set(gca,'TickLabelInterpreter','latex');
%xlim([0,1e5]);
%ylim([0,1.5]);
grid on;

%% Probabilities
close all
clear all
clc
load 'prob.mat' 
energies = [-800:4:800]/400;

figure();
bar(energies, counts1_r);
mean=(sum(energies'.*counts1_r))
grid on
xlim([-4.5/2 -3.5/2]);
xlabel('$E_i$','Interpreter','Latex','Fontsize',14);
ylabel('$p(E_i)$','Interpreter','Latex','Fontsize',14);
title('Energy per spin distribution - $T=1.0\,Jk_B^{-1}$','Interpreter','Latex','Fontsize',18);
set(gca,'TickLabelInterpreter','latex');

figure();
bar(energies, counts24_r)
mean=sum(energies'.*counts24_r)
grid on
xlim([-5/2 0])
xlabel('$E_i$','Interpreter','Latex','Fontsize',14);
ylabel('$p(E_i)$','Interpreter','Latex','Fontsize',14);
title('Energy per spin distribution - $T=2.4\,Jk_B^{-1}$','Interpreter','Latex','Fontsize',18);
set(gca,'TickLabelInterpreter','latex');

%% Phase transition
clear all
close all
clc

load 'Ising40.mat';
load 'Ising60.mat';
load 'Ising80.mat';
load 'Ising100.mat';

figure()
plot(T40,E40)
hold on
plot(T60,E60)
hold on
plot(T80,E80)
hold on
plot(T100,E100)
xlabel('T $[J\,k_B^{-1}]$','Interpreter','Latex','Fontsize',14);
ylabel('$<E>$ $[J]$','Interpreter','Latex','Fontsize',14);
title('Mean Energy - Temperature dependance','Interpreter','Latex','Fontsize',18);
ll=legend(gca,'show','$L=40$','$L=60$','$L=80$','$L=100$');
set(ll,'Interpreter','Latex','Location','northwest')
set(gca,'TickLabelInterpreter','latex');
grid on;

figure()
plot(T40,M40)
hold on
plot(T60,M60)
hold on
plot(T80,M80)
hold on
plot(T100,M100)
xlabel('T $[J\,k_B^{-1}]$','Interpreter','Latex','Fontsize',14);
ylabel('$<M>$','Interpreter','Latex','Fontsize',14);
title('Mean Magnetization - Temperature dependance','Interpreter','Latex','Fontsize',18);
ll=legend(gca,'show','$L=40$','$L=60$','$L=80$','$L=100$');
set(ll,'Interpreter','Latex','Location','northeast')
set(gca,'TickLabelInterpreter','latex');
grid on;

figure()
plot(T40,Cv40)
hold on
plot(T60,Cv60)
hold on
plot(T80,Cv80)
hold on
plot(T100,Cv100)
xlabel('T $[J\,k_B^{-1}]$','Interpreter','Latex','Fontsize',14);
ylabel('$c_V$ $[k_B]$','Interpreter','Latex','Fontsize',14);
title('Specific heat - Temperature dependance','Interpreter','Latex','Fontsize',18);
ll=legend(gca,'show','$L=40$','$L=60$','$L=80$','$L=100$');
set(ll,'Interpreter','Latex','Location','northeast')
set(gca,'TickLabelInterpreter','latex');
grid on;

figure()
plot(T40,chi_abs40)
hold on
plot(T60,chi_abs60)
hold on
plot(T80,chi_abs80)
hold on
plot(T100,chi_abs100)
xlabel('T $[J\,k_B^{-1}]$','Interpreter','Latex','Fontsize',14);
ylabel('$\chi$ $[J^{-1}]$','Interpreter','Latex','Fontsize',14);
title('Susceptibility - Temperature dependance','Interpreter','Latex','Fontsize',18);
ll=legend(gca,'show','$L=40$','$L=60$','$L=80$','$L=100$');
set(ll,'Interpreter','Latex','Location','northeast')
set(gca,'TickLabelInterpreter','latex');
grid on;
%% Find critical temp
[pks, loc] = findpeaks(chi_abs40);
Tc(1)=T40(loc);
[pks, loc] = findpeaks(chi_abs60);
Tc(2)=T60(loc);
[pks, loc] = findpeaks(chi_abs80);
Tc(3)=T80(loc(1));
[pks, loc] = findpeaks(Cv100);
Tc(4)=T100(loc(1));
L=[40,60,80,100];
Linv=1./L;
out = regressione_lineare(Linv,Tc)
x=[0:0.001:1]

plot(Linv,Tc,'linestyle','none','marker','x','MarkerSize',8,'color','r')
hold on
plot(x,out.b+out.m*x,'linestyle','--','color','b')
xlim([0,max(Linv)+0.005])
grid on
xlabel('$L^{-1}$','Interpreter','Latex','Fontsize',14);
ylabel('$T_c$ $[J\,k_B^{-1}]$','Interpreter','Latex','Fontsize',14);
title('Critical temperature - Linear Fit','Interpreter','Latex','Fontsize',18);
ll=legend(gca,'show','Data','Fit');
set(ll,'Interpreter','Latex','Location','northwest')
set(gca,'TickLabelInterpreter','latex');
