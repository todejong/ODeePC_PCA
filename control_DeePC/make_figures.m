close all; clear all; clc; 

load('plot_data_NMPC.mat')
load('plot_data.mat')

curr_fig = figure;
curr_axes1=axes('Parent',curr_fig,'FontSize',11,'FontName','Times New Roman');
box(curr_axes1,'on');
hold(curr_axes1,'all');
%your plots
subplot(2,1,1)
hold on
plot(r_plot(1:end),LineWidth=1.5,LineStyle="--",Color=[0 0 0])
plot(yk,LineWidth=1.5,Color="#0072BD")
plot(yk_NMPC,LineWidth=1.5,Color="#D95319")
legend('$r$','$y_{DeePC}$','$y_{NMPC}$',Interpreter='latex')
axis tight
grid on
ylabel('$y$',Interpreter='latex')
set(gca,'TickLabelInterpreter','Latex')
subplot(2,1,2)
hold on
plot(uk,LineWidth=1.5)
plot(uk_NMPC,LineWidth=1.5)
legend('$y_{DeePC}$','$y_{NMPC}$',Interpreter='latex')
ylabel('$u$',Interpreter='latex')
xlabel('$t [s]$',Interpreter='latex')
grid on
axis tight
%your plots
set(gca,'TickLabelInterpreter','Latex')
set(curr_fig,'Units','centimeters','PaperSize',[20.98 29.68],'PaperUnits','centimeters','PaperPosition',[0 0 12 8])
% you can change 9 and 6.3 to change the ratios of the plot...
savefig('product_kernel_control_comp.fig') %change it with the name you want to give to the .fig plot
print -depsc product_kernel_control_comp %change it with the name you want to give to the .eps figure
