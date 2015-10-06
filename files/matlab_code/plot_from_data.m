

Tskip = 100;
Tgen  = 750000;
dt = 0.1;
tvec = (1:Tskip:Tgen)'*dt;
% plot(tvec, R_percap, 'b')

% EN__ = sprintf('Het%g-%g', 0.15,50);
EN__ = 'Hom'
KI__ = 'G';
Nc__ = 100;
hc__ = .5;
RTc__ = 25;
Kd_sd__ = 8;
cDisp__ = 0.1;
bHarvest__ = .2;
rImit__ = 0.02;


dir = 'data_figs_RT25'; 

simname = [EN__, '_', KI__, sprintf('I_T%g_N%g_RT%g_b%5.4f_ir%3.2f_ch0.08h2', Tgen/1000, Nc__, RTc__, bHarvest__, rImit__)];

hfname  = [dir, '\\h_' , simname, '.mat'];
rcfname = [dir, '\\rc_', simname, '.mat'];
rtfname = [dir, '\\rt_', simname, '.mat'];
kdfname = [dir, '\\kd_', simname, '.mat'];
vfname  = [dir, '\\v_' , simname, '.mat'];
ndfname = [dir, '\\nd_', simname, '.mat'];


load(hfname);
load(rcfname);
% load(rtfname);
load(kdfname);
% load(vfname);
% load(ndfname);

f2=figure; colormap(hot);
set(f2, 'Position', [1024 48 560 512])

subplot(3,1,1);
set(gca,'YDir','reverse')
[h c] = hist(h_vec, linspace(0,2,50));
imagesc(tvec, c, log(h))
set(gca,'YDir','normal')
set(gca, 'XTick', [], 'FontName', 'Cambria', 'FontSize', 15);
ylabel(sprintf('Harvesting\nrate'))

subplot(3,1,2);
set(gca,'YDir','reverse')
[h c] = hist(Kd_vec, linspace(0,10,50));
%image(hist(Kd_vec))
imagesc(tvec, c, log(h))
set(gca,'YDir','normal')
set(gca, 'XTick', [], 'FontName', 'Cambria', 'FontSize', 15);
ylabel(sprintf('Dispersal\ndistance'))

subplot(3,1,3);
set(gca,'YDir','reverse')
[h c] = hist(Rc_vec, linspace(0,500,50));
imagesc(tvec, c, log(h))
set(gca,'YDir','normal')
set(gca, 'XTick', [], 'FontName', 'Cambria', 'FontSize', 15);
ylabel(sprintf('Resource\nconsumed'))
xlabel('Time')


% f = figure; colormap(hot);
% set(f, 'Position', [450 48 560 948])
% 
% subplot(6,1,1);
% set(gca,'YDir','reverse')
% [h c] = hist(h_vec, 50);
% imagesc(tvec, c, log(h))
% set(gca,'YDir','normal')
% 
% hold on;
% plot(tvec, mean(h_vec, 1), 'w')
% ylabel('h')
% %ylim([0 3]);
% 
% 
% subplot(6,1,2);
% set(gca,'YDir','reverse')
% [h c] = hist(Rc_vec, 50);
% imagesc(tvec, c, log(h))
% set(gca,'YDir','normal')
% 
% hold on;
% plot(tvec, mean(Rc_vec, 1), 'w')
% ylabel('Rc')
% 
% 
% subplot(6,1,3);
% set(gca,'YDir','reverse')
% [h c] = hist(RT_vec, 50);
% % image(hist(RT_vec))
% imagesc(tvec, c, log(h))
% set(gca,'YDir','normal')
% hold on;
% plot(tvec, mean(RT_vec, 1), 'w')
% ylabel('RT')
% %ylim([0,50])
% 
% 
% subplot(6,1,4);
% set(gca,'YDir','reverse')
% [h c] = hist(Kd_vec, 50);
% %image(hist(Kd_vec))
% imagesc(tvec, c, log(h))
% set(gca,'YDir','normal')
% 
% hold on;
% plot(tvec, mean(Kd_vec,1), 'w')
% ylabel('Kd_sd')
% 
% 
% subplot(6,1,5);
% set(gca,'YDir','reverse')
% [h c] = hist(V_vec, 50);
% imagesc(tvec, c, log(h))
% set(gca,'YDir','normal')
% 
% hold on;
% plot(tvec, mean(V_vec, 1), 'w')
% ylabel('V')
% 
% 
% subplot(6,1,6);
% set(gca,'YDir','reverse')
% [h c] = hist(nD_vec, 50);
% imagesc(tvec, c, log(h))
% set(gca,'YDir','normal')
% 
% hold on;
% plot(tvec, mean(nD_vec, 1), 'Color', 'w')
% ylabel('nD')
% ylim([0 1])
% % plot(tvec, mean(zD_vec, 1), 'Color', [0.2 0.6 0.2])
% hold off;
% 
% %export_fig (sprintf('%s.png', filename));
% 
% set(f,'PaperSize',[5.60, 9.48],'PaperPosition',[0 0 5.60 9.48])
% print('-dpng','-r100',[dir, '\\', simname, '.png'])
% close(f);

% %%
% 


% f3 = figure;
% hist(dV_cumm, 1000)

% 
% 
% % plot(hc__vec_, Rc_scan', '-ob', hc__vec_, nD_scan*250, '-or', hc__vec_, Vc_scan*5, '-og')
% % % legend('N = 10', 'N = 50', 'N = 100')
% % xlabel('h')
% % title('Kd = 4, RTc = 10')
% % ylabel('Resource consumed')
% 
% %%
% % plot the evoloved traits as function of cost of dispersal, for various
% % population sizes
% 
% figure();
% subplot(1,3,1);
% plot(cD__vec_, rt_evol_scan_hc0_25', '-o')
% xlabel('Cost of dispersal')
% ylabel('Evolved dispersal threshold')
% legend('4', '8', '16', '32', '64')
% title('hc = 0.25')
% ylim([0 40])
% 
% subplot(1,3,2);
% plot(cD__vec_, rt_evol_scan_hc0_5', '-o')
% xlabel('Cost of dispersal')
% ylabel('Evolved dispersal threshold')
% legend('4', '8', '16', '32', '64', '128', '256')
% title('hc = 0.5')
% ylim([0 40])
% xlim([0, 1])
% 
% subplot(1,3,3);
% plot(cD__vec_, rt_evol_scan_hc1', '-o')
% xlabel('Cost of dispersal')
% ylabel('Evolved dispersal threshold')
% legend('4', '8', '16', '32', '64')
% title('hc = 1.0')
% ylim([0 40])
% 
% %% 
% % plot the evoloved traits as function of cost of dispersal, for various
% % population sizes
% 
% figure();
% subplot(1,2,1);
% plot(cD__vec_, hc_evol_scan_rt10_t200k', '-o')
% xlabel('Cost of dispersal')
% ylabel('Evolved harvesting rate')
% legend('4', '8', '16', '32', '64')
% title('RT = 10')
% ylim([0, 2.4])
% 
% subplot(1,2,2);
% plot(cD__vec_, hc_evol_scan_rt30_t200', '-o')
% xlabel('Cost of dispersal')
% ylabel('Evolved harvesting rate')
% legend('4', '8', '16', '32', '64')
% title('RT = 30')
% ylim([0, 2.4])
% 
