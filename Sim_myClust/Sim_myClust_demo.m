clear all;
close all
clc
%%
set(0,'defaultAxesFontSize',8);
set(0,'defaultAxesFontName','Arial');
set(0,'defaultTextFontSize',8);
set(0,'defaultTextFontName','Arial');


scrsz = get(0,'ScreenSize');
X = 42.0;                  % two times A3 paper size
Y = 29.7;                  % two times A3 paper size
%%
cdir=cd;
addpath(genpath(cdir));

addpath('C:\toolbox\eeglab10.2.2.4b')
eeglab
close all
cd (cdir);

%%
fs = 1024;
t = 0:1/fs:4.6;
tst = 2.3;% task onset time in simulation
tspan = t;
wn=[1 100]./(fs/2); % cut-off frequency of bandpass filter in preprocessing 

save_flag=0;  % If set as save_flag = 1, simulation result is saved in mat-file. 
F = 30;
Trinum = 30;
Rwin=5;

savedir = [cdir, '\figures\'];
mkdir(savedir);

%%
if exist([cdir, '\mat\rnd_seed.mat'], 'file')==0
    data = zeros(6, length(t), Trinum);
    bindx = find(t >= -0.4 & t<= 0);
    
    [Tout,theta1, signal1, rnd_seed1]=Kuramoto_sim_v4(3, F, tspan, [.1, .4]+tst);
    [Tout,theta2, signal2, rnd_seed2]=Kuramoto_sim_v4(3, F, tspan, [.4, .7]+tst);
    
        %%
    parfor tri = 1:size(data,3)
        %%
        rnd_seed3(tri)=rng;
        data(:,:,tri) = [signal1;signal2]+normrnd(0, 1, size([signal1;signal2]));
        data(:,:,tri) = preprocessiong(data(:,:,tri),bindx,wn,'off');
        %%
        disp(tri)
    end
    t=tspan-tst;
    data = data(:, t>=-.5 & t<=1.5, :);
    t = t(t>=-0.5 & t<=1.5);
else
    load([cdir, '\mat\rnd_seed.mat']);    
    disp('data loading')
    t = 0:1/fs:4.6;
    data = zeros(6, length(t), Trinum);
    bindx = find(t >= -0.4 & t<= 0);
    
    %% Regenerate the template signals  
    [Tout,theta1, signal1, rnd_seed1]=Kuramoto_sim_v4(3, F, tspan, [.1, .4]+tst, rnd_seed1);
    [Tout,theta2, signal2, rnd_seed2]=Kuramoto_sim_v4(3, F, tspan, [.4, .7]+tst, rnd_seed2);
    
        %%
    parfor tri = 1:size(data,3)
        %%
        rng(rnd_seed3(tri));% reinitialize the random generator
        data(:,:,tri) = [signal1;signal2]+normrnd(0, 1, size([signal1;signal2]));
        data(:,:,tri) = preprocessiong(data(:,:,tri),bindx,wn,'off');
        %%
        disp(tri)
    end
    t=tspan-tst;
    data = data(:, t>=-.5 & t<=1.5, :);
    t = t(t>=-0.5 & t<=1.5);
end
%% plot template signal
fig = figure;
figure_setting(X*.8, Y*.7, fig);

t_tmp = tspan-tst;
subplot(2,2,1);
plot(t_tmp, theta1);
xlim([-.5, 1.5]); 
set(gca, 'Ytick', [-pi, 0, pi], 'Yticklabel', {'3.14', '0', '3.14'})
xlabel('time (s)');
ylabel('phase (rad)');

subplot(2,2,2);
plot(t_tmp, theta2);
xlim([-.5, 1.5]); %ylim([-3.14, 3.14]);
set(gca, 'Ytick', [-pi, 0, pi], 'Yticklabel', {'3.14', '0', '3.14'})
xlabel('time (s)');
ylabel('phase (rad)');

subplot(2,2,3);
plot(t_tmp, signal1);
xlim([-.5, 1.5]);ylim([-3, 3]);
xlabel('time (s)');
ylabel('amplitude (a.u.)');

subplot(2,2,4);
plot(t_tmp, signal2);
xlim([-.5, 1.5]);ylim([-3, 3]);
xlabel('time (s)');
ylabel('amplitude (a.u.)');

fname = [savedir, 'template_signal_Kuramoto'];
figure_save(fig, fname);
%% time-frequency analysis
tlimit = [t(1), t(end)];
twin=0.27;
tstep=0.01;
ntimesout=round(((length(t)-twin*fs)/fs)/tstep);% # of time samples for tf analysis

winsize=twin*fs;% sample number
tlimits=tlimit.*1000;
baseline=[-400 0];
flimits=[0, 40];
padratio=4;%% mod: padradio=4 -> padradio=2

parfor ch=1:size(data,1)
    tmp=squeeze(data(ch,:,:));
    [tf, f, t] = timefreq(tmp, fs, 'tlimits', tlimits, 'winsize', winsize, 'ntimesout', ntimesout, 'freqs', flimits, 'padratio', padratio, 'cycles', [1, .1]);
    alltf(:,:,:,ch)     =tf; % tf:[freqs x t x trial]
    freqs(ch,:) = f;
    t_tf(ch,:)  = t;
end
freqs = freqs(1,:);
t_tf = t_tf(1,:)./1000;
%% Calculate wPLI
mode  = 0;
wpli_tmp = cal_wpli(alltf, mode);% Set "mode = 0" -> Calculate normal w-PLI
                                 % Set "mode = 1" -> Calculate debiased w-PLI
bindx=find(t_tf>=-.4 & t_tf<=0); % Extract baseline index
fig=figure;
figure_setting(X, Y, fig); % Set figure configuration
%% Plot z-wPLI
count = 1;
for roi1=1:size(data,1)
    for roi2=1:size(data,1)
        if roi1~=roi2
            subplot(size(data,1), size(data,1), count)
            
            % z-scored wPLI
            val  = z_score(wpli_tmp(:,:,roi2,roi1), bindx); 
            % plot data in each channel-pairs
            imagesc(t_tf, freqs, val, [-6 6]);
            xlim([-.2, 1.1])
            ylim([20, 36])

            if roi1==3 && roi2==1
                text(-1, 40, 'frequency (Hz)', 'HorizontalAlignment', 'center', 'fontsize', 12, 'Rotation', 90)
            elseif roi1==size(data,1) && roi2==3
                text(1.25, 48, 'duration (s)', 'HorizontalAlignment', 'center', 'fontsize', 12)
            elseif roi1==1 && roi2==3
                text(1.25, 12, 'z-scored wPLI', 'HorizontalAlignment', 'center', 'fontsize', 15)
            elseif roi1==size(data,1) && roi2==size(data,1)
                cbar('vert', 33:64);
            end
        elseif roi1==roi2
            subplot(size(data,1), size(data,1), count)
            xlim([0, 1])
            ylim([0, 1])
            text(0.5, 0.5, ['Signal ', num2str(roi1)], 'HorizontalAlignment', 'center', 'fontsize', 12)
            axis off
            
            if roi1==size(data,1) && roi2==size(data,1)
                cbar('vert', 1:64);
                if mode ==1
                    ylabel('squared-wPLI (zscore)')
                else
                    ylabel('wPLI (zscore)')
                end
                set(gca, 'Ytick', [0, 1], 'YtickLabel', {'-6','6'});
            end
        end
        count = count +1;
    end
end
%%
fname = [savedir, 'wpli_all'];
figure_save(fig, fname); % Save figure
%% Data extraction of wPLI corresponding the relevant frequency
frq = find(freqs >= 28 & freqs<= 36); % Extract index
wpli_tmp = wpli_tmp(frq,:,:,:);
wpli_tmp = mean(z_score(wpli_tmp, bindx), 1);% z-scored and averaging over frequency

bindx=find(t_tf>=-.4 & t_tf<=0);
count=1;
for roi1=1:size(data,1)
    for roi2=roi1:size(data,1)
        if roi1~=roi2
            wpli(count,:) = wpli_tmp(:, :, roi2, roi1);
            count = count+1;

            disp(['pair-', num2str(count),': ROI', num2str(roi2), ' vs ROI', num2str(roi1)])
        end
    end
end
%% calculate DTW distance    
count=0;
for i=1:size(wpli,1)
    for j=1:size(wpli,1)
        Dist_dtw(j,i)=dtw_c(wpli(j, t_tf>=-0.0 & t_tf<=1.2)', wpli(i,t_tf>=-0.0 & t_tf<=1.2)', Rwin);
        % Compute DTW distance in each pair of z-wPLI with window size Rwin
    end
end
%% calculate average linkage
Z = linkage(Dist_dtw, 'average');
for num=1:10
    clust.DTW(:,num) = cluster(Z, 'maxclust', num);
end
%% evaluate cluster 
eva_DTW = evalclusters(Dist_dtw, clust.DTW,'CalinskiHarabasz')
% evaluate optimal cluster number with Pseudo-F index
%% Plot cluser averaged z-wPLI
clst_lab = clust.DTW(:,eva_DTW.OptimalK);
color = jet(eva_DTW.OptimalK+1);
    % hold on
fig=figure;
figure_setting(X, Y, fig);
for num = 1:eva_DTW.OptimalK
    
    subplot(2,5,num)
    hold on 
    plot(t_tf(t_tf>=-0.4 & t_tf<=1.2), wpli(clst_lab==num,t_tf>=-0.4 & t_tf<=1.2), 'color', [.6, .6, .6])
    plot(t_tf(t_tf>=-0.4 & t_tf<=1.2), mean(wpli(clst_lab==num,t_tf>=-0.4 & t_tf<=1.2),1), 'color', color(num+1,:), 'linewidth',2)
    
    
    hold off
    xlim([-.1, 1.2]);
    ylim([-4 12]);
    xlabel('duration (s)')
    ylabel('cluster-averaged wPLI (z-scored)')
    
    title(['cluster ', num2str(num)])
    
end
%%
fname = [savedir, 'cluster_ave_wpli'];
figure_setting(X, Y, fig);
figure_save(fig, fname);

%% Plot cluster label matrix
fig=figure;
figure_setting(X*.6, Y, fig);
imagesc(squareform(clust.DTW(:,eva_DTW.OptimalK)))

set(gca, 'XTick', 0.5:1:6.5, 'YTick', 0.5:1:6.5, 'XTicklabel', {[], '1', '2', '3', '4', '5', '6'}, 'YTicklabel', {[], '1', '2', '3', '4', '5', '6'})
xlabel('referece ROI');
ylabel('seed ROI');
colormap(jet(eva_DTW.OptimalK+1));
c=colorbar('location', 'SouthOutside');

pos = get(c, 'Position');

set(c, 'XTick', 0:3, 'XTicklabel', {[], '1', '2', '3'} ,'Ticklength',[0 0], 'Xlim', [1 eva_DTW.OptimalK], 'Position', [pos(1), pos(2)-0.117, pos(3), pos(4)]);
axis square
title('DTW')
grid on
%%
fname = [savedir, 'DTW_clust_label'];
figure_save(fig, fname);

if save_flag==1
    %%
    save_mat_dir = [cdir, '\mat\'];
    mkdir(save_mat_dir);
    time_str=datestr(now);
    time_str=[time_str(end-7:end-6), time_str(end-4:end-3)];
    fname = [save_mat_dir, 'Sim_myClust_dataset',date, '_', time_str,'.mat'];
    save(fname, 'data','t',...
                'Tout',...
                'rnd_seed1', 'rnd_seed2', ...
                'rnd_seed3',...
                'signal1', 'signal2',...
                'theta1', 'theta2',...
                'alltf', 'freqs', 't_tf', ...
                'wpli', ...
                'clust', 'eva_DTW', '-v7.3')
end