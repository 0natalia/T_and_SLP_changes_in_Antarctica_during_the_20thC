%_________________________________________________________________________
% Southern Annular Mode (SAM) index and its correlation with T/SLP changes                             
% SAM index from Gong & Wang (1999)
% NOAA20CR and ERA20C reanalysis datasets

% Natália Silva; natalia3.silva@usp.br
% (2020)
%________________________________________________________________________

clear all; close all

% load data
load('home/natalia/rotinas_MATLAB/TG/load_T_PNM_NOAA_ERA.mat')

slp = {slp1, slp2}; clear slp1; clear slp2; 
temp = {temp1, temp2}; clear temp2; clear temp1
var_lab = {'NOAACR', 'ERA20C'};

slp40a = zeros(193, 2, 111, 2); slp65a = zeros(193, 2, 111, 2); 
for k = 1:2
    [slp40a(:,:,:,k),tanos] = downsample_ts(slp{k}(:, lat>-42 & lat<-39, :), tempo, 'year');
    [slp65a(:,:,:,k),~] = downsample_ts(slp{k}(:, lat<-64 & lat> -66.5, :), tempo, 'year');
end; clear k

sam_anual = zeros(111,2); sam = zeros(1332,2);
slp40a = squeeze(mean(slp40a,2)); slp65a = squeeze(mean(slp65a,2));
slp40 = zeros(193, 1332, 2); slp65 = zeros(193, 1332, 2); 
meses = linspace(1,1332,1332); month = zeros(193,95,12); 
anom = zeros(193, 95, 1332, 2); latc = lat(lat<-54); v=(-1:0.01:1); 
corr = zeros(193,20,2); p = zeros(193,20,2); 

%%
for i = 1:2 
    %% SAM
    % Anual
    sam_anual(:,i)=(squeeze(mean(slp40a(:,:,i)))-mean2(slp40a(:,:,i)))/std(mean(slp40a(:,:,i)))...
        - (squeeze(mean(slp65a(:,:,i))) - mean2(slp65a(:,:,i)))/std(mean(slp65a(:,:,i))); 
    % Mensal
    slp40(:,:,i) = squeeze(mean(slp{i}(:, lat>-42 & lat<-39, :),2)); 
    slp65(:,:,i) = squeeze(mean(slp{i}(:, lat<-64 & lat> -66.5, :),2));
    sam(:,i) = (squeeze(mean(slp40(:,:,i)))-mean2(slp40(:,:,i)))/std(mean(slp40(:,:,i)))...
        - (squeeze(mean(slp65(:,:,i))) - mean2(slp65(:,:,i)))/std(mean(slp65(:,:,i)));
    
    lsq = fitlm(meses, sam(:,i)); p = lsq.Coefficients.pValue(2);
    line = lsq.Coefficients.Estimate(1)+lsq.Coefficients.Estimate(2)*meses; 
    if p < 0.01; SIG = ('p < 0.01');
    elseif p < 0.05; SIG = ('p < 0.05');
    elseif p < 0.1; SIG = ('p < 0.1');
    else SIG = ('p > 0.1');
    end
    t = (['Tend. SAM_{', var_lab{i}, '}: ', num2str(lsq.Coefficients.Estimate(2)*120,'% .3e'), ...
        ' \pm ', num2str(lsq.Coefficients.SE(2)*120,'% .3e'), '; ', SIG]);
    clear p; clear lsq; clear SIG
    
    figure('color',[1 1 1],'position',[108 305 850 700]); 
    subplot(2,1,i); plot(tempo, sam(:,i), 'color', [0.87 0.87 0.87],'linewidth',1);
    hold on; plot(tanos, sam_anual(:,i),'k','linewidth',1.3); legend(var_lab{i})
    plot(tempo,line,'r','linewidth',1.7); text(tanos(2),4,t); 
    datetick('x',10,'keepticks'); xlim([tanos(1) tanos(end)]); ylim([-5 5])
    ylabel('I_{SAM}'); hline(0, '--', 'color', [0.7 0.7 0.7]);
    if i == 1; title('Indice SAM'); 
    end
    
    %% T anomaly
    for w = 1:12
        month(:, :, w) = mean(temp{i}(:, :, w:12:end), 3);
        anom(:, :, w:12:end, i) = bsxfun(@minus, squeeze(temp{i}(:, :, w:12:end)), ...
            squeeze(month(: ,:, w)));
    end; clear w; clear month
    
    %% correlation SAM x T
    for k = 1:length(lon)
        for j = 1:length(latc)
            [corr_c, p_c] = corrcoef(sam(:,i), anom(k,j,:,i),'alpha',0.01);
            corr(k,j,i) = corr_c(2,1); p(k,j,i) = p_c(2,1);
        end; clear j; clear corr_c; clear p_c
    end; clear k; pp = p < 0.01; p(~pp) = NaN; 
    
    figure('color',[1 1 1],'position',[10 805 900 800]);
    colormap(flipud(cbrewer('div','Spectral',80)))
    m_proj('stereographic','lat',-90,'long',0,'radius',35);
    m_contourf(lon, latc, corr(:,:,i)',v); hold on
    m_contour(lon, latc, corr(:,:,i)',v)
    m_contour(lon, latc, p(:,:,i)',':k')
    m_grid('xtick',[0 45 90 135 180 -45 -90 -135],'tickdir','out','ytick',...
        [-60 -75 -90],'tickdir','out', 'xaxisLocation', 'top', ...
        'yaxisLocation', 'middle','fontsize',10,'linestyle',':','linewidth',0.5);
    m_coast('color','k','linewidth',1.5); caxis([-0.3 0.3])
    title(['Corr SAM x T_{anom}', var_lab{i}],'Position',[0 0.7 1]);
    t = colorbar('Position',[0.91 0.1 0.03 0.8]); tt = get(t,'title'); set(tt,'string','r')
    clear t; clear ttt; clear tt; clear ans;
    
end; clear i; clear tanos; clear t; clear meses; clear lat; clear pp
clear line; clear sam_anual; clear slp40a; clear slp65a; clear slp40; clear slp65
clear v; clear tempo

