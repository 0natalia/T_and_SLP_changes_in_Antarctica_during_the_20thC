%________________________________________________________________________
%           Zonal Wave 3 index and its correlation with T changes                             
% ZW3 index from Raphael (2004)

% Natália Silva; natalia3.silva@usp.br
% (2020)
%________________________________________________________________________

clear all; close all

%% SLP data
% ERA
pera = ('home/natalia/reanalises/ERA20C/era_slp.nc');
slpera = ncread(pera,'msl'); slpera = slpera/100; % hPa
slpera = cat(1,slpera, slpera(end,:,:)); slpera = [slpera(:,1,:) slpera];

lat = ncread(pera,'lat'); lat = [-90;lat]; 
lon = ncread(pera,'lon');  lon(end+1) = 180;
time2 = ncread(pera,'time'); tempo = datenum(datetime(1900,1,1)+hours(time2)); 
clear time2; clear pera; clear pnoaa

% selecionar os 3 pts
lat1 = find(lat<-47 & lat>-51); 
long = {find(lon>49 & lon<51), find(lon>165 & lon <167), find(lon<-75 & lon>-78)}; 

%% plot mean slp field
v = [960:1:990, 995:5:1000];
figure('color',[1 1 1],'position',[10 805 900 800]); 
m_proj('stereographic','lat',-90,'long',0,'radius',35);
[c] = m_contour(lon,lat(lat<-54),mean(slpera(:,lat<-54,:),3)',v, ...
    'linewidth', 2.7, 'color', 'k');
m_grid('xtick',[0 45 90 135 180 -45 -90 -135],'tickdir','out','ytick',...
    [-60 -75 -90],'tickdir','out', 'xaxisLocation', 'top', 'yaxisLocation', ...
    'middle','fontsize',10,'linestyle',':','linewidth',0.5);
m_coast('patch',[.7 .7 .7],'edgecolor','none');
title('PNM med (ERA20C)','Position',[0 0.7 1]);
clabel(c, 'manual','color', [0.5 0.5 0.5]); colormap(flipud(cmocean('balance')))
clear t; clear ttt; clear tt; clear ans;

%% ZW3 index

% climatologia
[~, sas, ~] = datevec(tempo);
iw3_era = []; iw3_noaa = [];
for mm = [1, 2, 3]
    
    djf = find(sas == 12 | sas == 1 | sas == 2); 
    mam = find(sas == 3 | sas == 4 | sas == 5);
    jja = find(sas == 6 | sas == 7 | sas == 8); 
    son = find(sas == 9 | sas == 10 | sas == 11);
   
    slp_era = squeeze(mean(slpera(long{mm},lat1,:),2));
    clmt_era = [mean(slp_era(djf)) mean(slp_era(mam))...
        mean(slp_era(jja)) mean(slp_era(son))];
    slp_era = [slp_era; 0; 0]; 
     
    % seasonal mean
    mntly_era = zeros(1,444); std_era = zeros(1,444); posi = 1;
    mntly_noaa = zeros(1,444); std_noaa = zeros(1,444); 
   
    for i = 1:3:length(mam)
        % ERA %
        mntly_era(posi) = mean(slp_era(mam(i:i+2)));
        mntly_era(posi+1) = mean(slp_era(jja(i:i+2)));
        mntly_era(posi+2) = mean(slp_era(son(i:i+2)));
        mntly_era(posi+3) = mean(slp_era(djf(i+2:i+4)));
        
        std_era(posi) = std(slp_era(mam(i:i+2)));
        std_era(posi+1) = std(slp_era(jja(i:i+2)));
        std_era(posi+2) = std(slp_era(son(i:i+2)));
        std_era(posi+3) = std(slp_era(djf(i+2:i+4)));

        posi = posi+4;
    end
    % adjusts
    mntly_era(end)=[]; mntly_era = [mean(slp_era([1,2])) mntly_era];
    std_era(end)=[]; std_era = [std(slp_era([1,2])) std_era];
    mntly_noaa(end)=[]; mntly_noaa = [mean(slp_noaa([1,2])) mntly_noaa];
    std_noaa(end)=[]; std_noaa = [std(slp_noaa([1,2])) std_noaa];
    clear i; clear posi;
    
    % standard units
    k = 1; izw3_era = zeros(size(mntly_era)); izw3_noaa = zeros(size(mntly_noaa));
    
    for i = 1:length(mntly_noaa)
    
        izw3_era(i) = ((mntly_era(i) - clmt_era(k)) / std_era(i));
    
        if k<4
            k = k+1;
        else
            k = 1;
        end
    end
    clear clmt_era; clear mntly_era; clear std_era; clear i; clear k
    clear clmt_noaa; clear mntly_noaa; clear std_noaa
    
    iw3_era = [iw3_era; izw3_era];
    clear izw3
    
end
clear djf; clear mam; clear jja; clear son; clear i; clear mm; clear lat1
clear sas; clear izw3_era; clear izw3_noaa
clear long; clear slp_era; clear slp_noaa;

% mean index 3 points
idx_era = mean(iw3_era); 
clear iw3_era; clear iw3_noaa

% outlier
condia = (idx_era > mean(idx_era)+4*std(idx_era) | idx_era < mean(idx_era)-4*std(idx_era)); 
condiv = (idx_noaa > mean(idx_noaa)+1*std(idx_noaa) | idx_noaa < mean(idx_noaa)-1*std(idx_noaa)); 
idx_noaa(condiv) = NaN; idx_era(condia) = NaN;
clear condia; clear condiv
    
x = 1:(length(idx_era)); inta = isnan(idx_era);
idx_era(inta) = interp1(x(~inta),idx_era(~inta),x(inta));
idx_era(1) = 0;

ve = 1:(length(idx_noaa)); int = isnan(idx_noaa);
idx_noaa(int) = interp1(ve(~int),idx_noaa(~int),ve(int));   
clear x; clear inta; clear ve; clear int

% Normalize
idx_era = idx_era/max(abs(idx_era)); idx_noaa = idx_noaa/max(abs(idx_noaa));

%% Calculate trends
lsq_e = fitlm(tempo(1:3:end),idx_era);
yo_e = lsq_e.Coefficients.Estimate(1);
trend_e = lsq_e.Coefficients.Estimate(2); erro_e = lsq_e.Coefficients.SE(2);
line_e = yo_e+trend_e*tempo(1:3:end); trend_e = (trend_e*10*4); erro_e = erro_e*10*4;
p_value_e = lsq_e.Coefficients.pValue(2); 

if p_value_e < 0.01
    SIG_e = ('p < 0.01');
elseif p_value_e < 0.05
    SIG_e = ('p < 0.05');
elseif p_value_e < 0.1
    SIG_e = ('p < 0.1');
else
    SIG_e = ('p > 0.1');
end

t = (['Trend_{ERA}: ', num2str(trend_e,'% .3e'),...
    ' \pm ', num2str(erro_e,'% .3e'),'; ', SIG_e]);

clear erro_e; clear erro_n; clear lsq_e; clear lsq_n; clear p_value_e; 
clear p_value_n; clear SIG_e; clear SIG_n; clear trend_e; clear trend_n
clear yo_e; clear yo_n

%% after 1950
tempo_novo = tempo(1:3:end); [ano,~,~] = datevec(tempo_novo);
ini = find(ano==1950); tempo_novo50 = tempo_novo(ini(1):end);
cinq_era = idx_era(ini(1):end); cinq_noaa = idx_noaa(ini(1):end);
clear ini; clear tempo_novo; clear ano

lsq_e_50 = fitlm(tempo_novo50,cinq_era); clear cinq_era
yo_e_50 = lsq_e_50.Coefficients.Estimate(1);
trend_e_50 = lsq_e_50.Coefficients.Estimate(2); erro_e_50 = lsq_e_50.Coefficients.SE(2);
line_e_50 = yo_e_50+trend_e_50*tempo_novo50; 
trend_e_50 = (trend_e_50*10*4); erro_e_50 = erro_e_50*10*4;
p_value_e_50 = lsq_e_50.Coefficients.pValue(2); 

if p_value_e_50 < 0.01
    SIG_e_50 = ('p < 0.01');
elseif p_value_e_50 < 0.05
    SIG_e_50 = ('p < 0.05');
elseif p_value_e_50 < 0.1
    SIG_e_50 = ('p < 0.1');
else
    SIG_e_50 = ('p > 0.1');
end

t50 = (['Trend_{ERA50}: ', num2str(trend_e_50, '% .3e'),...
    ' \pm ', num2str(erro_e_50,'% .3e'),'; ', SIG_e_50]);

%% PLOT

subplot(2,1,2); 
plot(tempo(1:3:end), idx_era, 'color', 'k', 'linewidth',1);
hold on; plot(tempo(1:3:end), line_e, 'color', 'r', 'linewidth',1.7)
plot(tempo_novo50, line_e_50, '--', 'color', 'g', 'linewidth',2.2); 
ylabel('I_{O3}'); datetick('x',10,'keepticks')
legend('ERA20C'); ylim([-1.1 1.1]); xlim([tempo(1) tempo(end)])
vline(tempo(601), ':', 'color', [0.7 0.7 0.7])
hline(0, '--', 'color', [0.7 0.7 0.7]); xlabel('Tempo'); 
text(tempo(6), -0.7, t, 'color', 'k', 'FontWeight', 'bold')
text(tempo(6), -0.9, t50, 'color', [0.4 0.4 0.4], 'FontSize', 8)

clear line_e; clear line_n;
clear t; clear tt; clear tempo; clear tempo_novo; clear tempo_novo50
clear line_e_50; clear line_n_50; clear t50; clear tt50

%% ZW3 x T correlation

% temperature data
% ERA
era20c = ('home/natalia/reanalises/ERA20C/era_temp.nc');
era20o = ('home/natalia/reanalises/ERA20C/era_sst.nc');
cont = ncread(era20c,'t'); cont = squeeze(cont(:,:,5,:));
sst = ncread(era20o,'sst');
cont = cont-273.15; sst = sst-273.15; % dgC
cont = cat(1,cont,cont(end,:,:)); cont = [cont(:,1,:) cont];
sst = cat(1,sst,sst(end,:,:)); sst = [sst(:,1,:) sst];
land = isnan(sst); temp2 = ones(size(sst));
temp2(land) = cont(land); temp2(~land) = sst(~land);
clear cont; clear sst; clear era20c; clear era20o; clear land

% remove T trend
detr_era = bsxfun(@minus,temp2(:,:,1:3:end),trend(temp2(:,:,1:3:end),[],3)); %total menos tendencia
% T anomaly
anom_era = bsxfun(@minus,detr_era,mean(detr_era,3));

% CORRELATION
corr_e_detr = zeros(193,20); p_e_detr = zeros(193,20); % zeros(lon,lat)
corr_n_detr = zeros(193,20); p_n_detr = zeros(193,20);

for i = 1:length(lon)
    for j = 1:length(latc)
        [ce_detr, pe_detr] = corrcoef(idx_era, anom_era(i,j,:),'alpha',0.01);
        [cn_detr, pn_detr] = corrcoef(idx_noaa, anom_noaa(i,j,:),'alpha',0.01);
        
        corr_e_detr(i,j) = ce_detr(2,1); p_e_detr(i,j) = pe_detr(2,1); 
        corr_n_detr(i,j) = cn_detr(2,1); p_n_detr(i,j) = pn_detr(2,1); 
    end
end
clear ce_detr; clear pe_detr; clear cn_detr; clear pn_detr
clear i; clear j; clear anom_era; clear anom_noaa

% significative corr: alpha = 1%
pp_detr = p_e_detr < 0.01; p_e_detr(~pp_detr) = NaN; 
ppp_detr = p_n_detr < 0.01; p_n_detr(~ppp_detr) = NaN;
clear pp_detr; clear ppp_detr

% CORR ERA 20C
v=(-1:0.01:1);
figure('color',[1 1 1],'position',[10 805 900 800]); 
colormap(flipud(cbrewer('div','Spectral',80)))
m_proj('stereographic','lat',-90,'long',0,'radius',35);
m_contourf(lon,latc,corr_e_detr',v); hold on
m_contour(lon,latc,corr_e_detr',v)
m_contour(lon,latc,p_e_detr',':k')
m_grid('xtick',[0 45 90 135 180 -45 -90 -135],'tickdir','out','ytick',...
    [-60 -75 -90],'tickdir','out', 'xaxisLocation', 'top', ...
    'yaxisLocation', 'middle','fontsize',10,'linestyle',':','linewidth',0.5);
m_coast('color','k','linewidth',1.5); caxis([-0.25 0.25])
title('Corr I_{O3} x T_{detrend} (ERA20C)','Position',[0 0.7 1]);
t = colorbar('Position',[0.91 0.1 0.03 0.8]); tt = get(t,'title'); set(tt,'string','r')
clear t; clear ttt; clear tt; clear ans;


%% CORRELACAO ZW3 x PNM
% ERA
% remove SLP trend
detr_era = bsxfun(@minus,slpera(:,:,1:3:end),trend(slpera(:,:,1:3:end),[],3)); %total menos tendencia
% SLP anomaly
anom_era = bsxfun(@minus,slpera(:,:,1:3:end),mean(slpera(:,:,1:3:end),3));

% CORRELATION SLP detrend
corr_e_detr = zeros(193,20); p_e_detr = zeros(193,20); % zeros(lon,lat)
corr_n_detr = zeros(193,20); p_n_detr = zeros(193,20);

for i = 1:length(lon)
    for j = 1:length(latc)
        [ce_detr, pe_detr] = corrcoef(idx_era, detr_era(i,j,:),'alpha',0.01);
        [cn_detr, pn_detr] = corrcoef(idx_noaa, detr_noaa(i,j,:),'alpha',0.01);
        
        corr_e_detr(i,j) = ce_detr(2,1); p_e_detr(i,j) = pe_detr(2,1); 
        corr_n_detr(i,j) = cn_detr(2,1); p_n_detr(i,j) = pn_detr(2,1); 
    end
end
clear ce_detr; clear pe_detr; clear cn_detr; clear pn_detr
clear i; clear j

% corr significativa ao nivel de 1%
pp_detr = p_e_detr < 0.01; p_e_detr(~pp_detr) = NaN; 
ppp_detr = p_n_detr < 0.01; p_n_detr(~ppp_detr) = NaN;
clear pp_detr; clear ppp_detr

% CORR ERA 20C
v=(-1:0.01:0);
figure('color',[1 1 1],'position',[10 805 900 800]); 
colormap(flipud(cmocean('curl')))
m_proj('stereographic','lat',-90,'long',0,'radius',35);
m_contourf(lon,latc,corr_e_detr',v); hold on
m_contour(lon,latc,corr_e_detr',v)
m_contour(lon,latc,p_e_detr',':k')
m_grid('xtick',[0 45 90 135 180 -45 -90 -135],'tickdir','out','ytick',...
    [-60 -75 -90],'tickdir','out', 'xaxisLocation', 'top', ...
    'yaxisLocation', 'middle','fontsize',10,'linestyle',':','linewidth',0.5);
m_coast('color','k','linewidth',1.5); caxis([-0.25 0.25])
title('Corr I_{O3} x PNM_{detrend} (ERA20C)','Position',[0 0.7 1]);
t = colorbar('Position',[0.91 0.1 0.03 0.8]); tt = get(t,'title'); set(tt,'string','r')
clear t; clear ttt; clear tt; clear ans;
