%_________________________________________________________________________
%                     Temperature time series and treds
% Calculates trends for the entire antactic regions;
% Compare trends in the continent vs. in the ocean
% Compare trends in the Antarctic Peninsula vs. in East Antarctica
 
% Natalia Silva; natalia3.silva@usp.br
% (2020)
%_________________________________________________________________________

clear all; close all
gap=[.1,.05]; marg_h=[.1,.05]; marg_w=[.1,.05]; % subtightplot

%% Read data: ECMWF- ERA20C
era20c = ('home/natalia/reanalises/ERA20C/era_temp.nc');
era20o = ('home/natalia/reanalises/ERA20C/era_sst.nc');
lat = ncread(era20o,'lat'); lon = ncread(era20o,'lon'); 
time = ncread(era20c,'time'); tempo = datenum(datetime(1900,1,1)+hours(time)); 
clear time
cont = ncread(era20c,'t'); cont = squeeze(cont(:,:,5,:));
sst = ncread(era20o,'sst');
cont = cont-273.15; sst = sst-273.15; % dgC
cont = cat(1,cont,cont(end,:,:)); lon(end+1) = 180; % fechar 360dg
cont = [cont(:,1,:) cont]; lat = [-90;lat]; % fechar -90dg
sst = cat(1,sst,sst(end,:,:)); sst = [sst(:,1,:) sst];

land = isnan(sst); temp2 = ones(size(sst));
temp2(land) = cont(land); temp2(~land) = sst(~land);
cont(~land) = nan; 

% Read data: NOAA - CIRES 20CR
noaa = ('home/natalia/reanalises/NOAA20CR/temp.mon.mean.nc');
temp1 = ncread(noaa,'air'); temp1 = temp1(:,:,349:1680); 
temp1 = temp1-273; % dgC
temp1 = cat(1,temp1,temp1(end,:,:)); % fechar em 360lon
temp1 = cat(2,temp1,temp1(:,1,:)); % fechar em 90lat

[tnoaa_ann,tanos] = downsample_ts(temp1,tempo,'year');
[tera_ann,~] = downsample_ts(temp2,tempo,'year');

clear noaa; clear era20c; clear era20o; clear land

anos = linspace(1900,2010,length(tanos));
austral = lat<-54; clear lat

%% Time series
% NOAA
tnoaa_anual_smean = squeeze(mean(mean(tnoaa_ann(:,austral,:))));
clear tnoaa_ann
% Trend
lsqma1 = fitlm(anos,tnoaa_anual_smean);
yoa1 = lsqma1.Coefficients.Estimate(1);
trenda1 = lsqma1.Coefficients.Estimate(2); erroa1 = lsqma1.Coefficients.SE(2);
linea1 = yoa1+trenda1*anos; trenda1 = (trenda1*10); erroa1 = erroa1*10;
p_valuea1 = lsqma1.Coefficients.pValue(2); 

if p_valuea1 < 0.01
    SIGa1 = ('p < 0.01');
elseif p_valuea1 < 0.05
    SIGa1 = ('p < 0.05');
elseif p_valuea1 < 0.1
    SIGa1 = ('p < 0.1');
else
    SIGa1 = ('p > 0.1');
end
t = (['Tend. Anual_{NOAA}: ', num2str(trenda1,'% .3e'), ...
    ' \pm ', num2str(erroa1,'% .3e'), '^oC/dec; ', SIGa1]);
clear stra1; clear strEra1; clear SIGa1; clear p_valuea1; 
clear trenda1; clear erroa1; clear lsqma1;clear yoa1;

% T anomaly in Antarctica
t_mean1 = mean(tnoaa_anual_smean); anom1 = [];
for i = 1:length(tanos)
anom1 = [anom1 tnoaa_anual_smean(i) - t_mean1];
end
clear i; clear t_mean1;

% ERA20C
tera_anual_smean = squeeze(mean(mean(tera_ann(:,austral,:))));
clear tera_ann

lsqma = fitlm(anos,tera_anual_smean);
yoa = lsqma.Coefficients.Estimate(1);
trenda = lsqma.Coefficients.Estimate(2); erroa = lsqma.Coefficients.SE(2);
linea = yoa+trenda*anos; trenda = (trenda*10); erroa = erroa*10;
p_valuea = lsqma.Coefficients.pValue(2); 
 
if p_valuea < 0.01
    SIGa = ('p < 0.01');
elseif p_valuea < 0.05
    SIGa = ('p < 0.05');
elseif p_valuea < 0.1
    SIGa = ('p < 0.1');
else
    SIGa = ('p > 0.1');
end
tt = (['Tend Anual_{ERA}: ', num2str(trenda,'% .3e'), ...
    ' \pm ', num2str(erroa,'% .3e'), '^oC/dec; ', SIGa]);
clear stra; clear strEra; clear SIGa; clear p_valuea; 
clear trenda; clear erroa; clear lsqma;clear yoa;

% T anomaly in Antarctica
t_mean = mean(tera_anual_smean); anom = [];

for i = 1:length(tanos)
    anom = [anom tera_anual_smean(i) - t_mean];
end
clear i; clear t_mean;

%% PLOT ANNUAL + TREND
figure('color',[1 1 1],'position',[108 305 850 700]);
subtightplot(2,1,1,gap,marg_h,marg_w); plot(tanos,tnoaa_anual_smean,'color', [0.6 0.6 1],'linewidth',1);
hold on; plot(tanos,tera_anual_smean,'color', [0.8706 0.8157 0.4471], 'linewidth',1);
title('Tendencia_T Anual');
xlabel('Anos');ylabel('^oC');datetick('x',10,'keepticks')
plot(tanos, linea1, 'color', [0.55 0.55 0.9],'linewidth',1.7)
plot(tanos, linea, 'color', [ 0.7804    0.6745    0.1216], 'linewidth',1.7) 
legend('NOAA20CR','ERA20C')
xlim([tanos(1) tanos(end)]); ylim([-19 -7])
text(tanos(2),-7.8, t,'color', [0.55 0.55 0.9])
text(tanos(2),-8.8, tt,'color', [ 0.7804    0.6745    0.1216])

subtightplot(2,1,2,gap,marg_h,marg_w);
plot(tanos,anom1','color', [0.6 0.6 1], 'linewidth',1.5); hold on
plot(tanos,anom','color', [0.8706 0.8157 0.4471], 'linewidth',1.5);
title('Anomalia');datetick('x',10,'keepticks')
xlabel('Anos');ylabel('^oC');ylim([-1.5 1.5]);xlim([tanos(1) tanos(end)])
plot([tanos(1) tanos(end)],[0 0],'Color',[0.8 0.8 0.8],'LineStyle','--')
legend('NOAA20CR','ERA20C')

clear tnoaa_anual_smean; clear anom1; clear linea1; clear t; clear tt
clear tera_anual_smean; clear anom; clear linea;


%% CONTINENT vs OCEAN

% NOAA
noaa_mask = ('home/natalia/reanalises/NOAA20CR/NOAA_landmask.nc');
noaa_m = ncread(noaa_mask, 'land');
noaa_m = cat(1,noaa_m,noaa_m(end,:,:)); noaa_m = [noaa_m(:,1,:) noaa_m];

noaa_sst = temp1; noaa_land = temp1; oc = noaa_m == 0; 
noaa_oc = ones(size(noaa_sst));
clear noaa_m

for i = 1:length(tempo)
    noaa_oc(:,:,i)=oc;
end
clear i; clear oc

noaa_sst(noaa_oc==0) = NaN; noaa_land(noaa_oc~=0) = NaN;
clear noaa_mask; clear noaa_oc;

%************************ OCEAN ****************************
noaa_sst_anual = downsample_ts(noaa_sst(:,austral, :),tempo,'year');
noaa_sst_anual(noaa_sst_anual < -5) = NaN;

tnoaa_oc = squeeze(nanmean(nanmean(noaa_sst_anual))); clear noaa_sst_anual
lsq_noaa_oc = fitlm(anos,tnoaa_oc);
yo_noaa_oc = lsq_noaa_oc.Coefficients.Estimate(1);
trend_noaa_oc = lsq_noaa_oc.Coefficients.Estimate(2); 
erro_noaa_oc = lsq_noaa_oc.Coefficients.SE(2);
line_noaa_oc = yo_noaa_oc+trend_noaa_oc*anos; 
trend_noaa_oc = (trend_noaa_oc*10); erro_noaa_oc = erro_noaa_oc*10;
p_noaa_oc = lsq_noaa_oc.Coefficients.pValue(2); 

if p_noaa_oc < 0.01
    SIG_noaa_oc = ('p < 0.01');
elseif p_noaa_oc < 0.05
    SIG_noaa_oc = ('p < 0.05');
elseif p_noaa_oc < 0.1
    SIG_noaa_oc = ('p < 0.1');
else
    SIG_noaa_oc = ('p > 0.1');
end

tt_noaa_oc = (['Tend Oceano_{NOAA}: ', num2str(trend_noaa_oc,'% .3e'), ...
    ' \pm ', num2str(erro_noaa_oc,'% .3e'), '^oC/dec; ', SIG_noaa_oc]);
clear SIG_noaa_oc; clear p_noaa_oc; clear yo_noaa_oc;
clear trend_noaa_oc; clear erro_noaa_oc; clear lsq_noaa_oc;

anom_noaa_oc = []; 
for i = 1:length(tanos)
    anom_noaa_oc(i) = (tnoaa_oc(i) - mean(tnoaa_oc));
end
clear i

%******************** CONTINENT ************************
noaa_land_anual = downsample_ts(noaa_land(:,austral,:),tempo,'year');
noaa_land_anual(noaa_land_anual > 0) = NaN;

tnoaa_lnd = squeeze(nanmean(nanmean(noaa_land_anual))); clear noaa_land_anual
lsq_noaa_lnd = fitlm(anos,tnoaa_lnd);
yo_noaa_lnd = lsq_noaa_lnd.Coefficients.Estimate(1);
trend_noaa_lnd = lsq_noaa_lnd.Coefficients.Estimate(2); 
erro_noaa_lnd = lsq_noaa_lnd.Coefficients.SE(2);
line_noaa_lnd = yo_noaa_lnd+trend_noaa_lnd*anos; 
trend_noaa_lnd = (trend_noaa_lnd*10); erro_noaa_lnd = erro_noaa_lnd*10;
p_noaa_lnd = lsq_noaa_lnd.Coefficients.pValue(2); 

if p_noaa_lnd < 0.01
    SIG_noaa_lnd = ('p < 0.01');
elseif p_noaa_lnd < 0.05
    SIG_noaa_lnd = ('p < 0.05');
elseif p_noaa_lnd < 0.1
    SIG_noaa_lnd = ('p < 0.1');
else
    SIG_noaa_lnd = ('p > 0.1');
end
tt_noaa_lnd = (['Tend Continente_{NOAA}: ', num2str(trend_noaa_lnd,'% .3e'), ...
    ' \pm ', num2str(erro_noaa_lnd,'% .3e'), '^oC/dec; ', SIG_noaa_lnd]);
clear str_noaa_lnd; clear strEr_noaa_lnd; clear SIG_noaa_lnd; clear p_noaa_lnd; 
clear trend_noaa_lnd; clear erro_noaa_lnd; clear lsq_noaa_lnd; clear yo_noaa_lnd;

anom_noaa_lnd = []; 
for i = 1:length(tanos)
    anom_noaa_lnd(i) = (tnoaa_lnd(i) - mean(tnoaa_lnd));
end
clear i

%%% ERA 
%************************** OCEANO *********************************
era_sst_anual = downsample_ts(sst(:,austral, :),tempo,'year'); clear sst

tera_oc = squeeze(nanmean(nanmean(era_sst_anual))); clear era_sst_anual
lsq_era_oc = fitlm(anos,tera_oc);
yo_era_oc = lsq_era_oc.Coefficients.Estimate(1);
trend_era_oc = lsq_era_oc.Coefficients.Estimate(2); 
erro_era_oc = lsq_era_oc.Coefficients.SE(2);
line_era_oc = yo_era_oc+trend_era_oc*anos; 
trend_era_oc = (trend_era_oc*10); erro_era_oc = erro_era_oc*10;
p_era_oc = lsq_era_oc.Coefficients.pValue(2); 

if p_era_oc < 0.01
    SIG_era_oc = ('p < 0.01');
elseif p_era_oc < 0.05
    SIG_era_oc = ('p < 0.05');
elseif p_era_oc < 0.1
    SIG_era_oc = ('p < 0.1');
else
    SIG_era_oc = ('p > 0.1');
end
tt_era_oc = (['Tend Oceano_{ERA}: ', num2str(trend_era_oc,'% .3e'), ...
    ' \pm ', num2str(erro_era_oc,'% .3e'), '^oC/dec; ', SIG_era_oc]);
clear SIG_era_oc; clear p_era_oc; clear yo_era_oc;
clear trend_era_oc; clear erro_era_oc; clear lsq_era_oc;

anom_era_oc = []; 
for i = 1:length(tanos)
    anom_era_oc(i) = (tera_oc(i) - mean(tera_oc));
end
clear i

%*********************** CONTINENTE ***************************
era_land_anual = downsample_ts(cont(:,austral,:), tempo,'year');
cont(cont>0) = NaN; clear cont

tera_lnd = squeeze(nanmean(nanmean(era_land_anual))); clear era_land_anual
lsq_era_lnd = fitlm(anos,tera_lnd);
yo_era_lnd = lsq_era_lnd.Coefficients.Estimate(1);
trend_era_lnd = lsq_era_lnd.Coefficients.Estimate(2); 
erro_era_lnd = lsq_era_lnd.Coefficients.SE(2);
line_era_lnd = yo_era_lnd+trend_era_lnd*anos; 
trend_era_lnd = (trend_era_lnd*10); erro_era_lnd = erro_era_lnd*10;
p_era_lnd = lsq_era_lnd.Coefficients.pValue(2); 

if p_era_lnd < 0.01
    SIG_era_lnd = ('p < 0.01');
elseif p_era_lnd < 0.05
    SIG_era_lnd = ('p < 0.05');
elseif p_era_lnd < 0.1
    SIG_era_lnd = ('p < 0.1');
else
    SIG_era_lnd = ('p > 0.1');
end
tt_era_lnd = (['Tend Continente_{ERA}: ', num2str(trend_era_lnd,'% .3e'), ...
    ' \pm ', num2str(erro_era_lnd,'% .3e'), '^oC/dec; ', SIG_era_lnd]);
clear SIG_era_lnd; clear p_era_lnd; clear yo_era_lnd;
clear trend_era_lnd; clear erro_era_lnd; clear lsq_era_lnd; 

anom_era_lnd = []; 
for i = 1:length(tanos)
    anom_era_lnd(i) = (tera_lnd(i) - mean(tera_lnd));
end
clear i

%% PLOT CONTvsOC
figure('color',[1 1 1],'position',[108 305 1400 700])
subtightplot(2,2,1,gap,marg_h,marg_w);  
plot(tanos,tnoaa_lnd,'color',[0.7451 0.7098 0.5608],'linewidth',0.8);
hold on; plot(tanos,tera_lnd,'color', [0.7235 0.8363 0.9565],'linewidth',0.8);
title('Continente');xlabel('Anos');ylabel('^oC');
datetick('x',10,'keepticks'); xlim([tanos(1) tanos(end)]); ylim([-31 -10])
plot(tanos, line_noaa_lnd, 'color', [0.6 0.5608 0.1647], 'linewidth',2) 
plot(tanos, line_era_lnd, 'color', [0 0.4471 0.7412], 'linewidth',2) 
legend('NOAA20CR', 'ERA20C')
text(tanos(2),-17, tt_noaa_lnd, 'color', [0.6000 0.5608 0.1647])
text(tanos(2),-20, tt_era_lnd, 'color', [0 0.4471 0.7412])

subtightplot(2,2,2,gap,marg_h,marg_w);
plot(tanos,anom_noaa_lnd,'color',[0.7451 0.7098 0.5608],'linewidth',1.5);
hold on; plot(tanos,anom_era_lnd,'color',[0.7235 0.8363 0.9565],'linewidth',1.5);
title('Anom_T Continente');xlabel('Anos');ylabel('^oC');
datetick('x',10,'keepticks');xlim([tempo(1) 734473]);ylim([-3 3])
legend('NOAA20CR', 'ERA20C')
plot([tanos(1) tanos(end)],[0 0],'Color',[0.85 0.85 0.85],'LineStyle','--')

subtightplot(2,2,3,gap,marg_h,marg_w);
plot(tanos,tnoaa_oc,'color',[0.7451 0.7098 0.5608],'linewidth',0.8);
hold on; plot(tanos,tera_oc,'color',[0.7235 0.8363 0.9565],'linewidth',0.8);
title('Oceano Austral');xlabel('Anos');ylabel('^oC');
datetick('x',10,'keepticks');xlim([tempo(1) 734473]); ylim([-0.6 0.7])
plot(tanos, line_noaa_oc, 'color',[0.6 0.5608 0.1647], 'linewidth',2)
plot(tanos, line_era_oc, 'color',[0 0.4471 0.7412], 'linewidth',2)
legend('NOAA20CR', 'ERA20C')
text(tanos(2),0.6, tt_noaa_oc, 'color', [0.6 0.5608 0.1647])
text(tanos(2),0.4, tt_era_oc, 'color', [0 0.4471 0.7412])

subtightplot(2,2,4,gap,marg_h,marg_w);
plot(tanos,anom_noaa_oc,'color',[0.7451 0.7098 0.5608],'linewidth',1.5);
hold on; plot(tanos,anom_era_oc,'color',[0.7235 0.8363 0.9565],'linewidth',1.5);
title('Anom_T Oc Austral');xlabel('Anos');ylabel('^oC');
datetick('x',10,'keepticks');xlim([tempo(1) 734473]);ylim([-3 3])
plot([tanos(1) tanos(end)],[0 0],'Color',[0.85 0.85 0.85],'LineStyle','--')
legend('NOAA20CR', 'ERA20C')

clear anom_noaa_oc; clear anom_era_oc; clear line_noaa_oc; clear line_era_oc;
clear anom_noaa_lnd; clear anom_era_lnd; clear line_noaa_lnd; clear line_era_lnd;
clear tnoaa_lnd; clear tera_lnd; clear tnoaa_oc; clear tera_oc
clear noaa_land; clear noaa_sst; clear tt_noaa_lnd; clear tt_era_lnd
clear tt_noaa_oc; clear tt_era_oc


%% EASTE vs WESTE

% NOAA 
%*************************** EASTE *****************************
leste = lon>44 & lon<136;
noaa_leste_anual = downsample_ts(temp1(leste,austral,:),tempo,'year');

tnoaa_leste = squeeze(nanmean(nanmean(noaa_leste_anual))); clear noaa_leste_anual
lsq_noaa_leste = fitlm(anos,tnoaa_leste);
yo_noaa_leste = lsq_noaa_leste.Coefficients.Estimate(1);
trend_noaa_leste = lsq_noaa_leste.Coefficients.Estimate(2); 
erro_noaa_leste = lsq_noaa_leste.Coefficients.SE(2);
line_noaa_leste = yo_noaa_leste+trend_noaa_leste*anos; 
trend_noaa_leste = (trend_noaa_leste*10); erro_noaa_leste = erro_noaa_leste*10;
p_noaa_leste = lsq_noaa_leste.Coefficients.pValue(2); 

if p_noaa_leste < 0.01
    SIG_noaa_leste = ('p < 0.01');
elseif p_noaa_leste < 0.05
    SIG_noaa_leste = ('p < 0.05');
elseif p_noaa_leste < 0.1
    SIG_noaa_leste = ('p < 0.1');
else
    SIG_noaa_leste = ('p > 0.1');
end
tt_noaa_leste = (['Tend Leste_{NOAA}: ', num2str(trend_noaa_leste,'% .3e'), ...
    ' \pm ', num2str(erro_noaa_leste,'% .3e'), '^oC/dec; ', SIG_noaa_leste]);
clear SIG_noaa_leste; clear p_noaa_leste; clear yo_noaa_leste;
clear trend_noaa_leste; clear erro_noaa_leste; clear lsq_noaa_leste; 

anom_noaa_leste = []; 
for i = 1:length(tanos)
    anom_noaa_leste(i) = (tnoaa_leste(i) - mean(tnoaa_leste));
end
clear i

%************************ PENINSULA ***************************
pen = lon>-91 & lon<-29; clear lon
noaa_pen_anual = downsample_ts(temp1(pen,austral,:),tempo,'year');
clear temp1

tnoaa_pa = squeeze(nanmean(nanmean(noaa_pen_anual))); clear noaa_pen_anual
lsq_noaa_pa = fitlm(anos,tnoaa_pa);
yo_noaa_pa = lsq_noaa_pa.Coefficients.Estimate(1);
trend_noaa_pa = lsq_noaa_pa.Coefficients.Estimate(2); 
erro_noaa_pa = lsq_noaa_pa.Coefficients.SE(2);
line_noaa_pa = yo_noaa_pa+trend_noaa_pa*anos; 
trend_noaa_pa = (trend_noaa_pa*10); erro_noaa_pa = erro_noaa_pa*10;
p_noaa_pa = lsq_noaa_pa.Coefficients.pValue(2); 

if p_noaa_pa < 0.01
    SIG_noaa_pa = ('p < 0.01');
elseif p_noaa_pa < 0.05
    SIG_noaa_pa = ('p < 0.05');
elseif p_noaa_pa < 0.1
    SIG_noaa_pa = ('p < 0.1');
else
    SIG_noaa_pa = ('p > 0.1');
end
tt_noaa_pa = (['Tend Peninsula_{NOAA}: ', num2str(trend_noaa_pa,'% .3e'), ...
    ' \pm ', num2str(erro_noaa_pa,'% .3e'), '^oC/dec; ', SIG_noaa_pa]);
clear SIG_noaa_pa; clear p_noaa_pa; clear yo_noaa_pa;
clear trend_noaa_pa; clear erro_noaa_pa; clear lsq_noaa_pa; 

anom_noaa_lnd = []; 
for i = 1:length(tanos)
    anom_noaa_lnd(i) = (tnoaa_pa(i) - mean(tnoaa_pa));
end
clear i

%%% ERA 
%*************************** EASTE **********************************
era_leste_anual = downsample_ts(temp2(leste,austral, :),tempo,'year');

tera_leste = squeeze(nanmean(nanmean(era_leste_anual))); clear era_leste_anual
lsq_era_leste = fitlm(anos,tera_leste);
yo_era_leste = lsq_era_leste.Coefficients.Estimate(1);
trend_era_leste = lsq_era_leste.Coefficients.Estimate(2); 
erro_era_leste = lsq_era_leste.Coefficients.SE(2);
line_era_leste = yo_era_leste+trend_era_leste*anos; 
trend_era_leste = (trend_era_leste*10); erro_era_leste = erro_era_leste*10;
p_era_leste = lsq_era_leste.Coefficients.pValue(2); 

if p_era_leste < 0.01
    SIG_era_leste = ('p < 0.01');
elseif p_era_leste < 0.05
    SIG_era_leste = ('p < 0.05');
elseif p_era_leste < 0.1
    SIG_era_leste = ('p < 0.1');
else
    SIG_era_leste = ('p > 0.1');
end
tt_era_leste = (['Tend Leste_{ERA}: ', num2str(trend_era_leste,'% .3e'), ...
    ' \pm ', num2str(erro_era_leste,'% .3e'), '^oC/dec; ', SIG_era_leste]);
clear SIG_era_leste; clear p_era_leste; clear yo_era_leste;
clear trend_era_leste; clear erro_era_leste; clear lsq_era_leste; 

anom_era_leste = []; 
for i = 1:length(tanos)
    anom_era_leste(i) = (tera_leste(i) - mean(tera_leste));
end
clear i

%************************* PENINSULA ****************************
era_pen_anual = downsample_ts(temp2(pen,austral,:), tempo,'year');
clear temp2

tera_pa = squeeze(nanmean(nanmean(era_pen_anual))); clear era_pen_anual
lsq_era_pa = fitlm(anos,tera_pa);
yo_era_pa = lsq_era_pa.Coefficients.Estimate(1);
trend_era_pa = lsq_era_pa.Coefficients.Estimate(2); 
erro_era_pa = lsq_era_pa.Coefficients.SE(2);
line_era_pa = yo_era_pa+trend_era_pa*anos; 
trend_era_pa = (trend_era_pa*10); erro_era_pa = erro_era_pa*10;
p_era_pa = lsq_era_pa.Coefficients.pValue(2); 

if p_era_pa < 0.01
    SIG_era_pa = ('p < 0.01');
elseif p_era_pa < 0.05
    SIG_era_pa = ('p < 0.05');
elseif p_era_pa < 0.1
    SIG_era_pa = ('p < 0.1');
else
    SIG_era_pa = ('p > 0.1');
end
tt_era_pa = (['Tend Peninsula_{ERA}: ', num2str(trend_era_pa,'% .3e'), ...
    ' \pm ', num2str(erro_era_pa,'% .3e'), '^oC/dec; ', SIG_era_pa]);
clear SIG_era_pa; clear p_era_pa; clear yo_era_pa;
clear trend_era_pa; clear erro_era_pa; clear lsq_era_pa;

anom_era_pa = []; 
for i = 1:length(tanos)
    anom_era_pa(i) = (tera_pa(i) - mean(tera_pa));
end
clear i

%% PLOT LESTE/WEDDELL
figure('color',[1 1 1],'position',[108 305 1400 700])
subtightplot(2,2,1,gap,marg_h,marg_w);
plot(tanos,tnoaa_pa,'color',[1 0.7608 0.9647],'linewidth',0.8);
hold on; plot(tanos,tera_pa,'color', [1 0.8784 0.4314],'linewidth',0.8);
title('Peninsula Antartica');xlabel('Anos');ylabel('^oC');
datetick('x',10,'keepticks'); xlim([tempo(1) 734473]); ylim([-15 -7])
plot(tanos, line_noaa_pa, 'color', [1 0 0.6], 'linewidth',2) 
plot(tanos, line_era_pa, 'color', [0.9804 0.7569 0.0706], 'linewidth',2) 
legend('NOAA20CR', 'ERA20C')
text(tanos(2),-7.5, tt_noaa_pa, 'color', [1 0 0.6])
text(tanos(2),-8.5, tt_era_pa, 'color', [0.9804 0.7569 0.0706])

subtightplot(2,2,2,gap,marg_h,marg_w);
plot(tanos,anom_noaa_lnd,'color',[1 0.7608 0.9647],'linewidth',1.5);
hold on; plot(tanos,anom_era_pa,'color',[1 0.8784 0.4314],'linewidth',1.5);
title('Anom_T Peninsula Antartica');xlabel('Anos');ylabel('^oC');
datetick('x',10,'keepticks');xlim([tempo(1) 734473]);ylim([-3 3])
legend('NOAA20CR', 'ERA20C')
plot([tanos(1) tanos(end)],[0 0],'Color',[0.85 0.85 0.85],'LineStyle','--')

subtightplot(2,2,3,gap,marg_h,marg_w);
plot(tanos,tnoaa_leste,'color',[1 0.7608 0.9647],'linewidth',0.8);
hold on; plot(tanos,tera_leste,'color',[1 0.8784 0.4314],'linewidth',0.8);
title('Antartica Leste'); xlabel('Anos');ylabel('^oC');
datetick('x',10,'keepticks');xlim([tempo(1) 734473]); ylim([-28 -10])
plot(tanos, line_noaa_leste, 'color',[1 0 0.6], 'linewidth',2)
plot(tanos, line_era_leste, 'color',[0.9804 0.7569 0.0706], 'linewidth',2)
legend('NOAA20CR', 'ERA20C')
text(tanos(2),-16, tt_noaa_leste, 'color', [1 0 0.6])
text(tanos(2),-18, tt_era_leste, 'color', [0.9804 0.7569 0.0706])

subtightplot(2,2,4,gap,marg_h,marg_w);
plot(tanos,anom_noaa_leste,'color',[1 0.7608 0.9647],'linewidth',1.5);
hold on; plot(tanos,anom_era_leste,'color',[1 0.8784 0.4314],'linewidth',1.5);
title('Anom_T Leste');xlabel('Anos');ylabel('^oC');
datetick('x',10,'keepticks');xlim([tempo(1) 734473]);ylim([-3 3])
plot([tanos(1) tanos(end)],[0 0],'Color',[0.85 0.85 0.85],'LineStyle','--')
legend('NOAA20CR', 'ERA20C')

clear anom_noaa_leste; clear anom_era_leste; clear line_noaa_leste; clear line_era_leste;
clear anom_noaa_pa; clear anom_era_pa; clear line_noaa_pa; clear line_era_pa;
clear tnoaa_pa; clear tera_pa; clear tnoaa_leste; clear tera_leste
clear tanos; clear pen; clear leste; clear austral; clear tempo; clear anos
clear tt_noaa_leste; clear tt_noaa_pa; clear tt_era_leste; clear tt_era_pa
clear gap; clear marg_h; clear marg_w



