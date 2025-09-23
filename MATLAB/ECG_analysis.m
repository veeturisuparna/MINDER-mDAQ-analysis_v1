clc; clear; close all; warning('off','all')
path = 'G:\My Drive\BSN internal testing\ECG+EDA testing\Biopac\4-11-24_electrode_study_ECG';
MatFiles = dir([path,'\**\*.mat']);
Names = {MatFiles.name};
Folders = {MatFiles.folder};
for n = 1:length(Names)
    % for n = 1:1
    fontsize = 40;
    name = Names(n);
    name = name{1}(1:end-4);
    folder = Folders(n);

    datafile = load(fullfile(Folders{n},Names{n}));
    srate = 125;
    drylead = normalize(datafile.data(:,2));
    drylead = drylead - mean(drylead);
    wetlead = normalize(datafile.data(:,1));
    wetlead = wetlead - mean(wetlead);

    dryhilbert = HilbertTransform(drylead);
    wethilbert = HilbertTransform(wetlead);

    [dryheartbeat,dryqrs,drytime,drylocs] = IdentifyQRS(dryhilbert,srate);
    [wetheartbeat,wetqrs,wettime,wetlocs] = IdentifyQRS(wethilbert,srate);

    dryskewness = QRSSkewness(dryqrs(:,13:39));
    wetskewness = QRSSkewness(wetqrs(:,13:39));

    drykurtosis = QRSKurtosis(dryqrs(:,13:39));
    wetkurtosis = QRSKurtosis(wetqrs(:,13:39));

    drypsd = QRSPSD(dryqrs(:,13:39),srate);
    wetpsd = QRSPSD(wetqrs(:,13:39),srate);


    figuretxt = ['Subject ',num2str(n),' MEAMA Data'];
    [F,T] = InitializeFigure(figuretxt,2,2,fontsize);

    t = tiledlayout(T,2,1,'TileSpacing','tight');
    t.Layout.Tile = 1;
    nexttile(t)
    PlotData(real(dryhilbert),srate,drylocs)
    title('Dry Lead')
    nexttile(t)
    PlotData(real(wethilbert),srate,wetlocs)
    title('Wet Lead')

    t = tiledlayout(T,2,1,'TileSpacing','tight');
    t.Layout.Tile = 2;
    nexttile(t)
    PlotHeartbeat(dryheartbeat,dryqrs,drytime,[-4 8],25,'Dry')
    nexttile(t)
    PlotHeartbeat(wetheartbeat,wetqrs,wettime,[-4 8],25,'Wet')

    t = tiledlayout(T,1,3,'TileSpacing','tight');
    t.Layout.Tile = 3;
    nexttile(t)
    BP(dryskewness,wetskewness,25,'Skewness',[-3 3]);
    yline(2,'k','Normal QRS Skewness','LabelHorizontalAlignment','center','linewidth',5)
    yline(-2,'k','Normal QRS Skewness','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center','linewidth',5)
    nexttile(t)
    BP(drykurtosis,wetkurtosis,25,'Kurtosis',[0 10]);
    yline(5,'k','Normal QRS Kurtosis','LabelHorizontalAlignment','center','linewidth',5)
    nexttile(t)
    BP(drypsd,wetpsd,25,'5-15hz PSD',[0 1]);
    yline(.5,'k','Normal QRS Kurtosis','LabelHorizontalAlignment','center','linewidth',5)


    t = tiledlayout(T,1,1,'tilespacing','tight');
    t.Layout.Tile = 4;
    nexttile(t)
    ShowHistogram(dryqrs(:,13:39),wetqrs(:,13:39),fontsize)

    FinalizeFigure(F,T,figuretxt,path,fontsize);

    dryresultsskew{n} = dryskewness;
    dryresultskur{n} = drykurtosis;
    dryresultspsd{n} = drypsd;

    wetresultsskew{n} = wetskewness;
    wetresultskur{n} = wetkurtosis;
    wetresultspsd{n} = wetpsd;
end
%%
figuretxt = ['Sorted SQI'];
[F,T] = InitializeFigure(figuretxt,2,3,fontsize/2);
nexttile(T)
SortSQI(dryresultsskew,'Dry Skewness',fontsize,[-3 3]);
yline(2,'k','Normal QRS Skewness','LabelHorizontalAlignment','center','linewidth',5)
yline(-2,'k','Normal QRS Skewness','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center','linewidth',5)

nexttile(T)
SortSQI(dryresultskur,'Dry Kurtosis',fontsize,[0 10]);
yline(5,'k','Normal QRS Kurtosis','LabelHorizontalAlignment','center','linewidth',5)

nexttile(T)
SortSQI(dryresultspsd,'Dry PSD',fontsize,[0 1]);
yline(.5,'k','Normal QRS Kurtosis','LabelHorizontalAlignment','center','linewidth',5)

nexttile(T)
SortSQI(wetresultsskew,'Wet Skewness',fontsize,[-3 3]);
yline(2,'k','Normal QRS Skewness','LabelHorizontalAlignment','center','linewidth',5)
yline(-2,'k','Normal QRS Skewness','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center','linewidth',5)

nexttile(T)
SortSQI(wetresultskur,'Wet Kurtosis',fontsize,[0 10]);
yline(5,'k','Normal QRS Kurtosis','LabelHorizontalAlignment','center','linewidth',5)

nexttile(T)
SortSQI(wetresultspsd,'Wet PSD',fontsize,[0 1]);
yline(.5,'k','Normal QRS Kurtosis','LabelHorizontalAlignment','center','linewidth',5)

FinalizeFigure(F,T,figuretxt,path,fontsize)

function SortSQI(sqi,figuretxt,fontsize,yrange)
collected_area = sqi;
for i = 1:size(collected_area,2)
    medianavg(i) = median(collected_area{i});
end
[a,n] = sort(medianavg,'descend');
% n = 1:length(n);
clear display
clear label
display_data = 0;
label = 0;
Names = {'S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11','S12','S13','S14','S15','S16','S17','S18'};
for m = n
    dataTextile = collected_area{m}';
    display_data = [display_data;dataTextile];
    name = Names{m};
    if length(name)>20
        name = name(1:20);
    else
        name = name;
    end
    txt = {name};
    label = [label;repmat(txt,size(dataTextile))];
end


display_data(1)=[];
label(1)=[];

bp = boxplot(display_data,label,'notch','on','Symbol','');
for count = 1:size(bp,2)
    set(bp(:,count),'linewidth',5);
end
yline(median(medianavg),'r','Median Value','linewidth',5)
ylim(yrange)
% ylim([0.35 .55])
% ylabel('Skewness','Fontsize',fontsize)
% FinalizeFigure(F,T,figuretxt,path,fontsize)
title(figuretxt,'fontsize',fontsize)
end

function RelPSD = QRSPSD(qrs,srate)
PSD_5to15 = sum(pwelch(real(qrs),[],[],5:15,srate).^2);
PSD_0to45 = sum(pwelch(real(qrs),[],[],0:45,srate).^2);
RelPSD = PSD_5to15./PSD_0to45;
end

function ShowHistogram(data1,data2,fontsize)
c = 100;
n = 4;
Period = reshape(data1,1,[]);
[N] = histcounts((real(Period)),linspace(-n,n ,c),'Normalization','pdf');
t = linspace(-n,n,c-1);
normal = 1/sqrt(2*pi) * exp(-(t.^2)./2);
plot(t,normal,'k','linewidth',5,'DisplayName','Normal Distribution')
hold on
plot(t,N,'linewidth',5,'DisplayName','Dry Lead')

Period = reshape(data2,1,[]);
[N] = histcounts((real(Period)),linspace(-n,n ,c),'Normalization','pdf');

plot(t,N,'linewidth',5,'DisplayName','Wet Lead')
legend
xlim tight
hold off
txt = {['Histogram of QRS Region']};
title(txt,'Fontsize',fontsize)
ylim([0 1])
end

function BP(data1SQI,data2SQI,fontsize,type,yrange)
area_int = [data1SQI,data2SQI];
label = [repmat({'Dry'},size(data1SQI)),repmat({'Wet'},size(data2SQI))];
bp = boxplot(area_int,label,'notch','on','Symbol','');
for count = 1:size(bp,2)
    set(bp(:,count),'linewidth',5);
end
ylim(yrange)
title(type,'Fontsize',fontsize)
% ylabel('Area of Joint PDF','Fontsize',fontsize)
end

function [skewresult] = QRSSkewness(qrs)
skewresult = skewness(real(qrs)');
end

function [kurresult] = QRSKurtosis(qrs)
kurresult = kurtosis(real(qrs)');
end

function PlotHeartbeat(heartbeat,qrs,time,yrange,fontsize,type)
plot(time,real(qrs'),'color',[0 0 0 2/size(qrs,1)],'linewidth',2)
hold on
plot(time,real(heartbeat),'color',[1 0 0 1],'linewidth',2)
xlim tight
ylim(yrange)
ylabel('Amplitude (uV)','Fontsize',fontsize)
xlabel('Time (ms)','Fontsize',fontsize)
txt = [type,' ',num2str(size(qrs,1)),' Individual Heartbeats'];
title(txt,'Fontsize',fontsize)
xlabel([])
ylabel([])
xticks([])
end

function [heartbeat,qrs,time,locs] = IdentifyQRS(data,srate)
minwidth = 0.01;
maxwidth = 0.1;
minheight = 0.5;
minprom = 0.5;
mindist = 0.5;
pre_ms = -0.2*srate;
post_ms = 0.8*srate;
time = [pre_ms:1:post_ms]/srate;
[~,locs,w,p] = findpeaks(abs(normalize(data)),'MinPeakWidth',minwidth*srate,'MaxPeakWidth',maxwidth*srate,'MinPeakProminence',minprom,'MinPeakDistance',mindist*srate);

epoch_begin = round(locs+pre_ms);
epoch_end = round(locs+post_ms);
epoch_end = epoch_end(epoch_begin>0);
epoch_begin = epoch_begin(epoch_begin>0);

epoch_begin = epoch_begin(epoch_end<size(data,1));
epoch_end = epoch_end(epoch_end<size(data,1));
qrs = zeros(length(epoch_begin),srate+1);
for j = 1:length(epoch_begin)
    epoch = data(epoch_begin(j):epoch_end(j));
    epoch = epoch-mean(epoch);
    qrs(j,:) = normalize(epoch);
end
heartbeat = mean(qrs);
end

function [data_hilbert] = HilbertTransform(data)
data_fft = fft(data);
data_fft_complex = 1i*data_fft;
freq_pos = 2:floor(size(data_fft,1)/2)+mod(size(data_fft,1),2);
freq_neg = ceil(size(data_fft,1)/2)+1+~mod(size(data_fft,1),2):size(data_fft,1);

data_fft(freq_pos) = data_fft(freq_pos) + -1i*data_fft_complex(freq_pos);
data_fft(freq_neg) = data_fft(freq_neg) + 1i*data_fft_complex(freq_neg);

data_hilbert = ifft(data_fft);
end

function PlotData(data,srate,loc)
time = linspace(0,length(data)/srate,length(data));
plot(time,data)
hold on
plot(time(loc),data(loc),'.','markersize',10)
xlim tight
ylim([-4 8])
end

function [F,T] = InitializeFigure(figuretxt,a,b,fontsize)
F = figure('Name',figuretxt,'NumberTitle','off');
set(F,'WindowState','Fullscreen','Color','White','DefaultAxesFontsize',25)
T=tiledlayout(F,a,b,'TileSpacing','loose');
end

function FinalizeFigure(F,T,figuretxt,path,fontsize)
title(T,figuretxt,'fontsize',fontsize)
export_fig(F,figuretxt)
close
end