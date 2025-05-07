function ephysBrowser

clear
display('Enter data file')

%%%%%%%%%IMPORT%%%%%%%%%%%%%
[fname, pname] = uigetfile({'*.atf'}, 'Select trace file to open');

pathFile = [pname fname];
a = importdata(pathFile, '\t', 11 );
b = a.data(:,2:end);
[r c] = size(b);
numTraces = c;
sr=10000;  %Sampling rate in Hz

%%%%%%%FILTER%%%%%%%%%%%%%%%%%
% for i =1:numTraces
% int = (1:1:r)/sr;
% p = b(:,i);
% d = timeseries(p, int');
% q = idealfilter(d, [10, 10000], 'pass');
% filtData(:,i) = q.Data;
% end
% b = filtData;
t = (1:r)/sr;
figure; plot(t,b); title([fname(1:end-4) ': raw'])  %plot figs for inspecting data


%%%%baseline subtract%%
base1= .001; base2= 5;                               %%%%%%%%%%%USER%%%%%%%%%%%%%
base = mean(b(base1*sr:base2*sr,:));
baseSub = repmat(base, r,1);
bs = b-baseSub;

figure; plot(t,bs);  title([fname(1:end-4) ': baseline subtracted']) 


%find peak and steady state minima values
ind1p=5 ; ind2p= 20  ;  % peak indices in seconds         %%%%%%%%%%%USER%%%%%%%%%%%%%
ind1ss= 15 ; ind2ss= 20 ;  % ss indices in seconds         %%%%%%%%%%%USER%%%%%%%%%%%%%

for i=1:numTraces
    %first, the peak min
    [bMin,I] = min(bs(ind1p*sr:ind2ss*sr, i));
    medMin = median(bs((I-5)+ind1p*sr:I+ind1p*sr,i));
    minB(i)=medMin;   %%% minB contains the median value w/in 1 msec of true min
    %then, the steady state min
    [ssMin,Iss] = min(bs(ind1ss*sr:ind2ss*sr, i));
     medMin = median(bs((Iss-5)+ind1ss*sr:Iss+ind1ss*sr,i));
    minBss(i)=medMin;   %%% minB contains the median value w/in 50 msec of end of pulse
end

% drugApp = input('time frame of drug application in sec:  ');
% drugApp = sr*[4.125 24.125];                              %%%%%%%%%%%USER%%%%%%%%%%%%%
drugApp = [5 20];                              %%%%%%%%%%%USER%%%%%%%%%%%%%
lengthApp = size(drugApp);
yPos = 0.1*max(abs(minB));

%%%MAKE CONC vector
highConc = 30; %highest drug conc               %%%%%%%%%%%USER%%%%%%%%%%%%%
conc=highConc; 
for i=1:numTraces-1
    conc(i+1) = conc(i)/3;
end
conc=fliplr(conc);

for i=1:numTraces
    if conc(i)>.999
       M{i} = [num2str(round(conc(i))) ' mM'];
    elseif conc(i)>.0999
        z = num2str(conc(i));
        M{i} = [z(1:3) ' mM'];
    else
        z = num2str(conc(i));
        M{i} = [z(1:4) ' mM'];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

response = input('Additional title information:  ');

%%%%%%%%%%PLOTTING%%%%%%%%%%%%%%%%

hfig=figure('Position', [10,10 900, 1200]);


%%ephys traces%%%
hsub1= subplot(3,2,1, 'replace');   %raw traces
    plot(t,b);
%     set(gca, 'XLim', [0, r])
    set(gca, 'YLim', [1.1*(min(minB)+min(base)), 0.2*abs(min(minB))])
    set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'LineWidth', 2)
%     xTicks = get(gca, 'XTick'); set(gca, 'XTickLabel', xTicks/sr);
    ylabel('pA', 'FontSize', 10, 'FontWeight', 'bold')
    xlabel('sec', 'FontSize', 10, 'FontWeight', 'bold')
    box off
    hold on
    hline = plot(drugApp ,repmat(yPos, lengthApp), '-k');
    set(hline, 'LineWidth', 2); box off
    hleg = legend(hsub1, M, 'Location', 'SouthEast'); set(hleg, 'FontSize', 7)
    title([fname(1:end-4) ':' response]);

 hsub2= subplot(3,2,2, 'replace');    %baseline subtracted traces
    plot(t,bs);
%     set(gca, 'XLim', [0, r])
    set(gca, 'YLim', [1.1*min(minB), 0.2*abs(min(minB))])
    set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'LineWidth', 2)
%     xTicks = get(gca, 'XTick'); set(gca, 'XTickLabel', xTicks/sr);
    ylabel('pA', 'FontSize', 10, 'FontWeight', 'bold')
    xlabel('sec', 'FontSize', 10, 'FontWeight', 'bold')
    hold on
    hline = plot(drugApp ,repmat(yPos, lengthApp), '-k');
    set(hline, 'LineWidth', 2); box off
    hleg = legend(hsub2, M, 'Location', 'SouthEast'); set(hleg, 'FontSize', 7)
    title([fname(1:end-4) ':' response ' BS']);

  %%%DR plots%%%%  
    
Z = abs(minB);
Zss= abs(minBss);

hsub3= subplot(3,2,3, 'replace');         %DR peak
    plot(log10(conc(1:end)),Z(1:end), 'ok')
    box off
    set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'LineWidth', 2)
    ylabel('pA', 'FontSize', 10, 'FontWeight', 'bold')
    xlabel('conc', 'FontSize', 10, 'FontWeight', 'bold')
    xLim1 = log10(str2num(M{1}(1:end-3))); xLim2 = log10(str2num(M{end}(1:end-3)));
    xLimits = xLim1:0.5:xLim2;
    set(gca, 'XTick', xLimits)
    for i =1:length(M)
        Mnum{i} =M{i}(1:end-3);
    end
    set(gca, 'XTickLabel',Mnum);
    results = ec50(conc', Z');
    h= legend(['EC50 = ' num2str(results(3)) ' mM']);
    set(h, 'Location', 'NorthWest', 'FontSize', 7)
    title([fname(1:end-4) ':' response ' peakDR']);
     
hsub4= subplot(3,2,4, 'replace');    %DR peak fitted
try
    plotSigmoid(results, conc);
    hold on
    plot(conc,Z, 'ok', 'MarkerFaceColor', 'k','MarkerSize', 9)
    box off
    set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'LineWidth', 1)
    set(gca, 'XScale', 'log', 'XLim', [0.01 10], 'XTickLabel',Mnum)

    xLim1 = log10(str2num(M{1}(1:end-3))); xLim2 = log10(str2num(M{end}(1:end-3)));
    xLimits = xLim1:0.5:xLim2;
    set(gca, 'XTick', 10.^xLimits)
    for i =1:length(M)
        Mnum{i} =M{i}(1:end-3);
    end
    set(gca, 'XTickLabel',Mnum);

    ylabel('pA', 'FontSize', 10, 'FontWeight', 'bold')
    xlabel('conc', 'FontSize', 10, 'FontWeight', 'bold')
    
    h= legend(['EC50 = ' num2str(results(3)) ' mM']);
    set(h, 'Location', 'NorthWest', 'FontSize', 7)
    title([fname(1:end-4) ':' response ' peakDR']);    
catch 
end

hsub5= subplot(3,2,5, 'replace');   %DR steady state
    plot(log10(conc(1:end)),Zss(1:end), 'ok')
    box off
    set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'LineWidth', 2)
    ylabel('pA', 'FontSize', 10, 'FontWeight', 'bold')
    xlabel('conc', 'FontSize', 10, 'FontWeight', 'bold')
    xLim1 = log10(str2num(M{1}(1:end-3))); xLim2 = log10(str2num(M{end}(1:end-3)));
    xLimits = xLim1:0.5:xLim2;
    set(gca, 'XTick', xLimits)
    for i =1:length(M)
        Mnum{i} =M{i}(1:end-3);
    end
    set(gca, 'XTickLabel',Mnum);
    resultsSS = ec50(conc', Zss');
    h= legend(['EC50 = ' num2str(resultsSS(3)) ' mM']);
    set(h, 'Location', 'NorthWest', 'FontSize', 7)
    title([fname(1:end-4) ':' response ' ssDR']);
     
hsub6= subplot(3,2,6, 'replace'); %DR steady state fitted
try
    plotSigmoid(resultsSS, conc);
    hold on
    plot(conc,Zss, 'ok', 'MarkerFaceColor', 'k','MarkerSize', 9)
    box off
    set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'LineWidth', 1)
    set(gca, 'XScale', 'log', 'XLim', [0.01 10], 'XTickLabel',Mnum)

    xLim1 = log10(str2num(M{1}(1:end-3))); xLim2 = log10(str2num(M{end}(1:end-3)));
    xLimits = xLim1:0.5:xLim2;
    set(gca, 'XTick', 10.^xLimits)
    for i =1:length(M)
        Mnum{i} =M{i}(1:end-3);
    end
    set(gca, 'XTickLabel',Mnum);

    ylabel('pA', 'FontSize', 10, 'FontWeight', 'bold')
    xlabel('conc', 'FontSize', 10, 'FontWeight', 'bold')
    
    h= legend(['EC50 = ' num2str(resultsSS(3)) ' mM']);
    set(h, 'Location', 'NorthWest', 'FontSize', 7)
    title([fname(1:end-4) ':' response ' ssDR']);    
catch 
end

gcf
hAx = axes('Position', [0 0 1 1], 'Visible', 'off');
set(gcf, 'CurrentAxes', hAx)

str = {['base1= ' num2str(base1) ',  base2= ' num2str(base2) ',  ind1p= ' num2str(ind1p) ',  ind2p= ' ...
    num2str(ind2p) ',  ind1ss= ' num2str(ind1ss) ',  ind2ss= ]' num2str(ind2ss) ']' ',  drugApp= ' ...
    num2str(drugApp/sr) ',  highConc= ' num2str(highConc)]};
text(.1,.06,str, 'FontSize', 9)

    orient tall   %forces figure to fill the page
    
    saveas(hfig, [num2str(fname(1:end-4)) '.jpg'])
    saveas(hfig, [num2str(fname(1:end-4)) '.m'])
    
    format short g

    T = horzcat(conc', round(Z'), round(Zss'));
    U = horzcat(results', resultsSS');
    openvar T
    openvar U%(1:4,:)
    
    
    



