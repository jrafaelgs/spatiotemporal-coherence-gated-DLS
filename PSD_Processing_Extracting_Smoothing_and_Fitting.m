% Extraction, Smoothing & Fitting of experimental LC-DLS spectra
% Function calls:
%      extractPSD - Extract and average PSD's from a file.
%      Sectioned_Smoothing - Smoothing is performed in 2 sections.
%      Lorentzian_Fit - Fit PSDn to Several Lorentzians
%      Lorentzian_Reconstruction - Reconstruct PSDn based on Fit Params
%      Exponential_Reconstruction - Reconstruct Sum of decaying exponentials based on Fit Params
% Rafael Guzman, 2015

%% Global variables
clear all
close all
clc
format long

%% Physical parameters [SI]
kb = 1.3806488e-23;
lambda = 0.680e-6;
refracInd = 1.37; % Refractive index of suspending medium
theta = pi;
particleRadius = 2.78e-6;   % RADIUS!!!
T = 298; % Temperature in Kelvin
q = (4*pi*refracInd*sin(theta/2)/(lambda));    % Magnitude of the Scattering Vector [um^-1]
% visc = 0.902176e-3;     % [Pa-s]
visc = (2.26509e-4) + 6.41005*exp(-0.030715*T);
% D0_Th = (1e12)*kb*T/(6*pi*visc*particleRadius);    % [um^2/s]

%% Configuration Parameters of Acquisition
xFreq = (1:1:1e4)';     % Frequencies that the PSD are plotted against
sampPerFile = 1500;

%% Filename for the PSD file
InputFileFolder1 = 'H:\Home\Random\Shared\CURRENT PROJECTS\LC-DLS\2014\Orlando Health - Arnold Palmer Hospital\Data and Analysis\';
InputFileFolder3 = 'Measurements\';
CaseNo = 10;
MeasNo = 5;

%%
switch CaseNo
    case 1
        InputFileFolder2 = '2014-06-23 - OH Case 1\';
        InputFileDate = '23-Jun-2014';
    case 2
        InputFileFolder2 = '2014-09-25 - OH Case 2\';
        InputFileDate = '25-Sep-2014';
    case 3
        InputFileFolder2 = '2014-10-02 - OH Case 3\';
        InputFileDate = '02-Oct-2014';
    case 4
        InputFileFolder2 = '2014-12-03 - OH Case 4\';
        InputFileDate = '03-Dec-2014';
    case 5
        InputFileFolder2 = '2015-01-07 - OH Case 5\';
        InputFileDate = '07-Jan-2015';
    case 6
        InputFileFolder2 = '2015-02-11 - OH Case 6\';
        InputFileDate = '11-Feb-2015';
    case 7
        InputFileFolder2 = '2015-02-16 - OH Case 7\';
        InputFileDate = '16-Feb-2015';
    case 8
        InputFileFolder2 = '2015-11-10 - OH Case 8\';
        InputFileDate = '10-Nov-2015';
    case 9
        InputFileFolder2 = '2015-11-18 - OH Case 9\';
        InputFileDate = '18-Nov-2015';
    case 10
        InputFileFolder2 = '2016-03-09 - OH Case 10\';
        InputFileDate = '09-Mar-2016';
end
switch MeasNo
    case 1
        InputFileName = 'Meas1_Baseline_';
    case 2
        InputFileName = 'Meas2_Heparin_';
    case 3
        InputFileName = 'Meas3_CPB_';
    case 4
        InputFileName = 'Meas4_Post_CPB_';
    case 5
        InputFileName = 'Meas5_Post_Reversal_';
end
%%
InputFileFolder = [InputFileFolder1 InputFileFolder2 InputFileFolder3];
OutputFileName = ['Case_' num2str(CaseNo) '_' InputFileName 'Processed.xls'];
Sheet_Raw = 'Raw';
Sheet_Summary = 'Summary';
WriteOutputFile = 0;

%% Definition of Smoothing parameters
lowFreqCutoff_Sm = min(xFreq);
highFreqCutoff_Sm = 1e4;
Breaking_freq = 1e3;   % Frequency to separate the smoothing...
Span_a = 0.01;
Span_b = 0.01;
SM_Method_a = 'lowess';
SM_Method_b = 'rlowess';

%% Definition of Fitting Parameters
lowFreqCutoff = 5e0;
highFreqCutoff = 1e4;
lorent_options = 5;
random_fit = 0;
Ntimes = 1;

% x0 = [0.25 0.25 0.25 0.25 0.5 0.02 0.002 0.0002];
% x0 = [0.6 0.2 0.2 0.05 0.001 0.0001];
% x0 = [0.6 0.2 0.2 0.01 0.0005 0.00001];
% x0 = [0.6 0.4 0.02 0.0002];
% x0 = [1 0.005];

% x0 = [0.6 0.2 0.2 10 0.01 0.0001];
% x0 = [0.7 0.1 0.1 0.1 5 0.02 0.001 0.0001];

% x0 = [0.2 0.2 0.2 0.2 0.2 0.5 0.1 0.02 0.001 0.0001];
x0 = [0.2 0.2 0.2 0.2 0.2 5 0.5 0.02 0.001 0.0001];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extraction and Normalization of the PSD
% [psdOut_raw, psdAvg_raw] = extractPSD([InputFileFolder InputFileName InputFileDate],xFreq,sampPerFile);
% [psdOut_raw0, psdAvg_raw0] = extractPSD([InputFileFolder InputFileName InputFileDate],xFreq,sampPerFile);
% [psdOut_raw0_bg, psdAvg_raw0_bg] = extractPSD([InputFileFolder InputFileNameBg InputFileDate],xFreq,sampPerFile);
[psdOut_raw, psdAvg_raw] = extractPSD([InputFileFolder InputFileName InputFileDate],xFreq,sampPerFile);

%Extract extra PSDs
% psdOut_raw = psdOut_raw(:,[1:7,9,11,13:22,25:end]);  
% psdOut_raw = psdOut_raw(:,[1:7,9,11,13:22,25:end]); % C1-M1
% psdOut_raw = psdOut_raw(:,[1:7,9,11,13:22]); % C1-M1
% psdOut_raw = psdOut_raw(:,[5:end]); % C4-M2
% psdOut_raw = psdOut_raw(:,[3:end]); % C4-M4
% psdOut_raw = psdOut_raw(:,[11:end]); % C6-M1
% psdOut_raw = psdOut_raw(:,[5:end]); % C6-M4
% psdOut_raw = psdOut_raw(:,[10:end]); % C7-M1, M4; C2-M1
% psdOut_raw = psdOut_raw(:,[10:end]); % C8-M3
% psdAvg_raw = mean(psdOut_raw,2);

%% Undersamplig set of PSD
subavg_block = 1;
psdOutUnder = subaveragePSD(psdOut_raw,subavg_block);
% % psdOut = psdOut_raw;
% psdOut = psdAvg_raw;
% psdOutUnder = psdOut_raw;
% break
%% Removing spectra with contrast less than 2 decades ("problematic" for fitting)
[M0 N0] = size(psdOutUnder);
time_interval0 = 30;   % Time between consecutive samples (s), from config file
time_interval = subavg_block*time_interval0;     % Time between consecutive acquisitions [s]
time = (time_interval.*(1:N0)/60)';

psdOut_raw2 = zeros(M0,1);
ContrastThres = 1;
for i = 1:N0
    psdContrast(i,1) = log10(mean(psdOutUnder(1:10,i))/mean(psdOutUnder(end-50:end,i)));
    if psdContrast(i,1)>=ContrastThres
        psdOut_raw2 = horzcat(psdOut_raw2,psdOutUnder(:,i));
    end
end
time = time(psdContrast>=ContrastThres);
psdOut_raw2 = psdOut_raw2(:,2:end);

psdOut = psdOut_raw2;
% psdOut = mean(psdOut_raw2,2);   % Work only with the average PSD
[M N] = size(psdOut);

%% DC Voltage
DCOut = extractDC([InputFileFolder InputFileName InputFileDate]);
DC_level = DCOut(:,2);
DC_time = time_interval0*(1:length(DC_level(:,1)))./60;   

%% Area of PSDs
for i = 1:N
    PSD_Int(i,1) = trapz(xFreq,psdOut(:,i));  % Total PSD power (same as beta, just a different variable name)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                     FIRST INSPECTION                  %%%%%%%%%%
%% PSD evolution colormap
h = figure;
AX = axes;
surf(time,xFreq,log10(psdOut));
set(AX,'YScale','log'), colormap(jet), caxis([-11 -5])
shading interp, view([0 0 1]), colorbar()
xlim([min(time) max(time)]), ylim([min(xFreq) max(xFreq)])
set(gca,'FontSize',24)
ylabel('Frequency (Hz)','FontSize',26)
xlabel('Time (min)','FontSize',26)
zlabel('PSD (a.u.)','FontSize',26)
%% Collection of (raw and undersampled) PSDs in log-log plot
figure,
subplot(1,2,1), loglog(xFreq,psdOut_raw,'-','LineWidth',1,'Color',[0.8 0.8 0.8])
hold on, loglog(xFreq,psdAvg_raw,'-k','LineWidth',4), hold off
xlim([min(xFreq) max(xFreq)]), ylim([1e-11 1e-5]), grid on, set(gca,'FontSize',24)
xlabel('Frequency, f (Hz)','FontSize',26), ylabel('PSD, P(f)','FontSize',26), title('Raw data','FontSize',16)
subplot(1,2,2), loglog(xFreq,psdOut,'-','LineWidth',1,'Color',[0.8 0.8 0.8])
hold on, loglog(xFreq,mean(psdOut,2),'-k','LineWidth',4), hold off
xlim([min(xFreq) max(xFreq)]), ylim([1e-11 1e-5]), grid on, set(gca,'FontSize',24)
xlabel('Frequency, f (Hz)','FontSize',26), ylabel('PSD, P(f)','FontSize',26), title('Undersampled data','FontSize',16)

%%
figure,
loglog(xFreq,psdOut_raw,'-','LineWidth',1,'Color',[0.8 0.8 0.8])
hold on, loglog(xFreq,psdAvg_raw,'-k','LineWidth',4), hold off
xlim([min(xFreq) max(xFreq)]), ylim([1e-11 1e-5]), grid on, set(gca,'FontSize',26)
xlabel('Frequency, f (Hz)','FontSize',28), ylabel('PSD, P(f)','FontSize',28)

%% PSD power metrics: contrast, total power (integral), and DC voltage
figure,
subplot(1,3,1), plot(time,psdContrast(1:length(time),1),'.b','MarkerSize',20)
hold on, plot(time,ContrastThres.*ones(length(time)),'--r'), hold off
xlabel('time, t (min)'), ylabel('PSD contrast (decades)')
subplot(1,3,2), plot(time,PSD_Int*1e6,'.b','MarkerSize',20)
xlabel('time, t (min)'), ylabel('Total PSD power, \intP(f)df  x10^-^6')
subplot(1,3,3), plot(DC_time,DC_level,'.b','MarkerSize',20)
xlabel('time, t (min)'), ylabel('Total scattering strength, DC (V)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local slope of RAW PSD using several smoothing algorithms
break
% xFreqLog = (xFreq);
% psdLog = log10(psdAvg_raw);
% 
% % Method 1 - Spline by cubic interpolation using length(xx)<length(xFreq)
% xx = linspace(min(xFreqLog),max(xFreqLog),1e3);
% psdLog1 = spline(xFreqLog,psdLog,xx);   % Spline - Cubic interpolation
% slope1 = diff(psdLog1)./diff(xx);
% 
% % Method 2 - Robust linear smoothing with weighted linear least squares and a 1st degree polynomial model
% psdLog2 = smooth(psdLog,0.1,'rloess');  % Smooth
% slope2 = diff(psdLog2)./diff(xFreqLog);
% 
% % Method 3 - Level-dependent denoizing assuming Non-Gausian noise
% deb = psdLog(1);
% scal = 'mln';   % Use a level-dependent estimation of the level noise
% psdLog3 = wden(psdLog-deb,'sqtwolog','s',scal,5,'db3') + deb;
% slope3 = diff(psdLog3)./diff(xFreqLog);

%%
% figure,
% plot(xFreqLog,psdLog,'Color',[0.8 0.8 0.8],'LineWidth',1), hold on
% plot(xx,psdLog1,'-b'),
% plot(xFreqLog,psdLog2,'-k'),
% plot(xFreqLog,psdLog3,'-r'),
% hold off
% 
% figure,
% plot(xx(1:end-1),slope1,'-b'), hold on
% plot(xFreqLog(1:end-1),slope2,'-k'),
% plot(xFreqLog(1:end-1),slope3,'-r'),
% hold off
% ylim([-2.5 0.5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% break
Fit_Params_Evol = zeros(N,2*lorent_options);
Fit_Statistics = zeros(N,5);
PSD_Area = zeros(N,1);
Dh = zeros(N,lorent_options);

Lorent_Reconst_Evol = zeros(length(xFreq(find(xFreq == lowFreqCutoff):find(xFreq == highFreqCutoff))),1);
MSD_Evol = zeros(length(xFreq(find(xFreq == lowFreqCutoff):find(xFreq == highFreqCutoff))),1);
G_Evol = zeros(length(xFreq(find(xFreq == lowFreqCutoff):find(xFreq == highFreqCutoff)))-1,1);
LossTan_Evol = zeros(length(xFreq(find(xFreq == lowFreqCutoff):find(xFreq == highFreqCutoff)))-1,1);
LossTanAvg = zeros(N,1);
entropy_param = zeros(N,1);
entropy_param_fit = zeros(N,1);

%%
% for i = 1:N
%     beta = trapz(xFreq,psdOut(:,i));
%     psdAvg = (psdOut(:,i))/beta;
%     % Smoothing
%     xFreq0 = xFreq;
%     psdAvg = psdOut(:,i);
%     [psd_Smooth0 psdAvg_Smooth_error] = Sectioned_Smoothing(xFreq0,psdAvg,Breaking_freq,SM_Method_a,Span_a,SM_Method_b,Span_b);
%     psd_Smooth1 = smooth(xFreq0,psdAvg,Span_b,'loess');
%     psd_Smooth2 = smooth(xFreq0,psdAvg,Span_b,'rloess');
%     psd_Smooth3 = smooth(xFreq0,psdAvg,Span_b,'lowess');
%     psd_Smooth4 = smooth(xFreq0,psdAvg,Span_b,'rlowess');
% 
%     % Display    
%     figure(4),
%     loglog(xFreq,psdOut(:,i),'-b','LineWidth',2)
%     hold on
%     loglog(xFreq0,psd_Smooth0,'-k','LineWidth',1)
%     loglog(xFreq0,psd_Smooth1,'-g','LineWidth',1)
%     loglog(xFreq0,psd_Smooth2,'-m','LineWidth',1)
%     loglog(xFreq0,psd_Smooth3,'-r','LineWidth',1)
%     loglog(xFreq0,psd_Smooth4,'-c','LineWidth',1)
%     hold off
% end


%%
NoMoments = 10;
psdn_moms = [];

% break
for i = 1:N
    beta = trapz(xFreq,psdOut(:,i));
    psdAvg = (psdOut(:,i))/beta;
    % Entropy of normalized, raw PSD
    entropy_param(i,1) = -trapz(xFreq,(psdAvg.*(log(psdAvg))));
    entropy_param(i,2) = -trapz(xFreq,(psdAvg.*(log2(psdAvg))));
    entropy_param(i,3) = -trapz(xFreq,(psdAvg.*(log10(psdAvg))));

    % Computing moments 'by hand' (not using "moment" command)
    for ii = 1:NoMoments
        aux2 = ((log10(xFreq)).^(ii-1)).*(-log10(psdAvg));    % Raw moments
        mom(ii,1) = trapz(xFreq,aux2);
        aux3 = ((log10(xFreq))-mean((-log10(psdAvg))).^(ii-1)).*(-log10(psdAvg));    % Central moments
        mom_cen(ii,1) = trapz(xFreq,aux3);
    end
    psdn_moms = horzcat(psdn_moms,[mom mom_cen]);
    display(mom)
    display(mom_cen)

    % Smoothing
    xFreq0 = xFreq(find(xFreq == lowFreqCutoff_Sm):find(xFreq == highFreqCutoff_Sm));
    psdAvg = psdAvg(find(xFreq == lowFreqCutoff_Sm):find(xFreq == highFreqCutoff_Sm),:);
    [psdAvg_Smooth psdAvg_Smooth_error] = Sectioned_Smoothing(xFreq0,psdAvg,Breaking_freq,SM_Method_a,Span_a,SM_Method_b,Span_b);

    % Fitting
%     x0 = [rand(1) rand(1) rand(1) rand(1)];
    xFreqn = xFreq(find(xFreq == lowFreqCutoff):find(xFreq == highFreqCutoff));
    PSDn = psdAvg_Smooth(find(xFreq == lowFreqCutoff):find(xFreq == highFreqCutoff));
    Fit_Params = Lorentzian_Fit(xFreqn,PSDn,lorent_options,x0,random_fit,Ntimes);
    Lorent_Reconst = Lorentzian_Reconstruction(xFreqn,lorent_options,Fit_Params);

    % Entropy of normalized, fit PSD
    entropy_param_fit(i,1) = -trapz(xFreqn,(Lorent_Reconst(:,end).*(log(Lorent_Reconst(:,end)))));
    entropy_param_fit(i,2) = -trapz(xFreqn,(Lorent_Reconst(:,end).*(log2(Lorent_Reconst(:,end)))));
    entropy_param_fit(i,3) = -trapz(xFreqn,(Lorent_Reconst(:,end).*(log10(Lorent_Reconst(:,end)))));

    % Average corner frequency from weighted average of Taus
    fc_avg(i,1) = 1/sum(Fit_Params(1:end/2).*Fit_Params(end/2+1:end));
    
    Lorent_Reconst_Evol = horzcat(Lorent_Reconst_Evol,[psdAvg(find(xFreq == lowFreqCutoff):find(xFreq == highFreqCutoff)) Lorent_Reconst]);

        % Fit Statistics
    Fit_error = PSDn - Lorent_Reconst(:,end);
    SS_res = sum(Fit_error.^2);
    SS_tot = sum((PSDn - mean(PSDn)).^2);
    R2 = SS_res/SS_tot;
    
    % Reconstruction of Decaying Exponentials
    t = flipud(1./(xFreqn));
    Exp_Reconst = Exponential_Reconstruction(t,lorent_options,Fit_Params);

    % Calculation of the MSD
%     MSD_Th = 6*D0_Th*t;     % Theoretical MSD (using theoretical Diff. coeff.)
    MSD = -(1e12)*(6/(q^2))*log(Exp_Reconst(:,end));   % [um^2-s^-1]
    MSD_Evol = horzcat(MSD_Evol,MSD);
    
    % G-Moduli
    omega = 1./t(1:end-1);
    alpha = diff(log(MSD))./diff(log(t));
    ampG = (kb*T)./(pi*particleRadius.*((1e-12)*MSD(1:end-1)).*gamma(1+alpha));
    Gp = ampG.*cos(pi*alpha/2);
    Gpp = ampG.*sin(pi*alpha/2);
    G_Evol = horzcat(G_Evol,[Gp Gpp]);

    % Loss Tangent
    LossTan = Gpp./Gp;
    LossTanAvg(i,1) = mean(LossTan);
    LossTan_Evol = horzcat(LossTan_Evol,LossTan);
    
    % Evolution of Fit Params, Area of PSD & Hydrodynamic Diameter
    Fit_Statistics(i,:) = [mean(Fit_error) mean(abs(Fit_error).^2) SS_res SS_tot R2];
    Fit_Params_Evol(i,:) = Fit_Params';
    PSD_Area(i,1) = beta;
    aux = (kb*q^2)*T./(6*(pi^2)*visc);
%     Taus = Fit_Params(2,1);   % 1-Lorentzian Fit
%     Taus = sort(Fit_Params(3:4,1),'descend');     % 2-Lorentzian Fit
%     Taus = sort(Fit_Params(4:6,1),'descend');     % 3-Lorentzian Fit
%     Taus = sort(Fit_Params(5:8,1),'descend');     % 4-Lorentzian Fit
    switch lorent_options
    case 1
        Taus = Fit_Params(2,1);   % 1-Lorentzian Fit
    case 2
        Taus = sort(Fit_Params(3:4,1),'descend');     % 2-Lorentzian Fit
    case 3
        Taus = sort(Fit_Params(4:6,1),'descend');     % 3-Lorentzian Fit
    case 4
        Taus = sort(Fit_Params(5:8,1),'descend');     % 4-Lorentzian Fit
    case 5
        Taus = sort(Fit_Params(6:10,1),'descend');    % 5-Lorentzian Fit
    end
    Dh(i,:) = (1e9)*aux.*(Taus');     % Hydrodynamic diameter in nm
    visc_exp(i,:) = ((q^2)*kb*T./(6*(pi^2)*(2*particleRadius)))*(Taus');
%     display(1e3.*[visc visc_exp])
%     display('Hydrodynamic size(s) in nm :')
%     display(Dh')


    % Display    
    figure(4),
    subplot(2,2,1), loglog(xFreqn,PSDn,'Color',[0.333333 0 0.666667],...
                'LineStyle','none', 'LineWidth',1,'Marker','.', 'MarkerSize',12)
    hold on
    loglog(xFreqn,Lorent_Reconst(:,1:(end-1)),'--k','LineWidth',2)
    loglog(xFreqn,Lorent_Reconst(:,end),'-r','LineWidth',2)
    hold off
    set(gca,'FontSize',20)
    title(['Processing PSD  ' num2str(i) '  of  ' num2str(N)],'FontSize',20)
%     cla
    subplot(2,2,2), semilogx(xFreqn(1:end-1),diff(log10(Lorent_Reconst(:,end)))./diff(log10(xFreqn)),'-r','LineWidth',2)
    ylim([-2.1 0.1])
    grid on
    set(gca,'FontSize',20)
    subplot(2,2,3), loglog(t,MSD,'-b','LineWidth',2)
    xlim([min(t) max(t)])
    xlabel('t [s]','FontSize',22)
    ylabel('MSD [um^2]','FontSize',22)
    set(gca,'FontSize',20)
    subplot(2,2,4), loglog(omega,Gp,'-r',omega,Gpp,'-b','LineWidth',2)
    xlabel('\omega [rad/s]','FontSize',22)
    ylabel('G-Moduli [Pa]','FontSize',22)
    xlim([min(omega) max(omega)])
    legend(' Gp',' Gpp','Location','NorthWest')
    set(gca,'FontSize',20)
%     pause(0.25)


%     % Fit of the MSD curve
%     fit_options = fitoptions('Method','NonlinearLeastSquares',...
%                'Lower',[0 0 0],...
%                'Upper',[10 5 1],...
%                'Startpoint',[0.01 0.5 0.5]);
% %     ft = fittype('6*(a^2)*((1 - exp(-(4.8221*x/(a^2))^(c)))^(1/c))*(1 + (b/(a^2))*x)');
%     ft = fittype('6*(a^2)*((1 - exp(-(4.8221*x/(a^2))^(c)))^(1/c))*(1 + (b/(a^2))*x)','options',fit_options);
%     
%     [curve_fit,goodness_fit] = fit(t,MSD_Evol(:,i),ft);
%     MSD_Fit = feval(curve_fit,t);
%     MSD_coeffs_evol(i,:) = coeffvalues(curve_fit);
%     display(MSD_coeffs_evol(i,:))
%     
%     figure(5),
%     loglog(t,MSD_Evol(:,i),'-b')
%     hold on
%     loglog(t,MSD_Fit,'-r')
   
end

%%

figure,
loglog(xFreqn,PSDn,'Color',[0.333333 0 0.666667],...
            'LineStyle','none', 'LineWidth',1,'Marker','.', 'MarkerSize',12)
hold on
loglog(xFreqn,Lorent_Reconst(:,1:(end-1)),'--k','LineWidth',2)
loglog(xFreqn,Lorent_Reconst(:,end),'-r','LineWidth',2)
hold off
set(gca,'FontSize',20)
xlabel('f (Hz)','FontSize',22), ylabel('P(f) (Hz)','FontSize',22)

figure, semilogx(xFreqn(1:end-1),diff(log10(Lorent_Reconst(:,end)))./diff(log10(xFreqn)),'-r','LineWidth',3)
ylim([-2.1 0.1])
grid on
set(gca,'FontSize',20)
xlabel('f (Hz)','FontSize',22), ylabel('dP_f_i_t(f)/df')

figure, loglog(omega,Gp,'-r',omega,Gpp,'-b','LineWidth',3)
xlim([min(omega) max(omega)])
xlabel('\omega (rad/s)','FontSize',22), ylabel('G-Moduli (Pa)','FontSize',22)
legend(' Gp',' Gpp','Location','NorthWest')
set(gca,'FontSize',20)

figure, loglog(omega,Gpp./Gp,'-k','LineWidth',3)
xlim([min(omega) max(omega)])
xlabel('\omega (rad/s)','FontSize',22), ylabel('Loss Tan.','FontSize',22)
set(gca,'FontSize',20)

% break
% figure,
% loglog(xFreq,psdOut_raw,'-b','LineWidth',1)
% xlim([1e0 1e4]), ylim([1e-12 1e-3]), grid on, set(gca,'FontSize',24)
% xlabel('Frequency, f (Hz)','FontSize',26), ylabel('PSD, P(f)','FontSize',26), title('Raw data','FontSize',16)
% 
%%
Lorent_Reconst_Evol = Lorent_Reconst_Evol(:,2:end);
MSD_Evol = MSD_Evol(:,2:end);
G_Evol = G_Evol(:,2:end);
LossTan_Evol = LossTan_Evol(:,2:end);

%%
% figure,
% subplot(1,3,1), plot(time, Fit_Params_Evol(:,1:lorent_options),'o')
% title('Fit Parameters  -  Relative Amplitude')
% xlabel('Time, t [min]','FontSize',26)
% ylabel('a_i [A.U.]','FontSize',26)
% set(gca,'FontSize',17)
% grid on
% subplot(1,3,2), semilogy(time, Fit_Params_Evol(:,lorent_options+1:end),'o')
% title('Fit Parameters  -  Relaxation Time')
% xlabel('Time, t [min]','FontSize',26)
% ylabel('\tau_i [s]','FontSize',26)
% set(gca,'FontSize',17)
% grid on
% subplot(1,3,3), semilogy(time, Dh,'o')
% title('Retrieved Hydrodynamic Size')
% xlabel('Time, t [min]','FontSize',26)
% ylabel('D_h [nm]','FontSize',26)
% set(gca,'FontSize',17)
% grid on
% 
% %%
% figure,
% semilogy(time, Dh,'o')
% title('Retrieved Hydrodynamic Size')
% xlabel('Time, t [min]','FontSize',26)
% ylabel('D_h [nm]','FontSize',26)
% set(gca,'FontSize',17)
% grid on
% 
% %%
% % break
% [maux naux] = size(MSD_Evol);
% h = figure;
% AX = axes;
% surf(time(1:naux),t,log10(MSD_Evol));
% set(AX,'YScale','log','XGrid','off','YGrid','off','ZGrid','off'),
% % set(gca,'YScale','log','ZScale','log','XGrid','off','YGrid','off','ZGrid','off')
% colormap jet,
% caxis([-4 1])
% shading interp, view([0 0 1]), colorbar()
% xlim([min(time) max(time)]), ylim([min(t) max(t)]),
% % zlim([1e-3 1e1])
% set(gca,'FontSize',24)
% ylabel('Lag time, \tau (s)','FontSize',26)
% xlabel('Time (min)','FontSize',26)
% title('MSD(\tau,t) (\mum^2)','FontSize',22)
% 
% %%
% figure,
% subplot(1,2,1), surf(time(1:naux),flipud(omega),G_Evol(:,1:2:end))
% set(gca,'YScale','log','ZScale','log','XGrid','off','YGrid','off','ZGrid','off')
% xlim([min(time) max(time)])
% ylim([min(omega) max(omega)])
% % zlim([1e-3 1e1])
% colormap jet, shading interp
% view([1 1 1])
% subplot(1,2,2), surf(time(1:naux),flipud(omega),G_Evol(:,2:2:end))
% set(gca,'YScale','log','ZScale','log','XGrid','off','YGrid','off','ZGrid','off')
% xlim([min(time) max(time)])
% ylim([min(omega) max(omega)])
% % zlim([1e-3 1e1])
% colormap jet, shading interp
% view([1 1 1])
% 
% %%
% figure,
% subplot(1,2,1), surf(time(1:naux),flipud(omega),LossTan_Evol)
% set(gca,'YScale','log','ZScale','log','XGrid','off','YGrid','off','ZGrid','off')
% xlim([min(time) max(time)])
% ylim([min(omega) max(omega)])
% % zlim([1e-3 1e1])
% colormap jet, shading interp
% view([1 1 1])
% subplot(1,2,2), plot(time,LossTanAvg,'.-b','LineWidth',1,'MarkerSize',10)
% % hold on
% % plot(min(time):1e-2:max(time),spline(time,LossTanAvg,min(time):1e-2:max(time)),'-r')
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit of the MSD curve & average cage size
% break
% fit_options = fitoptions('Method','NonlinearLeastSquares',...
%                'Startpoint',[0.01 0.5 0.01]);
[maux naux] = size(MSD_Evol);
MSD_coeffs_evol = zeros(naux,4);
MSD_fit_goodness_evol = zeros(naux,5);

%%
for i = 1:naux
    fit_options = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0 0 0 0],...
               'Upper',[10 5 1 50],...
               'Startpoint',[0.05 0.5 0.5 20]);
%                'Startpoint',[0.01 0.5 0.5]);
%     ft = fittype('6*(a^2)*((1 - exp(-(4.8221*x/(a^2))^(c)))^(1/c))*(1 + (b/(a^2))*x)');
%     ft = fittype('6*(a^2)*((1 - exp(-(4.8221*x/(a^2))^(c)))^(1/c))*(1 + (b/(a^2))*x)','options',fit_options);
    ft = fittype('6*(a^2)*((1 - exp(-(d*x/(a^2))^(c)))^(1/c))*(1 + (b/(a^2))*x)','options',fit_options);
    
    [curve_fit,goodness_fit] = fit(t,MSD_Evol(:,i),ft);
    MSD_Fit(:,i) = feval(curve_fit,t);
    MSD_coeffs_evol(i,:) = coeffvalues(curve_fit);
    MSD_fit_goodness_evol(i,:) = [goodness_fit.sse goodness_fit.rsquare goodness_fit.dfe goodness_fit.adjrsquare goodness_fit.rmse];

    figure(12),
    loglog(t,MSD_Evol(:,i),'-b','LineWidth',2)
    hold on
    loglog(t,MSD_Fit(:,i),'-r','LineWidth',2)
%     title(['Fitting MSD  ' num2str(i) '  of  ' num2str(naux)],'FontSize',20)
    hold off
    legend(' Exp',' Fit','Location','NorthWest')
    set(gca,'FontSize',20)
    xlabel('Time (s)'), ylabel('MSD (\mum^2)')
    xlim([1e-4 1]), ylim([1e-4 1])
%     pause(0.1)
end

%%
six_delta_squared = MSD_coeffs_evol(:,1);
six_delta_squared_sm = smooth(six_delta_squared,0.1,'rloess');

%%
figure,
plot(time(1:naux),MSD_fit_goodness_evol(:,2))
ylim([0.95 1.05]), ylabel('R^2')
figure,
plot(time(1:naux),six_delta_squared,'.b',time(1:naux),six_delta_squared_sm,'-r','MarkerSize',12,'LineWidth',2)
ylim([0 0.1])

%% Time evolution of the cage size
AvgCageSize_Evol = 2e3*(particleRadius*1e6 + (MSD_coeffs_evol(:,1)/6).^(1/2));
figure,
plot(time(1:naux),AvgCageSize_Evol)
xlabel('Time (min)'),ylabel('Cage size (nm)')
ylim([100 350])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Writing Reports
WriteOutputFile = 0;
if WriteOutputFile == 1
    [SUCCESS,MESSAGE] = xlswrite(OutputFileName,[xFreq psdOut(:,1:end)],'Raw1');
%     [SUCCESS,MESSAGE] = xlswrite(OutputFileName,[xFreq psdOut(:,101:200)],'Raw2');
%     [SUCCESS,MESSAGE] = xlswrite(OutputFileName,[xFreq psdOut(:,201:end)],'Raw3');
%     [SUCCESS,MESSAGE] = xlswrite(OutputFileName,[xFreq psdOut(:,301:400)],'Raw4');
%     [SUCCESS,MESSAGE] = xlswrite(OutputFileName,[xFreq psdOut(:,401:500)],'Raw5');
%     [SUCCESS,MESSAGE] = xlswrite(OutputFileName,[xFreq psdOut(:,501:end)],'Raw6');
    [SUCCESS,MESSAGE] = xlswrite(OutputFileName,[time Fit_Params_Evol Fit_Statistics PSD_Area Dh],Sheet_Summary);
end
