% extractPSD Extract the power spectral density data from the DLS setup
%
% extractPSD(psdFileName) Opens the .psd file associated with a particular
% data run from the DLS setup. The PSD at the optional timestamps are saved
% to the Matlab workspace.
%
% Kyle M. Douglass, 2010

function [psdOut, avgPSD, varargout] = extractPSD(psdFileName, freqX, sampPerFile, varargin)
%% Open the .psd file
% Set the file parameters
% psdSize = 5e3; % Number of points in each psd
% sampPerFile = 61; % The number of PSD's stored in each file
% freqX = 5:5:2.5e4; % Frequencies recorded by the DAQ

psdSize = length(freqX);

% Establish the filename structure
fileSuffix = 'psd';
files = dir([psdFileName, '*.', fileSuffix]);

% Loop over each psd file
psdOut = [];
for m = 1:length(files)
    fName = [psdFileName, '_', num2str(m), '.', fileSuffix];

    % Open the .psd file and read the data. Data is read into an array
    % where each column contains the PSD data from 10 Hz to 100 kHz at 10
    % Hz resolution.
    psdFile = fopen(fName, 'rb');
    psdTemp = fread(psdFile, [psdSize,sampPerFile], 'double');
    fclose(psdFile);
    
    psdOut = horzcat(psdOut,psdTemp);
end
%% Average the PSD's over the desired ranges
if ~isempty(varargin)
    timeInt = varargin{1}; % Time intervals to average over

    if ~isempty(timeInt)
        avgPSD = zeros(psdSize, length(timeInt)-1); % Initialize average array
        stdPSD = zeros(psdSize, length(timeInt)-1); % Initialize standard dev. array

        % Loop over each time interval and average the PSD's in those intervals
        for m = 1:(length(timeInt) - 1)
            avgPSD(:,m) = 10.^(mean(log10(psdOut(:,timeInt(m):timeInt(m+1))),2));
            stdPSD(:,m) = std(log10(psdOut(:,timeInt(m):timeInt(m+1))),0,2);
        end
    end
else
    avgPSD = mean(psdOut,2);
    stdPSD = 0;
end

varargout = {stdPSD}; % Error bars for the PSD's. Must take logarithm of PSD data first to use.
%% Saved parameters for specific data sets
% Nov. 18, 2010
% fName = 'H:\Home\Random\Shared\CURRENT PROJECTS\3d_speckle\data\nov_18_2010\test4\test4_18-Nov-2010';
% timeInt = [642 740 902 1003 1200 1403 1618 2590 2770];

% Nov. 21, 2010
% fName = 'H:\Home\Random\Shared\CURRENT PROJECTS\3d_speckle\data\nov_21_2010\1umSpheres_2_21-Nov-2010';
% timeInt = [1 208 305 408 507 606 711 807 910 1018];

% Nov. 29, 2010
% fName = 'H:\Home\Random\Shared\CURRENT PROJECTS\Radiation Forces\3d_speckle\data\nov_29_2010\0p33umSpheres_29-Nov-2010';

% Nov. 30, 2010
% fName = 'H:\Home\Random\Shared\CURRENT PROJECTS\Radiation Forces\3d_speckle\data\nov_30_2010\0.33spheres_30-Nov-2010';
% fName = 'H:\Home\Random\Shared\CURRENT PROJECTS\Radiation Forces\3d_speckle\data\nov_30_2010\For Kyle\nolaserforkyle2_30-Nov-2010';

% Jan. 21, 2010
% fName = 'H:\Random\Shared\CURRENT PROJECTS\Radiation
% Forces\3d_speckle\data\jan_21_2011\testDLS330nm_21-Jan-2011';
end
