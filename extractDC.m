% extractDC Extracts the DC voltage read from DAQ card in the DLS setup.
%
% extractDC(daqFileName) Opens the .par2 file associated with a particular
% data run from the DLS setup. 
%
% Rafael Guzman, 2015

function [DCOut, varargout] = extractDC(DCFileName, varargin)
%% Open the .par2 file
% Establish the filename structure
fileSuffix = 'par2';
files = dir([DCFileName, '*.', fileSuffix]);
header_lines = 19;

% Loop over each daq file (there is only one .par2 file at all times; code is generic)
DCOut = [];    DCTemp = [];
for m = 1:length(files)
    fName = [DCFileName, '_', num2str(m), '.', fileSuffix];
    % Open the .par2 file and read the data. Data is read into a 2-column array
    % where the first and second column contain the time (in seconds) and
    % the DC level, respectively.
    DCFile = fopen(fName, 'rb');

    res = {};
    while ~feof(DCFile)
      thisline = fgetl(DCFile);
      if ~ischar(thisline); break; end
      res{end+1,1} = thisline;
    end
    number_of_lines = numel(res);
    fclose(DCFile);

    DCFile = fopen(fName, 'rb');
    for i = 1:header_lines   % Ignore the header lines
        fgetl(DCFile);
    end
    for i = 1:(number_of_lines - header_lines)  % Read rest of file
        DCTemp_aux = str2num(fgetl(DCFile));
        DCTemp = vertcat(DCTemp,DCTemp_aux);
    end
    fclose(DCFile);
   
    DCOut = vertcat(DCOut,DCTemp);
end
