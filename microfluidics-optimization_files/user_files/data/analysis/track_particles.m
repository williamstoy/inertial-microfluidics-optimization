% Name of the file
folder = '../dpm/';

filename = 'dpm_data_cw120_ch60_cl500_nh30_nl60_q195_re36_step*.dpm';

strings = split(filename, '_');
ch_found = 0;
cw_found = 0;
for k = 1:length(strings)
    if startsWith(strings{k},'cw')
        cw = str2double(strings{k}(3:end));
        cw_found = 1;
    end

    if startsWith(strings{k},'ch')
        ch = str2double(strings{k}(3:end));
        ch_found = 1;
    end
end

if ~ch_found || ~cw_found
    error('Poorly formed filename string (does not contain cw and ch');
end

d = dir([folder, filename]);
for i = 1:length(d)
    % Open the file
    fid = fopen([folder, d(i).name], 'r');

    % get the step name from the file
    % Use regular expression to match the pattern and extract the step number
    pattern = 'step(\d+)\.dpm';
    tokens = regexp(d(i).name, pattern, 'tokens');
    
    % Convert the extracted number (stored in tokens{1}{1}) to double
    if ~isempty(tokens)
        stepNumber = str2double(tokens{1}{1});
        disp(['Step Number: ', num2str(stepNumber)]);
    else
        error('Pattern not found in filename.');
    end
    
    % Read the file line by line and store the data
    data = [];
    while ~feof(fid)
        line = fgetl(fid);
        
        % Check if line starts with '((' indicating data
        if startsWith(line, '((')
            
            % Split the line based on spaces
            value_texts = strsplit(line(4:end-1), ' ');
            value_texts{12} = value_texts{12}(1:end-1); % remove the trailing ')'
            values = str2double(value_texts);
            
            % Append to data
            data = [data; values];
        end
    end
    
    % Close the file
    fclose(fid);
    
    % Convert the array data to a table
    columnNames = {'x', 'y', 'z', 'u', 'v', 'w', 'diameter', 't', 'mass-flow', 'mass', 'frequency', 'time', 'name'};
    dataTable = array2table(data, 'VariableNames', columnNames);

    plotData{stepNumber} = dataTable;
end

%% display
close all; clc; figure;
for i = 1:length(plotData)
    % Display the table
    plot(plotData{i}.z * 1e6, plotData{i}.y * 1e6, 'ko');
    xlim([-cw/2, cw/2]);
    ylim([-ch/2, ch/2]);
    title(i);
    drawnow;
end
