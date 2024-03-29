clc
clear
close all
%% Import data from text file.
% Script for importing data from the following text file:
%
%    I:\Program\Data\He\Monte_Carlo\data_stats.txt
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2019/03/07 11:13:44

%% Initialize variables.
filename = 'C:\Users\minec\Desktop\Stuff\N2\data_stats.txt';
delimiter = {'       '};
startRow = 6;

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,5]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^[-/+]*\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
data_prog = table;
data_prog.VarName2 = cell2mat(raw(:, 1));
data_prog.VarName3 = cell2mat(raw(:, 2));
data_prog.VarName4 = cell2mat(raw(:, 3));
data_prog.VarName5 = cell2mat(raw(:, 4));
data_prog.VarName6 = cell2mat(raw(:, 5));

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp R;

%% Seperate data
data_prog = table2array(data_prog);

Energy = data_prog(1:50,1);

mob = data_prog(1:50,2);
mob_error = data_prog(1:50,4);

diffusion = data_prog(53:102,2);
diff_error = data_prog(53:102,4);

diff_mob = zeros(50,1);
diff_mob_error = zeros(50,1);

for i=1:50
    diff_mob(i,1) = diffusion(i,1)/mob(i,1);
    diff_mob_error(i,1) = diff_mob(i,1)*sqrt((mob_error(i,1)/mob(i,1))^2+(diff_error(i,1)/diffusion(i,1))^2);
    diff_mob_error(i,2) = (diff_mob_error(i,1)/diff_mob(i,1))*100;
end

town_ion = data_prog(365:414,2);
town_ion_error = data_prog(365:414,4);


%% Experimental data
% Script for importing data from the following text file:
%
%    F:\Program\Data\N2\Exp_data\combined_exp.txt
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2019/03/22 15:12:55

%% Initialize variables.
filename = 'C:\Users\minec\Desktop\Stuff\N2\combined_exp.txt';
delimiter = '\t';

%% Format for each line of text:
%   column1: text (%s)
%	column2: text (%s)
%   column3: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
exp_data = [dataArray{1:end-1}];

%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;

%% Seperate Experimental Data

mob_fisch_name = exp_data(3,1);
mob_fisch = str2double(exp_data(5:9,1:3));

mob_loke_name = exp_data(10,1);
mob_loke = str2double(exp_data(12:60,1:3));

mobility = [mob_fisch ;mob_loke];

diff_wagner_name = exp_data(63,1);
diff_wagner = str2double(exp_data(65:79,1:3));

diff_hae_name = exp_data(80,1);
diff_hae = str2double(exp_data(82:107,1:3));

diff_exp = [diff_wagner; diff_hae];

char_jory_name = exp_data(110,1);
char_jory = str2double(exp_data(112:142,1:3));

char_naidu_name = exp_data(143,1);
char_naidu = str2double(exp_data(145:152,1:3));

char_wagner_name = exp_data(153,1);
char_wagner = str2double(exp_data(155:168,1:3));

char_exp = [char_wagner ;char_naidu ;char_jory];

ion_bag_name = exp_data(171,1);
ion_bag = str2double(exp_data(173:190,1:3));

ion_cook_name = exp_data(191,1);
ion_cook = str2double(exp_data(193:210,1:3));

ion_hay_name = exp_data(211,1);
ion_hay = str2double(exp_data(213:238,1:3));

ion_exp = [ion_bag ;ion_cook ;ion_hay];

% 
 %% Plotting

h = figure;
errorbar(Energy(1:50), mob(1:50), mob_error(1:50), 'r')
hold on

for i=1:length(mobility)
    
    mobility(i,4) = mobility(i,3)*mobility(i,2);
    
end

errorbar(mobility(1:5,1), mobility(1:5,2), mobility(1:5,4),'b')
errorbar(mobility(6:54,1), mobility(6:54,2), mobility(6:54,4), 'k')
hold off
legend('Theoretical', mob_fisch_name, mob_loke_name, 'Position', [0.6,0.75,0.1,0.1])
xlabel('Reduced Electric Field (Td)')
ylabel('Mobility (V m s)^{-1}')
set(gca,'xscale','log','FontSize', 12);
set(gca,'yscale','log','FontSize', 12);
print(h,'C:\Users\minec\Desktop\Stuff\N2\mobility_comp','-djpeg','-r1500');

h = figure;
errorbar(Energy(1:50), diff_mob(1:50), diff_mob_error(1:50), 'r')
hold on

for i=1:length(char_exp)
    
    char_exp(i,4) = char_exp(i,3)*char_exp(i,2);
    
end

errorbar(char_exp(15:22,1), char_exp(15:22,2), char_exp(15:22,4), 'k')
errorbar(char_exp(23:53,1), char_exp(23:53,2), char_exp(23:53,4),'b')
hold off
legend('Theoretical', char_naidu_name, char_jory_name, 'Position', [0.6,0.2,0.1,0.1])
xlabel('Reduced Electric Field (Td)')
ylabel('Characteristic Energy (V)')
set(gca,'xscale','log','FontSize', 12);
set(gca,'yscale','log','FontSize', 12);
print(h,'C:\Users\minec\Desktop\Stuff\N2\diffmob_comp','-djpeg','-r1500');

h = figure;
errorbar(Energy(39:50), town_ion(39:50), town_ion_error(39:50), 'r')
hold on

for i=1:length(ion_exp)
    
    ion_exp(i,4) = ion_exp(i,3)*ion_exp(i,2);
    
end

errorbar(ion_exp(1:18,1), ion_exp(1:18,2), ion_exp(1:18,4), 'b')
errorbar(ion_exp(37:62,1), ion_exp(37:62,2), ion_exp(37:62,4), 'k')
hold off
legend('Theoretical', ion_bag_name ,ion_hay_name, 'Position', [0.6,0.2,0.1,0.1])
xlabel('Reduced Electric Field (Td)')
ylabel('Townsend Coef. (m^2)')
set(gca,'xscale','log','FontSize', 12);
set(gca,'yscale','log','FontSize', 12);
print(h,'C:\Users\minec\Desktop\Stuff\N2\town_ion_comp','-djpeg','-r1500');

%% Linear Plotting
% 
% h = figure;
% errorbar(Energy(26:43), mob(26:43), mob_error(26:43), 'r')
% hold on
% 
% for i=1:length(mobility)
%     
%     mobility(i,4) = mobility(i,3)*mobility(i,2);
%     
% end
% 
% errorbar(mobility(1:5,1), mobility(1:5,2), mobility(1:5,4), 'b')
% errorbar(mobility(32:54,1), mobility(32:54,2), mobility(32:54,4), 'k')
% hold off
% legend('Theoretical', mob_fisch_name, mob_loke_name, 'Position', [0.6,0.75,0.1,0.1])
% xlabel('Reduced Electric Field (Td)')
% ylabel('Mobility (V m s)^{-1}')
% set(gca,'FontSize', 12);
% set(gca,'FontSize', 12);
% print(h,'C:\Users\minec\Desktop\Stuff\N2\mobility_comp_lin','-djpeg','-r1500');

% h = figure;
% errorbar(Energy(26:43), diff_mob(26:43), diff_mob_error(26:43), 'r')
% hold on
% 
% for i=1:length(char_exp)
%     
%     char_exp(i,4) = char_exp(i,3)*char_exp(i,2);
%     
% end
% 
% errorbar(char_exp(15:22,1), char_exp(15:22,2), char_exp(15:22,4), 'b')
% errorbar(char_exp(23:53,1), char_exp(23:53,2), char_exp(23:53,4), 'k')
% hold off
% legend('Theoretical', char_naidu_name, char_jory_name, 'Position', [0.6,0.2,0.1,0.1])
% xlabel('Reduced Electric Field (Td)')
% ylabel('Characteristic Energy (V)')
% set(gca,'FontSize', 12);
% set(gca,'FontSize', 12);
% print(h,'C:\Users\minec\Desktop\Stuff\N2\diffmob_comp_lin','-djpeg','-r1500');
% 
% h = figure;
% errorbar(Energy(42:50), town_ion(42:50), town_ion_error(42:50), 'r')
% hold on
% 
% for i=1:length(ion_exp)
%     
%     ion_exp(i,4) = ion_exp(i,3)*ion_exp(i,2);
%     
% end
% 
% errorbar(ion_exp(1:9,1), ion_exp(1:9,2), ion_exp(1:9,4), 'b')
% errorbar(ion_exp(37:55,1), ion_exp(37:55,2), ion_exp(37:55,4), 'k')
% hold off
% legend('Theoretical', ion_bag_name ,ion_hay_name, 'Position', [0.6,0.2,0.1,0.1])
% xlabel('Reduced Electric Field (Td)')
% ylabel('Townsend Coef. (m^2)')
% set(gca,'FontSize', 12);
% set(gca,'FontSize', 12);
% print(h,'C:\Users\minec\Desktop\Stuff\N2\town_ion_comp_lin','-djpeg','-r1500');