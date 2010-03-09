function info = read_data_info(data_set_name)
% reads the CSV version of data_info for data_set_name
% if the data set is not foumd, it returns an empty array

FILENAME = fullfile('..', 'DATA_INFO.csv');
% assumes the first column is the name and the rest are
% numerical values

% count the number of columns by counting the number of commas in the first
% row and adding 1 (because # columns = # commas + 1)
fid = fopen(FILENAME);
count_commas = textscan(fid, '%s', 1);
count_commas = count_commas{1};
count_commas = count_commas{1};
num_cols = sum(count_commas == ',') + 1;
fclose(fid);

% get the labels and all the data
labels_format_string =              repmat('%s', 1, num_cols);
data_format_string   = strcat('%s', repmat('%n', 1, num_cols-1));

fid = fopen(FILENAME);
labels = textscan(fid, labels_format_string, 1, 'delimiter', ',');  % the 1 means just read this 1 time
data = textscan(fid, data_format_string, 'delimiter', ','); 
fclose(fid);


for i = 1:length(labels)
    labels{i} = char(labels{i});
end

all_data_set_names = data{1};

info = [];
for i = 1:length(all_data_set_names)
    if strcmp(all_data_set_names{i}, data_set_name)
        % create the "info" structure
%         info = struct(labels);
        info = cell2struct(cell(length(labels), 1), labels, 1);
        for j = 2:length(labels)  % start at 2 to skip the name
            if iscell(data{j})
                % this should never happen
                info.(labels{j}) = data{j}{i};
            else
                info.(labels{j}) = data{j}(i);
            end
        end
        break;
    end
end

      
dummy_val = NaN;
% if didn't find anything, fill the fields with a dummy val
if isempty(info)
    for j = 2:length(labels)
        info.(labels{j}) = dummy_val;
    end
    info.notfound = 1;
else
    info.notfound = 0;
end
    

