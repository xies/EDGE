function fields = changed_data_info(handles)

data_set_name = handles.data_set;
FILENAME = '../DATA_INFO.csv';

fid = fopen(FILENAME);
% count the number of columns
count_commas = textscan(fid, '%s', 1);
count_commas = count_commas{1};
count_commas = count_commas{1};
num_cols = sum(count_commas == ',') + 1;
fclose(fid);

labels_format_string =              repmat('%s', 1, num_cols);
data_format_string   = strcat('%s', repmat('%n', 1, num_cols-1));
% 
fid = fopen(FILENAME);
labels = textscan(fid, labels_format_string, 1, 'delimiter', ',');  % the 1 means just read this 1 time
data = textscan(fid, data_format_string, 'delimiter', ','); 
fclose(fid);

% get the lines literally
% fid = fopen(FILENAME);
% lines = textscan(fid, '%s', 'delimiter', '\n');
% lines = lines{1};
% fclose(fid);

for i = 1:length(labels)
    labels{i} = char(labels{i});
end

all_data_set_names = data{1};


%%%% because i don't know how to edit files, I just re-create
% the whole file (**sigh**) --> although this is not as silly if sorting
% alphabetically...
% lines_fs = '%s\n';
% fid = fopen(FILENAME, 'w');
% fprintf(fid, lines_fs, lines{1});
% fclose(fid);
fid = fopen(FILENAME);
% name_fs = '%s,'; % fs = 'format string'
% data_fs = repmat('%g,', 1, num_cols-1);
% data_fs = strcat(data_fs(1:end-1), '\n');

% success = 0;
fields = cell(0);  % tells you what fields were changed
for i = 1:length(all_data_set_names)
    if strcmp(all_data_set_names{i}, data_set_name)
        writedata = zeros(1, length(labels)-1);
        for j = 2:length(labels)  % start at 2 to skip the name
            writedata(j-1) = handles.info.(labels{j});
            if writedata(j-1) ~= data{j}(i) && ~(isnan(writedata(j-1)) && isnan(data{j}(i)))
                fields{length(fields)+1} = labels{j};
            end
        end
%         success = 1;
%         fprintf(fid, name_fs, data_set_name);
%         fprintf(fid, data_fs, writedata);
    else
%         fprintf(fid, lines_fs, lines{i+1});
    end
end

% if ~success  % if this is a new data set, we need to add an entry
%     writedata = zeros(1, length(labels)-1);
%     for j = 2:length(labels)  % start at 2 to skip the name
%         writedata(j-1) = handles.info.(labels{j});
%         fields{length(fields)+1} = labels{j};  % tells you what fields were changed (in this case, all)
%     end   
%     fprintf(fid, name_fs, data_set_name);
%     fprintf(fid, data_fs, writedata);
% end

fclose(fid);

