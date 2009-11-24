function [data_set] = import_get_dataset_name(src)

% get the name of the directory
slash_loc = strfind(src, filesep);
slash_loc = slash_loc(end);
actual_dirname = src(slash_loc+1:end);
prompt = {'Enter a name for this data set:'};
dlg_title = 'Data set name';
num_lines = 1;
defaultanswer = {actual_dirname};
options.Resize='on';
options.WindowStyle='normal';
    

while 1   % until you pick an ok name
    data_set = inputdlg(prompt, dlg_title, num_lines, defaultanswer, options);
    if isempty(data_set)
        return
    end
    data_set = data_set{1};

    foldername = fullfile('..', 'DATA_GUI', data_set);
    if exist(foldername, 'dir')
        res = questdlg('A data set with this name already exists. Do you want to replace it?', ...
            'Name already exists', 'Overwrite', 'Change name', 'Cancel', 'Overwrite');
        switch res
            case 'Overwrite'
                res2 = questdlg('Are you sure? All previous data will be deleted immediately.', ...
                    'Are you sure?', 'Yes', 'Cancel', 'Yes');
                switch res2
                    case 'Yes'
                        rmdir(foldername, 's'); % 's' removes subdirectories, contents
                        rmdir(fullfile('..', 'DATA_SEMIAUTO', data_set), 's');
                        % it will all be overwritten anyway
                        break;
                    case 'No'
                        continue;
                end
            case 'Change name'
                continue;
            case 'Cancel'
                data_set = [];
                return;
        end
    else
        break;
    end
end