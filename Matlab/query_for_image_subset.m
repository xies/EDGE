function q = query_for_image_subset()

q = [];
prompt = {'Enter Z-range (e.g., 5-10):', ...
            'Enter T-range (e.g., 0):'};
dlg_title = 'Input which images to process';
num_lines = 2;
defaultanswer = {'',''};
%     options.Resize='on';
options.WindowStyle='normal';
%         options.Interpreter='tex';
answer = inputdlg(prompt, dlg_title, num_lines, defaultanswer, options);

if isempty(answer) || isempty(answer{1})
    return;
end;

z_do   = answer{1}(1,:);        
if size(answer{1}, 1) == 2
    z_skip = answer{1}(2,:);
else
    z_skip = [];
end

t_do   = answer{2}(1,:);
if size(answer{2}, 1) == 2
    t_skip = answer{2}(2,:);
else
    t_skip = [];
end


q.z_do_range   = parse_range(z_do);
q.z_skip_range = parse_range(z_skip);
q.t_do_range   = parse_range(t_do);
q.t_skip_range = parse_range(t_skip);