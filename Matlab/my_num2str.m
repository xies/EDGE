function out = my_num2str(in)
% does num2str but if input is NaN outputs '-'

if ~isnan(in)
    out = num2str(in);
else
    out = '-';
end