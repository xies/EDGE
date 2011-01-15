% Outputs the date and time in a nice format to be printed to the screen.

function out = date_and_time()

x=clock;
y=x(1);m=x(2);d=x(3);h=x(4);mn=x(5);s=x(6); 
if mn < 10
    mn = strcat('0',num2str(mn));
else
    mn = num2str(mn);
end
if h >= 12
    if h > 12
        h = h-12;
    end
    ampm = 'PM';
else
    ampm = 'AM';
end
if h == 0
    h = 12;
end
out = sprintf('%d/%d/%d, %d:%s %s', m, d, y, h, mn, ampm);