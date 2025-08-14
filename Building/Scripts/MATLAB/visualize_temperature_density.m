clear; close all;

temp = ipt_xvg('temperature.xvg');
density = ipt_xvg('density.xvg');

subplot(1,2,1)
plot(temp(:,1),temp(:,2),'.-')
legend(['avg(the last half)= ', num2str(mean(temp(round(end/2)+1:end,2)))])
title('Temperature')
subplot(1,2,2)
plot(density(:,1),density(:,2),'.')
legend(['avg(the last half)= ', num2str(mean(density(round(end/2)+1:end,2)))])
title('Density')

function out1 = ipt_xvg(filename)
inputfile = fopen(filename, 'r');
C = textscan(inputfile, '%s', 'Delimiter', '\n');
fclose(inputfile);
D = C;
nondatarows1 = strfind(C{1}, '#');
nondatarows1 = find(~cellfun('isempty',nondatarows1));
nondatarows2 = strfind(C{1}, '@');
nondatarows2 = find(~cellfun('isempty',nondatarows2));
nondatarows3 = strfind(C{1}, '&');
nondatarows3 = find(~cellfun('isempty',nondatarows3));
D{1,1}([nondatarows1;nondatarows2;nondatarows3])=[];

D = regexp(D{1,1}, '\s+', 'split');
D = vertcat(D{:});
out1 = cellfun(@str2double,D);
end
