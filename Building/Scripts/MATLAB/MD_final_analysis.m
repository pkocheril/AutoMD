% cd '/Volumes/WeiLab/Phil/MD/Benzonitrile/DMSO-water 50-50/08_final_results/'
clear; close all; clc;


%% Load atomic charges from .itp file
uv1_idx = [7 8]; % Unit vector of nitrile C7-N8
uv2_idx = [8 7]; % Unit vector of nitrile N8-C7


workingdirectory = pwd; % get working directory path
itplist = dir(fullfile(workingdirectory,'*.itp')); itplist = {itplist(~[itplist.isdir]).name};

fid = fopen(string(fullfile(workingdirectory,itplist(1))));
i=1;
tline = fgetl(fid);
A{i} = tline;
while true
    tline = fgetl(fid);
    i=i+1;
    A{i} = tline;
    if ~ischar(tline)
        break
    end
    if contains(tline, '[ atoms ]')
        line_idx(1) = i+2;
    elseif contains(tline, '[ bonds ]')
        line_idx(2) = i-2;
    end
end
fclose(fid);
i=1;
chg=zeros(line_idx(2)-line_idx(1)+1,1);
atom=cell(line_idx(2)-line_idx(1)+1,1);
atom_idx=zeros(line_idx(2)-line_idx(1)+1,1);
for m=line_idx(1):line_idx(2)
    temp=strsplit(A{m});
    atom{i}=temp{6};
    atom_idx(i)=str2double(extract(temp{6}, digitsPattern));
    chg(i) = str2double(temp{8});
    i=i+1;
end
%% Load MD data
files = dir('*.xvg'); % Load all .xvg data
for i=1:numel(files)  % Read data into matrices
    if contains(files(i).name, 'xyz.xvg')
        coord = ipt_xvg(files(i).name);
        continue
    end
    if contains(files(i).name, 'f0.xvg')
        f0 = ipt_xvg(files(i).name);
        continue
    end
    if contains(files(i).name, 'f1.xvg')
        f1= ipt_xvg(files(i).name);
        continue
    end
end
coord(:, 1) = [];    % Remove the first column (time axis)
f0(:, 1) = [];
f1(:, 1) = [];

%% Obtain electric field from MD force
ef = f1-f0;     % Subtraction to exclude self-exerted force

uv1 = (coord(:, uv1_idx(1)*3-2:uv1_idx(1)*3) - coord(:, uv1_idx(2)*3-2:uv1_idx(2)*3)); % Obtain carbonyl vector within each frame
uv2 = (coord(:, uv2_idx(1)*3-2:uv2_idx(1)*3) - coord(:, uv2_idx(2)*3-2:uv2_idx(2)*3));
uv1 = uv1./vecnorm(uv1,2,2); % Normalize to unit vector
uv2 = uv2./vecnorm(uv2,2,2);

Fvib1 = zeros(length(coord), length(atom)); % To store 
Fvib2 = zeros(length(coord), length(atom));
for n = 1:length(atom) % Loop through each atom
    eforce = ef(:, atom_idx(n)*3-2:atom_idx(n)*3);
    efield = eforce./chg(atom_idx(n)).*0.1036427;    % Divide by atom charges and convert to MV/cm of field
    Fvib1(:,n) = dot(efield, uv1, 2); % Project electric filed along the carbonyl 
    Fvib2(:,n) = dot(efield, uv2, 2);
end

%% Average electric field between C and O
F1 = squeeze(mean(Fvib1(:, uv1_idx), 2));  
F2 = squeeze(mean(Fvib2(:, uv2_idx), 2));   

%% Histogram fitting and visualization
options=optimoptions(@lsqnonlin);
options.MaxFunctionEvaluations = 1e12; 
options.MaxIterations = 1e12; 
options.FunctionTolerance = 1e-12; 
options.OptimalityTolerance = 1e-12; 
options.StepTolerance = 1e-12; 
options.Display = 'off';

histN=128;
[N,edges] = histcounts(F1, histN);
edges = edges(2:end) - (edges(2)-edges(1))/2;
fun = @(r) r(1).*exp(-((edges-r(2))./r(3)).^2) - N; % Gaussian fitting
r0 = [median(N), 0, 10];
lb = [0, -inf, 0];
ub = [];
fitval = lsqnonlin(fun,r0,lb,ub,options);
fcurve = fitval(1).*exp(-((edges-fitval(2))./fitval(3)).^2);

figure
subplot(1,2,1)
hold on
plot(edges,N,'r.')
plot(edges, fcurve,'b')
formatSpec = ['Between %s and %s',newline,'%.2f ', char(177), ' %.2f MV/cm'];
title(sprintf(formatSpec, char(atom(uv1_idx(1))), char(atom(uv1_idx(2))), fitval(2), fitval(3)/sqrt(2)))
xlabel('Electric field (MV/cm)')
ylabel('Counts')

[N,edges] = histcounts(F2, histN);
edges = edges(2:end) - (edges(2)-edges(1))/2;
fun = @(r) r(1).*exp(-((edges-r(2))./r(3)).^2) - N;
r0(1) = median(N);
fitval = lsqnonlin(fun,r0,lb,ub,options);
fcurve = fitval(1).*exp(-((edges-fitval(2))./fitval(3)).^2);

subplot(1,2,2)
hold on
plot(edges,N,'r.')
plot(edges, fcurve,'b')
title(sprintf(formatSpec, char(atom(uv2_idx(1))), char(atom(uv2_idx(2))), fitval(2), fitval(3)/sqrt(2)))
xlabel('Electric field (MV/cm)')
ylabel('Counts')

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
