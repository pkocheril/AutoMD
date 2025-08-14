clear

inputfile = fopen('SP_VAC.chg', 'r');
C = textscan(inputfile, '%s', 'Delimiter', '\n');
fclose(inputfile);
D = regexp(C{1,1}, '\s+', 'split');
D = vertcat(D{:});

vac = str2double(D(:,5));

files = dir('*.chg*'); % Obtain all .chg files
totalFiles = numel(files);
fname = strings(totalFiles,1);

for n = 1:totalFiles
    fname(n)=extractBetween(string(files(n).name),'SP_','.chg');
end
varTypes = repelem("double", totalFiles);
resp = table('Size',[length(vac), totalFiles],...
    'VariableTypes',varTypes,'VariableNames',fname); 

for n =1:totalFiles % Loop through each solvent
    fid = fopen(files(n).name, 'r');
    C = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    D = regexp(C{1,1}, '\s+', 'split');
    D = vertcat(D{:});
    temp = str2double(D(:,5));
    temp = 0.5.*vac + 0.5.*temp; % RESP calculation
    resp{:,n}=temp;
    D(:,5)=num2cell(temp);
    writecell(D,append(fname(n),'.resp'),'FileType','text','Delimiter',' ');
end

writetable(resp,'all.resp','FileType','text','Delimiter',' ')
