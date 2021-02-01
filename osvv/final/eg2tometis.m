function [] = eg2tometis(inputDirectory)
%EG2TOMETIS Summary of this function goes here
%   Detailed explanation goes here

files = dir(fullfile(inputDirectory, '*.eg2'));
for f=1:length(files)
    [~, dataset, ~] = fileparts(files(f).name);
    fprintf('Processing %s...\n', dataset);
    outfile = fullfile(inputDirectory, sprintf('%s.metis', dataset));
    if isfile(outfile)
        continue;
    end
    [G, n, m] = loadeg2graph(fullfile(files(f).folder, files(f).name));
    fOut = fopen(outfile, 'w');
    weight = int64(full(sum(G)));
    [column, row, value] = find(G);
    fprintf(fOut, '%d %d 011\n', n, m/2);
    fprintf(fOut, '%d', weight(1));
    prev_row = 1;
    for i=1:length(value)
        % Print degree
        if row(i) ~= prev_row
            for k=prev_row+1:row(i)
                fprintf(fOut, '\n%d', weight(k));
            end
            prev_row = row(i);
        end
        
        % Print current edge
        fprintf(fOut, ' %d %d', column(i), value(i));
    end
end

end

