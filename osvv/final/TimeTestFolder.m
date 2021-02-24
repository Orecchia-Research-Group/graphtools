function[timeArr, graphNames] = TimeTestFolder(folderPath, runNumber)

if(~exist("runNumber", "var"))
        runNumber = 1;
end

eg2 = dir(fullfile(folderPath, '*eg2'));
ptn = dir(fullfile(folderPath, '*ptn'));

%Initilaization of array where 1st value of each row is the average runtime
%without matching, the 2nd value is the standard deviation without
%matching, the 3rd and 4th are the same, respectively, but with matching
timeArr = zeros(2*length(eg2), 4);
graphNames = strings([1, length(eg2)]);

for i = 1 : length(eg2)
    eg2path = fullfile(eg2(i).folder, eg2(i).name);
    [p, eg2name, ext] = fileparts(eg2path);
    for k = 1 : length(ptn)
        ptnpath = fullfile(ptn(k).folder, ptn(k).name);
        [p, ptnname, ext] = fileparts(ptnpath);
	    if strcmp(eg2name, ptnname)
		    [avgInit, avgS1, avgS2, avgMatch] = TimeTest(eg2path, ptnpath, int64(1), int64(1), runNumber, 'dinic');
		    timeArr(i, 1) = avgInit;
            timeArr(i, 2) = avgS1;
            timeArr(i, 3) = avgS2;
            timeArr(i, 4) = avgMatch;

		    [avgInit, avgS1, avgS2, avgMatch] = TimeTest(eg2path, ptnpath, int64(1), int64(1), runNumber, 'dynamic');
            timeArr(i+length(eg2), 1) = avgInit;
            timeArr(i+length(eg2), 2) = avgS1;
            timeArr(i+length(eg2), 3) = avgS2;
            timeArr(i+length(eg2), 4) = avgMatch;

            eg2name(1) = upper(eg2name(1));
            graphNames(i) = eg2name;
		    break;
	    end
    end
end
end
