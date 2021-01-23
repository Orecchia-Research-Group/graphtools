function[timeArr] = TimeTestFolder(folderPath, runNumber)

if(~exist("runNumber", "var"))
        runNumber = 1;
end

eg2 = dir(fullfile(folderPath, '*eg2'));
ptn = dir(fullfile(folderPath, '*ptn'));

%Initilaization of array where 1st value of each row is the average runtime
%without matching, the 2nd value is the standard deviation without
%matching, the 3rd and 4th are the same, respectively, but with matching
timeArr = zeros(length(eg2), 4);

for i = 1 : length(eg2)
    eg2name = eg2(i).name;
    eg2name = eg2name(1:end-3);
    for(k=1 : length(ptn))
	    ptnname = ptn(i).name;
	    ptnname = ptnname(1:end-3);
	    if(eg2name == ptnname)
            eg2path = fullfile(eg2(i).folder, eg2(i).name);
            ptnpath = fullfile(ptn(i).folder, ptn(i).name);
            [avgNM, avgM, stdNM, stdM] = TimeTest(eg2path, ptnpath, int64(1), int64(1), runNumber);
		    timeArr(i, 1) = avgNM;
            timeArr(i, 2) = stdNM;
            timeArr(i, 3) = avgM;
            timeArr(i, 4) = stdM;
		    break;
	    end
    end
end
end
