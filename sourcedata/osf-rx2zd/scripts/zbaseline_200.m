% subfunction to z-score time-series data on the baseline using statistical 
% bootstrapping
%
% Copyright (c) 2023
% Elizabeth Johnson, PhD
% eljohnson@northwestern.edu

function tfrz = zbaseline_200(data, baseline)

disp(' '); disp('performing statistical bootstrapping...'); disp(' ');

npermutes = 1000;
nrand = 200;

tfrz = data;

if length(size(baseline)) == 3
    tfrz.trial = zeros(size(tfrz.trial));

    ntrials = size(tfrz.trial, 1);
    ntimes = size(tfrz.trial, 3);

    for e = 1:size(baseline,2)
        base = baseline(:,e,:);
        base = base(:);

        basedist = zeros(1,npermutes);
        for z = 1:npermutes
            brand = randsample(base, nrand);
            basedist(z) = mean(brand); 
        end

        bmean = repmat(mean(basedist), [ntrials 1 ntimes]);
        bsd = repmat(std(basedist), [ntrials 1 ntimes]);

        tfrz.trial(:,e,:) = (data.trial(:,e,:)-bmean)./bsd;

        clear base basedist brand bmean bsd
    end
    
elseif length(size(baseline)) == 4
    tfrz.powspctrm = zeros(size(tfrz.powspctrm));

    ntrials = size(tfrz.powspctrm, 1);
    ntimes = size(tfrz.powspctrm, 4);

    for e = 1:size(baseline,2)
        for f = 1:size(baseline,3)
            base = baseline(:,e,f,:);
            base = base(:);

            basedist = zeros(1,npermutes);
            for z = 1:npermutes
                brand = randsample(base, nrand);
                basedist(z) = mean(brand); 
            end

            bmean = repmat(mean(basedist), [ntrials 1 1 ntimes]);
            bsd = repmat(std(basedist), [ntrials 1 1 ntimes]);

            tfrz.powspctrm(:,e,f,:) = (data.powspctrm(:,e,f,:)-bmean)./bsd;

            clear base basedist brand bmean bsd
        end
    end

clearvars -except tfrz

end