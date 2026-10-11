% subfunction to z-score continous data on itself
%
% Copyright (c) 2023
% Elizabeth Johnson, PhD
% eljohnson@northwestern.edu

function z_mtx = zcont_5(data_mtx)

disp(' '); disp('performing statistical bootstrapping...'); disp(' ');

npermutes = 1000;
nrand = 5;

z_mtx = nan(size(data_mtx));

ntimes = size(z_mtx,2);

for e = 1:size(z_mtx,1)
    tmp = data_mtx(e,:);

    tmp_dist = zeros(1,npermutes);
    for z = 1:npermutes
        tmp_rand = randsample(tmp, nrand);
        tmp_dist(z) = mean(tmp_rand); 
    end
    
    tmp_mean = repmat(mean(tmp_dist), [1 ntimes]);
    tmp_sd = repmat(std(tmp_dist), [1 ntimes]);

    z_mtx(e,:) = (data_mtx(e,:)-tmp_mean)./tmp_sd;
    clear tmp*
end

clearvars -except z_mtx

end