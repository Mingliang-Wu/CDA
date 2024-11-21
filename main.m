clc
clear
close

rand('state', sum(100*clock));
D = 30;
Xmin = -100;
Xmax = 100;
Xmin = repmat(Xmin,1,D);
Xmax = repmat(Xmax,1,D);
pop_size = 100;
FES_MAX = D*10000;
runs = 30;
fhd = str2func('cec14_func');
for i = 1:30
    func_num = i;
    for j = 1:runs
        i,j,
        tic
        [BF_CDA{i,j},BF_CDA_E(i,j),~,T_CDA(i,j)] = CDA(fhd,D,FES_MAX,Xmin,Xmax,func_num);
        t1 = toc; T_CDA(i,j) = t1;
        BF_CDA_E_mean = mean(BF_CDA_E,2); T_CDA_mean = mean(T_CDA,2);
        save ('.\CDA','BF_CDA_E_mean','BF_CDA','BF_CDA_E','T_CDA_mean','T_CDA');
    end
end


