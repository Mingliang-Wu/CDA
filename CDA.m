function  [best_f,best_ff,best_p,time] = CDA(fhd,D,FES_MAX,Xmin_vector,Xmax_vector,func_num) 
%  Community development algorithm (CDA) source codes demo V1.0           %
%                                                                         %
%  Developed in MATLAB R2022b                                             %
%                                                                         %
%  Author and programmer: Mingliang Wu, Date: 2024.05.08                                    %
%                                                                         %
%  e-Mail: 2210377@stu.neu.edu.cn                                         %
%          wuml_neu@163.com                                               %
%                                                                         %
%  Main paper: Wu M, Yang D, Wang Y, et al. A novel community             %
%  development algorithm and its application to optimize main steam       %
%  temperature of supercritical units[J]. Expert Systems with             %
%  Applications, 2024: 124190.                                            %

rand('state', sum(100*clock));
t1 = cputime;
%% Initialization parameters
ppp = 5; active_level = 0.7; aa = 0.25; bb = 0.5; cc = 0.75; N = 100; www = 0.9;
Xmin_matrix = repmat(Xmin_vector,N,1); Xmax_matrix = repmat(Xmax_vector,N,1);  
P = zeros(N,D);
for i = 1:N % Initializing populations
    P(i,:) = Xmin_vector + (Xmax_vector-Xmin_vector).*rand(1,D);
end
F = feval(fhd,P',func_num); % Calculate the fitness
FES = N;
t = 1;
EL = N/ppp;
%% Loop
while FES<FES_MAX
    lp = zeros(1,N);
    [~,index]= sort(F); 
    for iy = 1:ppp % Assign communities 
        take = index(ceil((iy-1)*EL+1):ceil(iy*EL));
        lp(take)=iy;
    end
    s = 1; [~,l] = find(lp==s); CD1 = P(l,:); AF = F(l); [~,FM1] = sort(AF);
    s = 2; [~,l] = find(lp==s); CD2 = P(l,:); AF = F(l); [~,FM2] = sort(AF);
    s = 3; [~,l] = find(lp==s); CD3 = P(l,:); AF = F(l); [~,FM3] = sort(AF);
    s = 4; [~,l] = find(lp==s); CD4 = P(l,:); AF = F(l); [~,FM4] = sort(AF);
    s = 5; [~,l] = find(lp==s); CD5 = P(l,:); AF = F(l); [~,FM5] = sort(AF);
    NP = zeros(N,D); 
    for i = 1:N
        inter = P(i,:); 
        if rand>active_level %SCAE
            if lp(i) == 1; NP(i,:) = inter; end
            if lp(i) > 1
                RR = rand;
                if RR<aa 
                    MP = FM1(1); NP(i,:) =  www*rand(1,D).*(CD1(MP,:)-inter)+inter; % M
                    FES = FES + 1;
                end
                if  aa<=RR && RR<bb 
                    NP(i,:) = inter +  www*rand(1,D).*(sum(CD1(FM1(randperm(7,3)),:)) - 3*inter); % E
                    FES = FES + 1;
                end   
                if bb<=RR && RR<cc
                    MP = FM1(7+1); NP(i,:) =  www*rand(1,D).*(CD1(MP,:)-inter)+inter; % W
                    FES = FES + 1;
                end
                if cc<=RR && RR<=1
                    NP(i,:) = inter +  www*rand(1,D).*(sum(CD1(FM1(randperm(7,3)+7),:)) - 3*inter); % U
                    FES = FES + 1;
                end
            end
         else %SCAE
                s = randi(ppp);
                if  s == 1; CD = CD1; FM = FM1; end
                if  s == 2; CD = CD2; FM = FM2; end
                if  s == 3; CD = CD3; FM = FM3; end
                if  s == 4; CD = CD4; FM = FM4; end
                if  s == 5; CD = CD5; FM = FM5; end
                NP(i,:) = inter + www*rand(1,D).*(sum(CD(FM(randperm(7,3)),:)) - 3*inter);
                FES = FES + 1;
         end
    end
    NP = max(NP, Xmin_matrix); % Repair of infeasible solutions
    NP = min(NP, Xmax_matrix);
    NF = feval(fhd,NP',func_num); % Calculating the fitness of new population
    for ii = 1:N % Retain if progressive
        if NF(ii)<F(ii)
           P(ii,:) = NP(ii,:);
           F(ii) = NF(ii);
        end
    end
%% Recording (in-loop)
    [best_f(t,1),bbb] = min(F);
    best_p(t,:) = P(bbb,:);
    t = t+1;
end
%% Output
best_ff = min(F);
t2 = cputime;
time = t2 - t1; 