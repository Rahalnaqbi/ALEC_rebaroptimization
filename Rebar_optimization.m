%This code is meant to be used for rebar optimization%
I = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]; % Set of beams
J = [8,10,12,14,16,18]; % Number of rebar types
T = 13; % Max number of rebars allowed in a beam
V = 6; % Max number of rebar sizes allowed in a beam
L = [1.55, 2.13, 2.87, 2.907, 3.01, 3.01, 3.36, 3.37, 3.37, 3.4, 3.98, 3.98,4.114, 4.61,5.2, 5.201, 5.23, 5.23, 5.231, 5.266,5.266, 5.4, 6.02, 6.59, 6.592, 6.741, 7.36, 7.9,7.9,7.9]; % Length of the bars for each beam
R = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]; % Number of identical rows in the rebar for each beam
As = [0.000105, 0.000105, 0.000105, 0.000105, 0.000105, 0.000105, 0.000105, 0.000105, 0.000105, 0.000105, 0.000105, 0.000105, 0.000105, 0.000105, 0.000105, 0.0021, 0.0021, 0.0021, 0.0021, 0.0021, 0.0021, 0.0021, 0.0021, 0.0021, 0.0021, 0.0021, 0.0021, 0.0021, 0.0021, 0.0021]; % Area of steel for reinforcing each beam
W = [0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]; % Width of each beam
b = [1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0]; % Indicator if the rebars need bending
D = [0.025; 0.032; 0.036; 0.043; 0.052; 0.057]; % Diameter of each rebar size
P = [2.765; 3.9178; 5.6477; 7.7138 ;9.99; 12.95]; % Purchase price of a standard length rebar
C = [0.0126; 0.0126; 0.0126; 0.0167; 0.0167; 0.0167]; % Cost of each time cutting the rebar
B = [0.013; 0.013; 0.013; 0.013; 0.0142; 0.0161]; % Cost of each time bending the rebar
CC = 0.03; % Concrete clear cover
DT = 0.03; % Diameter of the tie rebar assumed to be equal 
SL = 12; % Length of standard bars

%Finding the set of unique bar lengths
L_prime = unique(L);

%Generating length vector for unique lengths only
for i = 1:length(L_prime)
    for l = 1:length(I)
        if L(l) == L_prime(i)
            LL(i,l) = 1;
        else
            LL(i,l) = 0;
        end
    end
end

% Generate initial cutting patterns 
minL = min(L_prime);
Dmax = max(D);
i = 1;
k= 0;
M = zeros(length(L_prime), length(L_prime)); % Initialize matrix for initial patterns
N = zeros(1,length(L_prime));

%Create initial patterns, count number of patterns and pattern waste 
for l_prime = L_prime
    Mkl = floor(SL / l_prime);
    RL = SL - Mkl * l_prime;
    M(i, i) = Mkl;
    Mmin = floor(RL / minL);
    if l_prime == minL
        M(1,i) = Mkl+Mmin;
    else
        M(1, i) = Mmin;
    end
    N(1,i) = Mkl + Mmin - 1;
    k = k + 1;
    Patternswaste (i,1) = SL - sum(L_prime(1,:)*M(:,i));
    i = i + 1;
end

%Initalize parameters for loop of optimization using column-generation
T2 = T-1;
disp(N);
reducedCost = -Inf;
exitflag = 1;
reducedCostTolerance = -0.001;

while reducedCost < reducedCostTolerance && exitflag > 0 && k<=34 
%Decision variables define
rebarprob = optimproblem('ObjectiveSense','min');
ykj = optimvar('numberofJbarscutwithpatternk',k,length(J),'LowerBound',0);
xtij = optimvar('ifrebarissizej',T,length(J), length(I),'LowerBound',0,'UpperBound',1);
dit = optimvar('distancebetweenrebars',length(I),T2,'LowerBound',0);
fij = optimvar('ifbeamhassizej',length(I),length(J),'LowerBound',0,'UpperBound',1);

%Constraints:

%Constraint 1 ensures more than 2 rebars are used for each concrete element
constraint1 = optimconstr(length(I));
for i=1:length(I)
    constraint1(i) = sum(xtij(:,:,i),"all") >= 2;
end
rebarprob.Constraints.constraint1 = constraint1;

%Constraint 2 maintains the order of rebar 
constraint2 = optimconstr(length(I));
for i = 1:length(I)
    for t = 1:T2
        for t_prime = t+1:T
            sum_xtij_t = sum(xtij(t, :, i));
            sum_xtij_t_prime = sum(xtij(t_prime, :, i));
            constraint2(i) = sum_xtij_t >= sum_xtij_t_prime;
        end
    end
end
rebarprob.Constraints.constraint2 = constraint2;

%Constraint 3 guarantees rebar diameters, spacing, and CC/DT equals final width 
constraint3 = optimconstr(length(I));
sum_dit = optimexpr(length(I),1); 
sum_Dj_xtij = optimexpr(length(I),1); 
for i = 1:length(I)
    sum_dit(i,1) = sum(dit(i,:)); 
    sum_Dj_xtij(i,1) = sum(xtij(:, :, i)*D(:,1));
    constraint3(i) = sum_dit(i,1) + sum_Dj_xtij(i,1) +2*(CC+DT) == W(1,i);
end
rebarprob.Constraints.constraint3 = constraint3;

%Constraint 4 adheres to the set max variation for rebar size 
constraint4 = optimconstr(length(I));
for i = 1:length(I)
    sum_fij = sum(fij(i,:)); 
    constraint4(i) = sum_fij <= V;
end
rebarprob.Constraints.constraint4 = constraint4;

%Constraint 5 maintains uniformity among xtij and fij decision variables 
constraint5 = optimconstr(length(I),length(J));
sumxtij = optimexpr(length(I),length(J)); 
for i = 1:length(I)
    for j=1:length(J)
        sumxtij(i,j) = sum(xtij(:,j,i));
        constraint5(i,j) = sumxtij(i,j) == fij(i,j);
    end
end
rebarprob.Constraints.constraint5 = constraint5;

%Constraint 6 confirms that all (As) requirements are met 
constraint6 = optimconstr(length(I));
formula = zeros(T,1,length(I));
formula = optimexpr(T, 1, length(I)); 
finalarea = optimexpr(length(I),1); 
for i = 1:length(I)
    D2 = D.^2;
    formula(:,:,i) = xtij(:,:,i)*D2;
    finalarea(i,1) = sum(formula(:,1,i))*(pi()/4);
    constraint6(i) = finalarea(i,1) >= As(1,i)/R(1,i);
end
rebarprob.Constraints.constraint6 = constraint6;

%Constraint 7 verifies whether the demand for each rebar length is met 
constraint7 = optimconstr(length(L_prime),length(J));
sum_demand = optimexpr(T,length(I)); 
sum_demand_2 = optimexpr(T,1); 
provided = optimexpr(length(L_prime),length(J));
demand = optimexpr(length(L_prime),length(J)); 
sum_xtij = optimexpr(length(I),length(J)); 
for j = 1:length(J)
    for l=1:length(L_prime)
        provided (l,j)= sum(M(l,:)*ykj(:,j));
        
        for t=1:T
            for i=1:length(I)
            sum_demand(t,i)= LL(l,i)*xtij(t,j,i);
            end
            sum_demand_2(t,1) = sum(sum_demand(t,:));
        end
    
        demand(l,j) = sum(sum_demand_2(:,1));

        constraint7(l,j) = provided(l,j) >= demand(l,j);
    end
end
rebarprob.Constraints.constraint7 = constraint7;

%Objective Function (OF) 
%Define first term of OF as the total purchase cost of rebar
first_term = sum(ykj(:,:)*P(:,1));

%Define second term of OF as total cutting cost of all rebar patterns used
second_term_1 = optimexpr(k,length(J));
for j = 1:length(J)
    second_term_1 (:,j)= N(1,:)*ykj(:,j);
end
second_term = sum(second_term_1(:,:)*C(:,1));

%Define third term of OF as total bending costs 
third_term = optimexpr(length(I),1); 
for i = 1:length(I)
third_term(i,1) = b(1,i)*sum(xtij(:,:,i)*B(:,1));
end
third_termf = sum(third_term(:,1));

objective_function = first_term + second_term + third_termf;

rebarprob.Objective = objective_function;
options = optimoptions(rebarprob);
options.MaxTime = 200;

%Define column-generation process 
CGprocess = optimproblem('ObjectiveSense','min');

%Decision variables associated with column-generation process define 
Ncg = optimvar('numberofcutsinpattern','Type','integer','LowerBound',0);
Mcg = optimvar('numberofrebarsoflengthlinthepattern',length(L_prime),1,'Type','integer','LowerBound',0);

%Find for each rebar length max. amount of times it could be cut from
%standard-length rebar (SL) 
for l = 1:length(L_prime)
    Boundl(l,1)= floor(SL / L_prime(l));
end

%CG-process constraints: 
%Constraint 1 ensures that the proposed lengths in pattern don't exceed SL
constraintcg1 = optimconstr(length(J));
Lengths = sum(L_prime(1,:)*Mcg(:,1));
constraintcg1 = Lengths <= SL; 
CGprocess.Constraints.constraintcg1 = constraintcg1;

%Constraint 2 secures the number of cuts in pattern don't exceed limit 
constraintcg2 = optimconstr(length(J));
rebarcutcons = -1 + sum(Mcg(:,1)) + ((SL-sum(L_prime(1,:)*Mcg(:,1)))/minL);
constraintcg2 = Ncg >= rebarcutcons;
CGprocess.Constraints.constraintcg2 = constraintcg2;

%Constraint 3 verifies number of specific length in pattern isn't over max.
constraintcg3 = optimconstr(length(J));
    for l = 1:length(L_prime)
        constraintcg3 = Mcg(l,1) <= Boundl(l,1);
    end
CGprocess.Constraints.constraintcg3 = constraintcg3;

%Initalize terms for objective function for CG-process
First_termcg = zeros(length(J),1);
second_termcg = optimexpr(length(J),1);
third_termcg = optimexpr(length(J),1); 
Objective_CG = optimexpr(length(J),1);
initialptotalwaste=optimexpr(length(L_prime),1);

%Solve optimization using initial set of patterns 
    [sol,fval,exitflag,~,lambda] = solve(rebarprob,'Options',options);

    %Calculate the initial cost using initial set of patterns 
    if k == length(L_prime)
        Initialcost = fval;
        for f = 1:k
        initialptotalwaste(f,1)= sum(ykj(f,:))*Patternswaste(f,1); 
        end
        initialtotalwaste = sum(initialptotalwaste(:,1));
    end

    %If there is an optimal solution, extract the lagrange multipliers
    %necessary for (OF) of column-generation process from constraint 7
    if exitflag > 0
    f2 = transpose(lambda.Constraints.constraint7);

    %Generate improver pattern by solving CG-process for each rebar size 
for j = 1:length(J)
    First_termcg(j,1) = P(j,1);
    second_termcg(j) = C(j)*Ncg;
    third_termcg(j) = sum(f2(j,:)*Mcg(:,1));
    Objective_CG(j) = First_termcg(j) + second_termcg(j) - third_termcg(j);
    CGprocess.Objective = Objective_CG(j);
    [solt2,reducedCost,pexitflag] = solve(CGprocess);
    newpattern= solt2.numberofrebarsoflengthlinthepattern;
    newpatterncuts = solt2.numberofcutsinpattern;
    disp(reducedCost);
    disp(exitflag);
    disp(k);
    %Ensure the generated pattern cost is within the tolerance
    if pexitflag > 0 && reducedCost < reducedCostTolerance
        %Add the new pattern and cuts for it to initial matrices
        M = [M , newpattern];
        k = k +1;
        N = [N , newpatterncuts];
        %Compute the patteern waste and add it to initial matrix 
        Pi = SL - sum(L_prime(1,:)*newpattern(:,1));
        Patternswaste = [Patternswaste; Pi];
    end
    end
    end
end
disp(M);
disp(k);

% Final optimization with the most optimal cutting patterns generated 
%Decision variables define
rebarprob2 = optimproblem('ObjectiveSense','min');
ykj2 = optimvar('numberofJbarscutwithpatternk',k,length(J),'Type','integer','LowerBound',0);
xtij2 = optimvar('ifrebarissizej',T,length(J), length(I),'Type','integer','LowerBound',0,'UpperBound',1);
dit2 = optimvar('distancebetweenrebars',length(I),T2,'LowerBound',0);
fij2 = optimvar('ifbeamhassizej',length(I),length(J),'Type','integer','LowerBound',0,'UpperBound',1);

%Constraints 
constraint12 = optimconstr(length(I));
for i=1:length(I)
    constraint12(i) = sum(xtij2(:,:,i),"all") >= 2;
end
rebarprob2.Constraints.constraint12 = constraint12;

constraint22 = optimconstr(length(I));
for i = 1:length(I)
    for t = 1:T2
        for t_prime = t+1:T
            sum_xtij_t2 = sum(xtij2(t, :, i));
            sum_xtij_t_prime2 = sum(xtij2(t_prime, :, i));
            constraint22(i) = sum_xtij_t2 >= sum_xtij_t_prime2;
        end
    end
end
rebarprob2.Constraints.constraint22 = constraint22;

constraint32 = optimconstr(length(I));
sum_dit2 = optimexpr(length(I),1); 
sum_Dj_xtij2 = optimexpr(length(I),1); 
for i = 1:length(I)
    sum_dit2(i,1) = sum(dit2(i,:)); 
    sum_Dj_xtij2(i,1) = sum(xtij2(:, :, i)*D(:,1));
    constraint32(i) = sum_dit2(i,1) + sum_Dj_xtij2(i,1) +2*(CC+DT) == W(1,i);
end
rebarprob2.Constraints.constraint32 = constraint32;

constraint42 = optimconstr(length(I));
for i = 1:length(I)
    sum_fij2 = sum(fij2(i,:)); 
    constraint42(i) = sum_fij2 <= V;
end
rebarprob2.Constraints.constraint42 = constraint42;

constraint52 = optimconstr(length(I),length(J));
sumxtij2 = optimexpr(length(I),length(J)); 
for i = 1:length(I)
    for j=1:length(J)
        sumxtij2(i,j) = sum(xtij2(:,j,i));
        constraint52(i,j) = sumxtij2(i,j) == fij2(i,j);
    end
end
rebarprob2.Constraints.constraint52 = constraint52;

constraint62 = optimconstr(length(I));
formula2 = zeros(T,1,length(I));
formula2 = optimexpr(T, 1, length(I)); 
finalarea2 = optimexpr(length(I),1); 
for i = 1:length(I)
    D2 = D.^2;
    formula2(:,:,i) = xtij2(:,:,i)*D2;
    finalarea2(i,1) = sum(formula2(:,1,i))*(pi()/4);
    constraint62(i) = finalarea2(i,1) >= As(1,i)/R(1,i);
end
rebarprob2.Constraints.constraint62 = constraint62;

constraint72 = optimconstr(length(L_prime),length(J));
sum_demand2 = optimexpr(T,length(I)); 
sum_demand_22 = optimexpr(T,1); 
provided2 = optimexpr(length(L_prime),length(J));
demand2 = optimexpr(length(L_prime),length(J)); 
sum_xtij2 = optimexpr(length(I),length(J)); 
for j = 1:length(J)
    for l=1:length(L_prime)
        provided2 (l,j)= sum(M(l,:)*ykj2(:,j));
        
        for t=1:T
            for i=1:length(I)
            sum_demand2(t,i)= LL(l,i)*xtij2(t,j,i);
            end
            sum_demand_22(t,1) = sum(sum_demand2(t,:));
        end
    
        demand2(l,j) = sum(sum_demand_22(:,1));

        constraint72(l,j) = provided2(l,j) >= demand2(l,j);
    end
end
rebarprob2.Constraints.constraint72 = constraint72;

%Objective Function 
first_term2 = sum(ykj2(:,:)*P(:,1));
second_term_12 = optimexpr(k,length(J));
third_term2 = optimexpr(length(I),1); 

for j = 1:length(J)
    second_term_12 (:,j)= N(1,:)*ykj2(:,j);
end
second_term2 = sum(second_term_12(:,:)*C(:,1));

for i = 1:length(I)
third_term2(i,1) = b(1,i)*sum(xtij2(:,:,i)*B(:,1));
end
third_termf2 = sum(third_term2(:,1));

objective_function2 = first_term2 + second_term2 + third_termf2;
rebarprob2.Objective = objective_function2;
options = optimoptions(rebarprob2);
options.MaxTime = 200;

[sol2,fval2,exitflag2] = solve(rebarprob2,'Options',options);

%Show final solution and compute final total waste, cost, and cost saving
disp(sol);
disp(sol2.distancebetweenrebars);
disp(sol2.ifbeamhassizej);
disp(sol2.ifrebarissizej);
disp(sol2.numberofJbarscutwithpatternk);
disp(exitflag2);
finalptotalwaste = optimexpr(k,1);
for f = 1:k
        finalptotalwaste(f,1)= sum(ykj2(f,:))*Patternswaste(f,1); 
        column2 = finalptotalwaste(f,1);
 end
finaltotalwaste = sum(finalptotalwaste(:,1));
Finalcost = fval2;
Costsaving = fval - fval2;



