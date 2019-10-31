%% FANUC R-2000iC/165F

% lengthes are in mm
l = [670, 312, 1075, 225, 1280, 215];

T_base_real = Tx(1);
T_tool_real = Tz(5);

errors_real = [1 2 1 0 1 1.5 1 1 2 1 1 1.5 1 1 1.7 1 1 1 2 1];

%% Data generation
% TODO: 3 points instead of one

c = 30;
q_gen = [rand(1,c)*pi - pi/2; rand(1,c)*pi - pi/2; rand(1,c)*pi - pi/2; ...
        rand(1,c)*pi - pi/2; rand(1,c)*pi - pi/2; rand(1,c)*pi - pi/2];
disp('Generated joint angles');
disp(q_gen);

p_gen = zeros(30,3);
for i = 1:c
    T = T_base_real * FK_irreducible(q_gen(:,i),l,errors_real) * T_tool_real;
    p_gen(i,:) = T(1:3,4);
end
disp('Generated end effector positions')
disp(p_gen);

%% Identification

errors = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

for k = 1:20

    % Step 1

    %Sum_aTa = zeros(9);
    %Sum_aTdp = zeros(9,1);
    Sum_aTa = zeros(9);
    Sum_aTdp = zeros(9,1);
    for i = 1:c
        % p robot and R robot
        T_robot = FK_irreducible(q_gen(:,i),l,errors);
        p_robot = T(1:3,4);
        % A matrix and delta p
        tilda_p = Scew_Symm(p_robot);
        A(1:3,1:3) = eye(3);
        A(1:3,4:6) = tilda_p';
        A(1:3,7:9) = T_robot(1:3,1:3);
        dp = p_gen(i,:) - p_robot';
        % edit sums aTa and aTdp
        Sum_aTa = Sum_aTa + (A'*A);
        Sum_aTdp = Sum_aTdp + (A'*dp');
    end
    % T base and T tool
    solution = Sum_aTa\Sum_aTdp;
    %disp('p base, r base, u tool');
    %disp(solution);
    T_base = [eye(3), solution(1:3); 0, 0, 0, 1];
    T_tool = [eye(3), solution(7:9); 0, 0, 0, 1];

    % Step 2

    Observation = [];
    dp = [];

    for i = 1:c
        Observation = [Observation (J_p(q_gen(:,i), l, errors, T_base, T_tool))'];
        %disp('Observation');
        %disp(Observation);
        T_i = T_base * FK_irreducible(q_gen(:,i), l, errors) * T_tool;
        p = T_i(1:3,4);
        dp = [dp (p_gen(i,:)-p')];
    end

    J_arg = Observation';
    d_errors = ((J_arg'*J_arg)\(J_arg'*dp'))';

    % disp('Delta errors');
    % disp(d_errors); 

    errors = errors + d_errors;

end

disp('Errors');
disp(errors);





















