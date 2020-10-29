function lab1()
    %saveLinProg("Laby/TestFiles/linprog", 10, [2,2])
    
%     A = [1, 0, 1, 0; 0, 1, 0, 1]
%     b = [2; 2]
%     c = [1; 1]

    A = [2, 1, 1, 0, 0; 3, 3 ,0, 1, 0; 1.5, 0, 0, 0, 1];
    b = [10; 24; 6];
    c = [3; 2; 0; 0; 0];
    [result, status] = simplexSF(A, b, c);
end

function saveLinProg(path, max_it, demanded)
    if ~exist(path, 'dir')
       mkdir(path)
    end
    
    for it = 1:max_it
        var_num = 10
        cond_num = 5

        ALE = [];
        bLE = [];
        LB = [];
        c = [];

        ALE = randi([-5, 5], cond_num, var_num);
        ALE = [ALE, diag(ones(cond_num, 1))]
        c = randi([-5, 5], var_num, 1);
        c = [c; zeros(cond_num,1)]
        bLE = randi([1, 5], cond_num, 1)
        LB = zeros(var_num + cond_num, 1);

        options = optimset(@linprog);   
        options = optimset(options, 'Display', 'iter', 'Algorithm', 'dual-simplex');
        [x,fval,exitflag,output,lambda] = linprog(-c, [], [], ALE, bLE, LB, [], [], options)
        
        if exitflag == 1 & demanded(1) > 0
            demanded(1) = demanded(1) - 1
            fp = strcat(path, "/0_", int2str(demanded(1)), ".txt")
            fileID = fopen(fp,'w');
            writematrix(ALE(:, 1:var_num), fp,'Delimiter','tab')
            fprintf(fileID, "\n")
            writematrix(bLE, fp,'Delimiter','tab','WriteMode','append')
            fprintf(fileID, "\n")
            writematrix(c, fp,'Delimiter','tab','WriteMode','append')
            
        elseif exitflag == -3 & demanded(2) > 0
            demanded(2) = demanded(2) - 1  
            fp = strcat(path, "/-3_", int2str(demanded(2)), ".txt")
            fileID = fopen(fp,'a+');
            writematrix(ALE(:, 1:var_num), fp,'Delimiter','tab')
            fprintf(fileID, "\n")
            writematrix(bLE, fp,'Delimiter','tab','WriteMode','append')
            fprintf(fileID, "\n")
            writematrix(c, fp,'Delimiter','tab','WriteMode','append')
        end
        
        if demanded(1) == 0 & demanded(2) == 0 
            break
        end
        
    end
    fclose('all') 
end

function [result, status] = simplexSF(A, b, c) % algorytm simplex, który przyjmuje problem w postaci standardowej
    result = [];
    status = -1; %1 - ok, -1 - unlimited, 2 - more than one solution
    
    % krok 1
    
    [base, nonBase] = chooseBase(A);
    
    simplexTable = buildSimplexTable(A, b, c, base);
    
    while true
        % krok 2
        
        base = getBase(simplexTable);
        options = 1:length(c);
        nonBase = options(~ismember(options, base));
        
        OptRates = getOptRates(simplexTable);

        % krok 3 - potencjalny koniec

        isOptimal = all(OptRates>=0);
        
        if isOptimal
            break
        end

        % krok 4 - potencjalny koniec
        
        indWithNegativeOptR = find(OptRates < 0);
        
        A_wave = getAWave(simplexTable);
        
        for i = 1:length(indWithNegativeOptR)
            unlim = all(A_wave(:, i) <= 0);
            if unlim
                status = -1;
                result = [];
                break
            end
        end
        
        % krok 5
        
        [rows, ~] = size(OptRates);
        indicies = Inf(rows, 1);
        for i = nonBase
            indicies(i) = OptRates(i);
        end
        
        [~, enterIndex] = min(indicies);
        
        % krok 6
        
        b_wave = getBWave(simplexTable);
        
        [rows, columns] = size(A_wave);
        exitCritMatrix = Inf(rows, columns);
        
        base = getBase(simplexTable);
        
        
        % coś tu nie działa
        
        k = enterIndex;
        
        enteringColumn = A_wave(:,k)
        
        indices = find(enteringColumn > 0);
        
        leaveScores = Inf(1, length(b_wave));
        
        leaveScores(indices) = b_wave(indices) ./ enteringColumn(indices);
        
        [~, exitIndex] = min(leaveScores);
        
        %base(end + 1) = enterIndex;
        %base = base(base~=exitIndex);
        
        %nonBase = nonBase(nonBase~=enterIndex);
        %nonBase(end + 1) = exitIndex;
        
        %base = sort(base);
        %nonbase = sort(nonBase);
        
        % krok 7
        
        r = exitIndex;
        
        b_wave_p = b_wave;
        for i = 1:length(b_wave_p)
            if i == r
                b_wave_p(i) = b_wave(i) / A_wave(i, k);
            else
               b_wave_p(i) = b_wave(i) - A_wave(i, k) / A_wave(r, k) * b_wave(r);
            end
        end
        
        % krok 8
        
        A_wave_p = A_wave;
        
        for i = 1:size(A_wave_p, 1)
            for j = 1:size(A_wave_p, 2)
                if i == r
                    A_wave_p(i,j) = A_wave(i, j) / A_wave(i, k);
                else
                    A_wave_p(i,j) = A_wave(i, j) - A_wave(r, j) / A_wave(r, k) * A_wave(i, k);
                end
            end
        end
        
        
        
        
        
        % krok 9 - wróć do 2
        
        simplexTable = updateSimplexTable(simplexTable, A_wave_p, b_wave_p, k, r)
        
        
    end
    status = 0;
end
% 
% function exitInd = getIndToExit(simplexTable)
%     [r, c] = size(oldSimplexTable);
%     
%     
%     
% end

function simplexTable = updateSimplexTable(oldSimplexTable, A_wave, b_wave, enterInd, exitInd)
    [r, c] = size(oldSimplexTable);
    simplexTable = oldSimplexTable;
    simplexTable(2 + exitInd, 2) = enterInd;
    simplexTable(2 + exitInd, 1) = simplexTable(1, 2 + enterInd);
    simplexTable(3:r-2, 3:c-1) = A_wave;
    simplexTable(3:r-2, c) = b_wave;
    x = zeros(c-3, 1);
    x(simplexTable(3:r-2, 2)) = b_wave;
    base = simplexTable(3:r-2, 2);
    base = sort(base);
    c_b = simplexTable(3:r-2, 1);
    z = c_b.' * A_wave;
    simplexTable(r-1, 3:c-1) = z;
    simplexTable(r, 3:c-1) = z - simplexTable(1, 3:c-1);
end

function OptRates = getOptRates(simplexTable)
    [r, c] = size(simplexTable);
    simplexTable(r, 3:c-1) = simplexTable(r - 1, 3:c-1) - simplexTable(1, 3:c-1);
    OptRates = simplexTable(r, 3:c-1);
end

function A_wave = getAWave(simplexTable)
    [r, c] = size(simplexTable);
    A_wave = simplexTable(3:r-2, 3:c-1);
end

function b_wave = getBWave(simplexTable)
    [r, c] = size(simplexTable);
    b_wave = simplexTable(3:3+2, c);
end

function ind = getIndexOfExitingInBWave(simplexTable, variable)
    [r, c] = size(simplexTable);
    b = simplexTable(3:r-2, 2);
    ind = find(b==variable);
end

function base = getBase(simplexTable)
    [r, c] = size(simplexTable);
    base = simplexTable(3:r-2 ,2).';
end

function A_wave = calculateAWave(base, A)
    A_b = A(:,base);
    A_b_i = inv(A_b);
    A_wave = A_b_i * A;
end

function [base, nonBase] = chooseBase(A)
    base = [3,4,5];
    nonBase = [1,2];
end

function [Z, OptRates] = calculateOptRate(A_wave, c_b, c)
    Z = c_b.' * A_wave;
    Z = Z.';
    OptRates = Z - c;
end

function simplexTable = buildSimplexTable(A, b, c, base)
    simplexTable = NaN(size(b, 1) + 4, size(A,2) + 3);
    
    m = size(base, 2);
    n = size(c,1);
    mOffset = m - 1;
    nOffset = n - 1;
    
    A_b = A(:, base);
    c_b  = c(base,1);
    A_wave = calculateAWave(base, A);
    b_wave = inv(A_b) * b;
    x = zeros(size(c, 1),1);
    x(base) = b_wave;
    
    simplexTable(3:3 + mOffset,1) = c_b;
    simplexTable(3:3 + mOffset,2) = base.';
    simplexTable(1, 3:3 + nOffset) = c.';
    simplexTable(2, 3:3 + nOffset) = 1:size(c,1);
    simplexTable(3:3 + mOffset, 3:3 + nOffset) = A_wave;
    simplexTable(2 + m + 1, 3:3 + nOffset) = c.' * x;
    simplexTable(2 + m + 2, 3:3 + nOffset) = simplexTable(2 + m + 1, 3:3 + nOffset) - c.';
    simplexTable(3:3+mOffset, 2 + n + 1) = b_wave;
end
















