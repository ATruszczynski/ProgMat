function lab4()
%     linoprogTest(10)


    verbose = false;
    
%     %    single solution
%     
%     A = [1, 0; 0, 1];
%     b = [2; 2];
%     c = [1; 1];
%     ub = [10; 10];
%     [ZPx, ZDy, ZDflag, fval] = dual_symplex(A, b, c, ub, verbose)
%     
%     disp('----------------------------------------')


% %    two solutions
%     
%     A = [1, 0; 0, 1; 1, 1;];
%     b = [2; 2; 3;];
%     c = [1; 1;];
%     ub = [10; 10];
%     
%     A = [1,1;];
%     b = [3;];
%     c = [1; 1;];
%     ub = [2; 2];
%     
%     [ZPx, ZDy, ZDflag, fval] = dual_symplex(A, b, c, ub, verbose)
%     
%     disp('----------------------------------------')
% 
% %     three solutions
%    
%     A = [1, 1, 1];
%     b = [1;];
%     c = [1; 1; 1];
%     ub = [10; 10; 10;];
%     [ZPx, ZDy, ZDflag, fval] = dual_symplex(A, b, c, ub, verbose)
%     
%     disp('----------------------------------------')
%     
% %   pierwszy przykład z wykładu
%     A = [2, 1; 3, 3; 1.5, 0];
%     b = [10; 24; 6;];
%     c = [3; 2];
%     ub = [100; 100;];
%     [ZPx, ZDy, ZDflag, fval] = dual_symplex(A, b, c, ub, verbose)
    
    disp('----------------------------------------')
%     
% %   drugi przykład z wykładu
%     A = [1, 2, -1, -1; 2, -2, 3, -3;];
%     b = [4; 9;];
%     c = [3; 1; 3; -1;];
%     ub = [100; 100;];
%     [ZPx, ZDy, ZDflag, fval] = dual_symplex(A, b, c, ub, verbose)
%     
%     disp('----------------------------------------')

% %     unbound
%     
%     A = [-1, 1; 1, -1];
%     b = [2; 2];
%     c = [1; 1];
%     ub = [];
%     [ZPx, ZDy, ZDflag, fval] = dual_symplex(A, b, c, ub, verbose)
%     
%     disp('----------------------------------------')


     testDualSimplex(10000)
end

function linoprogTest(n)
    var_num = 5;
    cond_num = 10;
    rng(1001)
    
    for it = 1:n

        ALE = [];
        bLE = [];
        LB = [];
        c = [];

        ALE = randi([-5, 5], cond_num, var_num);
        ALE = [ALE, diag(ones(cond_num, 1))];
        c = randi([-5, 5], var_num, 1);
        c = [c; zeros(cond_num,1)];
        bLE = randi([1, 5], cond_num, 1);
        LB = zeros(var_num + cond_num, 1);
        UB = randi([1,30], var_num + cond_num, 1);

        options = optimset(@linprog);   
        options = optimset(options, 'Display', 'iter', 'Algorithm', 'dual-simplex');
        [x,fval,exitflag,output,lambda] = linprog(c, [], [], ALE, bLE, LB, UB, [], options)
    end
end

function testDualSimplex(n)
    rng(1001)
    
    confusion_matrix = zeros(2,2);
    error_list = [];
    
    for i = 1:n
        disp(i)
        var_num = randi([3,7]);
        cond_num = randi([5, 10]);
        A = randi([-5, 5], cond_num, var_num);
        c = randi([-5, 5], var_num, 1);
        b = randi([1, 5], cond_num, 1);
        lb = zeros(var_num, 1);
        ub = randi([1,30], var_num, 1);

        options = optimset(@linprog);   
        options = optimset(options, 'Display', 'off', 'Algorithm', 'dual-simplex');
        [x,fval,exitflag,output,lambda] = linprog(-c, A, b, [], [], lb, ub, [], options);
        fval = -fval;
        [ZPx, ZDy, ZDexitflag, funcVal] = dual_symplex(A, b, c, ub, false);
        
        prec = 1e-10;
        
        if isempty(fval)
            throw("Found inconsistent example")
        end
        
        if roughEquality(fval, funcVal, prec) % rozwiązania różnią się minimalnie, prawdopodobnie przez zaokrąglenia. roughEquality porównuje dwa wektory z pewną tolerancją dla różnic (ustaloną przez prec)
            cm_x = 1;
        else
            cm_x = 2;
        end
        
        if roughEquality(x, ZPx.', prec)
            cm_y = 1;
        else
            cm_y = 2;
        end
        
        if cm_x == 2
            error_list(end + 1) = i;
        end
        
        confusion_matrix(cm_x, cm_y) = confusion_matrix(cm_x, cm_y) + 1;
    end
    
    confusion_matrix
end

function isSame = roughEquality(x,y,e)
    isSame = true;
    for i = 1:length(x)
        if abs(x(i) - y(i)) > e
            isSame = false;
            break;
        end
    end
end

function [ZPx, ZDy, ZDexitflag, funcVal] = dual_symplex(A, b, c, ub, verbose)
    [row,col] = size(A);
    length(ub);
    A = [A; eye(length(ub))];
    b = [b; ub];
    
    [dA, db, dc] = transformProblemIntoDual(A, b, c);
    [dA, db, dc, dualAddedVar] = addComplimentaryVariables(dA, db, dc);
    
    [y, dFuncVal, status, message, y_d] = internal_dual_simplex(dA, db, dc, dualAddedVar, verbose);
    
    options = optimset(@linprog);   
    options = optimset(options, 'Algorithm', 'dual-simplex');
    lb = zeros(1, size(dc, 1));
%     [x,fval,exitflag,output,lambda] = linprog(-dc, [], [], -dA, -db, lb, [], options);
    
    dFuncVal = -dFuncVal;
    
    if status == 1
        [ZPx, ZDy] = calculateSolutions(A, b, y, y_d);
        ZDexitflag = 1;
        funcVal = dFuncVal;
    else
        ZPx = [];
        ZDy = [];
        ZDexitflag = 0;
        funcVal = [];
    end
end

function [ZPx, ZDy] = calculateSolutions(A, b, y, y_d)
    primalAddedVar = size(A, 1);
    Ap = [A, eye(primalAddedVar)];
    
    toRemove = [];
    
    for i = 1:length(y_d)
        if y_d(i) > 0
            toRemove(end+1) = i;
        end
    end
    
    for i = 1:length(y)
        if y(i) > 0
            toRemove(end + 1) = size(A,2) + i;
        end
    end
    
    kept = 1:size(Ap,2);
    kept = kept(~ismember(kept, toRemove));
    
    Ap = Ap(:, kept);
    
    sol = linsolve(Ap, b);
    
    ZPx = zeros(1,size(A,2)+primalAddedVar);
    ZPx(kept) = sol;
    ZDy = y;
    ZPx = ZPx(1:size(A,2));
end


function [dA, db, dc] = transformProblemIntoDual(A, b, c)
    dA = A.';
    db = c;
    dc = -b;
end

function [dA, db, dc, addedVars] = addComplimentaryVariables(A, b, c)
    [row, col] = size(A);
    addedVars = row;
    
    dA = [A, -eye(addedVars)];
    dc = [c; zeros(addedVars, 1)];
    db = b;
end

function [result, funcVal, status, message, x_d] = internal_dual_simplex(A, b, c, complVarCount, verbose) % algorytm simplex, który przyjmuje problem w postaci standardowej
    %status = 1 - ok; -1 - unlimited; 2 - more than one solution; -2 - inconsistent; 0 -
    %algorithm is ongoing (if this appears on output something went wrong)
    trueVariable = size(A,2) - complVarCount;
    
    % krok 1

    base = attemptChooseBase(A);
    
    base = trueVariable + 1: size(A,2);
    
    simplexTable = buildSimplexTable(A, b, c, base);
    
    % faza 2
    
    [simplexTable, status, negIndex] = internalDualSimplexLoop(simplexTable, verbose);
    
    % przypisz odpowiedni status, wyniki i wiadomość dla pojedynczego
    % rozwiązania
    if status == 0
        status = 1;
        message = 'Single solution found';
        result = getResult(simplexTable);
        x_d = result(trueVariable + 1:size(A,2));
        result = result(1:trueVariable);
        funcVal = getFuncVal(simplexTable);
        % sprawdź czy istnieje drugie rozwiązanie
        % jeśli tak, funkcja od razu zwraca indeksy zmiennej do dodadania do
        % bazy i zmiennej do usunięcia z bazy
%         [enterIndicies, exitIndcicies] = checkForAnotherSolution(simplexTable);
%         if length(enterIndicies) > 0
%             message = [];
%             for i = 1:length(enterIndicies)
%                 enterIndex = enterIndicies(i);
%                 exitIndex = exitIndcicies(i);
% 
%                 % przelicz drugie rozwiązanie powtarzajac kroki 7-8
%                 simplexTable = substituteBase(simplexTable, enterIndex, exitIndex);
% 
%                 % przypisz odpowiedni status, wyniki i wiadomość dla dwóch
%                 % rozwiązań
% 
%                 result2 = getResult(simplexTable);
%                 result2 = result2(1:trueVariable); 
%                 funcVal2 = getFuncVal(simplexTable);
%                 message = append(num2str(message), ['Solution ' num2str(1 + i) ' found at [' num2str(result2) ']. The value of function at this point is ' num2str(funcVal2) newline]);
%                 
%                 if verbose
%                     disp(['Simplex table for solution ' num2str(1 + i)]);
%                     disp(simplexTable);
%                 end
%             end
%         end
    elseif status == -1
        x_d=[];
        result = [];
        funcVal = [];
        message = calculateExtremeRadius(negIndex, simplexTable);
    else
        x_d=[];
        result = [];
        funcVal = [];
        message = 'Inconsistent problem';
    end
end

% algorithm operations

function base = attemptChooseBase(A)
    [r,c] = size(A);
    
%     base = (c-r+1:c);
    
    base = [];
    for j = 1:c
        nonZeroes = 0;
        varInd = -1;
        for i = 1:r
           a = A(i, j);
           if a ~= 0
               nonZeroes = nonZeroes + 1;
           end
           if a == 1
               varInd = j;
           end
        end
        if nonZeroes == 1 & varInd ~= -1
            base(end + 1) = varInd;
        end
    end
    
    base = unique(base);
    bl = size(base,2);
    if bl >= r
        base = base(1:r);
    else
        base = NaN;
    end
end

function [simplexTable, status, negInd] = internalDualSimplexLoop(simplexTable, verbose)
    negInd = -1;
    iteration = 0;
    status = 0;
    
    while true
        iteration = iteration + 1;
        if verbose
            printState(simplexTable, iteration);
        end
        
        % krok 2
        
        if all(getBWave(simplexTable) >= 0) % sprawdzenie warunku optymalności
            break
        end

        % krok 3
        
        exitIndex = chooseExitIndex(simplexTable);
        
        % krok 4
        
        b_wave = getBWave(simplexTable);
        if b_wave(exitIndex) >= 0
            throw("That should not happen")
        end
        
        A_wave = getAWave(simplexTable);
        rrow = A_wave(exitIndex,:);
        nonBase = getNonBase(simplexTable);
        rrow = rrow(nonBase); % Is that ok?
        
        if all(rrow >= 0)
            status = -2;
            break
        end
        
        % krok 5
        
        enterIndex = chooseEnterIndex(simplexTable, exitIndex);
        
        % krok 6
        
        simplexTable = substituteBase(simplexTable, enterIndex, exitIndex);
        
%         indiciesWithNegativeOptRates = find(optRates < 0); % indeksy, dla których należy sprawdzić warunek nieograniczności rozwiązania
%         
%         A_wave = getAWave(simplexTable);
%         
%         for i = indiciesWithNegativeOptRates
%             unlim = all(A_wave(:, i) <= 0);
%             if unlim
%                 status = -1;
%                 negInd = i;
%                 break
%             end
%         end
%            
%         if status ~= 0
%             break
%         end
        
    end
end

function [simplexTable, status] = phase1Computation(A, b, c)
    status = 0;
    [artVarCount, trueVarCount] = size(A)
    
    A = [A, eye(artVarCount)]
    c = [c; zeros(artVarCount,1)]
    
    artVars = trueVarCount + 1: trueVarCount + artVarCount
    
    simplexTable = buildSimplexTable(A, b, c, artVars)
    [simplexTable, status, negInd] = internalDualSimplexLoop(simplexTable, false)
    
    finished = all(ismember(artVars, getNonBase(simplexTable)))
    finished = false
    if finished
        return
    else
        result = getResult(simplexTable)
        toSearch = intersect(getBase(simplexTable), artVars)
        find(ismember())
    end
end

function enterIndex = chooseEnterIndex(simplexTable, exitIndex) % wyznacz indeks zmiennej wchodzącej do bazy
    optRates = getOptRates(simplexTable);
    A_wave = getAWave(simplexTable);
    nonBase = getNonBase(simplexTable);
    [~, cols] = size(optRates);
    indicies = -Inf(1, cols);
    for i = nonBase
        a = A_wave(exitIndex, i);
        if a < 0
            indicies(i) = optRates(i) / a;
        end
    end
        
    [~, enterIndex] = max(indicies);
end

function exitIndex = chooseExitIndex(simplexTable)
    b_wave = getBWave(simplexTable);
    [~, exitIndex] = min(b_wave);
end

function newSimplexTable = substituteBase(simplexTable, enterIndex, exitIndex)
    % krok 7

    A_wave = getAWave(simplexTable);
    b_wave = getBWave(simplexTable);

    k = enterIndex; % zmiana nazw zmiennych na te ze slajdów dla czytelności
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

    newSimplexTable = updateSimplexTable(simplexTable, A_wave_p, b_wave_p, k, r);
end

function [enterIndicies, exitIndicies] = checkForAnotherSolution(simplexTable)
    exitIndicies = [];
    optRates = getOptRates(simplexTable);
    
    nonBase = getNonBase(simplexTable);
    
    enterIndicies = nonBase(optRates(nonBase) == 0);
    
    for i = enterIndicies
        exitIndicies(end + 1) = chooseExitIndex(simplexTable, i);
    end
end

function message = calculateExtremeRadius(negInd, simplexTable)
    variableCount = size(simplexTable,2) - 3;
    base = getBase(simplexTable);
    A_wave = getAWave(simplexTable);
    
    x = getResult(simplexTable);
    
    % obliczenie wektora kiernkowego
    v = zeros(1, variableCount);
    v(base) = -A_wave(:,negInd);
    v(negInd) = 1;
    
    message = ['Solution is unbound. Extreme radius found: [' num2str(x) '] + a * [' num2str(v) ']'];
end

% gets

% funkcje wyciągające z tabelki sympleksu odpowiednie zmienne zapisane w
% komórkach. Ich użycie przypomina użycie akcesorów w programowaniu
% obiektowym

function base = getBase(simplexTable)
    [r, c] = size(simplexTable);
    base = simplexTable(3:r-2 ,2).';
end

function nonBase = getNonBase(simplexTable)
    [~,c] = size(simplexTable);
    
    % w nie-bazie jest dopełnienie bazy
    base = getBase(simplexTable);
    options = 1:c-3;
    nonBase = options(~ismember(options, base));
end

function A_wave = getAWave(simplexTable)
    [r, c] = size(simplexTable);
    A_wave = simplexTable(3:r-2, 3:c-1);
end

function b_wave = getBWave(simplexTable)
    [r, c] = size(simplexTable);
    b_wave = simplexTable(3:r-2, c);
end

function optRates = getOptRates(simplexTable)
    [r, c] = size(simplexTable);
    optRates = simplexTable(r, 3:c-1);
end

function result = getResult(simplexTable)
    [r,c] = size(simplexTable);
    base = getBase(simplexTable);
    b_wave = getBWave(simplexTable);
    result = zeros(1, c-3);
    result(base) = b_wave;
end

function funcVal = getFuncVal(simplexTable)
    [r,c] = size(simplexTable);
    funcVal = simplexTable(r-1, c);
end

% building simplex

function simplexTable = buildSimplexTable(A, b, c, base)
    % prawie wszystkie informacje w algorytmie są przechowywane w tabeli
    % sympleksowej, która ma taka samą strukturę jak na wykładzie, przy
    % czym ma dodane kilka NaNów, żeby mieć prostokątny kształt

    simplexTable = NaN(size(b, 1) + 4, size(A,2) + 3); % NaN są traktowane jak puste miejsca
    
    [row, col] = size(simplexTable);
    
    m = size(base, 2);
    n = size(c,1);
    mOffset = m - 1;
    nOffset = n - 1;
    
    % inicjalizacja elementów tabeli sympleksowej
    A_b = A(:, base);
    c_b  = c(base,1);
    A_wave = inv(A_b) * A;
    b_wave = inv(A_b) * b;
    x = zeros(size(c, 1),1);
    x(base) = b_wave;
    
    % przypisanie elementów tabeli sympleksowej w odpowiednich miejscach
    simplexTable(3:3 + mOffset,1) = c_b;
    simplexTable(3:3 + mOffset,2) = base.';
    simplexTable(1, 3:3 + nOffset) = c.';
    simplexTable(2, 3:3 + nOffset) = 1:size(c,1);
    simplexTable(3:3 + mOffset, 3:3 + nOffset) = A_wave;
    simplexTable(2 + m + 1, 3:3 + nOffset) = c.' * x;
    simplexTable(2 + m + 2, 3:3 + nOffset) = simplexTable(2 + m + 1, 3:3 + nOffset) - c.';
    simplexTable(3:3+mOffset, 2 + n + 1) = b_wave;
    simplexTable(row-1, col) = c.' * x;
end

function simplexTable = updateSimplexTable(oldSimplexTable, A_wave, b_wave, enterInd, exitInd)
    [r, c] = size(oldSimplexTable);
    
    simplexTable = oldSimplexTable;
    
    % podmiana zmiennej w bazie
    simplexTable(2 + exitInd, 2) = enterInd;
    simplexTable(2 + exitInd, 1) = simplexTable(1, 2 + enterInd);
    
    % aktualizacja tablic A_wave i b_wave
    simplexTable(3:r-2, 3:c-1) = A_wave;
    simplexTable(3:r-2, c) = b_wave;
    
    % aktualne rozwiązanie
    x = getResult(simplexTable);
    
    % aktualizacja z i z-c
    c_b = simplexTable(3:r-2, 1);
    z = c_b.' * A_wave;
    simplexTable(r-1, 3:c-1) = z;
    simplexTable(r, 3:c-1) = z - simplexTable(1, 3:c-1);
    
    % aktualizacja wartości funkcji
    cc = simplexTable(1, 3:c-1);
    simplexTable(r-1,c) = cc * x.';
end

% other

function printState(simplexTable, iteration)
    message = ['Iteration ' num2str(iteration) '. Current varaibles in base [' num2str(sort(getBase(simplexTable))) '].'];
    disp(message);
    disp(simplexTable);
end


