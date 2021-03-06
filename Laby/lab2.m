function lab1()
%     linoprogTest(5) % wywołuje linprog tyle razy, ile zostanie podane jako argument

    verbose = false;

% %    single solution
%     
%     A = [1, 0, 1, 0; 0, 1, 0, 1];
%     b = [2; 2];
%     c = [-1; -1; 0; 0;];
%     [result, funcVal, status, message] = simplex(A, b, -c, 2, verbose)
%     
%     disp('----------------------------------------')
% 
% %     unbound
%     
%     A = [-1, 1, 1, 0; 1, -1, 0, 1];
%     b = [2; 2];
%     c = [-1; -1; 0; 0;];
%     [result, funcVal, status, message] = simplex(A, b, -c, 2, verbose)
%     
%     disp('----------------------------------------')
% 
% %    two solutions
%     
%     A = [1, 0, 1, 0, 0; 0, 1, 0, 1, 0; 1, 1, 0, 0, 1;];
%     b = [2; 2; 3;];
%     c = [-1; -1; 0; 0; 0;];
%     [result, funcVal, status, message] = simplex(A, b, -c, 3, verbose)
%     
%     disp('----------------------------------------')

% %     three solutions
%    
%     A = [1, 1, 1, 1;];
%     b = [1;];
%     c = [-1; -1; -1; 0;];
%     [result, funcVal, status, message] = simplex(A, b, -c, 1, verbose)
%     
%     disp('----------------------------------------')
% %     
%   pierwszy przykład z wykładu
%     A = [2, 1, 1, 0, 0; 3, 3, 0, 1, 0; 1.5, 0, 0, 0, 1;];
%     b = [10; 24; 6;];
%     c = [3; 2; 0; 0; 0;];
%     [result, funcVal, status, message] = simplex(A, b, c, 3, verbose)
    
    disp('----------------------------------------')
%     
% %   drugi przykład z wykładu
%     A = [1, 2, -1, -1, 1, 0; 2, -2, 3, -3, 0, 1;];
%     b = [4; 9;];
%     c = [3; 1; 3; -1; 0; 0;];
%     [result, funcVal, status, message] = simplex(A, b, c, 2, verbose)
    
%     disp('----------------------------------------')
% 
% %   Beale
%     A = [1, 0, 0, 1/4, -8, -1, 9; 0, 1, 0, 1/2, -12, -1/2, 3; 0, 0, 1, 0, 0, 1, 0];
%     b = [0; 0; 1;];
%     c = [0; 0; 0; 3/4; -20; 1/2; -6;];
%     [result, funcVal, status, message] = simplex(A, b, c, 0, verbose)
%     
%     disp('----------------------------------------')

%   No easy base example
    A = [1, 1, 0; -1, 1, 0;];
    b = [1; 0;];
    c = [0; 1; 0;];
    [result, funcVal, status, message] = simplex(A, b, c, 0, verbose)
    
    disp('----------------------------------------')
        
    
end

function linoprogTest(max_it)
    var_num = 10;
    cond_num = 5;
    rng(1001)
    for it = 1:max_it

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

        options = optimset(@linprog);   
        options = optimset(options, 'Display', 'iter', 'Algorithm', 'dual-simplex');
        [x,fval,exitflag,output,lambda] = linprog(c, [], [], ALE, bLE, LB, [], [], options);
    end
end

function [result, funcVal, status, message] = simplex(A, b, c, complVarCount, verbose) % algorytm simplex, który przyjmuje problem w postaci standardowej
    %status = 1 - ok; -1 - unlimited; 2 - more than one solution; -2 - inconsistent; 0 -
    %algorithm is ongoing (if this appears on output something went wrong)
    trueVariable = size(A,2) - complVarCount;
    
    % krok 1 / faza 1

    base = attemptChooseBase(A);
    
    if isnan(base)
         disp("No base")
         [simplexTable, status] = phase1Computation(A, b, c);
         return
    else
        simplexTable = buildSimplexTable(A, b, c, base);
    end
    
    % faza 2
    
    [simplexTable, status, negIndex] = internalSimplexLoop(simplexTable, verbose);
    
    % przypisz odpowiedni status, wyniki i wiadomość dla pojedynczego
    % rozwiązania
    if status == 0
        status = 1;
        message = 'Single solution found';
        result = getResult(simplexTable);
        result = result(1:trueVariable);
        funcVal = getFuncVal(simplexTable);

        % sprawdź czy istnieje drugie rozwiązanie
        % jeśli tak, funkcja od razu zwraca indeksy zmiennej do dodadania do
        % bazy i zmiennej do usunięcia z bazy
        [enterIndicies, exitIndcicies] = checkForAnotherSolution(simplexTable);
        if length(enterIndicies) > 0
            message = [];
            for i = 1:length(enterIndicies)
                enterIndex = enterIndicies(i);
                exitIndex = exitIndcicies(i);

                % przelicz drugie rozwiązanie powtarzajac kroki 7-8
                simplexTable = substituteBase(simplexTable, enterIndex, exitIndex);

                % przypisz odpowiedni status, wyniki i wiadomość dla dwóch
                % rozwiązań

                status = 2;
                result2 = getResult(simplexTable);
                result2 = result2(1:trueVariable); 
                funcVal2 = getFuncVal(simplexTable);
                message = append(num2str(message), ['Solution ' num2str(1 + i) ' found at [' num2str(result2) ']. The value of function at this point is ' num2str(funcVal2) newline]);
                
                if verbose
                    disp(['Simplex table for solution ' num2str(1 + i)]);
                    disp(simplexTable);
                end
            end
        end
    elseif status == -1
        result = [];
        funcVal = [];
        message = calculateExtremeRadius(negIndex, simplexTable);
    end
end

% algorithm operations

function base = attemptChooseBase(A)
    [r,c] = size(A);
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

function [simplexTable, status, negInd] = internalSimplexLoop(simplexTable, verbose)
    negInd = -1;
    iteration = 0;
    status = 0;
    
    while true
        iteration = iteration + 1;
        if verbose
            printState(simplexTable, iteration)
        end
        
        % krok 2
        
        optRates = getOptRates(simplexTable); % optRates = z - c

        % krok 3
        
        if all(optRates>=0) % sprawdzenie warunku optymalności
            break
        end

        % krok 4 - potencjalny koniec
        
        indiciesWithNegativeOptRates = find(optRates < 0); % indeksy, dla których należy sprawdzić warunek nieograniczności rozwiązania
        
        A_wave = getAWave(simplexTable);
        
        for i = indiciesWithNegativeOptRates
            unlim = all(A_wave(:, i) <= 0);
            if unlim
                status = -1;
                negInd = i;
                break
            end
        end
        
        if status ~= 0
            break
        end
        
        % krok 5
        
        enterIndex = chooseEnterIndex(simplexTable);
        
        % krok 6
        
        exitIndex = chooseExitIndex(simplexTable, enterIndex);
        
        % krok 7,8
        
        simplexTable = substituteBase(simplexTable, enterIndex, exitIndex);
    end
end

function [simplexTable, status] = phase1Computation(A, b, c)
    status = 0;
    [artVarCount, trueVarCount] = size(A)
    
    A = [A, eye(artVarCount)]
    c = [c; zeros(artVarCount,1)]
    
    artVars = trueVarCount + 1: trueVarCount + artVarCount
    
    simplexTable = buildSimplexTable(A, b, c, artVars)
    [simplexTable, status, negInd] = internalSimplexLoop(simplexTable, false)
    
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

function enterIndex = chooseEnterIndex(simplexTable) % wyznacz indeks zmiennej wchodzącej do bazy
    optRates = getOptRates(simplexTable);
    nonBase = getNonBase(simplexTable);
    [~, rows] = size(optRates);
    indicies = Inf(rows, 1);
    for i = nonBase
        indicies(i) = optRates(i);
    end
        
    [~, enterIndex] = min(indicies);
end

function exitIndex = chooseExitIndex(simplexTable, enterIndex)
    A_wave = getAWave(simplexTable);
    b_wave = getBWave(simplexTable);
        
    enteringColumn = A_wave(:,enterIndex);
    indices = find(enteringColumn > 0);

    leaveScores = Inf(1, length(b_wave)); % wektor infów; jeśli inf nie zostanie zastąpiony inną wartością, to jest zasadniczo ignorowany przez min niżej
    leaveScores(indices) = b_wave(indices) ./ enteringColumn(indices);

    [~, exitIndex] = min(leaveScores);
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
















