clc
% Defina o diretório base para vetor1
baseDir1 = 'D:\ERICK\MD\extracted_voxels';

% Defina o diretório base para vetor2
baseDir2 = 'D:\ERICK\MD\Extracted salient events (ORIGINAL)\Extracted salient events';

% Liste os diretórios dos sujeitos em baseDir1
subjectDirs1 = dir(baseDir1);
subjectDirs1 = subjectDirs1([subjectDirs1.isdir]); % Seleciona apenas pastas

% Liste os diretórios dos sujeitos em baseDir2
subjectDirs2 = dir(baseDir2);
subjectDirs2 = subjectDirs2([subjectDirs2.isdir]); % Seleciona apenas pastas

% Verifique se o número de pastas em baseDir1 é igual ao número de pastas em baseDir2
if length(subjectDirs1) ~= length(subjectDirs2)
    error('O número de pastas em baseDir1 e baseDir2 não coincide.');
end

% Loop pelos sujeitos
for subjectIdx = 4:7 % Começa em 3 para ignorar '.' e '..'
    subjectDir1 = fullfile(baseDir1, num2str(subjectIdx));
    disp(subjectDir1)
    
    % Obtém o nome do sujeito atual em baseDir2
    subjectName2 = strcat('subj_', num2str(subjectIdx));
    disp(subjectName2)
    
    % Construa o caminho completo para o arquivo 'rBA1' em baseDir1
    rBA1File1 = fullfile(subjectDir1, 'rBA1');
    
    % Construa o caminho completo para o arquivo 'rBA1' em baseDir2
    rBA1File2 = fullfile(baseDir2, subjectName2, 'rBA1');
    
    % Carregue o arquivo 'rBA1' de baseDir1 como vetor1
    load(rBA1File1, 'subjdata');
    
    % Carregue o arquivo 'rBA1' de baseDir2 como vetor2
    load(rBA1File2, 'subjData');
    
    % Defina o vetor3 para o sujeito atual
    [numRuns, numTrials] = size(subjData);
    vetor3 = zeros(numRuns, numTrials); % Inicialize vetor3 com zeros
    
    % Loop através de todas as execuções (runs) e trials
    for run = 1:numRuns
        for trial = 1:numTrials
            % Compare vetor1 com vetor2 e atribua 0 ou 1 a vetor3
            if isequal(subjdata{run}.voxdiff{trial}, subjData{run, trial}.signedDifference)
                vetor3(run, trial) = 0;
            else
                vetor3(run, trial) = 1;
            end
        end
    end
    
    % Salve o vetor3 em um arquivo CSV
    csvFileName = fullfile(subjectDir1, 'time_corr.csv');
    csvwrite(csvFileName, vetor3);
    
    disp(['Vetor3 para o sujeito ' num2str(subjectIdx) ' salvo em ' csvFileName]);
end
