clear all;
close all;
clc;
% alz n = 40, controls = 18
IDlist = {'301', '302', '303', '305', '306', '307', '308', '309', '310',...
    '311', '312', '313', '314', '315', '316', '317', '318', '319', '320',...
    '322', '324', '325', '326', '327', '328', '329', '330', '331', '333',...
    '335', '336', '338', '341', '343', '345', '346', '347', '348', '349',...
    '351', '101', '102', '103', '104', '105', '107', '108', '109', '110',...
    '111', '112', '114', '115', '116', '117', '119', '120', '121'}; 


Channames = {'AF3', 'AF4', 'F7', 'F5', 'F1', 'FZ', 'F2', 'F4', 'F6', 'F8',...
    'FC5', 'FC3', 'FC1', 'FCZ', 'FC2', 'FC4', 'FC6', 'C5', 'C3', 'C1', 'CZ',...
    'C2', 'C4', 'C6', 'P7', 'P5', 'P3', 'P1', 'PZ', 'P2', 'P4', 'P6', 'P8',...
    'PO3', 'POZ', 'PO4', 'O1', 'OZ', 'O2'};
%% 
cd('E:\Bivariate_Granger_Loop_Trial\Granger_Causality_RSEEG\model_order_10\');
for fi = 1:length(IDlist)
    folto = ['E:\Bivariate_Granger_Loop_Trial\Granger_Causality_RSEEG\model_order_10\' IDlist{fi} '\'];
    cd(folto);
    dirList = dir;
    dirList(1:2) = [];
        
    GC.(strcat('A',IDlist{fi})) = [zeros(39)];
    
    for fi2 = 1:length(dirList);
        %the following lines look for the chan name and find it's file
        %index in dirList.name
        %reshape(strfind({dirList. name}, Channames{fi2}), size(dirList));
        Desc   = {dirList.name};
        Desc(~cellfun('isclass', Desc, 'char')) = {''};  % Care for non-strings
        matchC = reshape(strfind(Desc, Channames{fi2}), size(dirList));
        match  = ~cellfun('isempty', matchC);
        fto = dirList(match).name;
        T = readtable(fto);
        TT = removevars(T,{'Time'});
        TT = table2array(TT);
        TT(:,2:2:end) = [];
        TT = mean(TT); % mean of each column - now a row vector 
        if fi2 == 1
        GCMat(fi2,:) = [0 TT(1:end)];
        else GCMat(fi2,:) = [TT(1:fi2-1) 0 TT(fi2:end)];
        end       
        
        clear T;
        clear TT;
        clear match;
        clear matchC;
        clear Desc;
    end
   
GC.(strcat('A',IDlist{fi})) = GCMat;
GCNBS(:,:,fi) = GCMat;
clear GCMat;
end
%has to be channels x channels x subjects

cd('E:\Bivariate_Granger_Loop_Trial\Granger_Causality_RSEEG\Granger_NBS\Model_10');
save('Granger_RSEEG_DATASTRUCT_Model_10.mat','GC');
save('Granger_RSEEG_NBS_Model_10.mat','GCNBS');
