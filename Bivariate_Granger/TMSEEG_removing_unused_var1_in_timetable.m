clear all;
close all;
clc;

IDlist = {'301', '302', '305', '306', '307', '308', '309', '310',...
    '311', '312', '314', '315', '317', '318', '319', '320', '322',...
    '324', '325', '326', '327', '328', '329', '330', '333', '335',...
    '336', '338', '341', '343', '345', '346', '347', '348', '351',...
    '102', '104', '105', '106', '107', '108', '109', '110', '111',...
    '112', '114', '115', '116', '117', '118', '119', '120', '121'};  

cd('E:\Bivariate_Granger_Loop_Trial\Granger_Causality_TMS_EEG\Model_Order_10\');
for fi = 1:length(IDlist)
    folto = ['E:\Bivariate_Granger_Loop_Trial\Granger_Causality_TMS_EEG\Model_Order_10\' IDlist{fi} '\'];
    cd(folto);
    dirList = dir;
    dirList(1:2) = [];
        for fi2 = 1:length(dirList);
        fto = dirList(fi2).name;
        T = readtable(fto);
        T = removevars(T,{'Var1'});
        writetable(T,fto);
        clear T; clear newfn; 
        end
end
        
    