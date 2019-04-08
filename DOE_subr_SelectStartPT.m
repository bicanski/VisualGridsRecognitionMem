

% Code developed by Andrej Bicanski
% andrej.bicanski@gmail.com
% www.andrejbicanski.com
%
% model published in Current Biology
%
% Bicanski A, Burgess N. - A computational model of recognition 
% memory via grid cells. Current Biology, 2019, 29, 1–12. 
% DOI: 10.1016/j.cub.2019.01.077
%
% This script selects a fresh starting pt for the recognition process


function [new_startPT,prevSTpts] = DOE_subr_SelectStartPT(occ,startPT,IMPTS,xocc,yocc,prevSTpts)

count =1;

if occ == 0
    new_startPT = ceil(rand*9);
    while any(prevSTpts==new_startPT)
        new_startPT = ceil(rand*9);
    end
end

if occ == 1
    new_startPT     = ceil(rand*9);
    new_startCOORDs = IMPTS(new_startPT,:);
    while new_startCOORDs(1)>=xocc || any(prevSTpts==new_startPT)
        new_startPT     = ceil(rand*9);
        new_startCOORDs = IMPTS(new_startPT,:);      
        count=count+1;   
    end
end

if occ == 2
    new_startPT     = ceil(rand*9);
    new_startCOORDs = IMPTS(new_startPT,:);
    while new_startCOORDs(2)>=yocc || any(prevSTpts==new_startPT)
        new_startPT     = ceil(rand*9);
        new_startCOORDs = IMPTS(new_startPT,:);   
        count=count+1;
    end
end

prevSTpts = [prevSTpts new_startPT];
