function [val] = getparval(cfgcell,parstring,defval,indn)
%% getparval ==============================================================
% cfgcell: n-by-2 cell array containing parameters
% parstring: name of desired parameter (without colon)
% defval: default val if parameter is not found
% indn: index to use when multiple matches are found
% Andrew Watson @ Leeds, 15/06/2021
% Modified: make tolerant to trailing colon in parstring or cfgcell
% ========================================================================
    if nargin < 4
        indn = 1;
    end
    
    try
        % normalize both sides: strip trailing colon(s)
        keys = regexprep(cfgcell(:,1), ':$', '');
        query = regexprep(parstring, ':$', '');
        
        indx = find(strcmpi(keys, query));  % case-insensitive match
        if isempty(indx)
            val = defval;
        else
            val = cfgcell{indx(indn),2};
        end
    catch
        val = defval;
    end
end
