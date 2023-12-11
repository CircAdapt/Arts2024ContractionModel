function Value= Get(Type,Variable,Name)
% function Value= Get(Type,Variable,Name)
% Type=string, Variable=string, Name=string or {str1,str2}
% e.g. Get('Node','p',{'Ao','BrAr')
% inverse action of Put
% if Name=='All' -> all variables
% if Variable=='Index' -> index number(s) of Name
% Theo Arts, Maastricht University, Mar 20, 2013

global P;
Name=[Name,{}]; % convert Name to cell array
T=[Type,{}]; % convert Type to cell array
Names= P.(T{1}).Name; % all Name's in P.Type
n    = length(Names); % number of Name's in P.Type.Name
iN=(1:n)';

if strcmp(Name,'All')
    iName=1:n;
else
    iName=[];
    for i=1:length(Name)
        dName= strcmp(Name{i},Names)*iN;
        % if non-existent, value 0 is inserted
        iName= [iName,dName];
    end
end
% iName is array of indices

if strcmp(Variable,'Index')
    Value=iName; % zero if no name found
else
    iName=iName(logical(iName~=0)); % select only existing variable
    switch length(T)
        case 1
        Value= P.(T{1}).(Variable)(:,iName); % depth=1
        case 2
        Value= P.(T{1}).(T{2}).(Variable)(:,iName); % depth=2
    end
end
end

