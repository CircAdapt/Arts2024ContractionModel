function Put(Type,Variable,Name,Value)
% function Put(Type,Variable,Name,Value)
% Type=string, Variable=string, Name=string or {str1,str2}
% e.g. Put('ArtVen','Len',{'Ca','Br'},[0.11,0.12])
% inverse action of Get
% if Name=='All' -> all variables
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
        dName=strcmp(Name{i},Names)*iN;
        % if non-existent, value 0 is inserted
        iName=[iName,dName];
    end
end

iName=iName(logical(iName~=0)); % select only existing variable
% iName is array of indices

nCol= length(iName);
[nt,nVal]=size(Value);
if nVal==1
    Value=repmat(Value,[1,nCol]);
end

switch length(T)
    case 1 % depth=1
        if ~isfield(P.(T{1}),Variable) % if 'variable' does not exist yet
            P.(T{1}).(Variable)= zeros(nt,n);
        end
        ntMax= max([1,nt,size(P.(T{1}).(Variable),1)]); %number of rows
        P.(T{1}).(Variable)(ntMax-nt+1:ntMax,iName)= Value;
    case 2 % depth=2
        if ~isfield(P.(T{1}).(T{2}),Variable) % if 'variable' does not exist yet
            P.(T{1}).(T{2}).(Variable)= zeros(nt,n);
        end
        ntMax= max([1,nt,size(P.(T{1}).(T{2}).(Variable),1)]); %number of rows
        P.(T{1}).(T{2}).(Variable)(ntMax-nt+1:ntMax,iName)= Value;
end
end

