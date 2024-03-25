function MergemFiles
mFiles=dir('*.m');
nF=length(mFiles);
Text='';
for iF=1:nF
    Name=mFiles(iF).name;
    fid = fopen(Name,'r');
    c=fread(fid)';
    Text=[Text,[36,35,36],Name,[35,36,35],c];%[$ # $],[# $ #]
    fclose(fid);
end
Text=[Text,[36,35,36]];
fid = fopen('mFiles.txt','w');
fwrite(fid,Text);
fclose(fid);
end
