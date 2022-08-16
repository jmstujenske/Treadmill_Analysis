function out=dir2file(in)
if isstruct(in) && length(in)>0
out=fullfile(in.folder,in.name);
else
    out=[];
end
end
