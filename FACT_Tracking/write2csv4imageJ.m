function write2csv4imageJ(data, csvname) 
fid = fopen([csvname '.csv'],'wt'); 
fprintf(fid, '%s,%s,%s\n', ' ', 'X', 'Y'); 
for d = 1:size(data, 1)
     fprintf(fid, '%f,%f,%f\n', d, data(d, 1), data(d, 2));
end
fclose(fid);
end