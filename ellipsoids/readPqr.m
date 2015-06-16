function pqr = readPqr(filename)
fh = fopen(filename,'r');
C1 = textscan(fh, '%s%d%s%s%d%f%f%f%f%f');
frewind(fh);
C2 = textscan(fh, '%s%d%s%s%s%d%f%f%f%f%f');
fclose(fh);
if length(C1{1}) > length(C2{1})
  C = C1; C{10} = C{9}; C{9} = C{8}; C{8} = C{7}; C{7}=C{6}; C{6}= ...
      C{5}; C{5} = repmat({''},length(C{2}));
else
  C = C2; 
end

xyz = [C{7} C{8} C{9}];
q = C{10};
r = C{11};
resnum = C{6};
resid = C{4};
atomid = C{3};
atomnum = C{2};
atominfo = struct('resid',resid,'atomid',atomid);
pqr = struct('xyz', xyz, 'q', q, 'r', r, 'resnum',resnum, ...
	     'atomnum',atomnum,'atominfo',atominfo);
