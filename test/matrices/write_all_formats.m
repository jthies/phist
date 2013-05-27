function write_all_formats(name,description,A,C)

key='no key';

mmwrite(['s_',name,'.mm'],A,description,'real',8);
mmwrite(['d_',name,'.mm'],A,description,'real',16);

mmwrite(['c_',name,'.mm'],C,description,'complex',8);
mmwrite(['z_',name,'.mm'],C,description,'complex',16);

if (nnz(A)==0)
  A(1,1)=1.0;
  warning(['cannot write empty matrix in hb-format, you will have to \n fix the files ',...
           '*_',name,'.rua manually']);
end

if (nnz(C)==0)
  C(1,1)=1.0i;
  warning(['cannot write empty matrix in hb-format, you will have to \n fix the files ',...
           '*_',name,'.cua manually']);
end

dm2hb(['s_',name,'.rua'],A,[],description,key,'RUA',8,2);
dm2hb(['d_',name,'.rua'],A,[],description,key,'RUA',16,2);

zm2hb(['c_',name,'.cua'],C,[],description,key,'CUA',8,2);
zm2hb(['z_',name,'.cua'],C,[],description,key,'CUA',16,2);

end
