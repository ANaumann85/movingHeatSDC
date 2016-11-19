lapl=load('./laplace.mm');
lapls=sparse(lapl(2:end,1),lapl(2:end,2),lapl(2:end,3));
mass=load('./mass.mm');
masss=sparse(mass(2:end,1),mass(2:end,2),mass(2:end,3));
leigs=eigs(lapls,masss);

mvName=sprintf('mv-%0.2f.mm', 0.0);
mv=load(mvName);
mvs=sparse(mv(2:end,1), mv(2:end,2), mv(2:end,3));
mvEigs=eigs(mvs, masss);
leigs'
mvEigs'
