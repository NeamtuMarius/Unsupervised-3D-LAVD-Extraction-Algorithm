function [e1,e2,e3] = EigEval3D(Vmxx,Vmxy,Vmxz,Vmyx,Vmyy,Vmyz,Vmzx,Vmzy,Vmzz)
%EigEval3D evaluate the eigenvalues in each point of a 3d tensorial field
%   Vmij are the component of the tensor defined in space, so are 3D tensors

for i = 1:size(Vmxx,1)
    
    for j = 1:size(Vmxx,2)
    
        for k = 1:size(Vmxx,3)
        
            A = [Vmxx(i,j,k) Vmxy(i,j,k) Vmxz(i,j,k);
                 Vmyx(i,j,k) Vmyy(i,j,k) Vmyz(i,j,k);
                 Vmzx(i,j,k) Vmzy(i,j,k) Vmzz(i,j,k)];
            

            e = sort(eig(A));
            
            e3(i,j,k) = e(1);
            e2(i,j,k) = e(2);
            e1(i,j,k) = e(3);

            
        end
            
    end
    
end
end

