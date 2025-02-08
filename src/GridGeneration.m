function [physical_X,physical_Y]=GridGeneration(X,Y,angle)
computational_X = 0.5*(X(1:end-1,1) + X(2:end,1));
computational_Y = 0.5*(Y(1:end-1,1) + Y(2:end,1));
tan_ang = tan(pi/(180/angle));
Nx=length(computational_X);
Ny=length(computational_Y);
physical_X=zeros(Nx,Ny);
physical_Y=zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        physical_X(i,j)=computational_X(i);
        if computational_X(i)<1
            physical_Y(i,j)=computational_Y(j);
        else
            physical_Y(i,j)=computational_Y(j)+(1-computational_Y(j)/2.4)*tan_ang*(computational_X(i)-1);
        end
    end
end

end