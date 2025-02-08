function U_next = TimeIteration(dx,dy,dt,U,px1,py1,px2,py2,J,gamma,rec,epsilon, splitting_type)
% 采用三阶荣格库塔推进时间
Nx = size(U,1);    
Ny = size(U,2);

[F,G] = Flux_Splitting(U,px1,py1,px2,py2,J,gamma,rec,epsilon, splitting_type);

U1 = U - dt/dx*( F(2:Nx+1,:,:) - F(1:Nx,:,:) ) - dt/dy*( G(:,2:Ny+1,:) - G(:,1:Ny,:) );

[F1,G1] = Flux_Splitting(U1,px1,py1,px2,py2,J,gamma,rec,epsilon, splitting_type);

U2 = 0.75*U + 0.25*( U1 - dt/dx*( F1(2:Nx+1,:,:) - F1(1:Nx,:,:) ) - dt/dy*( G1(:,2:Ny+1,:) - G1(:,1:Ny,:) ) );


[F2,G2] = Flux_Splitting( U2,px1,py1,px2,py2,J,gamma,rec,epsilon, splitting_type);

U_next = 1/3*U + 2/3*( U2 - dt/dx*( F2(2:Nx+1,:,:) - F2(1:Nx,:,:) ) - dt/dy*( G2(:,2:Ny+1,:) - G2(:,1:Ny,:) ) );

end


