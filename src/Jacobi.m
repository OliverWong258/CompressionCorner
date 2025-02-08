function [py,Xp1,Xp2,Yp1,Yp2]=Jacobi(X,Y,dx,dy,angle,l)
computational_X = 0.5*(X(1:end-1,1) + X(2:end,1));
conputational_Y = 0.5*(Y(1:end-1,1) + Y(2:end,1));
computational_Xb=[computational_X(1)-3*dx;computational_X(1)-2*dx;computational_X(1)-dx;computational_X;computational_X(end)+dx;computational_X(end)+2*dx;computational_X(end)+3*dx];
computational_Yb=[conputational_Y(1)-3*dy;conputational_Y(1)-2*dy;conputational_Y(1)-dy;conputational_Y;conputational_Y(end)+dy;conputational_Y(end)+2*dy;conputational_Y(end)+3*dy];
k = tan(pi/(180/angle));
Nx=length(computational_Xb);
Ny=length(computational_Yb);
[Xm,Ym]=ndgrid(computational_Xb,computational_Yb);

py=1-k*(Xm-1)/2.4;

Xp1=ones(Nx,Ny); 
Xp2=zeros(Nx,Ny);
Yp1=-(k*(2.4-Ym)/2.4)./py; 
Yp2=1./py;

idx=find(computational_Xb<l, 1, 'last' );
py(1:idx,:)=1;

Xp1(1:idx,:)=1; Xp2(1:idx,:)=0;
Yp1(1:idx,:)=0; Yp2(1:idx,:)=1;

end