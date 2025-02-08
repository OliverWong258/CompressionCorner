function main(Nx,Ny, splitting_type)
%网格生成
corner_point=1;
x_l = 0; 
x_r = 4;
y_b = 0; 
y_t = 2.4;
angle = 15;
dx = (x_r - x_l) /Nx;
dy = (y_t - y_b) /Ny;
X = (x_l:dx:x_r)';
Y = (y_b:dy:y_t)';
[phisical_X,phisical_Y]=GridGeneration(X,Y,angle);
[J,px1,py1,px2,py2]=Jacobi(X,Y,dx,dy,angle,corner_point);
%设置初始条件
t = 0;
phisical_U = Initial(X,Y);
U=phisical_U.*J(4:3+Nx,4:3+Ny);
CFL = 0.5;
gamma = 1.4;
T = 5.0;
rec =[0.1;0.6;0.3];
epsilon = 1e-9;
tic;
while true
    lamdaX = MaxEigX(U,gamma);
    lamdaY  = MaxEigY(U,gamma);
    dt = CFL/(lamdaX/(dx)+lamdaY/(dy));
    if(t+dt>T)
        dt=T-t;
    end
    if dt<100*eps
        break;
    end
    U_next= TimeIteration(dx,dy,dt,U,px1,py1,px2,py2,J,gamma,rec,epsilon, splitting_type);
    %检测是否收敛
    if ~any(abs(U-U_next) > 1e-5)
        break;
    end
    % 更新U
    U = U_next;
    phisical_U=U./J(4:3+Nx,4:3+Ny);
    visualization(phisical_U, phisical_X, phisical_Y, gamma, t)
    t = t+dt;   
    disp(['Time: ', num2str(t), ' s.']); 
end
toc;
disp(['Cost: ', num2str(elapsedTime), ' seconds.']); 
end