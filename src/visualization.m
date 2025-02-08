function visualization(U_p, Xc_p, Yc_p, gamma, t)
    Ny = size(U_p,2);
    rho = U_p(:,:,1);
    u = U_p(:,:,2)./rho;
    v = U_p(:,:,3)./rho;
    E = U_p(:,:,4);
    e = E./rho - 0.5.*(u.^2+v.^2);
    P = (gamma-1).*rho.*e;
    c = sqrt(gamma*P./rho);
    M = sqrt(u.^2+v.^2)./c;

    persistent ax1 ax2 ax3 ax4

    if t < eps
        figure;

        % 初始化轴句柄
        ax1 = subplot(2, 2, 1);
        ax2 = subplot(2, 2, 2);
        ax3 = subplot(2, 2, 3);
        ax4 = subplot(2, 2, 4);
    else
        % 第一个子图
        contourf(ax1, Xc_p', Yc_p', rho', 30, 'LineColor', 'none');
        axis([0 4 0 2.4]);
        colormap('winter');
        colorbar(ax1);
        title(ax1, '密度等高线图');
        
        % 第二个子图
        contourf(ax2, Xc_p', Yc_p', P', 30, 'LineColor', 'none');
        axis([0 4 0 2.4]);
        colormap('winter');
        colorbar(ax2);
        title(ax2, '压强等高线图');
        
        % 第三个子图
        contourf(ax3, Xc_p', Yc_p', M', 30, 'LineColor', 'none');
        axis([0 4 0 2.4]);
        colormap('winter');
        colorbar(ax3);
        title(ax3, '马赫数等高线图');
        
        idx = floor(Ny/2);
        plot(ax4, Xc_p(:,idx)',u(:,idx)');
        title(ax4, '水平速度剖面图');
        
        % 使用 drawnow 来刷新图形窗口
        drawnow;
    end
end

