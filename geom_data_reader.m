clear;
% close all;
% movietitle=input('Name the movie file' , 's');
% 2-D Case
fid = fopen('two_D2.dat','r');
InputText=textscan(fid,'%s %f',1);
nt = InputText{2};nt=8197;
HeaderMain=sprintf('nt = %s',num2str(nt));
% disp(HeaderMain);
t = zeros(nt,1); nx = zeros(nt,1); ny = zeros(nt,1);dt = zeros(nt,1);
% M = zeros(
n = 0;
tag = 0;
while(~feof(fid))
    n = n + 1
    % for n = 1:nt
    t(n)=n*dt(n);
    sprintf('t = %s', num2str(t(n))); % Display block number
    % Read in t, nx, ny from block header
    InputText=textscan(fid,'%s %f %s %f %s %f %s %f',1);
    t(n) = InputText{2}; nx(n) = InputText{4};
    
    ny(n) = InputText{6};dt(n) = InputText{8};
    
    if n ~= 1
        t(n)=t(n-1)+dt(n);
    else
        t(n)=dt(n);
    end
    NumRows=nx(n)*ny(n);
    x=zeros(ny(n),nx(n));
    y=zeros(ny(n),nx(n));
    u=zeros(ny(n),nx(n));
    v=zeros(ny(n),nx(n));
    rho=zeros(ny(n),nx(n));
    p=zeros(ny(n),nx(n));
    vmag =zeros(ny(n),nx(n));
    theta=zeros(ny(n),nx(n));
    a=zeros(ny(n),nx(n));
    b=zeros(ny(n),nx(n));
    l=zeros(ny(n),nx(n));
    m=zeros(ny(n),nx(n));
    h=zeros(ny(n),nx(n));
    delta=zeros(ny(n),nx(n));
    
    data = zeros(NumRows,11);
    
    InputText=textscan(fid,'%s',6);
    InputText=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f');
    data=cell2mat(InputText);
    for i = 1:nx(n)
        for j = 1:ny(n)
            x(j,i) = data(j+(i-1)*ny(n),1);
            y(j,i) = data(j+(i-1)*ny(n),2);
            u(j,i) = data(j+(i-1)*ny(n),3);
            v(j,i) = data(j+(i-1)*ny(n),4);
            rho(j,i) = data(j+(i-1)*ny(n),5);
            p(j,i) = data(j+(i-1)*ny(n),6);
            vmag(j,i) = sqrt(u(j,i)^2+v(j,i)^2);
            theta(j,i)= atan2(v(j,i),u(j,i));
            
            a(j,i) = data(j+(i-1)*ny(n),7 );
            b(j,i) = data(j+(i-1)*ny(n),8 );
            l(j,i) = data(j+(i-1)*ny(n),9 );
            m(j,i) = data(j+(i-1)*ny(n),10);
            delta(j,i) = a(j,i)*m(j,i) - b(j,i)*l(j,i);
            
            h(j,i) = data(j+(i-1)*ny(n),11);
            
        end
    end
    skip = 1;
    if size(x,2)==1
%         plot(y,rho)
    elseif mod(n,skip) == 0
        figure(1)
        maximize
        mach=vmag./sqrt(1.4*p./rho);
        e=0.5*rho.*vmag.^2+1.0/0.4*p;
        surf(x,y,rho);
        axis equal;
%     figure(1);plot(x,y,'k')
        xlim([0. .9]);ylim([-.5 .5]);%zlim([-0.1 25]);
        colorbar;
        view([0,0,1]);
        xlabel('x');ylabel('y')
%         pause
        %         figure(2)
        %         surf(x,y,u);axis equal;xlim([0 0.55]);ylim([-0.5 0.5]);
        %         figure(2)
        %         contour(x,y,u(:,:));axis equal;%xlim([0 0.5]);ylim([0,1]);%zlim([-0.1 2.1]);colorbar;
        %         fprintf('Step')
        title('Mach number')
%         if(n>1)
%             if tag == 0;
%                 pause;tag = 1;
%             end
%             M(n-1) = getframe(gcf);
%         end
        pause(.01)
    end
    %     clear('x','y','u','v','rho','p')
end
                figure(3)
                plot(y(:,50)./x(:,50),rho(:,50),'.');ylim([.4 1.]);xlim([-0.6 0.4 ]);
        ylabel('\rho');xlabel('^y/_x')
        title('Density for riemann test problem')
% movie2avi(M,movietitle);
