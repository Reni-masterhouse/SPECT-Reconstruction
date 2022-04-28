%%
clear

%定义系统传输矩阵
CIJ = zeros(128*128,60*128);

%定义初值
grid_size = 4.0625;
D = 289.37/grid_size;
radius = 200/grid_size;

%每6度计算一次，重复该角度的128个探测器
for theta = 0:6:354
	for n = 1:128
		%1点随角度的变化关系
		X1=D*cos(pi*theta/180)+(n-64.5)*sin(pi*theta/180);
		Y1=D*sin(pi*theta/180)-(n-64.5)*cos(pi*theta/180);
        
        %生成向量s
        s=[-sin(pi*theta/180),cos(pi*theta/180)];
        
        [x, y] = meshgrid((1:128)-64.5,(1:128)-64.5);
        x = reshape(x', [128*128,1]);
        y = reshape(y', [128*128,1]);
        p = [x-X1, y-Y1];
        d = abs(p * s');
        CIJ(:, (theta/6)*128 + n) = (d < 1 / 2);
        
        %衰减校正
        c = x.^2 + y.^2 - radius^2;
        dir = [X1-x, Y1-y];
        dir = dir ./ vecnorm(dir, 2, 2);
        b = 2 * dot([x,y], dir, 2);
        
        %计算距离
        delta = b.^2 - 4*c;
        t1 = (-b + sqrt(delta)) / 2;
        t2 = (-b - sqrt(delta)) / 2;
        dt = t1 - t2;
        
        dis = zeros(128*128, 1);
        inter = delta > 0;
        m1 = inter & (t2 > 0);
        m2 = inter & (t1 > 0) & (t2 < 0);
        
        dis(m1) = dt(m1);
        dis(m2) = t1(m2);
        
        %得到矩阵
        AIJ = exp(- 0.15454 * dis * grid_size/10);
        CIJ(:, (theta/6)*128 + n) = CIJ(:, (theta/6)*128 + n) .* AIJ;
    end
end

%%
%CRC与CNR计算

[x_im, y_im] = meshgrid((1: 128)-64.5, (1: 128)-64.5);
x_im = reshape(x_im', [128*128,1]);
y_im = reshape(y_im', [128*128,1]);

%确定热区和背景区位置
hr_theta = 0:60:300;
hr_radius_mm = 4:4:24;
hr_radius = hr_radius_mm/grid_size;
back_radius = 200/grid_size;
hot_region = zeros([128*128, 6], 'logical');
back_region = (x_im.^2 + y_im.^2 - back_radius.^2) < 0;

for i = 1:6
    this_x = 100/grid_size*cos(hr_theta(i)*pi/180);
    this_y = 100/grid_size*sin(hr_theta(i)*pi/180);
    hot_region(:, i) = ((x_im - this_x).^2 + (y_im - this_y).^2) < hr_radius(i) .^2;
    back_region = back_region & ~hot_region(:, i);
end

%%
%读取数据
filename = 'Proj_1e5Counts';
pid = fopen (filename, 'r');
data = fread (pid, 128*60, 'float32');
fclose(pid);

%%
fi=ones(1,128*128);	%给一个最初的猜测值

%可视化waitbar
wb = waitbar(0, 'iteration...');
num_iter = 200;

%定频率保存
save_iter = [1, 2, 5, 10, 20, 50, 100, 200];

CRCs = zeros(4*num_iter, 6);
CNRs = zeros(4*num_iter, 6);

%做迭代
for iter=1:num_iter
    for sub=1:4
        sub_size = 15;
        proj_size = sub_size * 128;
        Cij = CIJ(:, (sub-1)*proj_size+1:sub*proj_size);
        prosub = data((sub-1)*proj_size+1:sub*proj_size);

        p = fi * Cij;
        ratio = prosub ./ transpose(p);
        eff = (Cij * ratio) ./ sum(Cij, 2); 
        fi = fi .* transpose(eff);
        
        %计算CRC/CNR
        mu_bak = mean(fi(back_region));
        sigma_bak = std(fi(back_region));
        for hr = 1:6
            mu_rod = mean(fi(hot_region(:, hr)));
            CRC = (mu_rod - mu_bak) / mu_bak / 4;
            CNR = abs(mu_rod - mu_bak) / sigma_bak;
            CRCs((iter-1)*4+sub, hr) = CRC;
            CNRs((iter-1)*4+sub, hr) = CNR;
        end

        waitbar(((iter-1)*4+sub)/(4*num_iter), wb, ['iter: ', num2str(iter), ' subset: ', num2str(sub)]);
    end
    
    %保存图像和dat文件
    if ismember(iter, save_iter)
        img_name = ['secimage_e5/', num2str(iter), '.png'];
        dat_name = ['secdat_e5/', num2str(iter), '.dat'];
        img = reshape(fi, [128,128]);
        fid = fopen(dat_name, 'wb');
        fwrite(fid, img, 'single');
        fclose(fid);
        mv = max(img, [], 'all');
        imwrite(img/mv, img_name); 
    end
end
delete(wb);

%%另一个可行的存储方式

% %写入raw文件
% %result='OSEM.raw';
% %fid=fopen(result,'w+');
% %cnt=fwrite(fid, c,'double');
% %fclose(fid);
% 
% result = 'image.dat';   %%% 定义路径下的文件
% fid = fopen(result, 'wb');
% cnt=fwrite(fid, c, 'single');    %%% data为需要写入的对象
% fclose(fid);
% 
% %%显示图像
% %fid=fopen(result,'r');
% %Pix=fread(fid,[128 128],'single');
% %imshow(Pix,[])
% %imshow(double(Pix));


%%
%绘制CRC和CNR的趋势图线
subplot(2,1,1);
for hr = 1:6
    plot(1:(4*num_iter), CRCs(:, hr));
    hold on;
end 
leg = arrayfun(@(x) [num2str(2*x), 'mm'], hr_radius_mm, 'UniformOutput', false);
legend(leg);
title('CRC');
subplot(2,1,2);
for hr = 1:6
    plot(1:(4*num_iter), CNRs(:, hr));
    hold on;
end
legend(leg);
title('CNR');
