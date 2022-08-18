clc; clear;

tic

fhn_known = load('mat_fhn_sp135.mat');
v_known = fhn_known.v_sp3;
w_known = fhn_known.w_sp3;

% fhn spline interpolaation calculation

syms tt
syms baslangic_v baslangic_w bitis_v bitis_w 

n = (length(v_known) - 1)/1; % eleman sayisi
t = (1:1:n);

a_v = sym('a_v', [1 n]);
b_v = sym('b_v', [1 n]);
c_v = sym('c_v', [1 n]);

a_w = sym('a_w', [1 n]);
b_w = sym('b_w', [1 n]);
c_w = sym('c_w', [1 n]);

int_1_v = sym('int_1_v', [1 n-2]);
int_2_v = sym('int_2_v', [1 n-2]);
der_int_v = sym('der_int_v', [1 n-2]);

int_1_w = sym('int_1_w', [1 n-2]);
int_2_w = sym('int_2_w', [1 n-2]);
der_int_w = sym('der_int_w', [1 n-2]);

%---------------------------------------------------------------

v_fhn_tek = zeros(2,n/2);
w_fhn_tek = zeros(2,n/2);
v_fhn_ori = zeros(2,n);
w_fhn_ori = zeros(2,n);

for i=1:n
    
    v_fhn_ori(1,i) = v_known(i);
    v_fhn_ori(2,i) = t(i);
    
    w_fhn_ori(1,i) = w_known(i);
    w_fhn_ori(2,i) = t(i);
    
end

for i=1:2:n
    
    v_fhn_tek(1,(i+1)/2) = v_known(i);
    v_fhn_tek(2,(i+1)/2) = t(i);
    
    w_fhn_tek(1,(i+1)/2) = w_known(i);
    w_fhn_tek(2,(i+1)/2) = t(i);
    
end

n = length(v_fhn_tek(1,:));
%----------------------------------------------------------------
a_v(1) = 0; a_w(1) = 0;
%----------------------------------------------------------------
% spline interpoaltion
% end points

baslangic_v = a_v(1)*v_fhn_tek(2,1)^2 + b_v(1)*v_fhn_tek(2,1) + c_v(1) - v_fhn_tek(1,1);
bitis_v = a_v(n-1)*v_fhn_tek(2,2)^2 + b_v(n-1)*v_fhn_tek(2,2) + c_v(n-1) - v_fhn_tek(1,n);

baslangic_w = a_w(1)*w_fhn_tek(2,1)^2 + b_w(1)*w_fhn_tek(2,1) + c_w(1) - w_fhn_tek(1,1);
bitis_w = a_w(n-1)*w_fhn_tek(2,2)^2 + b_w(n-1)*w_fhn_tek(2,2) + c_w(n-1) - w_fhn_tek(1,n);

%interior points
for i=1:n-2
    
   int_1_v(i) = a_v(i)*v_fhn_tek(2,i+1)^2 + b_v(i)*v_fhn_tek(2,i+1) + c_v(i) - v_fhn_tek(1,i+1);
   int_2_v(i) = a_v(i+1)*v_fhn_tek(2,i+1)^2 + b_v(i+1)*v_fhn_tek(2,i+1) + c_v(i+1) - v_fhn_tek(1,i+1);
   der_int_v(i) = 2*a_v(i)*v_fhn_tek(2,i+1) + b_v(i) - 2*a_v(i+1)*v_fhn_tek(2,i+1) - b_v(i+1);
   
   int_1_w(i) = a_w(i)*w_fhn_tek(2,i+1)^2 + b_w(i)*w_fhn_tek(2,i+1) + c_w(i) - w_fhn_tek(1,i+1);
   int_2_w(i) = a_w(i+1)*w_fhn_tek(2,i+1)^2 + b_w(i+1)*w_fhn_tek(2,i+1) + c_w(i+1) - w_fhn_tek(1,i+1);
   der_int_w(i) = 2*a_w(i)*w_fhn_tek(2,i+1) + b_w(i) - 2*a_w(i+1)*w_fhn_tek(2,i+1) - b_w(i+1);
   
end

S_v = solve([int_1_v(:);int_2_v(:);der_int_v(:);baslangic_v(:);bitis_v(:)]); % tum a b c degerlerini elde eder
C_v = struct2cell(S_v); A_v = cell2sym(C_v); B_v = double(A_v);

S_w = solve([int_1_w(:);int_2_w(:);der_int_w(:);baslangic_w(:);bitis_w(:)]); % tum a b c degerlerini elde eder
C_w = struct2cell(S_w); A_w = cell2sym(C_w); B_w = double(A_w);

% bulunan degerler a1,a2, b1,b2 c1, c2 matrislerde yerine yazilir
aa_v = B_v(1:1:n-2); aa_v = reshape(aa_v,1,length(aa_v)); aa_v = [0 aa_v];
bb_v = B_v(n-1:1:2*n-3); bb_v = reshape(bb_v,1,length(bb_v));
cc_v = B_v(2*n-2:1:length(B_v)); cc_v = reshape(cc_v,1,length(cc_v));

aa_w = B_w(1:1:n-2); aa_w = reshape(aa_w,1,length(aa_w)); aa_w = [0 aa_w];
bb_w = B_w(n-1:1:2*n-3); bb_w = reshape(bb_w,1,length(bb_w));
cc_w = B_w(2*n-2:1:length(B_w)); cc_w = reshape(cc_w,1,length(cc_w));

%spline matrisi elde edilir
for i=1:n-1
    
    int_1_numerical_v(i) = aa_v(i)*tt^2 + bb_v(i)*tt + cc_v(i);
    int_1_numerical_w(i) = aa_w(i)*tt^2 + bb_w(i)*tt + cc_w(i);
    
end

for i=1:n-1
    
    sonuc_v(i) = symfun(int_1_numerical_v(i),tt);
    asd_v = sonuc_v(i);
    asd1_v = symfun(asd_v,tt);
    vv(i) = double(asd1_v(i));
    
    sonuc_w(i) = symfun(int_1_numerical_w(i),tt);
    asd_w = sonuc_w(i);
    asd1_w = symfun(asd_w,tt);
    ww(i) = double(asd1_w(i));
    
end

sp_int_time = toc

save spline_interpolation_data.mat

%% 

figure('units','normalized','outerposition',[0 0 1 1])
figure(1); clf(1);
plot(vv,'b','linewidth',4); grid on;
set(gca,'Fontsize',30)
xlabel(' Adým Sayýsý '); ylabel(' v ');
title('Spline Interpolation - v '); grid on;
% saveas(gcf,'Spline Interpolation - v.jpg');

figure('units','normalized','outerposition',[0 0 1 1])
figure(2); clf(2);
plot(ww,'b','linewidth',4); grid on;
set(gca,'Fontsize',30)
xlabel(' Adým Sayýsý '); ylabel(' w ');
title('Spline Interpolation - w '); grid on;
% saveas(gcf,'Spline Interpolation - w.jpg');

figure('units','normalized','outerposition',[0 0 1 1])
figure(3); clf(3);
plot(vv,ww,'b','linewidth',4); grid on;
set(gca,'Fontsize',30)
xlabel(' v '); ylabel(' w ');
title('Spline Interpolation - v/w '); grid on;
% saveas(gcf,'Spline Interpolation - vw.jpg');

%%

%istenilen zamanda fonk degeri tt araligi secilerek yeerine yazilir

clc; clear; 

n = 299+2; syms tt

fhn_known = load('mat_fhn_sp135.mat');
intpol_spline = load('fhn_sp1_spline_interpolation_data');

v_known = fhn_known.v_sp1;
w_known = fhn_known.w_sp1;

for i = 2:2:2*n-4

    sonuc_v_1 = symfun(intpol_spline.int_1_numerical_v(i/2),tt);
    fonk_degeri_v(1,i/2) = double(sonuc_v_1(i));
    fonk_degeri_v(2,i/2) = intpol_spline.t(i);

    sonuc_w = symfun(intpol_spline.int_1_numerical_w(i/2),tt);
    fonk_degeri_w(1,i/2) = double(sonuc_w(i));
    fonk_degeri_w(2,i/2) = intpol_spline.t(i);

    v_fhn_known(1,i/2) = v_known(i);
    v_fhn_known(2,i/2) = intpol_spline.t(i);
    
    w_fhn_known(1,i/2) = w_known(i);
    w_fhn_known(2,i/2) = intpol_spline.t(i);

end

m_s_e_v = immse(v_fhn_known(1,:),fonk_degeri_v(1,:))
m_a_e_v = mae(v_fhn_known(1,:),fonk_degeri_v(1,:))
r_m_s_e_v = sqrt(m_s_e_v)
n_r_m_s_e_v = (r_m_s_e_v/(max(v_fhn_known(1,:)) - min(v_fhn_known(1,:))))*100
% R = corrcoef(vvv,fonk_degeri_v(1,:),'Rows','complete')

m_s_e_w = immse(w_fhn_known(1,:),fonk_degeri_w(1,:))
m_a_e_w = mae(w_fhn_known(1,:),fonk_degeri_w(1,:))
r_m_s_e_w = sqrt(m_s_e_w)
n_r_m_s_e_w = (r_m_s_e_w/(max(w_fhn_known(1,:)) - min(w_fhn_known(1,:))))*100
% R = corrcoef(w_fhn_known,fonk_degeri_w,'Rows','complete')

%%

fig1 = figure('Position',get(0,'Screensize'));
plot(v_fhn_known(1,:),'LineStyle','-','Marker','*','Color','r','MarkerSize',25,'linewidth',2)
grid on; hold on;
plot(fonk_degeri_v(1,:),'LineStyle','-','Marker','o','Color','k','MarkerSize',25,'linewidth',2)
ylabel('v, [V]')
xlabel('Number of Samples');
set(gca,'Fontsize',60);
legend({'Original','Calculated by Spline-Interpolation'},'Location','southwest');
saveas(fig1, 'fhn_sp3_v_spinter.jpg');

fig2 = figure('Position',get(0,'Screensize'));
plot(w_fhn_known(1,:),'LineStyle','-','Marker','*','Color','r','MarkerSize',25,'linewidth',2)
grid on; hold on;
plot(fonk_degeri_w(1,:),'LineStyle','-','Marker','o','Color','k','MarkerSize',25,'linewidth',2)
ylabel('u, [V]')
xlabel('Number of Samples');
set(gca,'Fontsize',60);
legend({'Original','Calculated by Spline-Interpolation'},'Location','southwest');
saveas(fig2, 'fhn_sp3_w_spinter.jpg');

fig3 = figure('Position',get(0,'Screensize'));
plot(v_fhn_known(1,:),w_fhn_known(1,:),'LineStyle','-','Marker','*','Color','r','MarkerSize',25,'linewidth',2)
grid on; hold on;
plot(fonk_degeri_v(1,:),fonk_degeri_w(1,:),'LineStyle','-','Marker','o','Color','k','MarkerSize',25,'linewidth',2)
ylabel('u')
xlabel('v');
set(gca,'Fontsize',60);
legend({'Original','Calculated by Spline-Interpolation'},'Location','southwest');
saveas(fig3, 'fhn_sp3_cw_known.jpg');
