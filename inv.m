load('forward.mat'); % résultats de la modélisation directe enregistrée dans le fichier forward.mat
%Ces résultats sont les variables [G, Gz, x_pr, y_pr, npro, nstn, rho]
%obtenues en sortie de la fonction new3dgm

figure(21); imagesc(G);colorbar;
title('G_{mn,k}'); saveas(gcf,'G.png');
d0=sum(G,2)*rho;
n_d=size(G,1);
ratio_br=[zeros(1,n_d); d0'/100; d0'/33.34; d0'/10];
n_br=size(ratio_br,1);

txt_br=[0 1 3 10];
figure(11);clf;
figure(12);clf;
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

for j=1:n_br
    d=d0+rand(length(d0),1).*ratio_br(j,:)';
    gz=reshape(d,[nstn,npro])';
    if j==1
    figure(20); subplot(122);imagesc(gz);colorbar;
    axis equal;xlabel('m profils'); ylabel('n stations');
    end
    figure(11);
    if j==1
        subplot(2,2,j); imagesc(x_pr,y_pr, Gz-gz); colorbar;
        title('gz-gz_G (mgals)');
    else
        gz2=reshape(d0,[nstn,npro])';
        subplot(2,2,j); imagesc(x_pr,y_pr, gz2-gz); colorbar;
        title(['gz-(gz+ ' num2str(txt_br(j)) '% de bruit)']);
    end
    
    xlabel('x(m)'); ylabel('y(m)');
    %m=inv(G'*G)*(G'*d);
    [U,s,V] = csvd(G);
    cond_numb=s(1)/s(end);%cond_numb2=cond(G)  1.8627e+13
    
    %Test G singular -> p.87 eq. 3.105 Aster
    test_sing=zeros(n_d-1,1);
    for i=2:n_d
        test_sing(i)=(G(i-1,:)*G(i,:)')/(norm(G(i-1,:))*norm(G(i,:)));
    end
    mean(test_sing); %0.9735
    
    figure(12);
    subplot(n_br,3,3*j-2); alpha = l_curve(U,s,d);
    M=size(G'*G);
    I=eye(M);
    disp('m_levenberg')
    tic; m_levenberg=inv(G'*G+I)*(G'*d); toc
    
    disp('m_tikh'),tic;
    x_tikh_l = tikhonov(U,s,V,d,alpha);
    m_tikh=V*x_tikh_l;toc
    m_true=0.2*ones(size(G,2),1);
    disp('m_cgls')
    tic;[m_cgls,flag,resne,iter] = cgls(G,d,0,1e-2,20); toc,flag,iter,resne
    
    figure(12);
    subplot(n_br,3,[3*j-1, 3*j]);plot(m_true);hold on;
    plot(m_levenberg,'ro-');
    plot(m_tikh,'gx-');
    plot(m_cgls,'k*-');hold off;
    if j==1
    legend(gca, 'Model réel', 'Levenberg \alpha=1', (['Tikhonov \alpha=corner']),...
        'CGLS', 'Location','NorthEastOutside');
    end
    ylabel('\rho (g/cm^3)'); xlabel('Faces voxels triangulaires'); title('');
    
    disp('erreurs ml, tikh, cgls')
    disp([norm(m_true-m_levenberg), norm(m_true-m_tikh),...
        norm(m_true-m_cgls)]/norm(m_true));  
end

figure(11); saveas(gcf,'residu_gz.png');
figure(12); saveas(gcf,'models.png');


