% McDermott
% 3-7-2016
% curlchem.m
%
% Apply Curl's model to fast chemistry in a hypothetical cell
%
% Consider the following simple reaction scheme:
%
% F1 + A1 --> 2*P1
% F2 + A2 --> 2*P2
% F1 + A2 --> 2*P3

close all
clear all

% initial volume fractions

X_F1_0 = 0.1
X_F2_0 = 0.4
X_A2_0 = 0.2
X_P1_0 = 0.0
X_P2_0 = 0.0
X_P3_0 = 0.0
X_A1_0 = 1 - X_F1_0 - X_F2_0 - X_A2_0 - X_P1_0 - X_P2_0 - X_P3_0

N = 100000;
N_F1_0 = X_F1_0*N; N1 = N_F1_0;
N_F2_0 = X_F2_0*N; N2 = N_F1_0+N_F2_0;
N_A1_0 = X_A1_0*N; N3 = N_F1_0+N_F2_0+N_A1_0;
N_A2_0 = X_A2_0*N;
N_P1_0 = X_P1_0*N;
N_P2_0 = X_P2_0*N;
N_P3_0 = X_P3_0*N;

N_F1 = N_F1_0; % p.ispec=1
N_F2 = N_F2_0; % p.ispec=2
N_A1 = N_A1_0; % p.ispec=3
N_A2 = N_A2_0; % p.ispec=4
N_P1 = N_P1_0; % p.ispec=5
N_P2 = N_P2_0; % p.ispec=6
N_P3 = N_P3_0; % p.ispec=7

% assign species to particles, add random position for visualization

%figure(1)
for i = 1:N
    %p(i).x = rand(1);
    %p(i).y = rand(1);
    if i<=N1
        p(i).spec = 'F1';
        p(i).ispec = 1;
        %plot(p(i).x,p(i).y,'bo'); hold on
    elseif i>N1 & i<=N2
        p(i).spec = 'F2';
        p(i).ispec = 2;
        %plot(p(i).x,p(i).y,'ko')
    elseif i>N2 & i<=N3
        p(i).spec = 'A1';
        p(i).ispec = 3;
        %plot(p(i).x,p(i).y,'g.')
    elseif i>N3
        p(i).spec = 'A2';
        p(i).ispec = 4;
        %plot(p(i).x,p(i).y,'b.')
    end
end
%title('initial condition')

iter=0;
n_iter=100000;

while iter<n_iter

    % compute new mass fractions

    X_F1 = N_F1/N;
    X_F2 = N_F2/N;
    X_A1 = N_A1/N;
    X_A2 = N_A2/N;
    X_P1 = N_P1/N;
    X_P2 = N_P2/N;
    X_P3 = N_P3/N;

    figure(2)
    plot(iter,X_F1,'bo'); hold on
    plot(iter,X_F2,'ko')
    plot(iter,X_A1,'g^')
    plot(iter,X_A2,'b^')
    plot(iter,X_P1,'c+')
    plot(iter,X_P2,'k+')
    plot(iter,X_P3,'m+')
    legend('F1','F2','A1','A2','P1','P2','P3','location','eastoutside')

    iter = iter + 1;
    N_iter = N_F1 + N_F2 + N_A1 + N_A2;

    % sort particles

    [dummy, order] = sort([p(:).ispec],'ascend');
    p = p(order);

    % create a random permutation vector

    %I = randperm(N);
    I = randi(N_iter,N_iter,1);

    % now react pairs of particles

    for i=2:2:N_iter
        % p(i-1) reacts (if possible) with p(i)
        if ( p(I(i-1)).spec=='F1' & p(I(i)).spec=='A1' ) | ( p(I(i-1)).spec=='A1' & p(I(i)).spec=='F1' )
            p(I(i-1)).spec='P1'; p(I(i-1)).ispec=5;
            p(I(i)).spec='P1';   p(I(i)).ispec=5;
            N_F1 = N_F1-1;
            N_A1 = N_A1-1;
            N_P1 = N_P1+2;
        elseif ( p(I(i-1)).spec=='F2' & p(I(i)).spec=='A2' ) | ( p(I(i-1)).spec=='A2' & p(I(i)).spec=='F2' )
            p(I(i-1)).spec='P2'; p(I(i-1)).ispec=6;
            p(I(i)).spec='P2';   p(I(i)).ispec=6;
            N_F2 = N_F2-1;
            N_A2 = N_A2-1;
            N_P2 = N_P2+2;
        elseif ( p(I(i-1)).spec=='F1' & p(I(i)).spec=='A2' ) | ( p(I(i-1)).spec=='A2' & p(I(i)).spec=='F1' )
            p(I(i-1)).spec='P3'; p(I(i-1)).ispec=7;
            p(I(i)).spec='P3';   p(I(i)).ispec=7;
            N_F1 = N_F1-1;
            N_A2 = N_A2-1;
            N_P3 = N_P3+2;
        end
    end

    % figure(3)
    % for i = 1:N
    %     if p(i).spec=='F1'
    %         plot(p(i).x,p(i).y,'bo'); hold on

    %     elseif p(i).spec=='F2'
    %         plot(p(i).x,p(i).y,'ko'); hold on

    %     elseif p(i).spec=='A1'
    %         plot(p(i).x,p(i).y,'g.'); hold on

    %     elseif p(i).spec=='A2'
    %         plot(p(i).x,p(i).y,'b.'); hold on

    %     elseif p(i).spec=='P1'
    %         plot(p(i).x,p(i).y,'c+'); hold on

    %     elseif p(i).spec=='P2'
    %         plot(p(i).x,p(i).y,'y+'); hold on

    %     elseif p(i).spec=='P3'
    %         plot(p(i).x,p(i).y,'m+'); hold on

    %     end
    % end

    pause(0.001)

    if (N_F1==0 | N_A1==0) & (N_F2==0 | N_A2==0) & (N_F1==0 | N_A2==0)
        disp('all reactions complete!')
        X_F1 = N_F1/N;
        X_F2 = N_F2/N;
        X_A1 = N_A1/N;
        X_A2 = N_A2/N;
        X_P1 = N_P1/N;
        X_P2 = N_P2/N;
        X_P3 = N_P3/N;
        X = [X_F1,X_F2,X_A1,X_A2,X_P1,X_P2,X_P3]
        X_Sum = sum(X)
        break
    end

end
























