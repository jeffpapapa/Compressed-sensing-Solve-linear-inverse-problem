% X = [5:35];
X = [1:35];
tic
run_times = 500;
for coin = 1:5
    coin
    time1=tic;
    Y = zeros([1,35]);
    for N = 1:35
        N
        A = round(rand(N,100));
        success_time = 0;
        for running_time =1:run_times
            running_time;
            answ = zeros([100,1]);
            fake_weight = random('uniform',-1,1,coin,1);
            fake = randsample(100,coin,false);
            for i = 1:coin
                answ(fake(i)) = fake_weight(i); 
            end
            b = A*answ;
%             (u,phi,tao,initial,threshold);
new_x=zeros([100,1]);

% [new_x, ERR] = gradient_decent(b,A,abs((max(abs(A'*b))/(10*(coin^1.05)))),new_x,10^-5,coin);
% [new_x, ERR] = gradient_decent(b,A,abs((max(abs(A'*b))/(10*(coin^0.775)))),new_x,10^-5,coin);
%[new_x, ERR] = gradient_decent(b,A,abs((max(abs(A'*b))/(10*(coin^0.7)))),new_x,10^-8,coin);
[new_x, ERR,tao_record,alpha_record] = gradient_decent(b,A,abs((max(abs(A'*b))/(9.8*(coin^(0.735))))),new_x,10^(-5),coin);
x1 = new_x;
err = norm(x1-answ,2);
[value,myans] = maxk(abs(x1),coin);
sucess = 0;
if(length(setdiff(fake,myans))<=0) 
    sucess = 1;
end
sucess;
% for i = 1:coin
% 
%          if(ismember(myans(i),fake) ~= 1)
%              sucess = 0;
%              i;
%          end;
% 
% end
% sucess;

            if(sucess>0)
%             if(abs(answ - x1(1:100))<(10^(-5)))
                 success_time = success_time + 1;
            end
             success_time;
        end
        
        percent = success_time/run_times;
        Y(N-0)=percent;
    end
    plot(X,Y);
    hold on;
    time(coin) = toc(time1);
end
time
hold off;
xlabel('number of weightings') 
ylabel('probability of success')
legend({'fake coin : 1','fake coin : 2','fake coin : 3','fake coin : 4','fake coin : 5'},'Location','northeast')
toc