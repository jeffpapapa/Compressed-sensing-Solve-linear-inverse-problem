clear all;
%A = round(rand(35,100));
x_size = 35;
y_size = 100;
A = round(rand(x_size,y_size));
coin=5
answ = zeros([y_size ,1]);
            fake_weight = random('uniform',-1,1,coin,1);
            fake = randsample(y_size ,coin,false);
            for i = 1:coin
                answ(fake(i)) = fake_weight(i); 
            end
            b = A*answ;
%             (u,phi,tao,initial,threshold);
new_x=ones([y_size ,1])*0;
[new_x, ERR,tao_record,alpha_record] = gradient_decent(b,A,abs((max(abs(A'*b))/(9.8*(coin^(0.735))))),new_x,10^(-5),coin);
%[new_x, ERR] = gradient_decent(b,A,abs((max(abs(A'*b))/(9.8*(coin^0.735)))),new_x,10^-8,coin);
%[new_x, ERR] = gradient_decent(b,A,abs((max(abs(A'*b))/(10*(coin^0.775)))),new_x,10^-8,coin);
% func = @(z)(0.5*((norm(A*z-b,2))^2));
%     options = optimset('Display','off');
%     opt_x = fminsearch(func,new_x,options);
max(abs(A'*b));
% x1 = new_x;
x1 = new_x;
err = norm(x1-answ,2);
sucess = 0;
[value,myans] = maxk(abs(x1),coin);
figure;
plot([1:length(tao_record)],tao_record);
title('tao record')
figure;
plot([1:length(alpha_record)],alpha_record);
title('alpha record')
figure;
plot([1:length(ERR)],ERR);
title('error record')
(coin-length(setdiff(fake,myans)))/coin;
if(length(setdiff(fake,myans))<=0) 
    sucess = 1;
end
sucess;
% for i = 1:coin
% 
%          if(ismember(myans(i),fake) ~= 1)
%              sucess = 0;
%              myans(i)
%          end;
% 
% end
sucess
if((norm(A*x1)>coin)&&sucess==0)
   norm(A*x1) ;
end
norm(A*x1-b) ;
fake;
myans;
org_mse = immse(zeros(100,1),answ);
mse = immse(x1,answ);
