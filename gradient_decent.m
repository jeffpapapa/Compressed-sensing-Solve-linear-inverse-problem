function [new_x, ERR,tao_record,alpha_record] = gradient_decent(u,phi,org_tao,initial,threshold,coin)
    tao_record = 0;
    alpha_record = 0;
    alpha_k = 1;
    alpha=(norm(phi'*phi)/2)+1;
    lower_alpha = (norm(phi'*phi)/2);
    upper_alpha = (norm(phi'*phi)/2)+301;
    (norm(phi'*phi)/2);
    gamma = 1;
    x = initial;
    k = 1;
    ERR=[1];
    little_sigma = 0.1;
    
    org_tao;
    
    %tasi = 0.375;
    tasi = 0.375;
    new_u = u;
    tao = max( tasi*(norm(phi'*new_u,inf)),org_tao);
    tao_record(k) = tao;
    %tao = org_tao;
    alpha_record(alpha_k) = alpha;
    alpha_k = alpha_k + 1;
    ut = (x-(phi')*(phi*x-u)/alpha);
    opt_x = sign(ut).*max(abs(ut) - (tao/alpha)/2,0);
    %opt_x = ut*( (max(norm(ut,1)-(tao/alpha),0))/ ((max(norm(ut,1)-(tao/alpha),0)) + (tao/alpha))  );
            
    %%%check if opt_x satisfied condition
    %%%if not, then we need to find new opt_x with new alpha
    history_phi_transfer(1) = 0.5*(norm(u-phi*x)^2)+tao*norm(x,1);
    phi_optx = 0.5*(norm(u-phi*opt_x)^2)+tao*norm(opt_x,1);
    check_optx = history_phi_transfer(1) - (little_sigma/2)*alpha*(norm(opt_x-x)^2);
    while(phi_optx > check_optx)
       alpha =  alpha*1.5;
       alpha_record(alpha_k) = alpha;
       alpha_k = alpha_k + 1;
       ut = (x-(phi')*(phi*x-u)/alpha);
       opt_x = sign(ut).*max(abs(ut) - (tao/alpha)/2,0);
       %opt_x = ut*( (max(norm(ut,1)-(tao/alpha),0))/ ((max(norm(ut,1)-(tao/alpha),0)) + (tao/alpha))  );
       phi_optx = 0.5*(norm(u-phi*opt_x)^2)+tao*norm(x,1);
       check_optx = history_phi_transfer(1) - (little_sigma/2)*alpha*(norm(opt_x-x)^2);
    end
    
        
    new_x = x + gamma*(opt_x-x);
    new_u = u-phi*new_x;
    phi_old = 0.5*(norm(u-phi*x)^2)+tao*norm(x,1);
    phi_new = 0.5*(norm(u-phi*new_x)^2)+tao*norm(new_x,1);
    error = abs(phi_new-phi_old)/phi_old;
    %error = norm(new_x-x)/norm(new_x);
    ERR(1)=error;
    set(1) = norm(x-new_x);
    max_k = coin*650;
%     max_k = 0.5*(10^coin);
%     if(max_k<100)
%        max_k = 100; 
%     end
    x_value = maxk(abs(new_x),coin+1);
    x_max = x_value(coin);
    x_topcoin = x_value(coin+1);
    stop_cond = 0;
    if((x_max-x_topcoin >0.001)&&(x_topcoin<0.001))
           stop_cond = 1; 
    end
%     sol_err = norm(phi*new_x-u);
%     
%     old_set = 0;
%     old_set_count = 0;
%     for i=1:length(x)
%        if(new_x(i)~=0)
%            old_set_count = old_set_count+1;
%            old_set( old_set_count ) = i;
%        end
%     end
%     
%    new_set = old_set;
%    new_set_count = old_set_count;
%    terminate_error = new_set_count/(10^-3)
    
%     v = new_x;
%     w = -new_x;
%     for i=1:length(new_x)
%        if(v(i)<0)
%            v(i)=0;
%        end
%        if(w(i)<0)
%            w(i)=0;
%        end
%     end
%     terminate_error = (norm(u-phi*x,2)^2)/(norm(u-phi*(v-w),2)^2)
    %terminate_error = norm(v)*norm(w)
    %&&(ter2_error > threshold)
    %ter2_error = norm(tao*(ones(length(x)))+phi'*(phi*(v-w)-u))*norm(tao*(ones(length(x)))-phi'*(phi*(v-w)-u))
    %while((terminate_error > threshold))
    %while((error>threshold)&&(stop_cond==0)&&(sol_err>0.25))
    while((error>threshold)&&(k<max_k)&&(stop_cond==0))
        k = k+1;
        s = new_x - x;
        alpha=( ( (norm((phi'*phi*s),2))^2) /  ((norm(phi*s,2))^2)  );
        
        if(alpha > upper_alpha)
            alpha = upper_alpha;
        end
        if(alpha < lower_alpha)
            alpha =lower_alpha;
        end
        alpha_record(alpha_k) = alpha;
        alpha_k = alpha_k + 1;
        x = new_x;
        if(tao~=org_tao)
            new_u = u - phi*x;
            %tao = org_tao;
            k;
            tao = max( tasi*(norm(phi'*new_u,inf)),org_tao);
        end
        tao_record(k) = tao;
        %func = @(z)( ((z-x)')*((phi')*(phi*x-u)) + 0.5*alpha*((norm(z-x)^2)) + tao*(norm(z,1))  );
        %func = @(z)(0.5*((norm(z-(x-(phi')*(phi*x-u)/alpha),2))^2) + (tao/alpha)*(norm(z,1))    );
        %opt_x = fminsearch(func,initial,options);
        ut = (x-(phi')*(phi*x-u)/alpha);
        opt_x = sign(ut).*max(abs(ut) - (tao/alpha)/2,0);
        %opt_x = ut*( (max(norm(ut,1)-(tao/alpha),0))/ ((max(norm(ut,1)-(tao/alpha),0)) + (tao/alpha))  );
            
        %%%we need to find if opt_x satisfies some conditions
        %%%if not, then we should choose new alpha and find a new opt_x
        history_phi_transfer(k) = 0.5*(norm(u-phi*x)^2)+tao*norm(x,1);
        phi_optx = 0.5*(norm(u-phi*opt_x)^2)+tao*norm(opt_x,1);
        check_optx = max(history_phi_transfer) - (little_sigma/2)*alpha*(norm(opt_x-x)^2);
        
        while(phi_optx > check_optx)
            alpha =  alpha*1.5
            alpha_record(alpha_k) = alpha;
            alpha_k = alpha_k + 1;
            ut = (x-(phi')*(phi*x-u)/alpha);
            opt_x = sign(ut).*max(norm(ut,1) - (tao/alpha)/2,0);
            %opt_x = ut*( (max(norm(ut,1)-(tao/alpha),0))/ ((max(norm(ut,1)-(tao/alpha),0)) + (tao/alpha))  );
            phi_optx = 0.5*(norm(u-phi*opt_x)^2)+tao*norm(opt_x,1);
            check_optx = history_phi_transfer(1) - (little_sigma/2)*alpha*(norm(opt_x-x)^2);
        end
        
        
        set(k) = max(abs(x-opt_x));
        new_x = x + gamma*(opt_x-x);
        phi_old = 0.5*(norm(u-phi*x)^2)+tao*norm(x,1);
        phi_new = 0.5*(norm(u-phi*new_x)^2)+tao*norm(new_x,1);
        %error = abs(phi_new-phi_old)/phi_old;
%         sol_err = norm(phi*new_x-u)
        error = norm(new_x-x)/(norm(new_x));
        ERR(k)=error;
        x_value = maxk(abs(new_x),coin+1);
        x_max = x_value(coin);
        x_topcoin = x_value(coin+1);
        if((x_max-x_topcoin >0.02)&&(x_topcoin<0.02))
           stop_cond = 1;
        end
%         if(sol_err<=0.3)
%            sol_err ;
%         end
%     old_set = new_set;
%     old_set_count = new_set_count;
%     new_set = 0;
%     new_set_count = 0;
%     for i=1:length(x)
%        if(new_x(i)~=0)
%            new_set_count = new_set_count +1;
%            new_set( new_set_count ) = i;
%        end
%     end
%     check_count = 0;
%     for i=1:length(x)
%         if( ((ismember(i,new_set)==1)&&(ismember(i,old_set)~=1)) || ((ismember(i,new_set)~=1)&&(ismember(i,old_set)==1)) )
%             check_count = check_count + 1;
%         end
%     end
%     check_count
%     new_set_count
%    terminate_error = check_count/new_set_count
        
        
        
        
%          v = new_x;
%          w = -new_x;
%         for i=1:length(new_x)
%         if(v(i)<0)
%              v(i)=0;
%          end
%         if(w(i)<0)
%             w(i)=0;
%         end
%         end
        
        %terminate_error = (norm(u-phi*x,2)^2)/(norm(u-phi*(v-w),2)^2)
        %terminate_error = norm(v)*norm(w)
        %ter2_error = norm(tao*(ones(length(x)))+phi'*(phi*(v-w)-u))*norm(tao*(ones(length(x)))-phi'*(phi*(v-w)-u))
        
    end
%     x_max;
%     x_topcoin;
%     sol_err;
    
end