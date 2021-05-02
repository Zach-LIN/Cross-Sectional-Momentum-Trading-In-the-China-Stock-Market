clear;
load('CloseData2.mat');
load('Top50.mat');
%Select Data from former result
for i = 1:500
    stock{i} = string([stock{i}(1:2),'.',stock{i}(3:end)]);
end

newData = [];
%for i = 1:50
    for j = 1:500
        %if(Top50{i,1}==stock{j})
            newData = [newData,data(:,j)];
        %end
    end
%end
num = 500;

res = zeros(40,7);


for lookback=[1 2 3 6 12]
    j = find([1 2 3 6 12]==lookback);
    for holdmons=[1 2 3 4 6 8 12]
        k = find([1 2 3 4 6 8 12]==holdmons);
        temp = [];
        t = [];
        for now=245:20*holdmons:245+23*20 %245 is start date of 2018 and we train on 24 mons
            rtn_lag_1 = newData(now-1,:)./newData(now-lookback*20-1,:);
            rtn_lag_2 = newData(now-1,:)./max(newData(now-lookback*20-1:now-1,:));
            rtn_fut = newData(now+holdmons*20,:)./newData(now,:);
            [rtn_lag_1_sorted,I1] = sort(rtn_lag_1);
            [rtn_lag_2_sorted,I2] = sort(rtn_lag_2);
            r1 = 1:length(rtn_lag_1);
            r1(I1) = r1;
            r2 = 1:length(rtn_lag_2);
            r2(I2) = r2;
            r = r1+r2;
            [~,I] = sort(r);
            rtn_fut_sorted = rtn_fut(I2);
            res(j*8-7:j*8-3,k) = res(j*8-7:j*8-3,k) + [ mean(rtn_fut_sorted(1:num/5));...
                                                        mean(rtn_fut_sorted(1*num/5+1:2*num/5));...
                                                        mean(rtn_fut_sorted(2*num/5+1:3*num/5));...
                                                        mean(rtn_fut_sorted(3*num/5+1:4*num/5));...
                                                        mean(rtn_fut_sorted(4*num/5+1:num))];
            temp = [temp,mean(rtn_fut_sorted(4*num/5+1:num)) - mean(rtn_fut_sorted(1:num/5))];
            t = [t;[mean(rtn_fut_sorted(4*num/5+1:num)) , mean(rtn_fut_sorted(1:num/5))]];
            [h,p] = ttest2(t(:,1),t(:,2),'Vartype','unequal');
            if(p>0.3)
                p = -p;
            end
            res(j*8-1,k) = p;
            
        end
        temp = (1+temp).^(1/holdmons)-1;
        res(j*8-7:j*8-3,k) = res(j*8-7:j*8-3,k)/length(temp);
        res(j*8-2,k) = mean(temp);
        res(j*8,k) = sum((temp-mean(temp)).^2);
        res(j*8,k) = sqrt(res(j*8,k)/(length(temp)-1));
        rf = data(now+holdmons*20,:)./data(now,:);
        rf = mean(rf.^(1/holdmons))-1;
        res(j*8,k) = (res(j*8-2,k)-rf)./res(j*8,k);
        
        
    end
end

    
