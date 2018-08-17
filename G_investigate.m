% test = zeros(125*9,125*9);
% for temp = 1:125*9
%     bavec = zeros(125*9,1);
%     bavec(temp) = 1;
%     test(:,temp) = G_dP(bavec);
% end
% test1 = test(1:375,1:375);
% test2 = test(376:750,376:750);
% test3 = test(751:1125,751:1125);

test_find = find(abs(test1-test1(26,1))<1e-10);
test_find(:,2) = 0;
test_find(:,3) = 0;
test_find(:,2) = mod(test_find(:,1),375);
test_find(:,3) = ( test_find(:,1) - test_find(:,2) ) / 375 + 1;

test_unique = uniquetol(test1);
test_unique(:,2) = 0;
for temp = 1:length(test_unique)
    temp2 = find( abs( test1-test_unique( temp,1 ) ) < 1e-10 );
    test_unique(temp,2) = length(temp2);
end

[test_C,test_ia,test_ic] = uniquetol(test1);
test11 = zeros(375,375);
for temp = 1:140625
    test11(temp) = test_C(test_ic(temp));
end