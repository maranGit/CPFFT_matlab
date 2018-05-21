fid = fopen('temp','wt');
i2v = [1,4,6; 4,2,5; 6,5,3];
i2f = [1,2,3; 4,5,6; 7,8,9];
for a=1:3
    for B=1:3
        for c=1:3
            for D=1:3
                temp = [i2f(a,B),i2f(c,D),B,1,i2v(a,1),i2v(c,1),D,1,B,1,i2v(a,1),i2v(c,2),D,2,B,1,i2v(a,1),i2v(c,3),D,3,B,2,i2v(a,2),i2v(c,1),D,1,B,2,i2v(a,2),i2v(c,2),D,2,B,2,i2v(a,2),i2v(c,3),D,3,B,3,i2v(a,3),i2v(c,1),D,1,B,3,i2v(a,3),i2v(c,2),D,2,B,3,i2v(a,3),i2v(c,3),D,3];
                fprintf(fid,'A(%i,%i) = fn1inv(%i,%i)*cep(%i,%i)*fn1inv(%i,%i) + fn1inv(%i,%i)*cep(%i,%i)*fn1inv(%i,%i) + fn1inv(%i,%i)*cep(%i,%i)*fn1inv(%i,%i) + fn1inv(%i,%i)*cep(%i,%i)*fn1inv(%i,%i) + fn1inv(%i,%i)*cep(%i,%i)*fn1inv(%i,%i) + fn1inv(%i,%i)*cep(%i,%i)*fn1inv(%i,%i) + fn1inv(%i,%i)*cep(%i,%i)*fn1inv(%i,%i) + fn1inv(%i,%i)*cep(%i,%i)*fn1inv(%i,%i) + fn1inv(%i,%i)*cep(%i,%i)*fn1inv(%i,%i);\n',temp);
            end
        end
    end
end
fclose(fid);
