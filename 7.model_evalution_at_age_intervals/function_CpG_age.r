age.predict.mn.f = function(Mat, ind, CpG_list) {
            X1v = Mat[ind, CpG_list[1]];
            X3v = Mat[ind, CpG_list[2]];
            X4v = Mat[ind, CpG_list[3]];
            X5v = Mat[ind, CpG_list[4]];
            X6v = Mat[ind, CpG_list[5]];
            X7v = Mat[ind, CpG_list[6]];
            X8v = Mat[ind, CpG_list[7]];
            X9v = Mat[ind, CpG_list[8]];
            Y5v = Mat[ind, CpG_list[9]];
            Y3v = Mat[ind, CpG_list[10]];
            Y2v = Mat[ind, CpG_list[11]];
            Y4v = Mat[ind, CpG_list[12]];
            age_p = 51.8299137765218-0.0894037801682898*X1v-0.387451058991341*X3v-3.38269035822311*X4v+2.11279716245427*X5v
            -0.118310126060979*(X6v+X7v)/2-1.90897713105373*X8v+2.11361675771719*X9v-0.662488450736703*Y5v
            +2.00567273735279*Y3v+1.86692955317212*Y2v-0.541876580154293*Y4v;
            c(age_p);
            }
            
age.predict.m1.f = function(Mat, ind, CpG_list) {
            X1v = Mat[ind, CpG_list[1]];
            X2v = Mat[ind, CpG_list[2]];
            X3v = Mat[ind, CpG_list[3]];
            X4v = Mat[ind, CpG_list[4]];
            X5v = Mat[ind, CpG_list[5]];
            X6v = Mat[ind, CpG_list[6]];
            X7v = Mat[ind, CpG_list[7]];
            X8v = Mat[ind, CpG_list[8]];
            X9v = Mat[ind, CpG_list[9]];
            age_p = 29.31+31.863*X1v-37.797*X2v-17.163*X3v-20.293*X4v+32.982*X5v+29.76*X6v+33.684*X7v-28.528*X8v+33.212*X9v;
            c(age_p);
            }
            
age.predict.m2.f = function(Mat, ind, CpG_list) {
            Y1v = Mat[ind, CpG_list[1]];
            Y2v = Mat[ind, CpG_list[2]];
            Y3v = Mat[ind, CpG_list[3]];
            Y4v = Mat[ind, CpG_list[4]];
            Y5v = Mat[ind, CpG_list[5]];
            Y6v = Mat[ind, CpG_list[6]];
            age_p = 72.8800-43.2336*Y1v+140.2192*Y2v+44.5862*Y3v+9.3521*Y4v-13.8511*Y5v-64.9604*Y6v
            c(age_p);
            }
            
normalize.rows.f = function (Mat) {
    if (length(nrow(Mat)) != 0) {
        t(apply(Mat, 1, function(y) (y-mean(y))/sd(y))) 
    } else {
        c(apply(t(Mat), 1, function(y) (y-mean(y))/sd(y))) 
    }
}
