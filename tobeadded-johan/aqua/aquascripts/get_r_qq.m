%% written by Norman Zacharias
% Evaluation of aqua.m output
N=42;
for ii=1:N
    r_qq(ii,1) = QA(1,ii).stat.r_qq;
    PSC(ii,1) = QA(1,ii).stat.psc;
    disp(QA(1,ii).stat.run)
end




