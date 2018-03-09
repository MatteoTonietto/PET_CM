function est = calcMacro(est)

model_name  = est.FUN;
modelFormat = est.fixed_par.modelFormat;  

par = est.par;
Cov = est.Cov;

switch model_name
    case {'model_2TCM','model_2TCM1K'}
        if strcmp(modelFormat,'K')
            K1 = par(1);
            k2 = par(2);
            k3 = par(3);
            k4 = par(4);
            
            %--------------------------------------------------------------
            Vt = (K1/k2)*(1 + k3/k4);
            
            D_VtvsK1 = 1/k2*(1 + k3/k4);
            D_Vtvsk2 = -K1/(k2^2)*(1 + k3/k4);
            D_Vtvsk3 = K1/(k2*k4);
            D_Vtvsk4 = -K1*k3/(k2*(k4^2));
            
            VAR_Vt = 2*D_VtvsK1*D_Vtvsk2*Cov(1,2) ...
                +2*D_Vtvsk3*D_Vtvsk4*Cov(3,4) ...
                +2*D_VtvsK1*D_Vtvsk3*Cov(1,3) ...
                +2*D_VtvsK1*D_Vtvsk4*Cov(1,4) ...
                +2*D_Vtvsk2*D_Vtvsk3*Cov(2,3) ...
                +2*D_Vtvsk2*D_Vtvsk4*Cov(2,4) ...
                +D_VtvsK1*D_VtvsK1*Cov(1,1)   ...
                +D_Vtvsk2*D_Vtvsk2*Cov(2,2)   ...
                +D_Vtvsk3*D_Vtvsk3*Cov(3,3)   ...
                +D_Vtvsk4*D_Vtvsk4*Cov(4,4)   ;
            
            Vt_cv = 100*(sqrt(VAR_Vt)/Vt);
            
            %--------------------------------------------------------------
            BP = k3/k4;
            
            D_BPvsk3 = 1/k4;
            D_BPvsk4 = -k3/(k4^2);
            
            VAR_BP = D_BPvsk3*D_BPvsk3*Cov(3,3)    ...
                +D_BPvsk4*D_BPvsk4*Cov(4,4)   ...
                +2*D_BPvsk3*D_BPvsk4*Cov(3,4) ;
            
            BP_cv = 100*(sqrt(VAR_BP)/BP);
            
            %--------------------------------------------------------------
            est.Vt    = Vt;
            est.Vt_cv = Vt_cv;
            est.BP    = BP;
            est.BP_cv = BP_cv;
            
        else
            
        end
    case 'model_1TCM'
        if strcmp(modelFormat,'K')
            K1 = par(1);
            k2 = par(2);
            
            %--------------------------------------------------------------
            Vt    = K1/k2;
            
            D_VtvsK1 = 1/k2;
            D_Vtvsk2 = -K1/(k2^2);
            
            VAR_Vt = 2*D_VtvsK1*D_Vtvsk2*Cov(1,2) ...
                +D_VtvsK1*D_VtvsK1*Cov(1,1)   ...
                +D_Vtvsk2*D_Vtvsk2*Cov(2,2) ;
            
            Vt_cv = 100*(sqrt(VAR_Vt)/Vt);
            
            est.Vt    = Vt;
            est.Vt_cv = Vt_cv;
            
        else
        end
end
        

