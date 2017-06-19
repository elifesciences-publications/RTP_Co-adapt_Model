%--------------------------------------------------------------------------
% This MATLAB code simulates topographic axonal mapping along the anterior/
% posterior axis of the retinotectal system. It is based on the model 
% published in Gebhardt, C., et al., Development, 139, 335 (2012) and 
% modified to include growth cone adaptation. It is associated with the 
% paper Fiedeling, F., et al., „Ephrin-A/EphA specific co-adaptation as a 
% novel mechanism in topographic axon guidance“.

% Copyright (C) 2017, Franco Weth (franco.weth@kit.edu)

% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU general Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version. This program is distributed in the hope that 
% it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
% GNU General Public License for more details.
%--------------------------------------------------------------------------

%%
clear all

% PROGRAM VERSION
%-------------------Version 1.2------------------------------------07.2016

% GENERAL PARAMETERS

NoGrowthCone    = 200;  
SizeGrowthCone  = 3;
offset          = (SizeGrowthCone-1)/2;
GCcutoff        = 0.01;
steps           = 30000;        % i
Qx              = 0;         
Qy              = 0;
sigma           = 0.12;          
mu              = 0.006;              
lambda          = 0.0045;    
knockIn         = 0;                  
cis_factor      = 1;
pre_adap        = 1;            
no_adap         = 10;
x_shift         = 5;            
C_dynamic       = 1;            
C               = 100;          % Fiber-Fiber factor C(i)

Pedestal_Receptor_Retina=0;     
Pedestal_Ligand_Retina=0;       
Pedestal_Receptor_Target=0;
Pedestal_Ligand_Target=0;

FieldSizeX      = 50;
FieldSizeXtd    = FieldSizeX + 2*offset;
o               = 100/FieldSizeX;
FieldSizeY      = 8;
FieldSizeYtd    = FieldSizeY + 2*offset;
kappa_retina    = o*0.025;       
omega_retina    = 0.4;           

ftw = 0;                         % 1: GC forced to take step in x
adap           = 1;              % 1:on
adapHistory     = 10;            % k           
               

%% DEFINING VECTORS AND MATRICES

YDrang = [(1-Qy)/3 1/3+(1-Qy)/3 1];

GrowthConeLigand   = zeros(FieldSizeXtd,FieldSizeYtd);
Lxy   = zeros(FieldSizeXtd,FieldSizeYtd);

GrowthConeReceptor = zeros(FieldSizeXtd,FieldSizeYtd);
Rxy = zeros(FieldSizeXtd,FieldSizeYtd);

AxonReceptor = zeros(1,NoGrowthCone);
AxonLigand   = zeros(1,NoGrowthCone);
AxonReceptor_REF = zeros(1,NoGrowthCone);
AxonLigand_REF   = zeros(1,NoGrowthCone);

xtHistory  = zeros(steps,NoGrowthCone);
ytHistory  = zeros(steps,NoGrowthCone);

DxHistory       = zeros(steps,NoGrowthCone);
AbsDxHistory    = zeros(steps,NoGrowthCone);
QxHistory       = zeros(steps,NoGrowthCone);
FactorHistory   = zeros(steps,1);

adapmeandenom   = sum(1:adapHistory);

ValAdapRec         = zeros(1,NoGrowthCone);
ValAdapLig         = zeros(1,NoGrowthCone);
ValAdapRecHistory  = zeros(steps,NoGrowthCone);
ValAdapLigHistory  = zeros(steps,NoGrowthCone);
AbsAdapRec         = zeros(1,NoGrowthCone);
AbsAdapLig         = zeros(1,NoGrowthCone);
AbsAdapRecHistory  = zeros(1,NoGrowthCone);
AbsAdapLigHistory  = zeros(1,NoGrowthCone);

ValResRec          = zeros(steps,NoGrowthCone); 
ValResLig          = zeros(steps,NoGrowthCone); 
ValResRecHistory   = zeros(steps,NoGrowthCone); 
ValResLigHistory   = zeros(steps,NoGrowthCone); 

AdapCoeff       = zeros(steps,NoGrowthCone);
AdapmeanCoeff   = ones (steps,NoGrowthCone);
ResRecCoeff     = ones (steps,NoGrowthCone);
ResLigCoeff     = ones (steps,NoGrowthCone);

ReceptorHistory = zeros(steps,NoGrowthCone);
LigandHistory   = zeros(steps,NoGrowthCone);

GC_GCfactor_History = zeros(1,steps);

FSfac= 50/FieldSizeX; 

%% SUBSTRATES

Substrate_Tectum;                          % n-t/a-p mapping (Figure 4B)

%Substrate_Gap_Assay;                      % Gap assays (Figure 4C, D)
%Substrate_Tectal_Innervation;             % Tectal innervation (Figure 10)


%% ALLOCATION OF STARTPOSITIONS

if NoGrowthCone == 1
     YStartPos = ceil(FieldSizeX/2);
else
      YStartPos =  round(linspace(offset,FieldSizeX,NoGrowthCone));    
end

a=0.5;c=0.5;
[W, V] = meshgrid(1:FieldSizeY, 1:FieldSizeX);
gaussian=@(x0,y0) exp( - (a*(V-x0).^2 + c*(W-y0).^2));  
weightxtd=zeros(FieldSizeXtd, FieldSizeYtd);

tic
datestr(now);

%% INITIALISATION

for n=1:NoGrowthCone
    xt  = 1+x_shift+offset;   
    
    if NoGrowthCone == 1
        yt  = ceil(FieldSizeY/2);
    else
        yt  = round(((FieldSizeY-1)/(NoGrowthCone-1))*...
              n+((NoGrowthCone-FieldSizeY)/(NoGrowthCone-1))) + offset;
    end
    

     AxonReceptor(n) = (exp(FSfac*0.05*(YStartPos(n)-FieldSizeX/2))+...
         Pedestal_Receptor_Retina)*pre_adap;
     AxonLigand(n)   = (exp(FSfac*-0.05*(YStartPos(n)-1-FieldSizeX/2))+...
         Pedestal_Receptor_Retina)*pre_adap;
     
     
   
if knockIn > 0
    for f=1:floor(NoGrowthCone/2)
        if n==2*f                                 %Isl EphA3 knockIn
            AxonReceptor(n) = AxonReceptor(n)+knockIn;
            AxonLigand(n)   = 1/AxonReceptor(n);
            break
        else                                      %wildtype
            AxonReceptor(n) = AxonReceptor(n);              
            AxonLigand(n)   = AxonLigand(n);
        end
    end
    
end


     ReceptorHistory(1,n)= AxonReceptor(n);
     LigandHistory(1,n)= AxonLigand(n);
     
     AxonReceptor_REF(n)= AxonReceptor(n);
     AxonLigand_REF(n)  = AxonLigand(n);
     
   
     
    
    xrandom = rand;
    if(xrandom < (1-Qx)/SizeGrowthCone)
        xtDirection = -1;
    elseif(xrandom < 1/SizeGrowthCone+(1-Qx)/SizeGrowthCone)
        xtDirection = 0;
    else
        xtDirection = 1;
    end

    yrandom = rand;
    if(yrandom < YDrang(1))
        ytDirection = -1;
    elseif(yrandom < YDrang(2))
        ytDirection = 0;
    else
        ytDirection = 1;
    end

    if xt+xtDirection<1+offset
        xtDirection=0;
    elseif xt+xtDirection>FieldSizeX+offset
        xtDirection=0;
    end
    if yt+ytDirection<1+offset
        ytDirection=0;
    elseif yt+ytDirection>FieldSizeY+offset
        ytDirection=0;
    end

    
    
    Lxy = SubstrateLigand;
    Rxy = SubstrateReceptor;
    
    Dx1=abs(log(((AxonReceptor(n)*(Lxy(xt,yt)+AxonLigand(n)))/...
        (AxonLigand(n)*(Rxy(xt,yt)+AxonReceptor(n))))));
    Dx2=abs(log(((AxonReceptor(n)*(Lxy(xt+xtDirection,yt+ytDirection)+...
        AxonLigand(n)))/(AxonLigand(n)*(Rxy(xt+xtDirection,yt+...
        ytDirection)+AxonReceptor(n))))));
     
    DxHistory(1,n) = Dx1;
   
    AdapCoeff(1,n) = 1;
                 
    WDx1 = wkeitpd(Dx1,sigma);
    WDx2 = wkeitpd(Dx2,sigma);

    if ftw == 1
    xt=xt+1;        
    yt=yt;            
    
    elseif ftw == 0
    
        if rand>wkeit01(WDx1,WDx2)
            xt = xt+xtDirection;
            yt = yt+ytDirection;
        end
    
    end
    
    xtHistory(1,n) = xt;
    ytHistory(1,n) = yt;
end


%% ITERATION

for ii=2:steps
    
    if C_dynamic == 1
        GC_GCfactor=C*(-exp(-log(2.^((ii./(steps./4)).^5)))+1);
    elseif C_dynamic == 0
        GC_GCfactor=C;
    end
        
    GC_SUBfactor=1;
    FactorHistory(ii)=GC_GCfactor;     

    allGrowthConeLigand   = zeros(FieldSizeXtd,FieldSizeYtd);
    allGrowthConeReceptor = zeros(FieldSizeXtd,FieldSizeYtd);
    

    for nn=1:NoGrowthCone 
                
        xn=xtHistory(ii-1,nn)-1;
        yn=ytHistory(ii-1,nn)-1;
        weight=gaussian(xn,yn);        
        weight(weight<GCcutoff)=0;
        weightxtd(2:FieldSizeXtd-1,2:FieldSizeYtd-1)= ...
                   weight(1:FieldSizeX,1:FieldSizeY);
        
        allGrowthConeLigand   = allGrowthConeLigand   + ...
            AxonLigand(nn)  .*weightxtd;
        allGrowthConeReceptor = allGrowthConeReceptor + ...
            AxonReceptor(nn).*weightxtd;
        
    end
    
    
    for n=1:NoGrowthCone  
     
        if ii<=no_adap                                                            
            ReceptorHistory(ii,n)=ReceptorHistory(1,n);                                
            LigandHistory(ii,n)=LigandHistory(1,n);                                        
        end                                                                                

        xt = xtHistory(ii-1,n);
        yt = ytHistory(ii-1,n);       
        
        weight=gaussian(xt-1,yt-1);
        weight(weight<GCcutoff)=0;
        weightxtd(2:FieldSizeXtd-1,2:FieldSizeYtd-1) = ...
                   weight(1:FieldSizeX,1:FieldSizeY);
        
        currentGrowthConeLigand   = AxonLigand(n)  .*weightxtd;
        currentGrowthConeReceptor = AxonReceptor(n).*weightxtd;     
        
        meancurrentGCLigand   = sum(sum(currentGrowthConeLigand))/...
                                nnz(currentGrowthConeLigand);
        meancurrentGCReceptor = sum(sum(currentGrowthConeReceptor))/...
                                nnz(currentGrowthConeReceptor);
        
            
        currentGCLigandRef   = AxonLigand_REF(n)  .*weightxtd;                               
        currentGCReceptorRef = AxonReceptor_REF(n).*weightxtd;                                  
            
        meancurrentGCLigandRef   = sum(sum(currentGCLigandRef))/...
                                   nnz(currentGCLigandRef);
        meancurrentGCReceptorRef = sum(sum(currentGCReceptorRef))/...
                                   nnz(currentGCReceptorRef);
       
        
        GrowthConeLigand   = allGrowthConeLigand   - ...
            currentGrowthConeLigand;
        GrowthConeReceptor = allGrowthConeReceptor - ...
            currentGrowthConeReceptor;
      

        xrandom = rand;
        if(xrandom < (1-Qx)/SizeGrowthCone)
            xtDirection = -1;
        elseif(xrandom < 1/SizeGrowthCone+(1-Qx)/SizeGrowthCone)
            xtDirection = 0;
        else
            xtDirection = 1;
        end

        yrandom = rand;
        if(yrandom < YDrang(1))
            ytDirection = -1;
        elseif(yrandom < YDrang(2))
            ytDirection = 0;
        else
            ytDirection = 1;
        end

        if xt+xtDirection<1+offset
            xtDirection=0;
        elseif xt+xtDirection>FieldSizeX+offset
            xtDirection=0;
        end
        if yt+ytDirection<1+offset
            ytDirection=0;
        elseif yt+ytDirection>FieldSizeY+offset
            ytDirection=0;
        end                 
        

          
        rev = GC_SUBfactor.*currentGrowthConeLigand.*SubstrateReceptor+...
            GC_GCfactor.*currentGrowthConeLigand.*GrowthConeReceptor+...
            cis_factor.*currentGrowthConeLigand.*...
             currentGrowthConeReceptor;
        
        
        fwd = GC_SUBfactor.*currentGrowthConeReceptor.*SubstrateLigand+...
             GC_GCfactor.*currentGrowthConeReceptor.*GrowthConeLigand+...
             cis_factor.*currentGrowthConeReceptor.*...
             currentGrowthConeLigand;
         
    
        fwdmean=sum(sum(fwd))/nnz(fwd);
        revmean=sum(sum(rev))/nnz(rev);
        
        revHistory(ii,n)=revmean;
        fwdHistory(ii,n)=fwdmean;
        
        Dx=abs(log(revmean/fwdmean));
          
        DxHistory(ii,n) = Dx;

        
        weight=gaussian(xt-1+xtDirection,yt-1+ytDirection);
        weight(weight<GCcutoff)=0;
        weightxtd(2:FieldSizeXtd-1,2:FieldSizeYtd-1) = ...
                   weight(1:FieldSizeX,1:FieldSizeY);
               
        currentGrowthConeLigand_target   = AxonLigand(n)  .*weightxtd;
        currentGrowthConeReceptor_target = AxonReceptor(n).*weightxtd;
        
         

        rev_target = GC_SUBfactor.*currentGrowthConeLigand_target.*...
            SubstrateReceptor+GC_GCfactor.*...
            currentGrowthConeLigand_target.*GrowthConeReceptor+...
            cis_factor.*currentGrowthConeLigand_target.*...
            currentGrowthConeReceptor_target;
        
        fwd_target = GC_SUBfactor.*currentGrowthConeReceptor_target.*...
            SubstrateLigand+GC_GCfactor.*...
            currentGrowthConeReceptor_target.*GrowthConeLigand+...
             cis_factor.*currentGrowthConeReceptor_target.*...
             currentGrowthConeLigand_target;
        
        fwdmean_target=sum(sum(fwd_target))/nnz(fwd_target);
        revmean_target=sum(sum(rev_target))/nnz(rev_target);
        
        Dx_target=abs(log(revmean_target/fwdmean_target));

        AbsDxHistory(ii,n) = log(revmean_target/fwdmean_target);
        
        
%% --------------------------Adaptation------------------------------------

        if adap==1 
            
            adaprev = currentGrowthConeLigand.*...
                      (currentGrowthConeReceptor+SubstrateReceptor);
            adapfwd = currentGrowthConeReceptor.*...
                      (currentGrowthConeLigand+SubstrateLigand);
            adaprevmean = sum(sum(adaprev))/nnz(adaprev);
            adapfwdmean = sum(sum(adapfwd))/nnz(adapfwd);
            
            AdapCoeff(ii,n)    = 1+log((1+(mu*(abs(log(adaprevmean/...
                adapfwdmean))))));    
         

            ResRecCoeff(ii,n) = lambda*(meancurrentGCReceptorRef-...
                                meancurrentGCReceptor);     
            ResLigCoeff(ii,n) = lambda*(meancurrentGCLigandRef-...
                                meancurrentGCLigand);

            ResRecCoeff(2,n) = 1;       
            ResLigCoeff(2,n) = 1;
                
                
            if ii>adapHistory && ii>no_adap                               
                
                AdapmeanCoeff(ii,n) = 0;
               
                for k=0:adapHistory
                    AdapmeanCoeff(ii,n) = AdapmeanCoeff(ii,n) + ...
                                          k*AdapCoeff(ii-adapHistory+k,n);    
                end
                
               AdapmeanCoeff(ii,n) = AdapmeanCoeff(ii,n)/adapmeandenom;  

               AxonReceptor(n) = AxonReceptor(n)*...
                                 AdapmeanCoeff(ii,n)+ResRecCoeff(ii,n);      
               AxonLigand(n)   = AxonLigand(n)*...
                                 AdapmeanCoeff(ii,n)+ResLigCoeff(ii,n);
               ReceptorHistory(ii,n) = AxonReceptor(n);
               LigandHistory(ii,n)   = AxonLigand(n);
                
            
            end
            
% canonical, unproportional adaptation of fwd and rev signals
 
%         AdapCoeff_rec(ii,n)=(log(1/ReceptorHistory(1,n)*...
%                             SubstrateLigand(xt,yt))); 
%         AdapCoeff_lig(ii,n)=(log(1/LigandHistory(1,n)*...
%                             SubstrateReceptor(xt,yt)));
%                         
%         AdapCoeff_rec(ii,n)=(AdapCoeff_rec(ii,n)-(2*AdapCoeff_rec(ii,n)));
%         AdapCoeff_lig(ii,n)=(AdapCoeff_lig(ii,n)-(2*AdapCoeff_lig(ii,n)));
%  
%         
%         AxonReceptor(n)=((AxonReceptor(n)-SubstrateReceptor(xt,yt))*...
%                         exp(-tau*AdapCoeff_rec(ii,n)^2))+...
%                         SubstrateReceptor(xt,yt); 
%         AxonLigand(n)=(AxonLigand(n)-SubstrateLigand(xt,yt))*...
%                        exp(-tau*AdapCoeff_lig(ii,n)^2)+...
%                        SubstrateLigand(xt,yt);   
%                    
%         ReceptorHistory(ii,n) = AxonReceptor(n);
%         LigandHistory(ii,n) = AxonLigand(n);
         
        else
            ReceptorHistory(ii,n) = AxonReceptor(n);
            LigandHistory(ii,n)   = AxonLigand(n);
        
        end
%%

        wDx = wkeitpd(Dx,sigma);
        wDx_target = wkeitpd(Dx_target,sigma);

        if ftw == 1
            xt=xt+1;   
            yt=yt;     
            
        elseif ftw == 0
        
            if rand>=wkeit01(wDx,wDx_target)
                xt = xt+xtDirection;
                yt = yt+ytDirection;
            end  
        end

  
        xtHistory(ii,n) = xt;
        ytHistory(ii,n) = yt;
        
    end

    clc;
    steps-ii
    
end


%% PLOTS

%--------------------------------------------------------------------------
%Map plot


if knockIn > 0      % Isl2-EphA3 constructs
    for i=1:round(NoGrowthCone/2)             
        xEndeWildTyp(i,1) = xtHistory(steps, 2*i-1)-offset;
        yStartWildTyp(i,1) = YStartPos(1, 2*i-1)-offset;
        
    end

    for i=1:floor(NoGrowthCone/2)
        xEndeKnockIn(i,1) = xtHistory(steps, 2*i)-offset;
        yStartKnockIn(i,1) = YStartPos(1, 2*i)-offset;
        
    end

    hold on
    scatter(xEndeWildTyp, yStartWildTyp, 'b', '*');
    scatter(xEndeKnockIn, yStartKnockIn, 'r', '*');
    axis([0 50 0 50]);
    hold off

elseif knockIn == 0
    
    xend=xtHistory(steps,:)-offset;
    xstart=xtHistory(1,:)-offset;

    yend=ytHistory(steps,:)-offset;
    ystart=round(YStartPos)-offset;
    
    scatter(xend, ystart, 'b', '*');
    axis([0 50 0 50]);
    
end

%--------------------------------------------------------------------------
% Gap Plot
% blue=[0 0.7 1];
% red =[1 0.3 0.4];

% figure
% hold on
% rectangle('Position',[0,0,unterkante,FieldSizeX],'FaceColor',...
%[1 0.3 0.4]);
% rectangle('Position',[oberkante,0,FieldSizeX-oberkante,FieldSizeX],...
%'FaceColor',[1 0.3 0.4]);
% scatter(xend,ystart,300,'MarkerEdgeColor','k','MarkerFaceColor',...
%[0.4 0.4 0.4],'LineWidth',1);
% alpha(0.5);
% axis([0 FieldSizeX 0 FieldSizeX]);
% hold off

