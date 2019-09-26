function [ A , S , runtime, itera , returnInfo ] = CNMF( HS , MS , opt )
returnInfo.status = 'null' ;
if( isfield( opt , 'A_init'       ) ) ; A            = opt.A_init       ; else error( 'ERROR @ HIS_CNMF(): A_init is not provided !' )              ; end ;
if( isfield( opt , 'S_init'       ) ) ; S            = opt.S_init       ; else error( 'ERROR @ HIS_CNMF(): S_init is not provided !' )              ; end ;
if( isfield( opt , 'F'            ) ) ; F            = opt.F            ; else error( 'ERROR @ HIS_CNMF(): Spectral Response F is not provided !' ) ; end ;
if( isfield( opt , 'G'            ) ) ; G            = opt.G            ; else error( 'ERROR @ HIS_CNMF(): Spatial Response G is not provided !' )  ; end ;
if( isfield( opt , 'maxIteraNum'  ) ) ; maxIteraNum  = opt.maxIteraNum  ; else maxIteraNum  = 1                                                     ; end ;
if( isfield( opt , 'AStep'        ) ) ; AStep        = opt.AStep        ; else AStep        = 130                                                   ; end ;
if( isfield( opt , 'SStep'        ) ) ; SStep        = opt.SStep        ; else SStep        = 130                                                   ; end ;
if( isfield( opt , 'convthresh'   ) ) ; convThresh   = opt.convthresh   ; else convThresh   = 1e-2                                                  ; end ;
if( isfield( opt , 'convTh_hyper' ) ) ; convTh_hyper = opt.convTh_hyper ; else convTh_hyper = 1e-8                                                  ; end ;
if( isfield( opt , 'convTh_multi' ) ) ; convTh_multi = opt.convTh_multi ; else convTh_multi = 1e-8                                                  ; end ;
HS_bandNum = size(HS,3)                     ;
MS_bandNum = size(MS,3)                     ;
HS(HS<0)   = 0                              ;
MS(MS<0)   = 0                              ;
Y_H        = reshape(HS,[],HS_bandNum)'     ;
Y_M        = reshape(MS,[],MS_bandNum)'     ;
obj_curr   = 0.5 * sum(sum((Y_H-A*S*G).^2)) + ...
             0.5 * sum(sum((Y_M-F*A*S).^2)) ;
obj_it     = zeros(maxIteraNum+1,1)         ;
time_it    = zeros(maxIteraNum+1,1)         ;
obj_it(1)  = obj_curr                       ;
% Iteration
% ============== %
% Initialization %
% ============== %
% --------------------------------------------------------------------------- %
% remark : from original implementation by Naoto YOKOYA we see that the order %
% of solving the CNMF subproblems in initialization stage is                  %
% 0th : vca to detect A , set SG as onesMatrix / modelOrder                   %
% 1st : solve hyper-problem to obtain A & SG                                  %
% 2nd : interpolate SG to obtain S                                            %
% 3rd : obtain FA by F*A                                                      %
% 4th : solve multi-problem to obtain FA & S                                  %
% we can see that the following approach is more meaningful                   %
% 0th : vca do detect A , set FA as F*A                                       %
% 1st : solve multi-problem to obtain FA & S                                  %
% 2nd : obtain SG by S*G                                                      %
% 3rd : solve hyper-problem to obtain A & SG                                  %
% here we don't need to interpolate from SG to S                              %
% --------------------------------------------------------------------------- %
FA = F * A ;
for i = 1:SStep
    if( i == 1 )
        AtFtFA = FA'*FA  ;
        S_n    = FA'*Y_M ;
        for q = 1:SStep
            S_d = AtFtFA *S   ;
            S   = max( S.*S_n./S_d , eps );
        end
    end
    % Update FA & S
    FA_n = Y_M *S'         ;
    FA_d = FA  *(S   *S')  ;
    FA   = max( FA .*FA_n./FA_d , eps ) ; % LeeSeung MU on FA update
    S_n  = FA' *Y_M        ;
    S_d  = (FA'*FA)  *S    ;
    S    = max( S  .*S_n ./S_d  , eps ) ; % LeeSeung MU on S  update
end
SG = S * G    ;
for i = 1 : AStep
    if( i == 1 )
        AtA  = A'*A   ;
        SG_n = A'*Y_H ;
        for q = 1 : AStep
            SG_d = AtA*SG         ;
            SG   = max( SG.*SG_n./SG_d , eps ) ; % LeeSeung MU on SG update
        end
    end
    % Update A & SG
    A_n  = Y_H*SG'                     ;
    A_d  = A  *(SG  *SG')              ;
    A    = max( A .*A_n ./A_d  , eps ) ; % LeeSeung MU on A  update
    SG_n = A' *Y_H                     ;
    SG_d = (A'*A)   *SG                ;
    SG   = max( SG.*SG_n./SG_d , eps ) ; % LeeSeung MU on SG update
end
% ===================== %
% End of Initialization %
% ===================== %
% ========= %
% Iteration %
% ========= %
tic;
for itera = 1:maxIteraNum
    FA = F * A ;
    for i = 1:SStep
        if( i == 1 )
            AtFtFA = FA'*FA  ;
            S_d    = FA'*Y_M ;
            for q = 1:SStep
                S_n = AtFtFA *S                ;
                S   = max( S.*S_d./S_n , eps ) ;
            end
            multiObjCurr = 0.5 * sum(sum( (Y_M-FA*S).^2 )) ;
        end
        % Update FA & S
        multiObjPrev = multiObjCurr ;
        FA_d = Y_M *S'                      ;
        FA_n = FA  *(S*S')                  ;
        FA   = max( FA .*FA_d./FA_n , eps ) ;
        S_d  = FA' *Y_M                     ;
        S_n  = (FA'*FA)  *S                 ;
        S    = max( S  .*S_d ./S_n  , eps ) ;
        multiObjCurr = 0.5 * sum(sum( (Y_M-FA*S).^2 )) ;
        if( abs(multiObjCurr-multiObjPrev)/multiObjCurr < convTh_multi )
            break ;
        end
    end
    SG = S * G ;
    for i = 1:AStep
        if( i == 1 )
            SGGtSt = SG *SG' ;
            A_d    = Y_H*SG' ;
            for q = 1:AStep
                A_n = A *SGGtSt   ;
                A   = max( A.*A_d./A_n , eps ) ;
            end
            hyperObjCurr = 0.5 * sum(sum( (Y_H-A*SG).^2 )) ;
        end
        % Update A & SG
        hyperObjPrev = hyperObjCurr ;
        SG_d = A' * Y_H       ;
        SG_n = (A'*A)   *SG   ;
        SG   = max( SG.*SG_d./SG_n , eps ) ;
        A_d  = Y_H*SG'        ;
        A_n  = A  *(SG  *SG') ;
        A    = max( A .*A_d ./A_n  , eps ) ;
        hyperObjCurr = 0.5 * sum(sum( (Y_H-A*SG).^2 )) ;
        if( abs(hyperObjCurr-hyperObjPrev)/hyperObjCurr < convTh_hyper )
            break ;
        end
    end
    obj_prev = obj_curr ;
    obj_curr = 0.5 * sum(sum((Y_H-A*S*G).^2)) + 0.5 * sum(sum((Y_M-F*A*S).^2)) ;
    obj_chng = abs(obj_curr-obj_prev)/obj_curr ;
    time_it(itera+1) = toc ;
    obj_it(itera+1)  = obj_curr ;
    if( obj_chng < convThresh )
        break ;
    end
end
% ================ %
% End of Iteration %
% ================ %
runtime = toc;
returnInfo.obj_it  = obj_it  ;
returnInfo.time_it = time_it ;
end
