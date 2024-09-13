function [Bound, Weyl] = LANbound(A,singval,tV,tU,varargin)
% LANbound computes a bound on the accuracy in extracting singular values, 
% given approximations tV (and tU) to the leading singular subspaces, for
% one of the methods: Generalized Nystrom, HMT, Rayleigh-Ritz, or SVD. 
% The bound is computed according to "Matrix Perturbation Analysis of 
% Methods for Extracting Singular Values from Approximate Singular 
% Subspaces" by L. Lazzarino, H. Al Daas, Y. Nakatsukasa.
%
%   Bound = LANbound(A,singval,tV,tU) computes a (forward) bound on the 
%   singular value extraction error of the Generalized Nystrom 
%   approximation on the matrix A, with exact singular values singval, 
%   and approximation to leading subspaces tV and tU.
%
%   Bound = LANbound(A,singval,tV) computes a (forward) bound on the 
%   singular value extraction error of the HMT method on the matrix A, 
%   with exact singular values singval, and approximation to leading 
%   subspaces tV and tU.
%
%   Bound = LANbound(A,singval,tV,tU,varargin) specifies integration
%   option values in the fields of a structure. It is possible to specify:
%   - 'Method', specifies the method for which the bound is computed.
%               Possibilities: 'GN' for Generalized Nystrom, 'HMT' for HMT 
%               method, 'RR' for Rayleigh-Ritz method, 'SVD' for SVD method. 
%               For the methods where only tV is needed, insert tU=[];
%   - 'Pseudoinverse', specifies how the pseudoinverse is computed.
%                       Possibilities: 'QR' (default) computes it by QR 
%                       factorization of the matrix to inverte,
%                       'lsqminnorm' uses the Matlab function lsqminnorm,
%                       'Backslash' uses the Matlab command \;
%   - 'ImproveBound', specifies if the heuristic improvement strategy has
%                     to be used. Possibilities: 'N' (default), 'Y'. Note,
%                     this strategy can be used only for the oversampled
%                     cases;
%   - 'Orientation', specifies which kind of bound has to be computed.
%                   Possibilities: 'Forward' (default), 'Backward',
%                   'ApproxBackward' (valid only with Method = GN) where 
%                   some quantities from the Backward bound are approximated
%                   to improve computability. In case of either Backward
%                   bounds, insert as singval the computed singular values.
%
%   [Bound, Weyl] = LANbound(A,singval,tV,...) gives also the bound from
%   Weyl's Theorem.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------- DEALING WITH OPTIONS

    % Define valid option names and their corresponding valid values
    validOptions = {'Method', 'Pseudoinverse', 'ImproveBound','Orientation'};
    validValues.Method = {'GN', 'HMT', 'SVD','RR'}; % Valid values for Method
    validValues.Pseudoinverse = {'QR', 'Backslash', 'lsqminnorm'}; % Valid values for Pseudoinverse
    validValues.ImproveBound = {'Y', 'N'}; % Valid values for ImproveBound
    validValues.Orientation = {'Forward', 'Backward','ApproxBackward'}; % Valid values for Orientation
    % Parse input options
    options = struct('Method', 'GN', 'Pseudoinverse', 'QR', 'ImproveBound', 'N','Orientation','Forward'); % Set default values
    
    nArgs = length(varargin);       
    if mod(nArgs,2) ~= 0     % Check options are given as name/value pairs
        error('Optional parameters should be provided as name/value pairs');
    end

    % Checking validity of options
    for pair = reshape(varargin,2,[])
        inputName = pair{1};
        if any(strcmp(inputName, validOptions))
            validVals = validValues.(inputName);
            if any(strcmp(pair{2}, validVals))
                options.(inputName) = pair{2};  % Storing if valid
            else
                error('%s is not a valid value for %s', pair{2}, inputName);
            end
        else
            error('%s is not a valid parameter name', inputName);
        end
    end

%---------------------------- 

    % Set default method when tU is not given
    if (nargin < 4)
        options.Method = 'HMT';    
        if (nargin < 3 || isempty(tV))   % Check enough input are given
            error('Not enough input')
        end
    end

    % Check if tU is given when needed from the method (for GN or RR)
    if (nargin < 4 || isempty(tU)) && any(strcmp(options.Method, ['GN','RR'])) 
        error('You need approximation to both left and right subspaces for Method: %s', options.Method) 
    end
    
    % Consider HMT as a GN with tU = A*tV 
    if strcmp(options.Method, 'HMT')
        tU = A*tV;
        options.Method = 'GN';
    end

    % Checking sizes
    [m,n] = size(A);
    
    [nV,r] = size(tV);  % obtain parameter r
        if nV ~= n 
            error('Wrong size of approximation')
        end 

    if ~strcmp(options.Method, 'SVD')
            
        % When meaningful, obtain parameter r+l=:rl (we can have rl=r)
        [mU,rl] = size(tU);     
            if mU ~= m 
                error('Wrong size of approximation')
            end
    end

    % Check suitable situation to use the ApproxBackward when requested
    if strcmp(options.Orientation, 'ApproxBackward') && ~strcmp(options.Method, 'GN')
        error('Approximated Backward Bound available only for Generalized Nystrom approximation')
    end

    % Check suitable situation to use the Improving strategy when requested
    if strcmp(options.ImproveBound, 'Y') && (strcmp(options.Method, 'SVD') || rl == r)
        error('The improving strategy can be used only for oversampled cases');
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BODY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------- Apply orthogonal Transformations

    % Compute \bar{A} = ([tU tU_\perp]*)A*[tV tV_\perp]
    [Q2,~] = qr(tV);
    if isempty(tU)
        bA = A*Q2;
    else
        [Q1,~] = qr(tU);
        bA =Q1'*A*Q2;
    end
    
    
%---------------------------- Define Blocks
    
    if strcmp(options.Method, 'SVD')
        tar_nsv = r;        % target number of sing.val. to extract
        bA2 = bA(:,r+1:end);
    else

        tar_nsv = min(r,rl); % target number of sing.val. to extract
        bA11 = bA(1:rl,1:r);
        
        % To improve the bound (when requested/possible), other
        % tranformation are needed
        if strcmp(options.ImproveBound, 'Y') 
            
            [U11,~,V11] = svd(bA11);    % SVD of the A11 block
    
            Up = blkdiag(U11,eye(m-rl));    % Transformations to 
            Vp = blkdiag(V11,eye(n-r));     % diagonalize A11 block
            
            bA = Up'*bA*Vp; % Redefine teh matrix

            rl = r;     % Set to have A11 block square
            bA11 = bA(1:rl,1:r);    % Redefine A11 block

        end
        
        bA12 = bA(1:rl, r+1:end);
        bA21 = bA(rl+1:end, 1:r);
        bA22 = bA(rl+1:end,r+1:end);
        svbA22 = svd(bA22);     % Singular Values of the  bA22 block
    end
    
    
%---------------------------- Compute Perturbation E
    
    % Compute perturbation for Generalized Nystrom
    if strcmp(options.Method, 'GN')
        
        
        if (rl> r) % Oversampling on the rows (as in the paper) 
            
            if strcmp(options.Pseudoinverse, 'QR')
                [Q,R] = qr(bA11);   % QR factorization of the matrix to invert
                P = Q'*bA12;
                B1 = (bA11/R) * P;
                B2 = (bA21/R) * P;
            elseif strcmp(options.Pseudoinverse,'lsqminnorm')
                B1 = bA11*(lsqminnorm(bA11,bA12));
                B2 = bA21*(lsqminnorm(bA11,bA12));
            else
                B1 = (bA11/bA11)*bA12; 
                B2 = (bA21/bA11)*bA12;
            end
            
            offdE = bA12 - B1;      % Off-diagonal block of Perturbation E
            brE = bA22 - B2;        % Bottom-rigth block of Perturbation E

            % Form Perturbation E
            E = [zeros(rl,r) offdE; zeros(m-rl,r) brE];
        
        elseif (rl == r)    % No Oversampling
            
            % There is only a Bottom-rigth Perturbation 
            if strcmp(options.Pseudoinverse, 'QR')
                [Q,R] = qr(bA11);   % QR factorization of the matrix to invert
                P = Q'*bA12;
                B2 = (bA21/R) * P;
            elseif strcmp(options.Pseudoinverse,'lsqminnorm')
                B2 = bA21*(lsqminnorm(bA11,bA12));
            else
                B2 = (bA21/bA11)*bA12;
            end
            
            brE = bA22 - B2;  
            offdE = 0;

            % Form Perturbation E
            E = [zeros(rl,n); zeros(m-rl,r) brE];

        else        % Oversampling on the columns
           
            if strcmp(options.Pseudoinverse, 'QR')
                [Q,R] = qr(bA11);   % QR factorization of the matrix to invert
                P = Q'*bA11;
                B1 = (bA21/R) * P;
                B2 = (bA21/R) * P;
            elseif strcmp(options.Pseudoinverse,'lsqminnorm')
                B1 = bA21*(lsqminnorm(bA11,bA11));
                B2 = bA21*(lsqminnorm(bA11,bA12));
            else
                B1 = (A21/bA11)*bA11; 
                B2 = (bA21/bA11)*bA12;
            end
            
            offdE = bA21 - B1;  % Off-diagonal block of Perturbation E
            brE = bA22 - B2;    % Bottom-rigth block of Perturbation E
            
            % Form Perturbation E
            E = [zeros(rl,n); offdE brE];
        end
            
        norm_offdE = norm(offdE);
        norm_brE = norm(brE);
        
    % Compute perturbation for Rayleigh-Ritz
    elseif strcmp(options.Method, 'RR')
    
        norm_offdE = max(norm(bA12),norm(bA21));
        norm_brE = norm(bA22);

        % Form Perturbation E
        E = [zeros(rl,r) bA12; bA21 bA22];
    
    % Compute perturbation for SVD
    else  
    
        norm_offdE = norm(bA2);
        norm_brE = 0;

        % Form Perturbation E
        E = [zeros(m,r) bA2];
    
    end

    % Compute Weyl's bound
    Weyl = norm(E);
    
%---------------------------- Compute Original Matrix-dependent terms

% delta := maximum norm between off-diagonal blocks of transformed matrix
% gap := quantity dependent on exact/computed singular values

    % Compute delta and gap for Forward bound
    if strcmp(options.Orientation, 'Forward')
            
        if strcmp(options.Method, 'SVD')
            gap = singval(1:tar_nsv);
            delta = norm(bA2);
        else
            for i=1:tar_nsv
                gap(i) = min(abs(singval(i)-svbA22));
            end
            gap = gap';
    
            delta = max(norm(bA12),norm(bA21));
        end
    
    % Compute delta and gap for Backward bound
    elseif strcmp(options.Orientation, 'Backward')                    
    
        if strcmp(options.Method, 'GN')
            svB2 = svd(B2);
            for i=1:tar_nsv
                gap(i) = min(abs(singval(i)-svB2));
            end
            gap = gap';
    
            if r == rl
                delta = max(norm(bA12),norm(bA21));
            else
                delta = max(norm(B1),norm(bA21));
            end

        elseif strcmp(options.Method, 'RR')
            gap = singval(1:tar_nsv);
            delta = 0;

        else  % i.e., if strcmp(options.Method, 'SVD')
            gap = singval(1:tar_nsv);
            delta = 0;
        end

    % Compute delta and gap for ApproxBackward bound
    else
        normB2 = norm(B2); 
        gap = singval(1:tar_nsv) - normB2;

        delta = max(norm(bA12),norm(bA21));
    end
    
%---------------------------- Compute Bound   
    
    % Form Numerator of tau 
    if strcmp(options.Orientation, 'ApproxBackward') && (rl ~= r)

       Num = delta + norm(bA12);

    else    

        Num = delta + norm_offdE;  % Note: norm_offdE = 0 if rl = r

    end

    for i=1:tar_nsv
        Den = gap(i) - 2*Weyl; % Form Denominator of tau
    
        % Compute bound if tau > 0, i.e, it has strictly positive denominator
        if Den > 0  
       
            tau(i) = Num/Den;
            
            % For GN and r==rl norm_offdE is zero and so it is term1
            term1(i) = 2*norm_offdE*tau(i); 
            
            term2(i) = norm_brE*(tau(i))^2;
            Bound(i) = term1(i)+term2(i);
       
        else    % Assign NaN when the bound does not exist
            
            Bound(i) = NaN;
        
        end

    end
end
    
