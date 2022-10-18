classdef ( ConstructOnLoad = false ) SSD
    
    % An Abstract Class to Define a Supersaturated Design
    
    properties ( Constant = true, Access = protected )
        ListOfProperties = ["Measure", "Es2", "MaxAbsCorr", "MeanAbsCorr", "C"];   % List of design evaluation measures
    end
    
    properties ( SetAccess = immutable )
        M               int16    {mustBeGreaterThanOrEqual(M, 2),...
                                  mustBeInteger(M), mustBeReal(M)}          % Number of factors
        CreationDate = date                                                 % Date obect was created
        UserName = getenv("username")                                       % Username for creator of class
        ComputerName = getenv("computername")                               % Computer name      
    end
    
    properties ( SetAccess = protected )
        Lo              double   {mustBeNumeric(Lo), mustBeFinite(Lo),...
                                  mustBeReal(Lo)}                           % Factor low levels
        Hi              double   {mustBeNumeric(Hi), mustBeFinite(Hi),...
                                  mustBeReal(Hi)}                           % Factor high levels
        Name            string                                              % Descriptive Name for Factors
        Abbreviation    string                                              % Factor abbreviation
        N               int16    {mustBeNumeric(N), mustBeFinite(N),...     % Design size
                                  mustBeReal(N), mustBeInteger(N)}          
    end
    
    properties ( Access = protected )
        S               double                                              % Design in scaled units
    end
    
    properties ( SetAccess = protected, Dependent )
        D               double                                              % Design in natural units
        X               double                                              % Regression matrix
        Es2             double                                              % E(s^2) criterion
        MaxAbsCorr      double                                              % Maximum absolute correlation
        MeanAbsCorr     double                                              % Mean absolute correlation
        No              int16                                               % Number of non-orthogonal column pairings  
        NumOrthCol      int16                                               % Number of orthogonal columns pairings
        C               double                                              % Jones' et al singular value criterion
    end
    
    properties ( Access = protected, Dependent )
        N_              double                                              % Number of runs
        M_              double                                              % Number of factors
    end
    
    properties (SetAccess = protected, Dependent, Abstract = true )
        Measure         double                                              % Design measure to optimise
    end
    
    methods ( Abstract = true )
        obj = designGenerator( obj, varargin );                             % Generate the design
        Ok = isBetter( obj, OldCost )                                       % Abstraction of interpretation of cost function
    end % end abstract methods
    
    methods
        function obj = SSD( M, varargin)
            %-----------------------------------------------------------------
            % Class constructor. Generate supersaturated designs.
            %
            % obj = SSD( M, 'PARAM1', VALUE1, ..., 'PARAM#', VALUE#);
            %
            % Where M is the number of factors, 'PARAM#' is a STRING defining
            % the nature of the following VALUE#. Valid forms for 'PARAM#'
            % are:
            %
            % Lo           --> Level corresponding to minus one setting
            %                  for each factor. Associated value is an
            %                  (Mx1) DOUBLE.
            % Hi           --> Level corresponding to plus one setting
            %                  for each factor. Associated value is an
            %                  (Mx1) DOUBLE.
            % Factor       --> Descriptive name of factor. Corresponding
            %                  VALUE is a (1xM) STRING array.
            % Abbreviation --> Abbreviation for each factor. Corresponding
            %                  VALUE is a (1xM) STRING array.
            %
            % Note if a 'PARAM#' is not assigned, then suitable defaults
            % are provided. For example, the defaults for "Lo" and "Hi" are
            % -1 and +1 respectively. Defaults for Factor and Abbreviation
            % are the letters of the alphabet.
            %
            % Examples:
            %
            % obj = SSD( 10, "HI", ones(10 , 1), "LO", -ones(10, 1));
            %-----------------------------------------------------------------
            obj.M = M;                                                      % Number of factors.
            %-----------------------------------------------------------------
            % Assign supplied factor information & set defaults if required
            %-----------------------------------------------------------------
            obj = obj.assignDefaults();                                     % Assign default values
            if ~isempty( varargin )
                obj = obj.assignFactorData( varargin{:} );                  % Overwrite with userdefined values if required
            end
        end
        
        function obj = importScaledDesign( obj, DesignScaled )
            %--------------------------------------------------------------
            % Import a scaled design. Useful for developing comparisons
            % between the current design and a reference, say from a paper.
            %
            % obj = obj.importScaledDesign( DesignScaled );
            %
            % Input Arguments:
            %
            % DesignScaled  --> (NxM) scaled design
            %--------------------------------------------------------------
            
            [Ok, ErrorMsg] = obj.checkDesign( DesignScaled );
            if Ok
                %----------------------------------------------------------
                % Import the design
                %----------------------------------------------------------
                R = size( DesignScaled, 1 );
                obj.S = DesignScaled;
                obj.N = R;
            else
                %----------------------------------------------------------
                % Display eror message
                %----------------------------------------------------------
                error(ErrorMsg);
            end
        end
    end % constructor and ordinary methods
    
    methods ( Hidden = true )
        function DM = getListOfProperties( obj )
            %--------------------------------------------------------------
            % Return the list of available design properties
            %
            % DM = obj.getListOfProperties();
            %--------------------------------------------------------------
            DM = obj.ListOfProperties();
        end
    end
    
    methods
        % SET/GET Methods
        function X = get.X( obj )
            % Return the regression matrix
            if ~isempty( obj.S )
                X = horzcat( ones(obj.N_, 1), obj.S );
            else
                X = [];
            end
        end
        
        function Rho = get.MeanAbsCorr( obj )
            % Return average absolute correlation
            Rho = abs( corrcoef( obj.S ) );
            A = true( size( Rho ) );
            A = tril(A, -1);
            Idx = find( A == true, numel( A ), 'first');
            Rho = reshape( Rho, numel( Rho ), 1);
            Rho = mean( Rho(Idx) );
        end
        
        function Rho = get.MaxAbsCorr( obj )
            % Return maximum absolute correlation
            Rho = corrcoef( obj.S );
            Rho = abs( tril( Rho, -1 ) );
            Rho = max( max( Rho ) );
        end
        
        function N = get.NumOrthCol( obj )
            % Return total number of orthogonal column pairs
            N = obj.M_;
            N = 0.5*N*( N + 1 );
            N = int16( N );
        end
        
        function N = get.No( obj )
            % Return number of nonorthogonal columns
            N = corrcoef( obj.S );
            N = abs( tril( N, -1 ) );
            N = reshape( N>0, numel( N ), 1 );
            N = sum( N );
            N = int16( N );
        end
        
        function C = get.C( obj )
            % Return singular value criterion
            [~, Sv, ~] = svd(obj.S,'econ');
            Sv = diag( Sv );
            Idx = 2/( obj.N_ - 1);
            C = (obj.N_ - 1)*prod( Sv(1:(obj.N_-1)).^(Idx) )/ obj.M_/ obj.N_;
        end
        
        function E = get.Es2( obj )
            % Return Es2 Criterion
            E = 0;
            for Q = 1:obj.M_
                for R = Q+1:obj.M_
                    E = E + (obj.S(:, Q).'*obj.S(:, R))^2;
                end
            end
            E = 2*E/(obj.M_ - 1)/obj.M_;
        end
            
        function D = get.D( obj )
            % Return the design matrix in natural or engineering units
            D = decode( obj );
        end
        
        function N = get.N_( obj )
            % Make copy of N for calcs
            N = double( obj.N );
        end
        
        function M = get.M_( obj )
            % Make copy of M for calcs
            M = double( obj.M );
        end
        
        function obj = set.N( obj, N )
            %----------------------------------------------------------------------
            % Ensure the design size is less than the number of factors
            %----------------------------------------------------------------------
            Ok = floor( N ) <= (obj.M_ - 1);                                          %#ok<MCSUP>
            if Ok
                obj.N = int16( floor( N ) );
            else
                error( '\nMaximum Design Size is %4.0f\n', obj.M_ - 1 );     %#ok<MCSUP>
            end
        end
    end % get/set methods
    
    methods ( Access = protected )
        function P = pointer2Columns( obj, J )
            %--------------------------------------------------------------
            % Return a logical pointer to included factors (columns of the
            % design matrix)
            %
            % P = obj.pointer2Columns( J );
            %
            % Input Arguments:
            %
            % J     --> Column to be deleted
            %--------------------------------------------------------------
            if nargin<2
                J = NaN;
            end
            P = ~ismember( 1:obj.M_, J );
        end
        
        function X = computeX( obj, J )
            %--------------------------------------------------------------
            % Generate the corresponding regression matrix.
            %
            % X = obj.computeX( J );
            %
            % Input Arguments:
            %
            % J     --> Deleted factor
            %--------------------------------------------------------------
            P = obj.pointer2Columns( J );
            X = obj.S( :, P );
            X = horzcat( ones( obj.N_, 1 ), X );
        end
    end % End protected methods
    
    methods ( Access = private )
        function [Ok, Msg] = checkDesign( obj, D )
            %--------------------------------------------------------------
            % Check the imported design is a SSD.
            %
            % [Ok, Msg] = obj.checkDesign( D );
            % 
            % Input Arguments:
            %
            % D     --> New design in scaled units
            % 
            % Output Arguments:
            %
            % Ok    --> Will be true if no errors are found
            % Msg   --> Frmatted error message
            %--------------------------------------------------------------
            Ok = true;
            Msg = '';
            [Rsize, Csize] = size( D );
            Plus1 = sum( ( D == 1) );
            Minus1 = sum( ( D == -1) );
            P = ( D == 1) | ( D == -1);
            if (obj.M_ ~= Csize)
                % Wrong number of factors supplied!
                Ok = false;
                Msg = sprintf('Number of factors supplied must be %4.0f. Imported design has %4.0f factors', obj.M_, Csize);
            elseif ( Rsize > (obj.M_ - 1) )
                % Too many runs supplied
                Ok = false;
                Msg = sprintf('Maximum allowable number of runs for this SSD is %4.0f', (obj.M_ - 1));
            elseif ~all( all( P ) )
                % Design must consist of -1 and +1 only
                Ok = false;
                Msg = sprintf('Scaled design must contain only "-1" or "+1"');
            elseif any( Plus1 ~= Minus1 )
                % Not a balanced design
                Ok = false;
                Msg = sprintf('Design must be balanced');
            end
        end
        
        function [LOW, MID, HIGH] = codingHelper( obj )
            %--------------------------------------------------------------
            % Return metrics to support coding and decoding of the design
            %
            % [LOW, MID, HIGH] = obj.codingHelper();
            %
            % Output Arguments:
            %
            % LOW   --> (1xM) row vector of min factor values
            % MID   --> (1XM) row vector of median factor values
            % HIGH  --> (1xM) row vector of max factor values
            %--------------------------------------------------------------
            HIGH = reshape( obj.Hi, 1, obj.M_);
            LOW = reshape( obj.Lo, 1, obj.M_);
            MID = median([ LOW; HIGH ]);
        end
        
        function X = decode( obj )
            %--------------------------------------------------------------
            % Decode design from the interval [-1, 1] to natural or
            % engineering units
            %
            % X = obj.decode();
            %
            % Output Arguments:
            %
            % X     --> Candidate list. The maximum allowable size for the
            %           candidate list is 2^17 rows. If 2^M is greater than
            %           this, then
            %--------------------------------------------------------------
            [LOW, MID, HIGH] = obj.codingHelper();
            X = 0.5*( HIGH - LOW ).*obj.S + MID;
        end
        
        function L = assignDefaultNames( obj )
            %--------------------------------------------------------------
            % Assign default names or abbreviations
            %
            % L = obj.assignDefaultNames();
            %
            % Output Arguments:
            %
            % L     --> Mx1 string array of names or abbreviations
            %--------------------------------------------------------------
            number2letter = @(X)( char(X-1+'A') );
            L = strings(1, obj.M_);
            for Q = 1:obj.M_
                K = floor( Q/26 );
                if K>0
                    LL = number2letter( K );
                else
                    LL = '';
                end
                Letter = rem( Q, 26 );
                if Letter == 0
                    Letter = 26;
                end
                L(Q)= string( horzcat( LL, number2letter( Letter ) ) );
            end
        end
        
        function obj = assignDefaults( obj )
            %--------------------------------------------------------------
            % Assign default values to properties
            %
            % obj = obj.assignDefaults();
            %--------------------------------------------------------------
            obj.Lo = -ones(obj.M_, 1);
            obj.Hi = ones(obj.M_, 1);
            obj.Name = obj.assignDefaultNames();
            obj.Abbreviation = obj.assignDefaultNames();
        end
        
        function obj = assignFactorData( obj, V )
            %-----------------------------------------------------------------
            % Assign initial factor data
            %
            % obj = obj.assignFactorData( V );
            %
            % Input Arguments:
            %
            % V is a cell array of form {'param1', value1, ... 'param#',
            % value#}
            %
            % Where 'param#' is a char defining the nature of the following
            % value. Valid  'param#' arguments are:
            %
            % Lo           --> Level corresponding to minus one setting (Mx1)
            % Hi           --> Level corresponding to plus one setting (Mx1)
            % Name         --> Full name of factor (Mx1) string
            % Abbreviation --> Abbreviation for factor (Mx1) string
            %
            % Note if a 'param#' is not assigned, then suitable defaults
            % are provided. For example, the defaults for 'Lo' and 'Hi' are
            % -1 and +1 respectively.
            %-----------------------------------------------------------------
            
            %-----------------------------------------------------------------
            % Overwrite defaults with user defined values if supplied
            %-----------------------------------------------------------------
            Props = V(1:2:end);
            Values = V(2:2:end);
            for Q = 1:numel(Props)
                switch lower( Props{ Q } )
                    case "lo"
                        LOW = Values{ Q };
                        if numel(LOW) == obj.M_
                            obj.Lo = LOW(:);
                        end
                    case "hi"
                        HIGH = Values{ Q };
                        if numel(HIGH) == obj.M_
                            obj.Hi = HIGH(:);
                        end
                    case "name"
                        NAMES = Values{ Q };
                        if numel(NAMES)==obj.M_ && isstring(NAMES)
                            obj.Name = NAMES;
                        end
                    case "abbreviation"
                        ABBREV = Values{ Q };
                        if numel(ABBREV)==obj.M_ && isstring(ABBREV)
                            obj.Abbreviation = ABBREV;
                        end
                end
            end
        end
    end % private/helper
    
    methods ( Static = true, Hidden = true )      
        function M = maxNumFactors( N )
            %--------------------------------------------------------------
            % Compute maximum number of allowable factors for a given
            % design size.
            %
            % M = SSD.maxNumFactors( N );
            % M = obj.maxNumFactors( N );
            %
            % Input Arguments:
            %
            % N     --> Number of runs (must be even)
            %--------------------------------------------------------------
            if mod( N, 2) == 1
                error('\nArgument "N" must be an even number\n');
            end
            % Compute maximum number of factors
            M = factorial( N - 1)/factorial( N/2 - 1 )/factorial( N/2 );
        end
    end
end

