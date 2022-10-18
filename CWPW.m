classdef CWPW < SSD
    % Li & Wu's Columnwise-Pairwise optimal design method for SSD
    properties
        NumInitDesigns  int16   { mustBeNumeric(NumInitDesigns), mustBeFinite(NumInitDesigns),...
            mustBeReal(NumInitDesigns), mustBeGreaterThan(NumInitDesigns, 0)} = int16( 200 )
    end
    
    properties ( SetAccess = protected )
        Kexch   int16  {mustBeGreaterThanOrEqual(Kexch, 1), mustBeFinite(Kexch),...
                        mustBeInteger(Kexch), mustBeReal(Kexch)} = int16( 5 )
    end
    
    properties ( Access = protected, Dependent )
        Kexch_                                                              % Number of columns to exchange
        NumInitDesigns_                                                     % Number of initial designs to choose from
    end
    
    methods ( Access = protected, Abstract = true )
        L = costFcn( obj, J )                                               % Change in the determinant for deleting column J
    end
    
    methods ( Access = protected, Static = true, Abstract = true )
        J = sortDesigns( D );
    end
    
    methods
        function obj = CWPW( M, varargin )
            %-----------------------------------------------------------------
            % Class constructor. Generate supersaturated designs.
            %
            % obj = CWPW( M, 'PARAM1', VALUE1, ..., 'PARAM#', VALUE#);
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
            % obj = CWPW( 10, "HI", ones(10 , 1), "LO", -ones(10, 1));
            %-----------------------------------------------------------------   
            if nargin<1
                % Apply default number of factors
                M = 10;
            end
            obj = obj@SSD( M, varargin{:});
        end
        
        function obj = designGenerator( obj, N, C, MaxIter, DispFlg )
            %--------------------------------------------------------------
            % Design generator method. Determines the column to delete and
            % then produces (N/2)^2 possible C1 replacements for evaluation
            %
            % X = obj.designGenerator( N, C, Iter, DispFlg );
            %
            % Input Arguments:
            %
            % N         --> Design size. Must be a multiple of 2. Default 
            %               is ceil( obj.M/2 )
            % C         --> Number of columns to exchange
            % MaxIter   --> Maximum number of iterations {100}
            % DispFlg   --> Set to true to output intermediate results to
            %               the screen {true}
            %--------------------------------------------------------------
            if nargin<3 || isempty( C )
                C = floor( obj.M_/2 );
            else
                C = floor( C );
            end
            obj.Kexch = C;
            if nargin<4 || isempty( MaxIter ) || ~isnumeric( MaxIter ) || ...
                           ~isreal( MaxIter ) || MaxIter<0 || ~isfinite( MaxIter )
                MaxIter = 100;
            end
            if nargin<5 || DispFlg~=false
                DispFlg = true;
            else
                DispFlg = false;
            end
            %--------------------------------------------------------------
            % Set the design size
            %--------------------------------------------------------------
            obj = obj.setDesignSize( N );
            %--------------------------------------------------------------
            % Generate the initial design
            %--------------------------------------------------------------
            obj = obj.initialDesign();
            %--------------------------------------------------------------
            % Generate the design
            %--------------------------------------------------------------
            obj = obj.optimiseDesign( C, MaxIter, DispFlg );
        end
    end % constructor and ordinary methods
    
    methods
        function K = get.Kexch_( obj )
            % Return number of columns to exchange as a double
            K = double( obj.Kexch );
        end
        
        function obj = set.Kexch( obj, Value )
            % Set then number of columns to be exchanged.
            Value = min( [obj.M_, Value] );
            Value = max( [ 1, Value ] );
            obj.Kexch = int16( Value );
        end
        
        function obj = set.NumInitDesigns( obj, Value )
            % Set the number of initial trials
            V = int16( Value );
            if ( V > 1 )
                obj.NumInitDesigns = V;
            end
        end
        
        function N = get.NumInitDesigns_( obj )
            % Return number of inital designs as a double
            N = double( obj.NumInitDesigns );
        end
    end % end Set/Get methods
    
    methods ( Access = protected )
    end % end protected methods
    
    methods ( Access = private )
        function obj = optimiseDesign( obj, C, MaxIter, DispFlg )
            %--------------------------------------------------------------
            % Optimise the design using the appropriate cost function:
            %
            % 
            % Input Arguments:
            %
            % C         --> 
            % MaxIter   --> Maximum number of iterations 
            % DispFlg   --> Set to true to output intermediate results to
            %               the screen {true}
            %--------------------------------------------------------------
            if nargin<4 || DispFlg~=false
                % Display intermediate results
                DispFlg = true;
            else
                % Don't display intermediate results
                DispFlg = false;
            end
            %--------------------------------------------------------------
            % Optimimise the design
            %--------------------------------------------------------------
            Stopflg = false;
            Counter = 0;
            if DispFlg
                fprintf('\n=====================================================');
                fprintf('\n             CWPW Exchange Algorithm                 ');
                fprintf('\n=====================================================\n');
            end
            J = zeros(1, obj.M_);
            while ~Stopflg
                Counter = Counter + 1;
                %--------------------------------------------------------------
                % Identify the columns contributing least to the determinant
                %--------------------------------------------------------------
                Jlast = J;
                J = obj.column2delete();
                %--------------------------------------------------------------
                % Determine the best pariwise exchange for each column to be
                % exchanged.
                %--------------------------------------------------------------
                for Q = 1:C
                    obj = obj.makeExchange( J( Q ) );
                    if DispFlg
                        fprintf('\nIteration %4.0f, Column Exchanged = %4.0f, New Design Measure = %10.4f',...
                            Counter, J( Q ), obj.Measure);
                    end
                end
                %----------------------------------------------------------
                % Evaluate convergence conditions
                %----------------------------------------------------------
                TerminateLoop = ismember( Jlast(1:C), J(1:C) );
                if (Counter >= MaxIter) || all( TerminateLoop )
                    Stopflg = true;
                    if DispFlg
                        fprintf('\n\n');
                    end
                end
            end
        end
        
        function [obj, NoUpdateFlg] = makeExchange( obj, J )
            %--------------------------------------------------------------
            % Perform the columnwise, pairwise exchange for column J
            %
            % obj = obj.makeExchange( J );
            %
            % Input Arguments:
            %
            % J             --> Column to exchange
            %--------------------------------------------------------------
            P1J = find( (obj.S(:, J) == 1), obj.N_, 'first' );              % Point to +1 in col J
            M1J = find( (obj.S(:, J) ~= 1), obj.N_, 'first' );              % Point to -1 in col J
            DesignScaled = obj.S;
            BestDesign = obj.S;
            OldCost = obj.Measure;                                          % Define the current best design measure
            %--------------------------------------------------------------
            % For each of (obj.N/2)^2 pairs in C1 evaluate the cost
            % function. Exchange the optimal pair in the column identified.
            %--------------------------------------------------------------
            for Q = 1:numel( P1J )
                for R = 1:numel( M1J )
                    %------------------------------------------------------
                    % Generate a new design using a first order adjustment
                    % of column J
                    %------------------------------------------------------
                    NewDesign = DesignScaled;
                    NewDesign( P1J( Q ), J ) = -1;
                    NewDesign( M1J( R ), J ) = 1;
                    obj.S = NewDesign;
                    %------------------------------------------------------
                    % Make the exchange if the measure is improved and the
                    % maximum absolute correlation is less than unity
                    %------------------------------------------------------
                    Exchange = obj.isBetter( OldCost );
                    if Exchange && ( obj.MaxAbsCorr < 1 )
                        OldCost = obj.Measure;
                        NoUpdateFlg = false;
                        BestDesign = NewDesign;
                    end
                end
            end
            %--------------------------------------------------------------
            % Perform the CWPW exchange
            %--------------------------------------------------------------
            obj.S = BestDesign;
        end
        
        function obj = initialDesign( obj )
            %--------------------------------------------------------------
            % Define the initial design. Keep selecting random balanced
            % designs of the correct size until one is feasible.
            %
            % obj = obj.initialDesign();
            %--------------------------------------------------------------
            D = zeros( obj.N_, obj.M_ );
            NumDesigns = obj.NumInitDesigns_;
            OldMeasure = 0;
            for R = 1:NumDesigns
                %----------------------------------------------------------
                % Generate a random design
                %----------------------------------------------------------
                for Q = 1:obj.M_
                    D(:, Q) = randperm( obj.N_ ).';
                end
                P = rem(D, 2) == 0;
                D( P ) = 1;
                D( ~P ) = -1;
                obj.S = D;
                %----------------------------------------------------------
                % Evaluate the cost function. And update the best design if
                % the cost function is improved
                %----------------------------------------------------------
                Exchange = obj.isBetter( OldMeasure );
                if Exchange || (R == 1)
                    OldMeasure = obj.Measure;
                    BestDesign = D;
                end
            end
            %--------------------------------------------------------------
            % Select the best design
            %--------------------------------------------------------------
            obj.S = BestDesign;
        end
        
        function J = column2delete( obj )
            %--------------------------------------------------------------
            % Define column to delete using largest reduction in the 
            % cost function
            %
            % J = obj.column2delete();
            %--------------------------------------------------------------
            D = zeros(1, obj.M_); % Define storage
            for Q = 1:obj.M_
                %----------------------------------------------------------
                % Compute determinant corresponding to deleting each column
                % in turn
                %----------------------------------------------------------
                D( Q ) = obj.costFcn( Q );
            end
            %--------------------------------------------------------------
            % Direction of sort depends on the cost function for the
            % algorithm specified. 
            %--------------------------------------------------------------
            J = obj.sortDesigns( D );
            %--------------------------------------------------------------
            % If maximum absolute correlation is unity then place these
            % columns first
            %--------------------------------------------------------------
            if obj.MaxAbsCorr == 1
                Rho = abs( corrcoef( obj.S ) );
                Rho = tril( Rho, -1 );
                [P, Q] = find( Rho == 1, numel(Rho), 'first' );
                Cols = unique( [P(:); Q(:)] );
                Cols = reshape( Cols, 1, numel( Cols ) );
                P = ~ismember( J, Cols );
                J = [Cols J( P )];
            end
        end           
                
        function obj = setDesignSize( obj, N )
            %--------------------------------------------------------------
            % Set the design size so it is the nearest pure multiple of 2.
            % That is generate balanced designs!
            %
            % obj = setDesignSize( N );
            %
            % Input Arguments:
            %
            % N     --> Desired design size
            %--------------------------------------------------------------
            if nargin<2
                N = ceil( obj.M_/2);
            else
                N = floor( N );
            end
            if rem(N, 2) == 1
                N = N + 1;
            end
            obj.N = N;
        end
    end % end private methods
    
    methods ( Static )
    end
end