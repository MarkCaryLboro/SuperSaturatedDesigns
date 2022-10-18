classdef LiSSD < CWPW
    % Generate an SSD using the algorithm developed by Li and Wu
    properties ( Access = protected )
        Active          activeFactorSSD                                     % Class to handle subset generation
    end
    
    properties ( Constant = true )
        Algorithm       SSDalgorithm = "Li"                                 % Generating algorithm
    end
    
    properties ( SetAccess = protected, Dependent )
        Measure     double                                                  % Design measure
    end
    
    methods
        function obj = LiSSD(M, F, varargin)
            %-----------------------------------------------------------------
            % Class constructor. Generate supersaturated designs.
            %
            % obj = LiSSD( M, F, 'PARAM1', VALUE1, ..., 'PARAM#', VALUE#);
            %
            % Where M is the number of factors, F is the number of assumed
            % active factors, 'PARAM#' is a STRING defining the nature of
            % the following VALUE#. Valid forms for 'PARAM#' are:
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
            if nargin<1
                % Apply default number of factors
                M = 10;
            end
            if nargin<2
                F = 2;
            end
            obj = obj@CWPW( M, varargin{:});
            obj.Active = activeFactorSSD( obj.M_, F );
        end
    end % Constructor and ordinary methods
    
    methods
        function M = get.Measure( obj )
            % Return average determinant for submatrices of f+1 columns
            try
                M = obj.costFcn();
            catch
                M = NaN;
            end
        end
    end % end Set/Get methods
    
    methods ( Access = protected )
        function L = costFcn( obj, J )
            %--------------------------------------------------------------
            % Evaluate the cost function for the design matrix provided
            %
            % L = obj.costFcn( J );
            %
            % Input Arguments:
            %
            % J         --> Column to be delted from the current design
            %--------------------------------------------------------------
            if nargin<2 || isempty( J ) || (J < 1) || (J > obj.M_)
                J = NaN;
            end
            %--------------------------------------------------------------
            % Compute all combinations of columns of size F involving
            % column J
            %--------------------------------------------------------------
            [Cmb, R] = obj.Active.calcCombs( J );
            Factors = 1:obj.M_;
            %--------------------------------------------------------------
            % Calculate the corresponding cost function
            %--------------------------------------------------------------
            L = 0;
            for Q = 1:R
                P = ~ismember( Factors, Cmb( Q,:) );
                X = obj.computeX( Factors( P ) );
                E = eig(X.'*X)/obj.N_;
                L = L + prod(E).^(1/obj.Active.F_);
            end
        end
    end % end protected methods
    
    methods ( Hidden = true )
        function Ok = isBetter( obj, OldCost )
            %--------------------------------------------------------------
            % Return a logical variable  to indicate if the new design is 
            % better than the current design with resepct to the design 
            % measure.
            %
            % Ok = obj.isBetter( OldCost );
            %
            % Input Arguments:
            %
            % OldCost   --> Cost function for the current design
            %--------------------------------------------------------------
            Ok = obj.Measure > OldCost;
        end
    end % End hidden methods
    
    methods ( Access = protected, Static = true )
        function J = sortDesigns( D )
            %--------------------------------------------------------------
            % Rank design measures from smallest to largest
            %
            % J = obj.sortDesigns( D );
            % J = LiSSD.sortDesigns( D );
            %
            % Input Arguments:
            %
            % D     --> Unsorted list of design measures
            %--------------------------------------------------------------
            [~, J] = sort( D, 'ascend' );
        end
    end % End protected static methods  
end % End classdef