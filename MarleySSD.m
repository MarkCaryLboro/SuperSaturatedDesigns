classdef MarleySSD < CWPW
    % Generate a super saturated design using the algorithm by Marley

    properties ( Access = protected )
        Active          activeFactorSSD                                     % Class to handle subset generation
    end
    
    properties ( SetAccess = protected, Dependent )
        Measure     double                                                  % Design measure to be optimised
    end
    
    properties ( Constant = true )
        Algorithm       SSDalgorithm = "Marley"                             % Generating algorithm
    end
    
    properties 
        W   double      {mustBeGreaterThan(W, 0), mustBeLessThan(W, 1),...
                         mustBeReal(W), mustBeNonempty(W), ...
                         mustBeFinite(W)} = 0.5                             % Criterion weight factor
    end    
    
    methods
        % CONSTRUCTOR AND ORDINARY METHODS
        function obj = MarleySSD(M, F, varargin)
            %-----------------------------------------------------------------
            % Class constructor. Generate supersaturated designs.
            %
            % obj = MarleySSD( M, F, 'PARAM1', VALUE1, ..., 'PARAM#', VALUE#);
            %
            % Where M is the number of factors, F is the anticipated number
            % of active factors and 'PARAM#' is a STRING defining the
            % nature of the following VALUE#. Valid forms for 'PARAM#' are:
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
            obj = obj@CWPW( M, varargin{:});
            obj.Active = activeFactorSSD( obj.M_, F );
        end
    end % Constructor & ordinary methods
    
    methods
        % GET/SET METHODS
        function M = get.Measure( obj )
            % Return average VIF for subsets comprised of f+1 columns
            try
                M = obj.costFcn();
            catch
                M = NaN;
            end
        end
    end % End get/set methods
    
    methods ( Access = protected )
        % PROTECTED METHODS
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
            %--------------------------------------------------------------
            % Calculate the corresponding cost function
            %--------------------------------------------------------------
            V = 0;                                                          % Average VIF criterion
            A = 0;                                                          % Average A-optimal criterion
            for Q = 1:R
                Ds = obj.S( :,Cmb( Q,: ) );
                E = eig( Ds.'*Ds);
                L = obj.N_*sum( 1./E );
                V = V + L;
                E1 = [obj.N_; E];
                L = sum( 1./E1 );
                A = A + L;
            end
            V = V/R/obj.Active.F_;
            A = A/(obj.Active.F_ + 1)/R;
            L = obj.W*log( A ) + (1 - obj.W)*log( V );
        end
    end % End protected methods
    
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
            Ok = obj.Measure < OldCost;
        end
    end % End hidden methods
    
    methods ( Access = protected, Static = true )
        function J = sortDesigns( D )
            %--------------------------------------------------------------
            % Rank design measures from largest to smallest
            %
            % J = obj.sortDesigns( D );
            % J = BayesianSSD.sortDesigns( D );
            %
            % Input Arguments:
            %
            % D     --> Unsorted list of design measures
            %--------------------------------------------------------------
            [~, J] = sort( D, 'descend' );
        end
    end % End protected, static methods
    
    methods ( Access = private )
    end % End private methods
end

