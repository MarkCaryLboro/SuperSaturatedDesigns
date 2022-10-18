classdef BayesianSSD < CWPW
    % Generate a super saturated design using the algorithm by Jones et al
    
    properties ( Access = public )
        Tau             double  {mustBeNumeric(Tau), mustBeFinite(Tau),...
                                 mustBeReal(Tau),...
                                 mustBeGreaterThan( Tau, 1e-6)} = sqrt(5)   % Prior distribution variance
    end
    
    properties ( SetAccess = immutable, GetAccess = protected )
        ActiveFacs   logical                                                % Set each element to true to indicate factor is anticipated to be active
    end
    
    properties ( Constant = true )
        Algorithm       SSDalgorithm = "Jones"                              % Generating algorithm
    end
        
    properties( SetAccess = protected, Dependent)
        K               double                                              % Prior distribution covariance matrix
        F               double                                              % Information matrix
        Measure         double                                              % Determinant of information matrix
        V               double                                              % Variance covariance matrix
        NumScreenFacs                                                       % Number of factors to be screened
        NumActiveFacs                                                       % Number of assumed active factors
        Active                                                              % Abbreviation for the active factors
        Potential                                                           % Abbreviation for the potential factors
    end
    
    methods
        % CONSTRUCTOR AND ORDINARY METHODS
        function obj = BayesianSSD(M, varargin)
            %-----------------------------------------------------------------
            % Class constructor. Generate supersaturated designs.
            %
            % obj = BayesianSSD( M, 'PARAM1', VALUE1, ..., 'PARAM#', VALUE#);
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
            % Active       --> (Mx1) logical array. Set each element to
            %                  true if corresponding factor is expected to
            %                  be active from physical considerations.
            %
            % Note if a 'PARAM#' is not assigned, then suitable defaults
            % are provided. For example, the defaults for "Lo" and "Hi" are
            % -1 and +1 respectively. Defaults for Factor and Abbreviation
            % are the letters of the alphabet.
            %
            % Examples:
            %
            % obj = BayesianSSD( 10, "Hi", ones(10 , 1), "Lo", -ones(10, 1));
            %-----------------------------------------------------------------   
            if nargin<1
                % Apply default number of factors
                M = 10;
            end
            %--------------------------------------------------------------
            % Define the parent object
            %--------------------------------------------------------------
            obj = obj@CWPW( M, varargin);
            %--------------------------------------------------------------
            % Set the active terms if any
            %--------------------------------------------------------------
            P = strcmpi("Active", [varargin{1:2:end}]);
            if any(P)
                % set the active factors
                P = find(P, numel(P), 'first');
                obj.ActiveFacs = varargin{ 2*P };
            else
                % define the default - only screening factors
                obj.ActiveFacs = false( 1, obj.M_ );
            end
        end
    end % constuctor and ordinary methods
    
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
    
    methods ( Access = protected )
        function L = costFcn( obj, J )
            %--------------------------------------------------------------
            % Evaluate the cost function for the design matrix provided
            %
            % L = obj.costFcn( J );
            %
            % Input Arguments:
            %
            % J     --> Column to be delted from the current design
            %--------------------------------------------------------------
            if nargin<2
                J = NaN;
            end
            %--------------------------------------------------------------
            % Recompute the Prior information matrix
            %--------------------------------------------------------------
            R = obj.priorInformationMatrix( J );
            %--------------------------------------------------------------
            % Define the regression matrix
            %--------------------------------------------------------------
            X = obj.computeX( J );
            %--------------------------------------------------------------
            % Calculate the corresponding cost function
            %--------------------------------------------------------------
            L = det( X.'*X + R).^(1/obj.M_);
        end
    end % End protected methods
    
    methods ( Access = private )
        function R = priorInformationMatrix( obj, J )
            %--------------------------------------------------------------
            % Compute the corresponding information matrix for a missing
            % column
            %
            % R = obj.priorInformationMatrix( J );
            %
            % Input Arguments:
            %
            % J     --> Column to be deleted
            %--------------------------------------------------------------
            if nargin<2
                J = NaN;
            end
            P = obj.pointer2Columns( J );
            R = zeros( 1, obj.M_ );
            R( ~obj.ActiveFacs ) = 1;
            R = R( P );
            R = horzcat( 0, R );
            R = diag( R );
            R = R / obj.Tau^2;
        end
    end % private and helper methods
    
    methods
        % SET/GET Methods
        function S = get.NumScreenFacs( obj )
            % Return the number of factors to be screened
            S = obj.M_ - obj.NumActiveFacs;
        end
        
        function A = get.NumActiveFacs( obj )
            % Return number of active factors
            A = sum( obj.ActiveFacs);
        end
        
        function A = get.Active( obj )
            % Display active factor abbreviations
            A = obj.Abbreviation( obj.ActiveFacs );
        end
        
        function P = get.Potential( obj )
            % Display potential factor abbreviations
            P = obj.Abbreviation( ~obj.ActiveFacs );
        end
        
        function F = get.F( obj )
            % Return information matrix
            F = obj.X.'*obj.X + obj.K;
        end
        
        function V = get.V( obj )
            % Calculate the variance matrix from the information matrix
            V = obj.F\eye( obj.M_ + 1);
        end
        
        function DF = get.Measure( obj )
            % Return determinant of information matrix
            try
                DF = det( obj.F )^(1/obj.M_);
            catch
                DF = NaN;
            end
        end
      
        function K = get.K( obj )
            % Return the information matrix
            K = obj.priorInformationMatrix();
        end
    end % get/set methods
    
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
end

