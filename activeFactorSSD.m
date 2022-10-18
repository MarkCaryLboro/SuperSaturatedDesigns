classdef activeFactorSSD
    % A class to handle active factors and calculating subset sizes
    
    properties ( SetAccess = protected )
        F   int16 {mustBeGreaterThan(F, 0), mustBeReal(F),...               % Assumed number of active factors
                   mustBeNumeric(F), mustBeFinite(F)} = int16(2)
        M   double                                                          % Number of factors       
    end

    properties ( SetAccess = protected, Dependent ) 
        F_  double                                                          % Assumed number of active factors as a double            
    end

    methods
        function obj = activeFactorSSD( M, F )
            %--------------------------------------------------------------
            % Class constructor
            %
            % obj = activeFactorSSD( M, F );
            %
            % Input Arguments:
            %
            % M     --> Total number of factors
            % F     --> Number of anticipated active factors
            %--------------------------------------------------------------
            obj.M = M;
            obj.F = F;
        end
        
        function [Cmb, NumCmb] = calcCombs( obj, J )
            %--------------------------------------------------------------
            % Return list of all combinations involving factor J of size
            % obj.F
            %
            % Cmb = obj.calcCombs( J );
            %
            % Input Arguments:
            %
            % J     --> Column to be deleted from the current design
            %
            % Output Arguments:
            %
            % Cmb       --> Full list of possible combinations
            % NumCmb    --> Number of combinations
            %--------------------------------------------------------------
            Cmb = combnk( (1:obj.M), obj.F_);
            if ~isnan( J )
                P = any( Cmb == J, 2);
                Cmb = Cmb( P,: );
            end
            if nargout == 2
                NumCmb = size( Cmb, 1 );
            end
        end
    end % End constructor & ordinary methods
    
    methods
        % GET/SET Methods
        function obj = set.F( obj, Value )
            % Set the number of active factors
            V = int16( Value );
            if V < 2
                V = 2;
            elseif V > floor( obj.M / 3 )                                   %#ok<MCSUP>
                V = floor( obj.M / 3 );                                     %#ok<MCSUP>
            end
            obj.F = V;
        end
        
        function F = get.F_( obj )
            % Return the number of active factors
            F = double( obj.F );
        end
     end % End get/set methods
end