classdef makeSSD < TrialSSD
    % Generate a Supersaturated Design
    
    methods
        function obj = makeSSD( SSDobj )
            %--------------------------------------------------------------
            % Class constructor
            %
            % obj = makeSSD( SSDobj );
            %
            % Input Arguments:
            %
            % SSDobj    --> A BayesianSSD, LiSSD or MarleySSD object.
            %--------------------------------------------------------------
            obj.SSDobj = SSDobj;
            obj.BestSSD = SSDobj;
        end
    end % End constructor and ordinary methods

    methods ( Access = protected )
        function Ok = isBest( obj )
            %--------------------------------------------------------------
            % Return true if the current meaure is better than the last
            % best known measure
            %
            % Ok = obj.isBest();
            %--------------------------------------------------------------
            Ok = obj.SSDobj.isBetter( obj.BestSSD.Measure );
        end
    end % End protected methods
end
