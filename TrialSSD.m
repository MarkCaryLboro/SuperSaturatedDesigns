classdef TrialSSD
    % Generate a Supersaturated Design
    
    properties ( Access = protected )
        SSDobj       {validateObj(SSDobj)}
    end
    
    properties ( SetAccess = protected)
        NumTrials   int16   {mustBePositive(NumTrials), mustBeLessThanOrEqual(NumTrials, 32767),...
                             mustBeReal(NumTrials), mustBeNonempty(NumTrials)} = 1000;
        DesignProperties  table
    end
    
    properties ( SetAccess = protected )
        BestSSD      {validateObj(BestSSD)}
    end
    
    properties( SetAccess = private, Dependent )
        M           int16                                                   % Number of factors
        Algorithm   SSDalgorithm                                            % Design algorithm
        N           int16                                                   % Design size
        Measure     double                                                  % Design measure for the best design
        Es2         double                                                  % E(s^2) criterion for the best design
        MaxAbsCorr  double                                                  % Max ||correlation|| for the best design
        MeanAbsCorr double                                                  % Average ||correlation|| for the best design
        No          int16                                                   % Number of non-orthogonal column pairings  
        NumOrthCol  int16                                                   % Number of orthogonal columns pairings
        C           double                                                  % "C" diagnostic
    end
    
    methods ( Abstract = true, Access = protected )
        Ok = isBest( obj )                                                   % Return true if new trial design is new best design
    end
    
    methods
        function obj = generateDesign( obj, NumTrials, varargin)
            %--------------------------------------------------------------
            % Generate a series of designs and select the one associated
            % with the best design measure.
            %
            % obj = obj.generateDesign( NumTrials, 'Param1', Value1,...
            %                            'Param#', Value#  );
            %
            % Where NumTrials is the number of design trials and must be in
            % the range: (1 <= NumTrials <= 32767). 
            %
            % 'PARAM?' is a STRING defining the nature of the following 
            % VALUE?. Valid forms for 'PARAM?' are:
            %
            % N             --> Design Size
            % Kexch         --> Number of columns to exchange per iteration
            % MaxIter       --> Maximum number of iterations per trial
            %--------------------------------------------------------------
            obj.NumTrials = int16(NumTrials);
            obj = obj.conductDesignTrials( varargin{:} );
        end
        
        function makeHistograms( obj, NumBins )
            %--------------------------------------------------------------
            % Make histograms of all the relevant design evaluation 
            % properties
            %
            % obj.makeHistograms( NumBins );
            %
            % Input Arguments:
            %
            % NumBins --> Number of bins for histogram. Default is 21.
            %--------------------------------------------------------------
            if nargin<2 || isempty( NumBins ) || ~isreal( NumBins) || isnan( NumBins ) || ~isfinite( NumBins )
                NumBins = 21;
            else
                NumBins = ceil( NumBins );
            end 
               
            Props = obj.DesignProperties.Properties.VariableNames;
            for Q = 1:numel( Props )
                figure;
                histogram( obj.DesignProperties.(Props{Q}), NumBins );
                grid on;
                xlabel( Props{ Q } );
                ylabel( 'Frequency [#]');
            end
        end
    end % End constructor and ordinary methods
    
    methods
        function N = get.N( obj )
            % retrieve the design size
            N = obj.BestSSD.N;
        end
        
        
        
        function Rho = get.MaxAbsCorr( obj )
            % retrieve max absolute correlation for the best design
            Rho = obj.BestSSD.MaxAbsCorr;
        end
        
        function M = get.M( obj )
            % Return number of factors
            M = obj.BestSSD.M;
        end
        
        function M = get.Measure( obj )
            % Return design measure for best design
            M = obj.BestSSD.Measure;
        end
        
        function Rho = get.MeanAbsCorr( obj )
            % return the average absolute correlation for the best design
            Rho = obj.BestSSD.MeanAbsCorr;
        end
        
        function E = get.Es2( obj )
            % Return the E(s^2) for the best design
            E = obj.BestSSD.Es2;
        end
        
        function No = get.No( obj )
            % Retrieve number of non-orthogonal columns for the best design
            No = obj.BestSSD.No;
        end
        
        function N = get.NumOrthCol( obj )
            % Return number of possible orthogonal combinations
            N = obj.BestSSD.NumOrthCol;
        end
        
        function C = get.C( obj )
            % Return the "c" diagnostic for the best design
            C = obj.BestSSD.C;
        end
        
        function A = get.Algorithm( obj )
            % Return generation algorithm
            A = obj.BestSSD.Algorithm;
        end
    end % End set/get methods
    
    methods ( Access = private )
        function obj = conductDesignTrials( obj, varargin )
            %--------------------------------------------------------------
            % Conduct the design trials and set the BestSSD property. Plots
            % histograms of all availabel design evaluation measures for
            % the "best" design.
            %
            % obj = obj.conductDesignTrials( obj, 'Param1', Value1,...
            %                            'Param#', Value# );
            %
            % 'PARAM?' is a STRING defining the nature of the following 
            % VALUE?. Valid forms for 'PARAM?' are:
            %
            % N             --> Design Size
            % Kexch         --> Number of columns to exchange per iteration
            % MaxIter       --> Maximum number of iterations per trial
            %--------------------------------------------------------------
            [DesignSize, NumColExch, MaxIterations] = ...
                obj.defineArgs( varargin{:} );
            %--------------------------------------------------------------
            % Perform the trials
            %--------------------------------------------------------------
            DesignMeasures = obj.makeTable();
            fprintf('\n====================================================');
            fprintf('\n              Making Trial Designs                  ');
            fprintf('\n====================================================');
            
            for Q = 1:obj.NumTrials
                %----------------------------------------------------------
                % Generate a trial design
                %----------------------------------------------------------
                obj.SSDobj = obj.SSDobj.designGenerator( DesignSize,...
                    NumColExch, MaxIterations, false );
                %----------------------------------------------------------
                % Update the table of evaluation measures
                %----------------------------------------------------------
                DesignMeasures = obj.updateTable( DesignMeasures, Q );
                %----------------------------------------------------------
                % Indicate whether the design is best
                %----------------------------------------------------------
                Exchange = obj.isBest();
                if isnan(obj.BestSSD.Measure) || Exchange
                    %------------------------------------------------------
                    % Update the best design if required
                    %------------------------------------------------------
                    obj.BestSSD = obj.SSDobj;
                    fprintf('\nIteration #%4.0f, Best Design Updated. Measure = %8.5f', Q, obj.BestSSD.Measure);
                else
                    fprintf('\nIteration #%4.0f, Best Design NOT Updated.', Q);
                end
            end
            fprintf('\n\n');
            obj.DesignProperties = DesignMeasures;
            %--------------------------------------------------------------
            % Make histograms of all the design evaluation properties
            %--------------------------------------------------------------
            obj.makeHistograms();
        end
        
        function DM = makeTable( obj, varargin )
            %--------------------------------------------------------------
            % Set up the table of design evaluation measures
            %
            % DM = obj.makeTable();
            %--------------------------------------------------------------
            Props = obj.SSDobj.getListOfProperties();
            DM = table( 'Size', [obj.NumTrials, numel(Props)], ...
                'VariableTypes', repmat("double", 1, numel(Props)));
            DM.Properties.RowNames = string( 1:obj.NumTrials );
            DM.Properties.VariableNames = Props;
        end
        
        function DM = updateTable( obj, DM, CurrentEntry )
            %--------------------------------------------------------------
            % Update the table of design measures for each trial
            %
            % DM = obj.updateTable( DM, EntryCounter );
            %
            % Input Arguments:
            %
            % DM            --> Table of design properties.
            % CurrentEntry  --> Entry to update
            %--------------------------------------------------------------
            Props = obj.SSDobj.getListOfProperties();
            for Q = 1:numel( Props )
                DM.( Props{ Q } )( CurrentEntry ) = obj.SSDobj.( Props{ Q } );
            end
        end
    end % end private and helper methods
    
    methods ( Static = true )
        function [DesignSize, NumColExch, MaxIterations] = defineArgs( varargin )
            %--------------------------------------------------------------
            % Parse the settings for the design generator
            %
            % [DesignSize, NumColExch, MaxIterations] = ...
            %       makeSSD.defineArgs( 'Param1', Value1,...
            %                           'Param#', Value# );
            %
            % [DesignSize, NumColExch, MaxIterations] = ...
            %       obj.defineArgs( 'Param1', Value1,...
            %                       'Param#', Value# );
            % 
            % 'PARAM?' is a STRING defining the nature of the following 
            % VALUE?. Valid forms for 'PARAM?' are:
            %
            % N             --> Design Size
            % Kexch         --> Number of columns to exchange per iteration
            % MaxIter       --> Maximum number of iterations per trial
            %--------------------------------------------------------------
            Props = varargin(1:2:end);
            Values = varargin(2:2:end);
            for Q = 1:numel(Props)
                switch lower( Props{ Q } )
                    case 'n'
                        DesignSize = Values{ Q };
                    case 'kexch'
                        NumColExch = Values{ Q };
                    case 'maxiter'
                        MaxIterations = Values{ Q };
                end
            end
        end
    end
end

%--------------------------------------------------------------------------
% Validation Function
%--------------------------------------------------------------------------
function validateObj( ObjSSD )
    %----------------------------------------------------------------------
    % Ensure the supplied object is of the correct type
    %----------------------------------------------------------------------
    if ~isempty(ObjSSD) && ~isa( ObjSSD, 'SSD')
        %------------------------------------------------------------------
        % Return an error
        %------------------------------------------------------------------
        error('Objects of Class "%s" are unsupported. Supported classes are "%s", "%s" or "%s"',...
            class( ObjSSD ), "BayesianSSD", "LiSSD", "MarleySSD");
    end
end