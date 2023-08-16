classdef GMM
    
    properties
        m = [0 0 0];
        id = 0;
        lineage = 0;
        parent = -1;
        
        area = 0;
        meanIns = 0;
        maxIns = 0;
        
        svIdx = 0;
        dims = 3;
        splitScore = 3;
        scale = [1 1 1.4];

        nu = 1; beta = 1; alpha = 1;
        nuPrior = 1; betaPrior = 1; alphaPrior = 1;
        W = [0.0230,0,0,0,0.0230,0,0,0,0.0452];
        WPrior = [0.05555, 0, 0, 0, 0.05555, 0, 0, 0, 0.10888];
        mPrior = [0 0 0];
        distMRFPrior = 1;        
    end
    
    methods
        
        function obj = initializer(obj, m, id, lineage, parent, area, meanIns, maxIns)
            obj.m = m;
            obj.id = id;
            obj.svIdx = id;
            obj.mPrior = m;
            obj.lineage = lineage;
            obj.parent = parent;
            obj.area = area;
            obj.meanIns = meanIns;
            obj.maxIns = maxIns;
        end
        
       
    end
    
    
end