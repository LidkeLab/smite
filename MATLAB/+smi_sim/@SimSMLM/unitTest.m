function unitTest()

   obj = smi_sim.SimSMLM();
   obj.SZ = 64;
   obj.Rho = 30;
   obj.NDatasets = 2;
   obj.NFrames = 10;
   obj.ZoomFactor = 2;
   obj.K_OnToOff = 1;
   obj.K_OffToOn = 0.005;
   obj.K_OnToBleach = 0.2;
   obj.EmissionRate = 1000;
   obj.Bg = 15;
   obj.PSFSigma = 1.3;
   obj.StartState = 'Equib';
   obj.LabelingEfficiency = 1;

   NWings = 16;
   obj.simStar(NWings);
   SMD_Data = obj.genNoisySMD(obj.SMD_Model);
   %[SMD_Model] = obj.genBlinks(SMD_True,'Equib'); 
   %[SMD_Data] = obj.genNoisySMD(SMD_Model);
   
   % To generate the blobs without noise, execute the following:
   [Model] = smi_sim.GaussBlobs.gaussBlobImage(obj.SZ,obj.NDatasets*obj.NFrames,obj.SMD_Model,0,0,0);
   dipshow(Model)
   
    % To generate the blobs having poisson noise, execute the following:
   %[Data] = obj.genNoisyData(Model);
   [Data] = smi_sim.GaussBlobs.gaussBlobImage(obj.SZ,obj.NDatasets*obj.NFrames,SMD_Data,obj.Bg,0,0);
   dipshow(Data)
    
end
