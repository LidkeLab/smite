function unitTest()

   obj = smi_sim.SimSMLM();
   obj.SZ = 256;
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
   
   [Model, Data] = obj.genImageStack();

   % Display the blobs without noise:
   imshow(Model)
   %dipshow(Model)
   
    % Display the blobs having Poisson noise:
   imshow(Data)
   %dipshow(Data)
    
end
