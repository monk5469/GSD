function [DCM] = DCM_Connection_Model(dcm,options)
%------------------------------------------------------------
%This function is used to edit the connection model.
%This simulation edited three connection models and focused on only
%three stimuli (full stimulus, immediate repetition, delayed repetition).
%-----------------------------------------------------------------------
DCM                = dcm; 
DCM.U.u            = [dcm.U.u(:,1) dcm.U.u(:,3) dcm.U.u(:,4)];
DCM.U.name         = [dcm.U.name(1) dcm.U.name(3) dcm.U.name(4)];
DCM.b(:, :, [2,5]) = [ ];     %This simulation focuses on only three effects: 
                               %total stimulation, immediate repetition, and delayed repetition 
DCM.c(:,[4,5])     = [ ];                           
if options.connection_model == 1
    
    DCM.b(:,1,2) = 0;
    DCM.b(1,:,3) = 0;
    
end

if options.connection_model == 2
    
     DCM.b(2,:,3) = 0;
     DCM.b(1,1,3) = 0;
     DCM.b(1,1,2) = 0;
     DCM.b(2,:,2) = 0;
    
end

if options.connection_model == 3
          
       DCM.b(:,:,2) = 0;
       DCM.b(:,:,3) = 0;
          
end


end

