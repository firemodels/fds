function Uface_new = project_faces(Uface, Uface1, dt, FA_face_n, FB_face_n, dHdn_f, stage)
% Inputs: all Nf x 1, owner→neigh sign convention
% Output: Uface_new (Nf x 1)
if(stage==1)
  % Predictor:
  Uface_new = Uface - dt*( FA_face_n + FB_face_n + dHdn_f );
else
  % Corrector:
  Uface_new = 0.5*(Uface+Uface1) - dt/2*( FA_face_n + FB_face_n + dHdn_f );
end
