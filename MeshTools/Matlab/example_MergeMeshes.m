% MergeMeshes.m
% merges two SU2 meshes
% beware of large meshes!!
%

%% SETUP

meshname1   = 'mesh_bipara_1.su2';
meshname2   = 'mesh_bipara_2.su2';
meshnameout = 'mesh_bipara_merged.su2';


%% READ

mesh1 = ReadSU2(meshname1);
mesh2 = ReadSU2(meshname2);

figure(1); clf; hold on;
plotMarkers(mesh1,{mesh1.mark(:).tag},'b-');
plotMarkers(mesh2,{mesh2.mark(:).tag},'r-');


%% MERGE

meshout = MergeSU2(mesh1,mesh2);

figure(2); clf; hold on;
plotMarkers(meshout,'k-');


%% WRITE

WriteSU2(meshout,'mesh_bipara_merged.su2');


%% DONE


