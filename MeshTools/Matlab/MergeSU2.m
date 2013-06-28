function meshout = MergeSU2(mesh1,mesh2)
% meshout = MergeSU2(mesh1,mesh2)
% concatenates mesh2 information below mesh1

%% SETUP

assert( mesh1.dim == mesh2.dim , 'mesh dimension mismatch' );

meshout = mesh1;


%% POINTS

mesh2.poin(:,end) = mesh2.poin(:,end) + mesh1.npoin;

meshout.poin  = [ meshout.poin ; mesh2.poin ];
meshout.npoin = size( meshout.poin, 1 );


%% ELEMENTS

for ielem = 1:mesh2.nelem
  elem_leng  = getElemTypeInfo( mesh2.elem(ielem,1) );
  elem_range = 2:(1+elem_leng);
  mesh2.elem( ielem, elem_range )  = mesh1.npoin + mesh2.elem( ielem, elem_range );
  mesh2.elem( ielem, elem_leng+2 ) = mesh1.nelem + mesh2.elem( ielem, elem_leng+2 );
end

meshout.elem  = [ meshout.elem ; mesh2.elem ];
meshout.nelem = size( meshout.elem, 1 );


%% MARKERS

mesh1_marklist = { mesh1.mark(:).tag };

for imark = 1:mesh2.nmark;
  
  % Elements
  for ielem = 1:mesh2.mark(imark).nelem
    elem_leng  = getElemTypeInfo( mesh2.mark(imark).elem(ielem,1) );
    elem_range = 2:(1+elem_leng);
    mesh2.mark(imark).elem( ielem, elem_range )  = mesh1.npoin + mesh2.mark(imark).elem( ielem, elem_range );
  end
  
  jmark = find( strcmp( mesh2.mark(imark).tag, mesh1_marklist ) );
  njmark = length(jmark);
  
  % Existing Marker
  if njmark == 1
    
    meshout.mark(jmark).elem  = [ meshout.mark(jmark).elem ; mesh2.mark(imark).elem ];
    meshout.mark(jmark).nelem = size( meshout.mark(jmark).elem, 1 );
    
  % New Marker
  elseif njmark == 0
    
    jmark = length(meshout.mark)+1;
    
    meshout.mark(jmark).tag   = mesh2.mark(imark).tag;
    meshout.mark(jmark).elem  = mesh2.mark(imark).elem;
    meshout.mark(jmark).nelem = size( meshout.mark(jmark).elem, 1 );
    
  else
    error('duplicate markers?')
    
  end
  
end

meshout.nmark = length(meshout.mark);


%% DONE



