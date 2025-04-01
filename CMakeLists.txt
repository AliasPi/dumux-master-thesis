dune_symlink_to_source_files(FILES "params.input")
dune_symlink_to_source_files(FILES "params3d.input")

dune_symlink_to_source_files(FILES grids)

dumux_add_test(NAME 2d
               LABELS custom co2 basic
               SOURCES main.cc boundarycellmanager.cc
               COMPILE_DEFINITIONS TYPETAG=Heterogeneous PUBLIC ISOTHERMAL=1)

dumux_add_test(NAME 2dni
               SOURCES main.cc boundarycellmanager.cc
               LABELS porousmediumflow co2 co2ni
               COMPILE_DEFINITIONS TYPETAG=HeterogeneousNI PUBLIC ISOTHERMAL=0)

dumux_add_test(NAME 3d
               LABELS custom co2 basic
               SOURCES main.cc boundarycellmanager.cc
               COMPILE_DEFINITIONS TYPETAG=Heterogeneous PUBLIC ISOTHERMAL=1 PUBLIC THREEDIMENSIONAL=1)

dumux_add_test(NAME 3dni
               SOURCES main.cc boundarycellmanager.cc
               LABELS porousmediumflow co2 co2ni
               COMPILE_DEFINITIONS TYPETAG=HeterogeneousNI PUBLIC ISOTHERMAL=0 PUBLIC THREEDIMENSIONAL=1)