if (NOT USE_SM)
    set (tm_hydram
        ../sm/hydram.C
        )
endif ()

if (USE_CEMHYD)
    set (tm_cemhyd
        cemhyd/cemhydmat.C
        )
endif ()

set (tm_emodel 
    stationarytransportproblem.C
    nonstationarytransportproblem.C
    nltransienttransportproblem.C
    staggeredproblem.C
    darcyflow.C
    )

set (tm_element
    transportelement.C
    quad1_ht.C
    tr1_ht.C
    quadaxisym1_ht.C
    traxisym1_ht.C
    brick1_ht.C
    qbrick1_ht.C
    tetrah1_ht.C
    tr1darcy.C
    )

set (tm_status
    transportmaterial.C
    )

set (tm_material
    hydratingisoheatmat.C
    hydratinghemomat.C
    hemotkmat.C
    isoheatmat.C
    isomoisturemat.C
    isolinmoisturemat.C
    bazantnajjarmat.C
    nlisomoisturemat.C
    hydratingconcretemat.C
    nonlinearmassmat.C
    anisomassmat.C
    rvestokesflow.C
    )

set (tm
    ${tm_emodel}
    ${tm_element}
    ${tm_material}
    ${tm_status}
    ${tm_new}
    ${tm_cemhyd}
    ${tm_hydram}
    )

add_library (tm
    ${LIB_TYPE}
    ${tm}
    )