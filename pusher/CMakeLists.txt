add_library(aphi-pusher SHARED
    pv.h
    pv.cc
    magfield.h
    magfield.cc
    efield.h
    efield.cc
    pusher.h
    pusher.cc
)

option(GAMMA_CORRECTION_APHI "Perform the gamma correction for the aphi pusher" ON)
if (GAMMA_CORRECTION_APHI)
target_compile_definitions(aphi-pusher PRIVATE GAMMA_CORRECTION_APHI)
endif (GAMMA_CORRECTION_APHI)

set_target_properties(aphi-pusher PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)