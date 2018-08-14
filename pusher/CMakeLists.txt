option(APHI_GAMMA_CORRECTION "Perform gamma correction for aphi pusher" ON)
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

if (APHI_GAMMA_CORRECTION)
target_compile_definitions(aphi-pusher PRIVATE APHI_GAMMA_CORRECTION)
endif (APHI_GAMMA_CORRECTION)
