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
set_target_properties(aphi-pusher PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)