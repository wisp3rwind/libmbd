if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS_INIT
        "-fall-intrinsics -std=f2008ts -pedantic -Wall -Wcharacter-truncation -Wimplicit-procedure -Wextra -Wno-maybe-uninitialized -fstack-usage")
    set(CMAKE_Fortran_FLAGS_DEBUG_INIT
        "-Og -fcheck=all -fno-check-array-temporaries -fstack-usage")
endif()
