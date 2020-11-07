cmake_minimum_required(VERSION 3.16.5)
project (Constraiint_branching_for_SCP)

# add path to FindCPLEX.cmake that contains the directives to add the CPLEX library
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# search in the previously specified path for the file findCPLEX.cmake
FIND_PACKAGE(CPLEX REQUIRED)

# bring the headers into the project
#include_directories(include;/home/enrico/concorde;${CPLEX_INCLUDE_DIRS})
include_directories(include;${CPLEX_INCLUDE_DIRS})

# add debug information
set(CMAKE_BUILD_TYPE Debug)

# can manually add the sources using the set command as follows:
set(SOURCES "include/chrono.h"
            "include/scp.h"
            "include/utilities.h"
            "src/scoreComparison.c"
            "src/main.c"
            "src/chrono.c"
            "src/scp.c"
            "src/utilities.c"
            "src/preprocessing.c"
            "src/genericCallbacks.c"
            "src/constraintIntersections.c"
            
            
            )

# specify the executable name, followed by all the source files
add_executable( "main.o"
                "src/main.c"
                "src/scoreComparison.c"
                "src/chrono.c"
                "src/scp.c"
                "src/utilities.c"
                "src/preprocessing.c"
                "src/genericCallbacks.c"
                "src/constraintIntersections.c"
            
              )

# link CPLEX libraries
#TARGET_LINK_LIBRARIES(main.o;${CPLEX_LIBRARIES};/home/enrico/concorde/concorde.a)
TARGET_LINK_LIBRARIES(main.o;${CPLEX_LIBRARIES};)
