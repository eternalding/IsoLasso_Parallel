project(SAM_functions)
aux_source_directory(. DIR_SAM_SRCS)
message("Compiling SAM-Module ${DIR_SAM_SRCS}")

find_package(Boost 1.71.0 REQUIRED COMPONENTS iostreams)
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )


add_library(SAM_functions ${DIR_SAM_SRCS})