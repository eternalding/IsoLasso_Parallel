project(SAM_functions)
aux_source_directory(. DIR_SAM_SRCS)
message("Compiling SAM-Module ${DIR_SAM_SRCS}")
add_library(SAM_functions ${DIR_SAM_SRCS})