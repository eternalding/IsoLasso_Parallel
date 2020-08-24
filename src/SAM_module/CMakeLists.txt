project(SAM_functions)
aux_source_directory(. DIR_SAM_SRCS)
message(${DIR_SAM_SRCS})
add_library(SAM_functions ${DIR_SAM_SRCS})