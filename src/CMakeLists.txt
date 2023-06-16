
include(GenerateExportHeader)

set(PLACENTA_EXPORTDEFINITIONS_H "${CMAKE_CURRENT_BINARY_DIR}/exportdefinitions.h")

set(exes combie continue pennati_02 plac_flow plot_combie_fetal plot_combie_placenta)
foreach(exe ${exes})
  add_library(${exe} SHARED ${exe}.c interface.h placentastdlib.c placentastdlib.h ${PLACENTA_EXPORTDEFINITIONS_H})
  add_executable(${exe}_exe main.c)
  target_link_libraries(${exe}_exe ${exe})
  target_include_directories(${exe} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
  generate_export_header(${exe} EXPORT_FILE_NAME ${PLACENTA_EXPORTDEFINITIONS_H} BASE_NAME PLACENTA)
endforeach()

add_subdirectory(bindings)
