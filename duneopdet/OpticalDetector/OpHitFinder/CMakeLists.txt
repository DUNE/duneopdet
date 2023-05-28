install_fhicl()
install_headers()
install_source()

# Also put a copy in the source tree.

FILE(GLOB fcl_files *.fcl)
FILE(GLOB text_files *.txt)
FILE(GLOB header_files *.h)
FILE(GLOB source_files *.cxx)
install_source( EXTRAS ${fcl_files} )
install_source( EXTRAS ${text_files} )
install_source( EXTRAS ${header_files} )
install_source( EXTRAS ${source_files} )